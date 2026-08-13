using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA.FeatureAlign
{
    /// <summary>How two proteoform identifications from different runs are decided to be "the same proteoform".</summary>
    public enum ProteoformMatchKey
    {
        /// <summary>Accession + the exact Full Sequence string (base sequence + every localized modification).
        /// Strictest: a PTM localized to a different residue in another run is a DIFFERENT proteoform.</summary>
        FullSequence,

        /// <summary>Accession + Base Sequence + the position-agnostic multiset of modification identities.
        /// Merges runs that localized the same PTM(s) to different residues (default).</summary>
        ModsIgnorePosition,

        /// <summary>Accession + Base Sequence only — ignores modifications entirely (collapses all proteoforms of a
        /// protein region to one row).</summary>
        BaseSequence
    }

    /// <summary>
    /// Aligns MetaMorpheus proteoform identifications ACROSS runs by IDENTITY (not by mass/RT tolerance the way
    /// <see cref="FeatureAligner"/> does). Each input is a MetaMorpheus proteoform/PSM .psmtsv
    /// (a per-run <c>*_Proteoforms.psmtsv</c>, or a pooled <c>AllProteoforms.psmtsv</c>/<c>AllPSMs.psmtsv</c> that
    /// carries a "File Name" column). Rows are grouped by a proteoform key (see <see cref="ProteoformMatchKey"/>)
    /// and collapsed to one row per proteoform with each run's retention time, precursor intensity and q-value —
    /// a proteoform × run presence/abundance matrix.
    ///
    /// One run = one input file by default. If a SINGLE file is passed it is instead split into runs by its
    /// "File Name" column (so a pooled AllPSMs.psmtsv works out of the box). Pass <c>splitByFileName</c> to force
    /// either behaviour. NOTE: a task-level AllProteoforms.psmtsv is deduplicated to one representative row per
    /// proteoform (tagged with its best-scoring run), so splitting one of those by File Name UNDER-counts sharing;
    /// use the per-run *_Proteoforms.psmtsv files (one per run) or AllPSMs.psmtsv for a true cross-run matrix.
    /// </summary>
    public static class ProteoformAligner
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

        // ---- column-name dictionaries (auto-detected, so column order/extra columns don't matter) --------------
        private static readonly string[] FileNameNames  = { "File Name" };
        private static readonly string[] RtNames        = { "Scan Retention Time", "Retention Time" };          // MINUTES
        private static readonly string[] MonoMassNames  = { "Monoisotopic Mass" };                              // theoretical
        private static readonly string[] IntensityNames = { "Precursor Intensity" };
        private static readonly string[] ChargeNames    = { "Precursor Charge" };
        private static readonly string[] ScoreNames     = { "Score" };
        private static readonly string[] BaseSeqNames   = { "Base Sequence" };
        private static readonly string[] FullSeqNames   = { "Full Sequence" };
        private static readonly string[] AccNames       = { "Accession", "Protein Accession" };
        private static readonly string[] GeneNames      = { "Gene Name" };
        private static readonly string[] NameNames      = { "Name" };
        private static readonly string[] QValNames      = { "QValue" };
        private static readonly string[] DctNames       = { "Decoy/Contaminant/Target" };   // single letter T/D/C
        private static readonly string[] DecoyNames     = { "Decoy" };
        private static readonly string[] ContamNames    = { "Contaminant" };

        // ====================================================================================================
        //  ENTRY POINTS
        // ====================================================================================================

        /// <summary>
        /// Reads the proteoform .psmtsv inputs, groups identifications across runs by <paramref name="key"/>,
        /// writes the proteoform × run matrix, and returns the output path.
        /// </summary>
        /// <param name="inputFiles">Proteoform/PSM .psmtsv files. One file = one run (unless a single file is
        /// split by its File Name column — see <paramref name="splitByFileName"/>).</param>
        /// <param name="outputPath">tsv to write.</param>
        /// <param name="key">What makes two IDs the same proteoform (default: mods ignoring position).</param>
        /// <param name="qValueThreshold">Drop identifications with QValue above this (default 0.01). Rows with no
        /// QValue column are kept.</param>
        /// <param name="includeDecoys">Keep decoy rows (default false).</param>
        /// <param name="includeContaminants">Keep contaminant rows (default false).</param>
        /// <param name="splitByFileName">null = auto (split a single input by its File Name column; treat each of
        /// several inputs as its own run). true/false forces the behaviour.</param>
        public static string Align(IEnumerable<string> inputFiles, string outputPath,
            ProteoformMatchKey key = ProteoformMatchKey.ModsIgnorePosition,
            double qValueThreshold = 0.01, bool includeDecoys = false, bool includeContaminants = false,
            bool? splitByFileName = null)
        {
            if (inputFiles == null) throw new ArgumentNullException(nameof(inputFiles));
            if (string.IsNullOrWhiteSpace(outputPath)) throw new ArgumentException("OutputPath is not set.");

            var files = inputFiles.Where(File.Exists).ToList();
            if (files.Count == 0) throw new ArgumentException("No existing input proteoform files.");

            // Decide run assignment. Default: each file is a run. If exactly one file is given, split it into runs
            // by its "File Name" column so a pooled AllPSMs/AllProteoforms works without extra arguments.
            bool split = splitByFileName ?? (files.Count == 1);

            // Read every identification, tagged with the run label it belongs to.
            var recs = new List<Rec>();
            var runOrder = new List<string>();                 // run labels in first-seen order (stable output columns)
            var seenRun = new HashSet<string>(StringComparer.Ordinal);
            void NoteRun(string r) { if (seenRun.Add(r)) runOrder.Add(r); }

            if (split)
            {
                // Split whatever files were given by their internal File Name column.
                foreach (var file in files)
                    ReadRows(file, useFileNameAsRun: true, fixedRunLabel: null, key,
                             qValueThreshold, includeDecoys, includeContaminants, recs, NoteRun);
            }
            else
            {
                // One run per file; label = unique base name.
                var labels = MakeUniqueLabels(files);
                for (int i = 0; i < files.Count; i++)
                    ReadRows(files[i], useFileNameAsRun: false, fixedRunLabel: labels[i], key,
                             qValueThreshold, includeDecoys, includeContaminants, recs, NoteRun);
            }

            if (recs.Count == 0)
                throw new InvalidOperationException("No identifications passed the filters (QValue/target). Nothing to align.");

            WriteMatrix(recs, runOrder, key, outputPath);
            return outputPath;
        }

        // ====================================================================================================
        //  READING
        // ====================================================================================================

        // Read one proteoform .psmtsv into Rec records (one per surviving identification), tagging each with its
        // run label and its precomputed proteoform group key. Rows failing the QValue / target filters are dropped.
        private static void ReadRows(string path, bool useFileNameAsRun, string fixedRunLabel, ProteoformMatchKey key,
            double qValueThreshold, bool includeDecoys, bool includeContaminants,
            List<Rec> outRecs, Action<string> noteRun)
        {
            var lines = ReadAllLinesShared(path);
            if (lines.Length < 2) return;

            var header = lines[0].Split('\t');
            int cFile   = FindColumn(header, FileNameNames);
            int cRt     = FindColumn(header, RtNames);
            int cMono   = FindColumn(header, MonoMassNames);
            int cInt    = FindColumn(header, IntensityNames);
            int cCharge = FindColumn(header, ChargeNames);
            int cScore  = FindColumn(header, ScoreNames);
            int cBase   = FindColumn(header, BaseSeqNames);
            int cFull   = FindColumn(header, FullSeqNames);
            int cAcc    = FindColumn(header, AccNames);
            int cGene   = FindColumn(header, GeneNames);
            int cName   = FindColumn(header, NameNames);
            int cQ      = FindColumn(header, QValNames);
            int cDct    = FindColumn(header, DctNames);
            int cDecoy  = FindColumn(header, DecoyNames);
            int cContam = FindColumn(header, ContamNames);

            if (cBase < 0 || cFull < 0 || cAcc < 0)
                throw new FormatException(
                    $"'{Path.GetFileName(path)}' is not a MetaMorpheus proteoform .psmtsv (missing " +
                    $"{(cBase < 0 ? "Base Sequence " : "")}{(cFull < 0 ? "Full Sequence " : "")}{(cAcc < 0 ? "Accession" : "")}). " +
                    $"Columns found: {string.Join(", ", header)}");
            if (useFileNameAsRun && cFile < 0)
                throw new FormatException(
                    $"'{Path.GetFileName(path)}' has no \"File Name\" column, so it cannot be split into runs. " +
                    $"Pass one file per run instead.");

            for (int r = 1; r < lines.Length; r++)
            {
                if (lines[r].Length == 0) continue;
                var t = lines[r].Split('\t');

                // target / decoy / contaminant filtering
                if (!includeDecoys && IsDecoy(t, cDct, cDecoy)) continue;
                if (!includeContaminants && IsContaminant(t, cDct, cContam)) continue;

                // FDR filter (only applied when a QValue column exists and parses)
                double q = double.NaN;
                if (cQ >= 0 && cQ < t.Length) double.TryParse(t[cQ], NumberStyles.Any, Ci, out q);
                if (!double.IsNaN(q) && q > qValueThreshold) continue;

                string baseSeq = Cell(t, cBase);
                string fullSeq = Cell(t, cFull);
                string acc     = Cell(t, cAcc);
                if (baseSeq.Length == 0 || acc.Length == 0) continue;

                string modsCanon = CanonicalMods(fullSeq);                 // sorted, position-stripped mod list
                string run = useFileNameAsRun ? Cell(t, cFile) : fixedRunLabel;
                if (string.IsNullOrEmpty(run)) run = "(unknown)";
                noteRun(run);

                outRecs.Add(new Rec
                {
                    Run       = run,
                    GroupKey  = BuildKey(key, acc, baseSeq, fullSeq, modsCanon),
                    Accession = acc,
                    Gene      = Cell(t, cGene),
                    Name      = Cell(t, cName),
                    BaseSeq   = baseSeq,
                    ModsCanon = modsCanon,
                    MonoMass  = ParseDouble(t, cMono),
                    Rt        = ParseDouble(t, cRt),
                    Intensity = ParseDouble(t, cInt),
                    Charge    = ParseInt(t, cCharge),
                    Score     = ParseDouble(t, cScore),
                    QValue    = double.IsNaN(q) ? 0 : q
                });
            }
        }

        // ====================================================================================================
        //  PROTEOFORM KEY
        // ====================================================================================================

        // Build the grouping key for the chosen strategy. Accession is always part of the key (a proteoform is
        // defined per protein). The U+0001 control char is the field delimiter: it cannot occur in a
        // sequence/accession/mod string, so keys can't collide across field boundaries.
        private static string BuildKey(ProteoformMatchKey key, string acc, string baseSeq, string fullSeq, string modsCanon)
        {
            switch (key)
            {
                case ProteoformMatchKey.FullSequence:      return acc + "\u0001" + fullSeq;
                case ProteoformMatchKey.BaseSequence:      return acc + "\u0001" + baseSeq;
                case ProteoformMatchKey.ModsIgnorePosition:
                default:                                   return acc + "\u0001" + baseSeq + "\u0001" + modsCanon;
            }
        }

        // Position-agnostic canonical modification string: extract every [Category:ModName on Residue] annotation
        // from the Full Sequence, drop the leading "Category:" prefix, sort, and join. Two runs that localize the
        // same PTM(s) to different residues therefore produce the SAME string (the position lives in WHERE the
        // bracket sat, which we discard). Returns "" for an unmodified proteoform.
        internal static string CanonicalMods(string fullSeq)
        {
            var mods = ExtractModAnnotations(fullSeq);
            if (mods.Count == 0) return "";
            mods.Sort(StringComparer.Ordinal);              // multiset order-independence
            return string.Join("; ", mods);
        }

        // Extract the text inside each top-level [...] of a MetaMorpheus Full Sequence, with the "Category:" prefix
        // removed (e.g. "[Common Fixed:Carbamidomethyl on C]" -> "Carbamidomethyl on C"). Bracket depth is tracked
        // so a mod annotation that itself contains brackets is captured whole rather than split.
        private static List<string> ExtractModAnnotations(string fullSeq)
        {
            var mods = new List<string>();
            if (string.IsNullOrEmpty(fullSeq)) return mods;
            int depth = 0, start = -1;
            for (int i = 0; i < fullSeq.Length; i++)
            {
                char c = fullSeq[i];
                if (c == '[') { if (depth == 0) start = i + 1; depth++; }
                else if (c == ']') { if (depth > 0) { depth--; if (depth == 0 && start >= 0) { mods.Add(StripCategory(fullSeq.Substring(start, i - start))); start = -1; } } }
            }
            return mods;
        }

        // Drop the leading "Category:" (everything up to and including the first colon) so the token matches the
        // "ModName on Residue" form the Mods column uses and stays stable across searches with different mod
        // category labels. If there is no colon the token is returned unchanged.
        private static string StripCategory(string inner)
        {
            int colon = inner.IndexOf(':');
            return (colon >= 0 ? inner.Substring(colon + 1) : inner).Trim();
        }

        // ====================================================================================================
        //  WRITING THE PROTEOFORM x RUN MATRIX
        // ====================================================================================================

        private static void WriteMatrix(List<Rec> recs, List<string> runOrder, ProteoformMatchKey key, string outputPath)
        {
            int runCount = runOrder.Count;
            var runIndex = new Dictionary<string, int>(StringComparer.Ordinal);
            for (int i = 0; i < runCount; i++) runIndex[runOrder[i]] = i;

            // Group all identifications by proteoform key.
            var rows = new List<OutRow>();
            foreach (var grp in recs.GroupBy(x => x.GroupKey))
            {
                var members = grp.ToList();
                // Representative = best identification overall (lowest QValue, then highest Score) for the summary cells.
                var rep = members.OrderBy(x => x.QValue).ThenByDescending(x => x.Score).First();

                var rtByRun    = new double?[runCount];
                var intByRun   = new double?[runCount];
                var qByRun     = new double?[runCount];
                int nRuns = 0;

                // Within each run, several rows can map to the same proteoform (e.g. two localization variants that
                // collapse under ModsIgnorePosition). Keep that run's best row (lowest QValue, then highest Score).
                foreach (var runGrp in members.GroupBy(x => x.Run))
                {
                    int ri = runIndex[runGrp.Key];
                    var best = runGrp.OrderBy(x => x.QValue).ThenByDescending(x => x.Score).First();
                    rtByRun[ri]  = best.Rt;
                    intByRun[ri] = best.Intensity;
                    qByRun[ri]   = best.QValue;
                    nRuns++;
                }

                rows.Add(new OutRow
                {
                    Accession = rep.Accession,
                    Gene      = rep.Gene,
                    Name      = rep.Name,
                    BaseSeq   = rep.BaseSeq,
                    ModsCanon = rep.ModsCanon,
                    MonoMass  = rep.MonoMass,
                    Length    = rep.BaseSeq.Length,
                    NRuns     = nRuns,
                    TotalIds  = members.Count,
                    BestScore = members.Max(x => x.Score),
                    BestQ     = members.Min(x => x.QValue),
                    Rt        = rtByRun,
                    Intensity = intByRun,
                    QValue    = qByRun
                });
            }

            // Most-reproducible proteoforms first, then most abundant (by best q, then accession) for stable ties.
            rows.Sort((a, b) =>
            {
                int c = b.NRuns.CompareTo(a.NRuns);
                if (c != 0) return c;
                c = a.BestQ.CompareTo(b.BestQ);
                if (c != 0) return c;
                return string.CompareOrdinal(a.Accession, b.Accession);
            });

            var dir = Path.GetDirectoryName(outputPath);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir)) Directory.CreateDirectory(dir);

            using var sw = new StreamWriter(outputPath, append: false);

            // Header: summary columns, then a per-run block for intensity, RT (minutes) and QValue.
            var header = new List<string>
            {
                "accession", "gene", "name", "base_sequence", "mods", "match_key",
                "mono_mass", "length", "n_runs", "total_ids", "best_score", "best_qvalue"
            };
            header.AddRange(runOrder.Select(l => "intensity_" + l));
            header.AddRange(runOrder.Select(l => "rt_" + l));
            header.AddRange(runOrder.Select(l => "qvalue_" + l));
            sw.WriteLine(string.Join("\t", header));

            string keyName = key.ToString();
            foreach (var row in rows)
            {
                var cells = new List<string>
                {
                    row.Accession, row.Gene, row.Name, row.BaseSeq, row.ModsCanon, keyName,
                    row.MonoMass.ToString(Ci), row.Length.ToString(Ci), row.NRuns.ToString(Ci),
                    row.TotalIds.ToString(Ci), row.BestScore.ToString(Ci), row.BestQ.ToString(Ci)
                };
                cells.AddRange(row.Intensity.Select(v => v.HasValue ? v.Value.ToString(Ci) : ""));
                cells.AddRange(row.Rt.Select(v => v.HasValue ? v.Value.ToString(Ci) : ""));
                cells.AddRange(row.QValue.Select(v => v.HasValue ? v.Value.ToString(Ci) : ""));
                sw.WriteLine(string.Join("\t", cells));
            }
        }

        // ====================================================================================================
        //  SHARED HELPERS
        // ====================================================================================================

        // True if the row is a decoy (via the Decoy/Contaminant/Target letter "D", or a "Decoy" == true column).
        private static bool IsDecoy(string[] t, int cDct, int cDecoy)
        {
            if (cDct >= 0 && cDct < t.Length && t[cDct].Trim().Equals("D", StringComparison.OrdinalIgnoreCase)) return true;
            if (cDecoy >= 0 && cDecoy < t.Length && t[cDecoy].Trim().Equals("Y", StringComparison.OrdinalIgnoreCase)) return true;
            if (cDecoy >= 0 && cDecoy < t.Length && bool.TryParse(t[cDecoy].Trim(), out bool d) && d) return true;
            return false;
        }

        // True if the row is a contaminant (via the Decoy/Contaminant/Target letter "C", or a "Contaminant" column).
        private static bool IsContaminant(string[] t, int cDct, int cContam)
        {
            if (cDct >= 0 && cDct < t.Length && t[cDct].Trim().Equals("C", StringComparison.OrdinalIgnoreCase)) return true;
            if (cContam >= 0 && cContam < t.Length && t[cContam].Trim().Equals("Y", StringComparison.OrdinalIgnoreCase)) return true;
            if (cContam >= 0 && cContam < t.Length && bool.TryParse(t[cContam].Trim(), out bool cc) && cc) return true;
            return false;
        }

        // Read all lines with a SHARED handle so a file open in Excel can still be read.
        private static string[] ReadAllLinesShared(string path)
        {
            var lines = new List<string>();
            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            using var reader = new StreamReader(fs);
            string line;
            while ((line = reader.ReadLine()) != null) lines.Add(line);
            return lines.ToArray();
        }

        // Index of the first candidate name present in the header (case-insensitive, trimmed); -1 if none.
        private static int FindColumn(string[] header, string[] candidates)
        {
            foreach (var c in candidates)
            {
                int idx = Array.FindIndex(header, h => string.Equals(h.Trim(), c, StringComparison.OrdinalIgnoreCase));
                if (idx >= 0) return idx;
            }
            return -1;
        }

        // A short, unique label per file from its base name; duplicates get a ".1", ".2"… suffix.
        private static List<string> MakeUniqueLabels(List<string> files)
        {
            var labels = new List<string>();
            var seen = new Dictionary<string, int>();
            foreach (var f in files)
            {
                string name = Path.GetFileNameWithoutExtension(f);
                if (seen.TryGetValue(name, out int count)) { seen[name] = count + 1; labels.Add($"{name}.{count}"); }
                else { seen[name] = 1; labels.Add(name); }
            }
            return labels;
        }

        private static string Cell(string[] t, int col) => (col >= 0 && col < t.Length) ? t[col].Trim() : "";

        private static double ParseDouble(string[] t, int col) =>
            (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? v : double.NaN;

        // Parses via double so charges written with decimals (MetaMorpheus "11.00000") still work.
        private static int ParseInt(string[] t, int col) =>
            (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? (int)Math.Round(v) : 0;

        // ---- record types (private) ------------------------------------------------------------------------

        private sealed class Rec
        {
            public string Run;         // run label this identification belongs to
            public string GroupKey;    // precomputed proteoform grouping key
            public string Accession;
            public string Gene;
            public string Name;
            public string BaseSeq;
            public string ModsCanon;   // sorted, position-stripped mod list ("" if unmodified)
            public double MonoMass;    // theoretical monoisotopic mass
            public double Rt;          // scan retention time (MINUTES, as MetaMorpheus reports)
            public double Intensity;   // precursor intensity
            public int    Charge;      // precursor charge
            public double Score;       // MetaMorpheus score (higher is better)
            public double QValue;      // lower is better
        }

        private sealed class OutRow
        {
            public string Accession, Gene, Name, BaseSeq, ModsCanon;
            public double MonoMass, BestScore, BestQ;
            public int Length, NRuns, TotalIds;
            public double?[] Rt, Intensity, QValue;   // one slot per run; null = that run did not identify it
        }
    }
}
