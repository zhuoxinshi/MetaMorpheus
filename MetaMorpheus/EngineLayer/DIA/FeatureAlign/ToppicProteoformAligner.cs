using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA.FeatureAlign
{
    /// <summary>How two TopPIC proteoform identifications from different runs are treated as the same proteoform.</summary>
    public enum ToppicMatchKey
    {
        /// <summary>Protein accession + proteoform monoisotopic mass within a ppm tolerance. Same protein + same total
        /// mass = same proteoform, regardless of where TopPIC localized the mass shifts (analytic analog of
        /// "mods ignore position"). Different modifications of a protein differ in mass and stay separate.</summary>
        MassPpm,

        /// <summary>Protein accession + TopPIC's exact "Proteoform" string (sequence + localized mass shifts).
        /// Strict, position-aware: the same shift localized to a different residue is a different proteoform.</summary>
        ExactProteoform
    }

    /// <summary>
    /// Aligns TopPIC proteoform identifications ACROSS runs and writes a proteoform × run matrix (per-run feature
    /// intensity, apex RT and E-value) so you can count how many proteoforms are shared. Reads TopPIC's
    /// <c>*_ms2_toppic_proteoform_single.tsv</c> files (one per run) — these have one representative protein per
    /// proteoform and are already filtered at the search's E-value cutoff, so every row is taken.
    ///
    /// TopPIC differs from MetaMorpheus in two ways this handles: the accession is wrapped as <c>db|ACC|NAME</c>
    /// (parsed down to ACC), and modifications are localized mass shifts inside the proteoform string
    /// (e.g. <c>(MPPKD)[-103.0109]DKKK…</c>) rather than named mods — so the two grouping options are mass-based
    /// (<see cref="ToppicMatchKey.MassPpm"/>) or exact-string (<see cref="ToppicMatchKey.ExactProteoform"/>).
    /// A TopPIC tsv begins with a parameter header block; the reader skips to the "Data file name" table header.
    /// </summary>
    public static class ToppicProteoformAligner
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

        // ---- TopPIC proteoform.tsv column names (auto-detected) --------------------------------------------
        private static readonly string[] DataFileNames = { "Data file name" };
        private static readonly string[] RtNames       = { "Feature apex time", "Retention time" };            // MINUTES
        private static readonly string[] MassNames     = { "Proteoform mass", "Adjusted precursor mass", "Precursor mass" };
        private static readonly string[] IntensityNames= { "Feature intensity", "Proteoform intensity" };
        private static readonly string[] AccNames      = { "Protein accession" };
        private static readonly string[] DescrNames    = { "Protein description" };
        private static readonly string[] ProteoformNames = { "Proteoform" };
        private static readonly string[] EvalueNames   = { "E-value" };
        private static readonly string[] ChargeNames   = { "Charge" };

        // ====================================================================================================
        //  ENTRY POINT
        // ====================================================================================================

        /// <summary>
        /// Reads the TopPIC proteoform_single tsv inputs (one run each), groups proteoforms across runs by
        /// <paramref name="key"/>, writes the proteoform × run matrix, and returns the output path.
        /// </summary>
        /// <param name="inputFiles">TopPIC <c>*_toppic_proteoform_single.tsv</c> files, one per run.</param>
        /// <param name="outputPath">tsv to write.</param>
        /// <param name="key">MassPpm (accession + mass within tolerance) or ExactProteoform (accession + string).</param>
        /// <param name="massTolerancePpm">ppm tolerance for the mass match (MassPpm only). Default 10 (the search's tol).</param>
        /// <param name="rtToleranceMinutes">If &gt; 0, also require apex RTs to agree within this many minutes
        /// (MassPpm only). Default 0 = ignore RT (TopPIC apex time is missing/0.001 for feature-less proteoforms).</param>
        public static string Align(IEnumerable<string> inputFiles, string outputPath,
            ToppicMatchKey key = ToppicMatchKey.MassPpm,
            double massTolerancePpm = 10, double rtToleranceMinutes = 0)
        {
            if (inputFiles == null) throw new ArgumentNullException(nameof(inputFiles));
            if (string.IsNullOrWhiteSpace(outputPath)) throw new ArgumentException("OutputPath is not set.");

            var files = inputFiles.Where(File.Exists).ToList();
            if (files.Count == 0) throw new ArgumentException("No existing input TopPIC proteoform files.");

            var labels = MakeUniqueLabels(files);
            var recs = new List<Rec>();
            for (int i = 0; i < files.Count; i++)
                recs.AddRange(ReadRows(files[i], i));

            if (recs.Count == 0)
                throw new InvalidOperationException("No proteoform rows read from the inputs. Nothing to align.");

            var groups = key == ToppicMatchKey.ExactProteoform
                ? GroupExact(recs)
                : GroupByMass(recs, massTolerancePpm, rtToleranceMinutes);

            WriteMatrix(groups, labels, key, outputPath);
            return outputPath;
        }

        // ====================================================================================================
        //  GROUPING
        // ====================================================================================================

        // Exact: one group per (accession + proteoform string). Position-aware strict identity.
        private static List<List<Rec>> GroupExact(List<Rec> recs)
        {
            var groups = new Dictionary<string, List<Rec>>(StringComparer.Ordinal);
            foreach (var r in recs)
            {
                string k = r.Accession + "\u0001" + r.Proteoform;
                if (!groups.TryGetValue(k, out var list)) { list = new List<Rec>(); groups[k] = list; }
                list.Add(r);
            }
            return groups.Values.ToList();
        }

        // MassPpm: within each accession, cluster proteoforms whose masses agree within the ppm tolerance (and,
        // if rtTol > 0, whose apex RTs also agree). Union-find over mass-sorted records gives transitive clusters.
        private static List<List<Rec>> GroupByMass(List<Rec> recs, double ppm, double rtTolMinutes)
        {
            var tol = new PpmTolerance(ppm);
            var groups = new List<List<Rec>>();

            foreach (var byAcc in recs.GroupBy(r => r.Accession))
            {
                var arr = byAcc.OrderBy(r => r.Mass).ToList();
                int n = arr.Count;
                var parent = new int[n];
                for (int i = 0; i < n; i++) parent[i] = i;

                for (int i = 0; i < n; i++)
                {
                    double upper = tol.GetMaximumValue(arr[i].Mass);
                    for (int j = i + 1; j < n; j++)
                    {
                        if (arr[j].Mass > upper) break;                       // sorted by mass -> stop early
                        bool rtOk = rtTolMinutes <= 0 || (HasRt(arr[i]) && HasRt(arr[j])
                                        ? Math.Abs(arr[i].Rt - arr[j].Rt) <= rtTolMinutes : true);
                        if (tol.Within(arr[j].Mass, arr[i].Mass) && rtOk)
                            Union(parent, i, j);
                    }
                }

                var clusters = new Dictionary<int, List<Rec>>();
                for (int i = 0; i < n; i++)
                {
                    int root = Find(parent, i);
                    if (!clusters.TryGetValue(root, out var list)) { list = new List<Rec>(); clusters[root] = list; }
                    list.Add(arr[i]);
                }
                groups.AddRange(clusters.Values);
            }
            return groups;
        }

        private static bool HasRt(Rec r) => !double.IsNaN(r.Rt) && r.Rt > 0.0011;   // 0.001 is TopPIC's "no feature" placeholder

        // ====================================================================================================
        //  WRITING
        // ====================================================================================================

        private static void WriteMatrix(List<List<Rec>> groups, List<string> labels, ToppicMatchKey key, string outputPath)
        {
            int runCount = labels.Count;

            var rows = new List<OutRow>();
            foreach (var g in groups)
            {
                var rep = g.OrderBy(x => x.Evalue).First();   // best (lowest E-value) representative for the summary cells

                var intByRun = new double?[runCount];
                var rtByRun  = new double?[runCount];
                var evByRun  = new double?[runCount];
                int nRuns = 0;
                foreach (var runGrp in g.GroupBy(x => x.FileIndex))
                {
                    var best = runGrp.OrderBy(x => x.Evalue).First();
                    intByRun[runGrp.Key] = best.Intensity;
                    rtByRun[runGrp.Key]  = HasRt(best) ? best.Rt : (double?)null;
                    evByRun[runGrp.Key]  = best.Evalue;
                    nRuns++;
                }

                rows.Add(new OutRow
                {
                    Accession   = rep.Accession,
                    Species     = rep.Species,
                    Gene        = rep.Gene,
                    Description = rep.Description,
                    Proteoform  = rep.Proteoform,
                    MassMean    = g.Average(x => x.Mass),
                    NRuns       = nRuns,
                    TotalIds    = g.Count,
                    BestEvalue  = g.Min(x => x.Evalue),
                    Intensity   = intByRun,
                    Rt          = rtByRun,
                    Evalue      = evByRun
                });
            }

            // Most-shared first, then lowest E-value, then accession — stable, informative ordering.
            rows.Sort((a, b) =>
            {
                int c = b.NRuns.CompareTo(a.NRuns);
                if (c != 0) return c;
                c = a.BestEvalue.CompareTo(b.BestEvalue);
                if (c != 0) return c;
                return string.CompareOrdinal(a.Accession, b.Accession);
            });

            var dir = Path.GetDirectoryName(outputPath);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir)) Directory.CreateDirectory(dir);

            using var sw = new StreamWriter(outputPath, append: false);
            var header = new List<string>
            {
                "accession", "species", "gene", "description", "proteoform_example", "match_key",
                "mass_mean", "n_runs", "total_ids", "best_evalue"
            };
            header.AddRange(labels.Select(l => "intensity_" + l));
            header.AddRange(labels.Select(l => "rt_" + l));
            header.AddRange(labels.Select(l => "evalue_" + l));
            sw.WriteLine(string.Join("\t", header));

            string keyName = key.ToString();
            foreach (var row in rows)
            {
                var cells = new List<string>
                {
                    row.Accession, row.Species, row.Gene, row.Description, row.Proteoform, keyName,
                    row.MassMean.ToString("F5", Ci), row.NRuns.ToString(Ci), row.TotalIds.ToString(Ci),
                    row.BestEvalue.ToString(Ci)
                };
                cells.AddRange(row.Intensity.Select(v => v.HasValue ? v.Value.ToString(Ci) : ""));
                cells.AddRange(row.Rt.Select(v => v.HasValue ? v.Value.ToString(Ci) : ""));
                cells.AddRange(row.Evalue.Select(v => v.HasValue ? v.Value.ToString(Ci) : ""));
                sw.WriteLine(string.Join("\t", cells));
            }
        }

        // ====================================================================================================
        //  READING
        // ====================================================================================================

        // Read one TopPIC proteoform_single tsv into Rec records. Skips the leading parameter block to the table
        // header ("Data file name"), then parses each proteoform row. rows with an unparseable mass are dropped.
        private static List<Rec> ReadRows(string path, int fileIndex)
        {
            var lines = ReadAllLinesShared(path);
            int headerIdx = Array.FindIndex(lines, l => l.StartsWith("Data file name", StringComparison.OrdinalIgnoreCase));
            if (headerIdx < 0)
                throw new FormatException($"'{Path.GetFileName(path)}' is not a TopPIC proteoform tsv (no 'Data file name' header row).");

            var header = lines[headerIdx].Split('\t');
            int cRt    = FindColumn(header, RtNames);
            int cMass  = FindColumn(header, MassNames);
            int cInt   = FindColumn(header, IntensityNames);
            int cAcc   = FindColumn(header, AccNames);
            int cDescr = FindColumn(header, DescrNames);
            int cProt  = FindColumn(header, ProteoformNames);
            int cEval  = FindColumn(header, EvalueNames);
            if (cMass < 0 || cAcc < 0 || cProt < 0)
                throw new FormatException($"'{Path.GetFileName(path)}' is missing a required TopPIC column " +
                    $"({(cMass < 0 ? "mass " : "")}{(cAcc < 0 ? "accession " : "")}{(cProt < 0 ? "proteoform" : "")}).");

            var recs = new List<Rec>(lines.Length - headerIdx);
            for (int r = headerIdx + 1; r < lines.Length; r++)
            {
                if (lines[r].Length == 0) continue;
                var t = lines[r].Split('\t');
                if (!double.TryParse(Cell(t, cMass), NumberStyles.Any, Ci, out double mass)) continue;

                string rawAcc = Cell(t, cAcc);
                string descr  = Cell(t, cDescr);
                recs.Add(new Rec
                {
                    FileIndex   = fileIndex,
                    Accession   = ParseAccession(rawAcc),
                    Species     = ParseSpecies(rawAcc, descr),
                    Gene        = ParseTag(descr, "GN="),
                    Description = descr,
                    Proteoform  = Cell(t, cProt),
                    Mass        = mass,
                    Rt          = ParseDouble(t, cRt),
                    Intensity   = ParseDouble(t, cInt),
                    Evalue      = ParseEvalue(t, cEval)
                });
            }
            return recs;
        }

        // "mz|P05109|S10A8_HUMAN" -> "P05109". If there are no pipes, return the string unchanged.
        private static string ParseAccession(string raw)
        {
            if (string.IsNullOrEmpty(raw)) return "";
            var parts = raw.Split('|');
            return parts.Length >= 2 ? parts[1] : raw;
        }

        // Species from the accession's entry name suffix ("..._HUMAN"/"..._YEAST"); falls back to the OS= tag.
        private static string ParseSpecies(string rawAcc, string descr)
        {
            if (!string.IsNullOrEmpty(rawAcc))
            {
                var parts = rawAcc.Split('|');
                string entry = parts.Length >= 3 ? parts[2] : "";
                int us = entry.LastIndexOf('_');
                if (us >= 0 && us < entry.Length - 1) return entry.Substring(us + 1);
            }
            return ParseTag(descr, "OS=");
        }

        // Value of a "TAG=..." field in a UniProt-style description, read up to the next " XX=" tag or end.
        private static string ParseTag(string descr, string tag)
        {
            if (string.IsNullOrEmpty(descr)) return "";
            int i = descr.IndexOf(tag, StringComparison.Ordinal);
            if (i < 0) return "";
            i += tag.Length;
            int end = i;
            while (end < descr.Length)
            {
                // stop at the next " <letters>=" tag boundary
                if (descr[end] == ' ' && end + 2 < descr.Length && descr[end + 2] == '=' ) break;
                if (descr[end] == ' ' && end + 3 < descr.Length && descr[end + 3] == '=' ) break;
                end++;
            }
            return descr.Substring(i, end - i).Trim();
        }

        // ====================================================================================================
        //  HELPERS (mirrors the other aligners so this file is self-contained)
        // ====================================================================================================

        private static string[] ReadAllLinesShared(string path)
        {
            var lines = new List<string>();
            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            using var reader = new StreamReader(fs);
            string line;
            while ((line = reader.ReadLine()) != null) lines.Add(line);
            return lines.ToArray();
        }

        private static int FindColumn(string[] header, string[] candidates)
        {
            foreach (var c in candidates)
            {
                int idx = Array.FindIndex(header, h => string.Equals(h.Trim(), c, StringComparison.OrdinalIgnoreCase));
                if (idx >= 0) return idx;
            }
            return -1;
        }

        private static List<string> MakeUniqueLabels(List<string> files)
        {
            var labels = new List<string>();
            var seen = new Dictionary<string, int>();
            foreach (var f in files)
            {
                // strip the TopPIC suffix so the label is just the run name
                string name = Path.GetFileName(f);
                foreach (var suf in new[] { "_ms2_toppic_proteoform_single.tsv", "_ms2_toppic_proteoform.tsv", ".tsv" })
                    if (name.EndsWith(suf, StringComparison.OrdinalIgnoreCase)) { name = name.Substring(0, name.Length - suf.Length); break; }
                if (seen.TryGetValue(name, out int count)) { seen[name] = count + 1; labels.Add($"{name}.{count}"); }
                else { seen[name] = 1; labels.Add(name); }
            }
            return labels;
        }

        private static string Cell(string[] t, int col) => (col >= 0 && col < t.Length) ? t[col].Trim() : "";

        private static double ParseDouble(string[] t, int col) =>
            (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? v : double.NaN;

        // E-value column; if absent/unparseable ("-") treat as worst (1.0) so it never beats a real value.
        private static double ParseEvalue(string[] t, int col) =>
            (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? v : 1.0;

        private static int Find(int[] parent, int x)
        {
            int root = x;
            while (parent[root] != root) root = parent[root];
            while (parent[x] != root) { int next = parent[x]; parent[x] = root; x = next; }
            return root;
        }

        private static void Union(int[] parent, int a, int b)
        {
            int ra = Find(parent, a), rb = Find(parent, b);
            if (ra != rb) parent[ra] = rb;
        }

        // ---- record types (private) ------------------------------------------------------------------------

        private sealed class Rec
        {
            public int FileIndex;
            public string Accession;
            public string Species;
            public string Gene;
            public string Description;
            public string Proteoform;   // TopPIC proteoform string (sequence + localized mass shifts)
            public double Mass;         // proteoform monoisotopic mass
            public double Rt;           // feature apex time, MINUTES (NaN/0.001 = no feature)
            public double Intensity;    // feature intensity
            public double Evalue;       // proteoform E-value (lower is better)
        }

        private sealed class OutRow
        {
            public string Accession, Species, Gene, Description, Proteoform;
            public double MassMean, BestEvalue;
            public int NRuns, TotalIds;
            public double?[] Intensity, Rt, Evalue;
        }
    }
}
