using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA.FeatureAlign
{
    /// <summary>
    /// Aligns intact-proteoform feature tables using ONLY monoisotopic mass and apex retention time.
    ///
    /// Two entry points:
    ///  • <see cref="Align(FeatureAlignParameters)"/> — matches feature tables ACROSS runs into a consensus
    ///    table (one row per proteoform, with each run's abundance).
    ///  • <see cref="AlignToIdentifications"/> — matches one feature table TO a MetaMorpheus proteoform/PSM
    ///    .psmtsv (one row per feature, annotated with its best-matching identification).
    ///
    /// A match requires the masses to agree within the mass tolerance (applied around monoisotope "notches",
    /// see <see cref="WithinNotch"/>) AND the apex RTs to agree within the RT tolerance. Column names are
    /// auto-detected so FLASHQuant, FLASHDeconv and TopFD tables all work. Port of the R script custom_consensus.R.
    /// </summary>
    public static class FeatureAligner
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

        // ---- column-name dictionaries ----------------------------------------------------------------------
        // For each field we list the header names used by the different tools; the first one present in a given
        // file's header wins. This is what lets the same code read FLASHQuant, FLASHDeconv and TopFD tables.
        private static readonly string[] MassNames = { "MonoisotopicMass", "Mass" };                                              // FLASH* / TopFD
        private static readonly string[] RtNames   = { "MedianApexRetentionTime", "ApexRetentionTime", "HighestApexRetentionTime", "Apex_time", "RetentionTime" };
        private static readonly string[] IdxNames  = { "FeatureGroupIndex", "FeatureIndex", "Feature_ID" };

        // RT columns reported in MINUTES (TopFD) instead of seconds. When the chosen RT column is one of these,
        // the value is multiplied by 60 on read so every RT in this class is in SECONDS.
        private static readonly HashSet<string> RtMinuteColumns = new(StringComparer.OrdinalIgnoreCase) { "Apex_time", "Min_time", "Max_time" };

        // Column names for a MetaMorpheus PSM / proteoform .psmtsv (the identification side).
        private static readonly string[] IdMassNames = { "Precursor Mass", "Monoisotopic Mass" };   // experimental precursor mass preferred (matches feature mass to <1 ppm)
        private static readonly string[] IdRtNames   = { "Scan Retention Time", "Retention Time" };  // MetaMorpheus RT is in MINUTES
        private static readonly string[] IdSeqNames  = { "Full Sequence", "Base Sequence" };
        private static readonly string[] IdAccNames  = { "Accession", "Protein Accession" };
        private static readonly string[] IdGeneNames = { "Gene Name" };
        private static readonly string[] IdNameNames = { "Name" };
        private static readonly string[] IdQValNames = { "QValue" };

        // Charge-state columns: features have a min/max range, an identification has a single precursor charge.
        private static readonly string[] ChargeMinNames = { "MinCharge", "Min_charge" };
        private static readonly string[] ChargeMaxNames = { "MaxCharge", "Max_charge" };
        private static readonly string[] IdChargeNames  = { "Precursor Charge" };

        // ====================================================================================================
        //  CROSS-RUN ALIGNMENT
        // ====================================================================================================

        /// <summary>Reads the inputs, aligns features across runs, writes the consensus tsv, and returns its path.</summary>
        public static string Align(FeatureAlignParameters parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

            // Keep only inputs that actually exist; a consensus needs at least two runs to compare.
            var files = parameters.InputFiles?.Where(File.Exists).ToList() ?? new List<string>();
            if (files.Count < 2)
                throw new ArgumentException("Need at least 2 existing input feature files.");
            if (string.IsNullOrWhiteSpace(parameters.OutputPath))
                throw new ArgumentException("OutputPath is not set.");

            // A short, unique label per file (its base name); used to name the per-run output columns.
            var labels = MakeUniqueLabels(files);

            // Pool every feature from every run into one flat list, then sort by mass. Sorting is what makes the
            // matching loop below fast: features that could match are adjacent, so we only compare near neighbours.
            var all = new List<Feature>();
            for (int f = 0; f < files.Count; f++)
                all.AddRange(ReadFeatures(files[f], f, parameters.PreferredQuantColumn));
            all.Sort((a, b) => a.Mass.CompareTo(b.Mass));
            // Notch-aware when notches are allowed (PpmToleranceWithNotch); the grouping then matches via massTol.Within().
            var massTol = MakeTolerance(parameters.MassTolerance, parameters.MassUnit, parameters.PositiveNotch, parameters.NegativeNotch);

            // Grouping strategies:
            //  - UnionFind (default): transitive connected components.
            //  - ReferenceAnchor: OpenMS-style greedy anchor (approximation — discards non-winning window members).
            //  - OpenMSConsensus: faithful port of OpenMS TopDownConsensusFeatureGroup (keeps leftovers; all-run only).
            var groups = parameters.Grouping switch
            {
                GroupingMethod.ReferenceAnchor => GroupReferenceAnchor(all, massTol, parameters),
                GroupingMethod.OpenMSConsensus => GroupOpenMSConsensus(all, parameters, files.Count),
                _ => GroupUnionFind(all, massTol, parameters),
            };

            WriteConsensusTsv(groups, labels, parameters.OutputPath, parameters);
            return parameters.OutputPath;
        }

        // ---- grouping strategy 1: union-find / connected components (original) -----------------------------
        // Links any two features within mass(+notch)+RT tolerance of EACH OTHER; links are transitive, so a
        // proteoform can span any number of runs (but a long ladder of near-mass features can chain into one blob).
        private static List<List<Feature>> GroupUnionFind(List<Feature> all, Tolerance massTol, FeatureAlignParameters p)
        {
            int n = all.Count;
            var parent = new int[n];
            for (int i = 0; i < n; i++) parent[i] = i;

            for (int i = 0; i < n; i++)   // all is sorted by mass, so the inner scan can stop early
            {
                double mi = all[i].Mass, ri = all[i].Rt;
                double upper = massTol.GetMaximumValue(mi);   // notch-aware tolerance already widens this to the top notch
                for (int j = i + 1; j < n; j++)
                {
                    if (all[j].Mass > upper) break;
                    if (Math.Abs(all[j].Rt - ri) <= p.RtTolerance && massTol.Within(all[j].Mass, mi))
                        Union(parent, i, j);
                }
            }

            var clusters = new Dictionary<int, List<Feature>>();
            for (int i = 0; i < n; i++)
            {
                int root = Find(parent, i);
                if (!clusters.TryGetValue(root, out var list)) { list = new List<Feature>(); clusters[root] = list; }
                list.Add(all[i]);
            }
            return clusters.Values.ToList();
        }

        // ---- grouping strategy 2: greedy reference-anchored (OpenMS ConsensusFeatureGroupDetector) ----------
        // Repeatedly take the most-abundant unused feature as a fixed reference/anchor, gather every unused
        // feature within mass(+notch)+RT tolerance OF THE ANCHOR, keep one per run (WhenDuplicate), and discard
        // the rest. Because candidates are compared to the anchor (not to each other) the group span is capped at
        // the tolerance width — it cannot chain an unbounded ladder the way union-find can.
        private static List<List<Feature>> GroupReferenceAnchor(List<Feature> all, Tolerance massTol, FeatureAlignParameters p)
        {
            int n = all.Count;
            var used = new bool[n];
            var byAbundance = Enumerable.Range(0, n).OrderByDescending(i => all[i].Quant).ToList();  // greedy order
            var groups = new List<List<Feature>>();

            foreach (int refIdx in byAbundance)
            {
                if (used[refIdx]) continue;
                var reference = all[refIdx];

                // gather all still-unused features within the anchor's window, grouped by run
                var perRun = new Dictionary<int, List<int>>();
                for (int j = 0; j < n; j++)
                {
                    if (used[j]) continue;
                    var c = all[j];
                    if (Math.Abs(c.Rt - reference.Rt) <= p.RtTolerance && massTol.Within(c.Mass, reference.Mass))
                    {
                        if (!perRun.TryGetValue(c.FileIndex, out var li)) { li = new List<int>(); perRun[c.FileIndex] = li; }
                        li.Add(j);
                    }
                }

                // keep one feature per run (WhenDuplicate), mark ALL window features used (duplicates discarded)
                var group = new List<Feature>();
                foreach (var li in perRun.Values)
                {
                    int chosen = p.WhenDuplicate == DuplicateResolution.NearestMass
                        ? li.OrderBy(j => Math.Abs(all[j].Mass - reference.Mass)).First()
                        : li.OrderByDescending(j => all[j].Quant).First();
                    group.Add(all[chosen]);
                    foreach (int j in li) used[j] = true;
                }
                groups.Add(group);
            }
            return groups;
        }

        // ---- grouping strategy 3: faithful port of OpenMS TopDownConsensusFeatureGroup -----------------------
        // (the ConsensusFeatureGroupDetector.exe). Equivalent of OpenMS's computeConsensusFeatureGroup:
        //   * greedy: the most-abundant remaining feature is the reference/anchor;
        //   * candidates are the still-unused features within an ABSOLUTE ppm/Da window (no notch) and RT of the
        //     anchor, taken from OTHER runs only (the anchor is its own run's representative);
        //   * a group is emitted ONLY if candidates cover every other run (present in ALL runs); otherwise just the
        //     reference is dropped and its candidates stay in the pool;
        //   * KEY difference from ReferenceAnchor: only the chosen members (reference + one winner per other run) are
        //     removed — the rest of the window stays in the pool, so two proteoforms sharing a window become two
        //     groups instead of one merged/discarded group.
        // repCount = number of input runs (OpenMS's REP_COUNT_).
        private static List<List<Feature>> GroupOpenMSConsensus(List<Feature> all, FeatureAlignParameters p, int repCount)
        {
            var pts = all.OrderBy(f => f.Mass).ToList();                       // mass-sorted working pool
            bool ppm = p.MassUnit == MassToleranceUnit.Ppm;
            double tolBase = ppm ? p.MassTolerance * 1e-6 : p.MassTolerance;   // OpenMS MASS_TOL_ (already *1e-6 for ppm)
            var groups = new List<List<Feature>>();

            while (pts.Count > 0)
            {
                // reference = the max-abundance remaining feature (first max on ties, like std::max_element)
                int refIdx = 0;
                for (int i = 1; i < pts.Count; i++)
                    if (pts[i].Quant > pts[refIdx].Quant) refIdx = i;
                var reference = pts[refIdx];
                double refMass = reference.Mass;
                double tol = ppm ? tolBase * refMass : tolBase;

                // candidates within mass+RT of the anchor, from OTHER runs only (scan right then left in mass order)
                var cand = new List<(Feature f, int idx)>();
                for (int j = refIdx + 1; j < pts.Count && pts[j].Mass - refMass <= tol; j++)
                    if (Math.Abs(pts[j].Rt - reference.Rt) <= p.RtTolerance && pts[j].FileIndex != reference.FileIndex)
                        cand.Add((pts[j], j));
                for (int j = refIdx - 1; j >= 0 && refMass - pts[j].Mass <= tol; j--)
                    if (Math.Abs(pts[j].Rt - reference.Rt) <= p.RtTolerance && pts[j].FileIndex != reference.FileIndex)
                        cand.Add((pts[j], j));

                // OpenMS keeps a group only if it spans every run (reference's run + all repCount-1 others)
                var repSet = cand.Select(c => c.f.FileIndex).Distinct().ToList();
                if (repSet.Count < repCount - 1)
                {
                    pts.RemoveAt(refIdx);   // drop only the reference; leave its candidates for later anchors
                    continue;
                }

                // one winner per other run (max-abundance or nearest-mass), plus the reference
                var chosenIdx = new List<int> { refIdx };
                var group = new List<Feature> { reference };
                foreach (int rep in repSet)
                {
                    var inRep = cand.Where(c => c.f.FileIndex == rep);
                    var pick = p.WhenDuplicate == DuplicateResolution.NearestMass
                        ? inRep.OrderBy(c => Math.Abs(c.f.Mass - refMass)).First()
                        : inRep.OrderByDescending(c => c.f.Quant).First();
                    group.Add(pick.f);
                    chosenIdx.Add(pick.idx);
                }
                groups.Add(group);

                // remove ONLY the chosen features (descending index so earlier removals don't shift later ones)
                foreach (int idx in chosenIdx.OrderByDescending(x => x)) pts.RemoveAt(idx);
            }
            return groups;
        }

        /// <summary>
        /// Convenience overload that builds the parameters from individual values.
        /// <paramref name="minNFiles"/>/<paramref name="maxNFiles"/> and <paramref name="minMass"/>/<paramref name="maxMass"/>
        /// are OUTPUT filters (inclusive; null = unbounded): a consensus row is written only if its n_files and mean
        /// mass fall inside those ranges. Leave them null to write every consensus row (the original behaviour).
        /// </summary>
        public static string Align(IEnumerable<string> inputFiles, string outputPath,
            double massTolerance = 20, MassToleranceUnit massUnit = MassToleranceUnit.Ppm,
            double rtTolerance = 150, string preferredQuantColumn = "FeatureGroupQuantity",
            int positiveNotch = 0, int negativeNotch = 0,
            GroupingMethod grouping = GroupingMethod.ReferenceAnchor, DuplicateResolution whenDuplicate = DuplicateResolution.MaxAbundance,
            int? minNFiles = null, int? maxNFiles = null, double? minMass = null, double? maxMass = null)
        {
            var p = new FeatureAlignParameters(inputFiles, outputPath, massTolerance, massUnit, rtTolerance, preferredQuantColumn);
            p.PositiveNotch = positiveNotch;
            p.NegativeNotch = negativeNotch;
            p.Grouping = grouping;
            p.WhenDuplicate = whenDuplicate;
            p.MinNFiles = minNFiles;
            p.MaxNFiles = maxNFiles;
            p.MinMass = minMass;
            p.MaxMass = maxMass;
            return Align(p);
        }

        // ====================================================================================================
        //  FEATURE -> IDENTIFICATION ALIGNMENT
        // ====================================================================================================

        /// <summary>
        /// Aligns features (one FLASHQuant/FLASHDeconv/TopFD table) to identifications (a MetaMorpheus
        /// PSM/proteoform .psmtsv) by mass and apex RT, writing one row per feature with its best-matching
        /// identification (blank if none).
        ///
        /// Notes that make the matches land:
        ///  • MetaMorpheus RT is in MINUTES; feature RT is in SECONDS. ID RT is converted ×60 (unless
        ///    <paramref name="idRetentionTimeInMinutes"/> is false).
        ///  • The feature mass is matched to the psmtsv's experimental "Precursor Mass" (not the theoretical
        ///    "Monoisotopic Mass"), so a real match is well under 1 ppm.
        ///  • The mass tolerance is applied around monoisotope notches (0, ±1..±k × <see cref="NotchStep"/> Da),
        ///    because the deconvolutor and the search engine often pick monoisotopic peaks 1-2 isotopes apart on
        ///    intact proteoforms. Set both notches to 0 to require an exact monoisotope match.
        ///  • Identifications with a QValue above <paramref name="qValueThreshold"/> are ignored.
        /// </summary>
        public static string AlignToIdentifications(string featureFile, string identificationFile, string outputPath,
            double massTolerance = 20, MassToleranceUnit massUnit = MassToleranceUnit.Ppm, double rtTolerance = 150,
            double qValueThreshold = 0.01, string preferredQuantColumn = "FeatureGroupQuantity",
            bool idRetentionTimeInMinutes = true, int positiveNotch = 2, int negativeNotch = 2,
            bool requireChargeInRange = false)
        {
            if (!File.Exists(featureFile)) throw new ArgumentException("Feature file not found: " + featureFile);
            if (!File.Exists(identificationFile)) throw new ArgumentException("Identification file not found: " + identificationFile);
            if (string.IsNullOrWhiteSpace(outputPath)) throw new ArgumentException("OutputPath is not set.");

            var features = ReadFeatures(featureFile, 0, preferredQuantColumn);
            var ids = ReadIdentifications(identificationFile, qValueThreshold, idRetentionTimeInMinutes);

            var dir = Path.GetDirectoryName(outputPath);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir)) Directory.CreateDirectory(dir);

            var massTol = MakeTolerance(massTolerance, massUnit, positiveNotch, negativeNotch);   // notch-aware match via massTol.Within
            var shifts  = NotchShifts(positiveNotch, negativeNotch);                              // only to report the best isotope offset below

            using var sw = new StreamWriter(outputPath, append: false);
            sw.WriteLine(string.Join("\t",
                "feature_mass", "feature_rt", "feature_quant", "feature_idx", "feature_charge_range", "num_id_matches",
                "accession", "gene", "name", "full_sequence", "id_mass", "id_rt", "isotope_shift", "mass_error_ppm", "rt_diff", "qvalue",
                "id_charge", "charge_in_range"));

            foreach (var f in features.OrderBy(x => x.Mass))
            {
                // Candidate IDs: within RT tolerance AND within mass tolerance of some notch offset. Optionally the
                // identified precursor charge must also fall inside the feature's observed charge range.
                var matches = ids.Where(id =>
                    Math.Abs(f.Rt - id.RtSeconds) <= rtTolerance &&
                    massTol.Within(f.Mass, id.Mass) &&
                    (!requireChargeInRange || ChargeInRange(f, id))).ToList();

                // Best candidate: the one whose mass is closest after allowing the nearest notch, RT breaks ties.
                var best = matches.OrderBy(id => Math.Abs(f.Mass - (id.Mass + BestNotch(f.Mass, id.Mass, shifts) * NotchStep)))
                                  .ThenBy(id => Math.Abs(f.Rt - id.RtSeconds)).FirstOrDefault();

                string fRange = f.MaxCharge > 0 ? $"{f.MinCharge}-{f.MaxCharge}" : "";
                var cells = new List<string>
                {
                    f.Mass.ToString(Ci), f.Rt.ToString(Ci), f.Quant.ToString(Ci), f.Index ?? "", fRange, matches.Count.ToString(Ci)
                };
                if (best != null)
                {
                    int shift = BestNotch(f.Mass, best.Mass, shifts);                          // which isotope offset matched (0, ±1, ±2…)
                    double residualPpm = (f.Mass - (best.Mass + shift * NotchStep)) / best.Mass * 1e6;   // residual error AFTER removing that offset
                    cells.Add(best.Accession);
                    cells.Add(best.Gene);
                    cells.Add(best.Name);
                    cells.Add(best.Sequence);
                    cells.Add(best.Mass.ToString(Ci));
                    cells.Add(best.RtSeconds.ToString(Ci));
                    cells.Add(shift.ToString(Ci));
                    cells.Add(residualPpm.ToString("F2", Ci));
                    cells.Add((f.Rt - best.RtSeconds).ToString("F2", Ci));
                    cells.Add(best.QValue.ToString(Ci));
                    cells.Add(best.Charge > 0 ? best.Charge.ToString(Ci) : "");                // id_charge
                    cells.Add(ChargeInRange(f, best) ? "1" : "0");                             // charge_in_range (1 = ID charge is inside feature's range)
                }
                else
                {
                    cells.AddRange(Enumerable.Repeat("", 12));   // no ID: leave the 12 identification columns blank
                }
                sw.WriteLine(string.Join("\t", cells));
            }
            return outputPath;
        }

        /// <summary>
        /// Rewrites a FLASHDeconv feature tsv into the FLASHQuant 17-column layout that OpenMS's
        /// ConsensusFeatureGroupDetector.exe expects (it reads columns by fixed position and CRASHES on any other
        /// layout). If the file is already FLASHQuant (has a FeatureGroupIndex column) it is copied unchanged.
        /// Returns the output path.
        /// </summary>
        public static string ConvertFlashDeconvToFlashQuant(string inPath, string outPath)
        {
            var lines = ReadAllLinesShared(inPath);
            if (lines.Length < 1) throw new FormatException("Empty file: " + inPath);
            var header = lines[0].Split('\t');
            int Col(string name) => Array.FindIndex(header, h => string.Equals(h.Trim(), name, StringComparison.OrdinalIgnoreCase));

            var dir = Path.GetDirectoryName(outPath);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir)) Directory.CreateDirectory(dir);

            if (Col("FeatureGroupIndex") >= 0) { File.Copy(inPath, outPath, true); return outPath; }   // already FLASHQuant

            int fi = Col("FeatureIndex"), fn = Col("FileName"), mm = Col("MonoisotopicMass"), am = Col("AverageMass"),
                srt = Col("StartRetentionTime"), ert = Col("EndRetentionTime"), art = Col("ApexRetentionTime"),
                si = Col("SumIntensity"), fq = Col("FeatureQuantity"), mnc = Col("MinCharge"), mxc = Col("MaxCharge"),
                cc = Col("ChargeCount"), ics = Col("IsotopeCosineScore"), mq = Col("MaxQscore");
            if (mm < 0 || art < 0 || fq < 0)
                throw new FormatException("Not a FLASHDeconv feature table: " + Path.GetFileName(inPath));

            using var sw = new StreamWriter(outPath, append: false);
            sw.WriteLine(string.Join("\t", "FeatureGroupIndex", "FileName", "MonoisotopicMass", "AverageMass",
                "StartRetentionTime(FWHM)", "EndRetentionTime(FWHM)", "HighestApexRetentionTime", "MedianApexRetentionTime",
                "FeatureGroupQuantity", "AllAreaUnderTheCurve", "SumIntensity", "MinCharge", "MaxCharge", "ChargeCount",
                "MostAbundantFeatureCharge", "IsotopeCosineScore", "FeatureScore"));
            for (int r = 1; r < lines.Length; r++)
            {
                if (lines[r].Length == 0) continue;
                var t = lines[r].Split('\t');
                sw.WriteLine(string.Join("\t", Cell(t, fi), Cell(t, fn), Cell(t, mm), Cell(t, am), Cell(t, srt),
                    Cell(t, ert), Cell(t, art), Cell(t, art), Cell(t, fq), Cell(t, fq), Cell(t, si), Cell(t, mnc),
                    Cell(t, mxc), Cell(t, cc), Cell(t, mnc), Cell(t, ics), Cell(t, mq)));
            }
            return outPath;
        }

        // ====================================================================================================
        //  MASS TOLERANCE WITH NOTCHES
        //  A "notch" is an integer number of C13-C12 steps. Matching within a tight ppm AROUND each notch lets
        //  us bridge monoisotope-assignment disagreements (±1-2 Da) without the chain-merging that a single wide
        //  Da window causes. This mirrors MzLibUtil.PpmToleranceWithNotch, built on the real mzLib Tolerance.
        // ====================================================================================================

        public const double NotchStep = 1.0033548381; // C13 - C12 mass difference (Da)

        // ppm -> PpmTolerance, or mzLib's PpmToleranceWithNotch when any monoisotope notch is allowed; Da ->
        // AbsoluteTolerance. Notches are a ppm/isotope concept, so a Dalton tolerance ignores them. Callers do the
        // mass match through the returned Tolerance's own Within(), so the notch search is mzLib's
        // PpmToleranceWithNotch (whose Within loops the same {0,±1..} shifts), not a hand-rolled loop here.
        // internal so the sibling ConsensusFeatureAligner reuses the exact same tolerance construction.
        internal static Tolerance MakeTolerance(double value, MassToleranceUnit unit, int positiveNotch = 0, int negativeNotch = 0)
        {
            if (unit != MassToleranceUnit.Ppm) return new AbsoluteTolerance(value);
            if (positiveNotch > 0 || negativeNotch > 0) return new PpmToleranceWithNotch(value, positiveNotch, negativeNotch);
            return new PpmTolerance(value);
        }

        // The allowed offsets: {0, +1..+pos, -1..-neg} × NotchStep. pos = neg = 0 gives just {0} (exact match).
        internal static double[] NotchShifts(int positiveNotch, int negativeNotch)
        {
            var shifts = new List<double>(1 + positiveNotch + negativeNotch) { 0.0 };
            for (int i = 1; i <= positiveNotch; i++) shifts.Add(NotchStep * i);
            for (int i = 1; i <= negativeNotch; i++) shifts.Add(-NotchStep * i);
            return shifts.ToArray();
        }

        // True if 'experimental' is within tolerance of 'theoretical + shift' for ANY notch (PpmToleranceWithNotch.Within).
        internal static bool WithinNotch(Tolerance tol, double experimental, double theoretical, double[] shifts)
        {
            foreach (double s in shifts)
                if (tol.Within(experimental, theoretical + s)) return true;
            return false;
        }

        // The whole-isotope offset k that best explains experimental - theoretical (used for reporting/residual).
        private static int BestNotch(double experimental, double theoretical, double[] shifts)
        {
            int best = 0; double bestAbs = double.MaxValue;
            foreach (double s in shifts)
            {
                double d = Math.Abs(experimental - (theoretical + s));
                if (d < bestAbs) { bestAbs = d; best = (int)Math.Round(s / NotchStep); }
            }
            return best;
        }

        // ====================================================================================================
        //  UNION-FIND (disjoint set) — tracks which features are in the same cluster as links are added
        // ====================================================================================================

        // Return the representative (root) of x's set, compressing the path so later lookups are O(1)-ish.
        private static int Find(int[] parent, int x)
        {
            int root = x;
            while (parent[root] != root) root = parent[root];
            while (parent[x] != root) { int next = parent[x]; parent[x] = root; x = next; } // path compression
            return root;
        }

        // Merge the sets containing a and b (no-op if already together).
        private static void Union(int[] parent, int a, int b)
        {
            int ra = Find(parent, a), rb = Find(parent, b);
            if (ra != rb) parent[ra] = rb;
        }

        // ====================================================================================================
        //  READING
        // ====================================================================================================

        // Read one feature table into Feature records. Resolves the mass/RT/quant/index columns by name (so any
        // of the supported tools works) and normalises RT to seconds. Rows that don't parse are skipped.
        private static List<Feature> ReadFeatures(string path, int fileIndex, string preferredQuantColumn)
        {
            var lines = ReadAllLinesShared(path);
            if (lines.Length < 2) return new List<Feature>();

            var header = lines[0].Split('\t');
            int massCol  = FindColumn(header, MassNames);
            int rtCol    = FindColumn(header, RtNames);
            int quantCol = FindColumn(header, QuantNames(preferredQuantColumn));
            int idxCol   = FindColumn(header, IdxNames);        // -1 -> fall back to the row number
            int cMinCol  = FindColumn(header, ChargeMinNames);  // -1 -> no charge info (0)
            int cMaxCol  = FindColumn(header, ChargeMaxNames);

            // Mass, RT and a quantity are required; a clear error naming the file beats a downstream crash.
            if (massCol < 0 || rtCol < 0 || quantCol < 0)
                throw new FormatException(
                    $"'{Path.GetFileName(path)}' is not a FLASHQuant/FLASHDeconv/TopFD feature table (missing " +
                    $"{(massCol < 0 ? "mass " : "")}{(rtCol < 0 ? "rt " : "")}{(quantCol < 0 ? "quant" : "")}). " +
                    $"Columns found: {string.Join(", ", header)}");

            double rtScale = RtMinuteColumns.Contains(header[rtCol].Trim()) ? 60.0 : 1.0;   // TopFD minutes -> seconds
            int maxCol = Math.Max(massCol, Math.Max(rtCol, quantCol));                       // guard against short/ragged rows
            var features = new List<Feature>(lines.Length - 1);
            for (int r = 1; r < lines.Length; r++)
            {
                if (lines[r].Length == 0) continue;
                var t = lines[r].Split('\t');
                if (t.Length <= maxCol) continue;
                if (!double.TryParse(t[massCol], NumberStyles.Any, Ci, out double mass)) continue;
                if (!double.TryParse(t[rtCol],   NumberStyles.Any, Ci, out double rt))   continue;
                rt *= rtScale;
                double.TryParse(t[quantCol], NumberStyles.Any, Ci, out double quant);      // blank/parse-fail -> 0
                string idx = (idxCol >= 0 && idxCol < t.Length) ? t[idxCol] : r.ToString(Ci);
                features.Add(new Feature { FileIndex = fileIndex, Index = idx, Mass = mass, Rt = rt, Quant = quant,
                    MinCharge = ParseInt(t, cMinCol), MaxCharge = ParseInt(t, cMaxCol) });
            }
            return features;
        }

        // Candidate quantity columns, most-preferred first (the caller's choice, then the known tool columns).
        private static string[] QuantNames(string preferred)
        {
            var names = new List<string> { preferred, "FeatureGroupQuantity", "FeatureQuantity", "SumIntensity", "Intensity", "Apex_intensity" };
            return names.Where(s => !string.IsNullOrWhiteSpace(s)).Distinct().ToArray();
        }

        // Read all lines with a SHARED handle so a file that is open in Excel/another program can still be read.
        // internal so the sibling ConsensusFeatureAligner reuses the exact same tsv reading.
        internal static string[] ReadAllLinesShared(string path)
        {
            var lines = new List<string>();
            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            using var reader = new StreamReader(fs);
            string line;
            while ((line = reader.ReadLine()) != null) lines.Add(line);
            return lines.ToArray();
        }

        // Index of the first candidate name present in the header (case-insensitive, trimmed); -1 if none.
        internal static int FindColumn(string[] header, string[] candidates)
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

        // Read the identification rows (mass, RT-in-seconds, sequence/accession/gene/name, QValue) from a
        // MetaMorpheus .psmtsv, dropping anything above the QValue threshold.
        private static List<Identification> ReadIdentifications(string path, double qValueThreshold, bool rtInMinutes)
        {
            var lines = ReadAllLinesShared(path);
            if (lines.Length < 2) return new List<Identification>();

            var header = lines[0].Split('\t');
            int massCol = FindColumn(header, IdMassNames);
            int rtCol   = FindColumn(header, IdRtNames);
            if (massCol < 0 || rtCol < 0)
                throw new FormatException(
                    $"'{Path.GetFileName(path)}' is not a MetaMorpheus PSM/proteoform tsv (missing " +
                    $"{(massCol < 0 ? "mass " : "")}{(rtCol < 0 ? "retention time" : "")}). Columns found: {string.Join(", ", header)}");

            int seqCol  = FindColumn(header, IdSeqNames);
            int accCol  = FindColumn(header, IdAccNames);
            int geneCol = FindColumn(header, IdGeneNames);
            int nameCol = FindColumn(header, IdNameNames);
            int qCol    = FindColumn(header, IdQValNames);
            int zCol    = FindColumn(header, IdChargeNames);

            double rtScale = rtInMinutes ? 60.0 : 1.0;   // MetaMorpheus RT (min) -> seconds
            int maxCol = Math.Max(massCol, rtCol);

            var ids = new List<Identification>(lines.Length - 1);
            for (int r = 1; r < lines.Length; r++)
            {
                if (lines[r].Length == 0) continue;
                var t = lines[r].Split('\t');
                if (t.Length <= maxCol) continue;
                if (!double.TryParse(t[massCol], NumberStyles.Any, Ci, out double mass)) continue;
                if (!double.TryParse(t[rtCol],   NumberStyles.Any, Ci, out double rt))   continue;

                double q = double.NaN;
                if (qCol >= 0 && qCol < t.Length) double.TryParse(t[qCol], NumberStyles.Any, Ci, out q);
                if (!double.IsNaN(q) && q > qValueThreshold) continue;   // FDR filter (only when a QValue column exists)

                ids.Add(new Identification
                {
                    Mass = mass,
                    RtSeconds = rt * rtScale,
                    Sequence  = Cell(t, seqCol),
                    Accession = Cell(t, accCol),
                    Gene      = Cell(t, geneCol),
                    Name      = Cell(t, nameCol),
                    QValue    = double.IsNaN(q) ? 0 : q,
                    Charge    = ParseInt(t, zCol)
                });
            }
            return ids;
        }

        // Safe cell access: value at col, or "" if the column is absent/short.
        private static string Cell(string[] t, int col) => (col >= 0 && col < t.Length) ? t[col] : "";

        // Safe integer parse of a cell (charge states); 0 if the column is absent/short/unparseable.
        // Parses via double so values written with decimals (e.g. MetaMorpheus "11.00000") still work.
        internal static int ParseInt(string[] t, int col) =>
            (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? (int)Math.Round(v) : 0;

        // True if the identified precursor charge sits inside the feature's observed charge range.
        // Permissive when charge info is missing on either side (can't check -> treated as not violated).
        private static bool ChargeInRange(Feature f, Identification id) =>
            f.MaxCharge <= 0 || id.Charge <= 0 || (id.Charge >= f.MinCharge && id.Charge <= f.MaxCharge);

        // ====================================================================================================
        //  WRITING THE CONSENSUS (cross-run) TABLE
        // ====================================================================================================

        private static void WriteConsensusTsv(IEnumerable<List<Feature>> clusters, List<string> labels, string outputPath,
            FeatureAlignParameters parameters)
        {
            int fileCount = labels.Count;

            // Collapse each cluster (one proteoform) to a single row with one value per run.
            var rows = new List<(double MassMean, double RtMean, int NumFiles, int MaxPerRun, string ChargeRange, double?[] Quant, string[] Idx)>();
            foreach (var cluster in clusters)
            {
                var quant = new double?[fileCount];   // per-run abundance; null means that run has no feature here
                var idx = new string[fileCount];      // per-run original feature index (of the representative)
                int maxPerRun = 1;

                // Group the cluster's features by run. In a clean cluster each run contributes exactly ONE
                // feature. If a run contributes several, their masses are within tolerance in the SAME run —
                // which normally can't be two real proteoforms, so it is almost always a DECONVOLUTION ARTIFACT
                // (one species reported twice). We keep that run's most-abundant feature as its representative
                // and record maxPerRun = the largest per-run count so these clusters can be flagged/filtered.
                foreach (var g in cluster.GroupBy(x => x.FileIndex))
                {
                    if (g.Count() > maxPerRun) maxPerRun = g.Count();
                    var best = g.OrderByDescending(x => x.Quant).First();
                    quant[g.Key] = best.Quant;
                    idx[g.Key] = best.Index;
                }

                // Observed charge range across the whole cluster (blank if the files carry no charge columns).
                var withCharge = cluster.Where(x => x.MaxCharge > 0).ToList();
                string chargeRange = withCharge.Count > 0 ? $"{withCharge.Min(x => x.MinCharge)}-{withCharge.Max(x => x.MaxCharge)}" : "";

                // Consensus mass/RT are the simple means over the cluster; NumFiles = how many runs it spans.
                rows.Add((cluster.Average(x => x.Mass), cluster.Average(x => x.Rt),
                          quant.Count(q => q.HasValue), maxPerRun, chargeRange, quant, idx));
            }

            // Order the output: most-shared proteoforms first, then by mass.
            rows.Sort((a, b) =>
            {
                int c = b.NumFiles.CompareTo(a.NumFiles);
                return c != 0 ? c : a.MassMean.CompareTo(b.MassMean);
            });

            var dir = Path.GetDirectoryName(outputPath);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir)) Directory.CreateDirectory(dir);

            using var sw = new StreamWriter(outputPath, append: false);

            // Header: fixed columns, then quant_<run> for every run, then idx_<run> for every run.
            // max_features_per_run == 1 is clean; > 1 flags a same-run duplicate (suspected deconvolution error).
            var header = new List<string> { "mass_mean", "rt_mean", "n_files", "max_features_per_run", "charge_range" };
            header.AddRange(labels.Select(l => "quant_" + l));
            header.AddRange(labels.Select(l => "idx_" + l));
            sw.WriteLine(string.Join("\t", header));

            foreach (var row in rows)
            {
                // OUTPUT FILTER: skip rows whose n_files / mass_mean fall outside the requested ranges. Grouping is
                // already done, so this only shrinks the written table (leaving every bound null writes everything).
                if (!parameters.PassesOutputFilter(row.NumFiles, row.MassMean)) continue;

                var cells = new List<string>
                {
                    row.MassMean.ToString(Ci),
                    row.RtMean.ToString(Ci),
                    row.NumFiles.ToString(Ci),
                    row.MaxPerRun.ToString(Ci),
                    row.ChargeRange
                };
                cells.AddRange(row.Quant.Select(q => q.HasValue ? q.Value.ToString(Ci) : ""));   // blank where a run has no feature
                cells.AddRange(row.Idx.Select(s => s ?? ""));
                sw.WriteLine(string.Join("\t", cells));
            }
        }

        // ---- record types (private) ------------------------------------------------------------------------

        private sealed class Feature
        {
            public int FileIndex;   // which input run this feature came from (index into the files/labels list)
            public string Index;    // the feature's original id in its source file
            public double Mass;     // monoisotopic mass
            public double Rt;       // apex retention time, in SECONDS
            public double Quant;    // abundance
            public int MinCharge;   // lowest charge state observed (0 if the file has no charge columns)
            public int MaxCharge;   // highest charge state observed (0 if none)
        }

        private sealed class Identification
        {
            public double Mass;       // experimental precursor monoisotopic mass
            public double RtSeconds;  // scan retention time, in SECONDS
            public string Sequence;
            public string Accession;
            public string Gene;
            public string Name;
            public double QValue;
            public int Charge;        // precursor charge of the identified scan (0 if absent)
        }
    }
}
