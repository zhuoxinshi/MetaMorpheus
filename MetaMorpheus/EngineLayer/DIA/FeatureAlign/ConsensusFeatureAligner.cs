using MassSpectrometry.Deconvolution.Consensus;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA.FeatureAlign
{
    /// <summary>Which consensus-pipeline intensity to treat as a feature's abundance (drives the greedy
    /// reference-anchor order and same-run duplicate resolution). All three measures are written to the output
    /// regardless of which one is chosen here.</summary>
    public enum ConsensusQuantMeasure
    {
        /// <summary><see cref="MassFeature.SummedIntensity"/> — summed over every charge trace (the "Intensity" column).</summary>
        SummedIntensity,
        /// <summary><see cref="MassFeature.SummedIntensityTop3"/> — summed over the 3 most abundant charge traces.</summary>
        SummedIntensityTop3,
        /// <summary><see cref="MassFeature.ApexIntensityTop3"/> — the 3 most abundant charge traces at the shared apex scan.</summary>
        ApexIntensityTop3
    }

    /// <summary>
    /// Cross-run aligner for CONSENSUS-pipeline features. It consumes the `_ms1.feature` tables written by
    /// <see cref="IsdMsAlignExporter.WriteMs1FeatureFile"/> (one per run) and preserves ALL THREE consensus quant
    /// measures per run — <c>Intensity</c> (=<see cref="MassFeature.SummedIntensity"/>), <c>Sum_intensity_top3</c>
    /// (=<see cref="MassFeature.SummedIntensityTop3"/>) and <c>Apex_intensity_top3</c>
    /// (=<see cref="MassFeature.ApexIntensityTop3"/>) — instead of collapsing to a single quantity column the way
    /// <see cref="FeatureAligner"/> does. Column names are auto-detected, so plain TopFD/FLASHDeconv `_ms1.feature`
    /// tables also work (the top-3 columns simply fall back to the summed intensity when absent).
    ///
    /// Matching is identical to <see cref="FeatureAligner"/>: two features are the same proteoform when their
    /// monoisotopic masses agree within the mass tolerance (applied around ±k isotope "notches") AND their apex
    /// retention times agree within the RT tolerance. The notch/tolerance/reading helpers are shared with
    /// <see cref="FeatureAligner"/> so both aligners behave identically. Two grouping strategies are supported:
    /// <see cref="GroupingMethod.UnionFind"/> (transitive components) and <see cref="GroupingMethod.ReferenceAnchor"/>
    /// (greedy, anchored to the most-abundant remaining feature — bounded group width, one feature kept per run).
    ///
    /// Units: the apex-RT column (e.g. "Apex_time") is in MINUTES; it is converted to SECONDS on read, so
    /// <c>rtTolerance</c> is in seconds and means the same thing here as in <see cref="FeatureAligner"/>.
    /// </summary>
    public static class ConsensusFeatureAligner
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

        // ---- column-name dictionaries (first present in a file's header wins) ------------------------------
        private static readonly string[] FMass     = { "Mass", "MonoisotopicMass" };
        private static readonly string[] FRt       = { "Apex_time", "MedianApexRetentionTime", "ApexRetentionTime", "HighestApexRetentionTime", "RetentionTime" };
        private static readonly string[] FSum      = { "Intensity", "SumIntensity", "FeatureGroupQuantity", "FeatureQuantity" };
        private static readonly string[] FSumTop3  = { "Sum_intensity_top3" };
        private static readonly string[] FApexTop3 = { "Apex_intensity_top3" };
        private static readonly string[] FId       = { "Feature_ID", "FeatureIndex", "FeatureGroupIndex" };
        private static readonly string[] FCMin     = { "Min_charge", "MinCharge" };
        private static readonly string[] FCMax     = { "Max_charge", "MaxCharge" };

        // RT columns reported in MINUTES (converted ×60 to seconds on read). Mirrors FeatureAligner.
        private static readonly HashSet<string> RtMinuteColumns =
            new(StringComparer.OrdinalIgnoreCase) { "Apex_time", "Min_time", "Max_time" };

        // ====================================================================================================
        //  ENTRY POINT: align a set of `_ms1.feature` files
        // ====================================================================================================

        /// <summary>
        /// Align the `_ms1.feature` tables of two or more runs into one consensus table (one row per proteoform,
        /// with each run's three quant measures and the feature id it came from) and write it to
        /// <paramref name="outputPath"/>. Returns the output path.
        /// </summary>
        public static string AlignFeatureFiles(
            IEnumerable<string> featureFiles, string outputPath,
            double massTolerance = 20, MassToleranceUnit massUnit = MassToleranceUnit.Ppm,
            double rtTolerance = 150, int positiveNotch = 2, int negativeNotch = 2,
            GroupingMethod grouping = GroupingMethod.ReferenceAnchor,
            DuplicateResolution whenDuplicate = DuplicateResolution.MaxAbundance,
            ConsensusQuantMeasure quantMeasure = ConsensusQuantMeasure.SummedIntensityTop3)
        {
            if (featureFiles == null) throw new ArgumentNullException(nameof(featureFiles));

            var files = featureFiles.Where(File.Exists).ToList();
            if (files.Count < 2) throw new ArgumentException("Need at least 2 existing feature files to align.");
            if (string.IsNullOrWhiteSpace(outputPath)) throw new ArgumentException("OutputPath is not set.");

            var labels = MakeUniqueLabels(files.Select(Path.GetFileNameWithoutExtension).ToList());

            var all = new List<CFeature>();
            for (int r = 0; r < files.Count; r++)
                all.AddRange(ReadFeatureFile(files[r], r, quantMeasure));

            return AlignCore(all, labels, outputPath, massTolerance, massUnit, rtTolerance,
                positiveNotch, negativeNotch, grouping, whenDuplicate);
        }

        /// <summary>
        /// In-memory overload: align consensus features that are already traced (e.g. straight from
        /// <see cref="ConsensusMassXicConstructor.TraceFeatures"/>) without writing/reading the intermediate
        /// `_ms1.feature` files. One entry per run: a label (used to name the columns) and that run's MassFeatures.
        /// </summary>
        public static string Align(
            IReadOnlyList<(string label, IReadOnlyList<MassFeature> features)> runs,
            string outputPath,
            double massTolerance = 20, MassToleranceUnit massUnit = MassToleranceUnit.Ppm,
            double rtTolerance = 150, int positiveNotch = 2, int negativeNotch = 2,
            GroupingMethod grouping = GroupingMethod.ReferenceAnchor,
            DuplicateResolution whenDuplicate = DuplicateResolution.MaxAbundance,
            ConsensusQuantMeasure quantMeasure = ConsensusQuantMeasure.SummedIntensityTop3)
        {
            if (runs == null) throw new ArgumentNullException(nameof(runs));
            if (runs.Count < 2) throw new ArgumentException("Need at least 2 runs to build a consensus.");

            var labels = MakeUniqueLabels(runs.Select(r => r.label).ToList());
            var all = new List<CFeature>();
            for (int r = 0; r < runs.Count; r++)
                foreach (var mf in runs[r].features ?? Array.Empty<MassFeature>())
                    all.Add(ToCFeature(mf, r, quantMeasure));

            return AlignCore(all, labels, outputPath, massTolerance, massUnit, rtTolerance,
                positiveNotch, negativeNotch, grouping, whenDuplicate);
        }

        // Shared core: sort by mass, group, write. (Both entry points funnel through here.)
        private static string AlignCore(
            List<CFeature> all, List<string> labels, string outputPath,
            double massTolerance, MassToleranceUnit massUnit, double rtTolerance,
            int positiveNotch, int negativeNotch, GroupingMethod grouping, DuplicateResolution whenDuplicate)
        {
            if (string.IsNullOrWhiteSpace(outputPath)) throw new ArgumentException("OutputPath is not set.");

            // Sort by mass so the matching loop only compares near-mass neighbours (same trick as FeatureAligner).
            all.Sort((a, b) => a.Mass.CompareTo(b.Mass));

            var massTol = FeatureAligner.MakeTolerance(massTolerance, massUnit);
            var shifts  = FeatureAligner.NotchShifts(positiveNotch, negativeNotch);

            var groups = grouping == GroupingMethod.ReferenceAnchor
                ? GroupReferenceAnchor(all, massTol, shifts, rtTolerance, whenDuplicate)
                : GroupUnionFind(all, massTol, shifts, rtTolerance, positiveNotch);

            WriteConsensusTsv(groups, labels, outputPath);
            return outputPath;
        }

        // ====================================================================================================
        //  GROUPING (mirrors FeatureAligner)
        // ====================================================================================================

        // Transitive connected components: link any two features within mass(+notch)+RT tolerance of each other.
        private static List<List<CFeature>> GroupUnionFind(
            List<CFeature> all, Tolerance massTol, double[] shifts, double rtTol, int positiveNotch)
        {
            int n = all.Count;
            var parent = new int[n];
            for (int i = 0; i < n; i++) parent[i] = i;
            double maxPosShift = positiveNotch * FeatureAligner.NotchStep;

            for (int i = 0; i < n; i++)   // sorted by mass, so the inner scan can stop early
            {
                double mi = all[i].Mass, ri = all[i].RtSeconds;
                double upper = massTol.GetMaximumValue(mi + maxPosShift);
                for (int j = i + 1; j < n; j++)
                {
                    if (all[j].Mass > upper) break;
                    if (Math.Abs(all[j].RtSeconds - ri) <= rtTol &&
                        FeatureAligner.WithinNotch(massTol, all[j].Mass, mi, shifts))
                        Union(parent, i, j);
                }
            }

            var clusters = new Dictionary<int, List<CFeature>>();
            for (int i = 0; i < n; i++)
            {
                int root = Find(parent, i);
                if (!clusters.TryGetValue(root, out var list)) { list = new List<CFeature>(); clusters[root] = list; }
                list.Add(all[i]);
            }
            return clusters.Values.ToList();
        }

        // Greedy reference-anchored grouping (OpenMS ConsensusFeatureGroupDetector style): take the most-abundant
        // unused feature as a fixed anchor, gather every unused feature within tolerance OF THE ANCHOR, keep one
        // per run (WhenDuplicate), discard the rest. Group width is capped at the tolerance (no runaway chaining).
        private static List<List<CFeature>> GroupReferenceAnchor(
            List<CFeature> all, Tolerance massTol, double[] shifts, double rtTol, DuplicateResolution whenDuplicate)
        {
            int n = all.Count;
            var used = new bool[n];
            var byAbundance = Enumerable.Range(0, n).OrderByDescending(i => all[i].Primary).ToList();
            var groups = new List<List<CFeature>>();

            foreach (int refIdx in byAbundance)
            {
                if (used[refIdx]) continue;
                var reference = all[refIdx];

                var perRun = new Dictionary<int, List<int>>();
                for (int j = 0; j < n; j++)
                {
                    if (used[j]) continue;
                    var c = all[j];
                    if (Math.Abs(c.RtSeconds - reference.RtSeconds) <= rtTol &&
                        FeatureAligner.WithinNotch(massTol, c.Mass, reference.Mass, shifts))
                    {
                        if (!perRun.TryGetValue(c.RunIndex, out var li)) { li = new List<int>(); perRun[c.RunIndex] = li; }
                        li.Add(j);
                    }
                }

                var group = new List<CFeature>();
                foreach (var li in perRun.Values)
                {
                    int chosen = whenDuplicate == DuplicateResolution.NearestMass
                        ? li.OrderBy(j => Math.Abs(all[j].Mass - reference.Mass)).First()
                        : li.OrderByDescending(j => all[j].Primary).First();
                    group.Add(all[chosen]);
                    foreach (int j in li) used[j] = true;
                }
                groups.Add(group);
            }
            return groups;
        }

        // ---- union-find primitives ----
        private static int Find(int[] parent, int x)
        {
            int root = x;
            while (parent[root] != root) root = parent[root];
            while (parent[x] != root) { int next = parent[x]; parent[x] = root; x = next; } // path compression
            return root;
        }

        private static void Union(int[] parent, int a, int b)
        {
            int ra = Find(parent, a), rb = Find(parent, b);
            if (ra != rb) parent[ra] = rb;
        }

        // ====================================================================================================
        //  READING A FEATURE FILE
        // ====================================================================================================

        // Parse one `_ms1.feature` table into CFeatures. Mass, an apex RT, and an intensity are required; the two
        // top-3 columns are optional and fall back to the summed intensity when the file doesn't carry them (e.g.
        // a plain TopFD table), so the chosen quant measure is always populated.
        private static List<CFeature> ReadFeatureFile(string path, int runIndex, ConsensusQuantMeasure measure)
        {
            var lines = FeatureAligner.ReadAllLinesShared(path);
            if (lines.Length < 2) return new List<CFeature>();

            var header    = lines[0].Split('\t');
            int massCol   = FeatureAligner.FindColumn(header, FMass);
            int rtCol     = FeatureAligner.FindColumn(header, FRt);
            int sumCol    = FeatureAligner.FindColumn(header, FSum);
            int sumT3Col  = FeatureAligner.FindColumn(header, FSumTop3);
            int apexT3Col = FeatureAligner.FindColumn(header, FApexTop3);
            int idCol     = FeatureAligner.FindColumn(header, FId);
            int cMinCol   = FeatureAligner.FindColumn(header, FCMin);
            int cMaxCol   = FeatureAligner.FindColumn(header, FCMax);

            if (massCol < 0 || rtCol < 0 || sumCol < 0)
                throw new FormatException(
                    $"'{Path.GetFileName(path)}' is not a consensus/TopFD `_ms1.feature` table (missing " +
                    $"{(massCol < 0 ? "mass " : "")}{(rtCol < 0 ? "apex-rt " : "")}{(sumCol < 0 ? "intensity" : "")}). " +
                    $"Columns found: {string.Join(", ", header)}");

            double rtScale = RtMinuteColumns.Contains(header[rtCol].Trim()) ? 60.0 : 1.0;
            int maxCol = Math.Max(massCol, Math.Max(rtCol, sumCol));

            var result = new List<CFeature>(lines.Length - 1);
            for (int r = 1; r < lines.Length; r++)
            {
                if (lines[r].Length == 0) continue;
                var t = lines[r].Split('\t');
                if (t.Length <= maxCol) continue;
                if (!double.TryParse(t[massCol], NumberStyles.Any, Ci, out double mass)) continue;
                if (!double.TryParse(t[rtCol],   NumberStyles.Any, Ci, out double rt))   continue;
                double.TryParse(t[sumCol], NumberStyles.Any, Ci, out double sum);           // blank -> 0
                double sumTop3  = ParseOr(t, sumT3Col, sum);                                // absent -> summed intensity
                double apexTop3 = ParseOr(t, apexT3Col, sum);

                double primary = measure switch
                {
                    ConsensusQuantMeasure.SummedIntensity   => sum,
                    ConsensusQuantMeasure.ApexIntensityTop3 => apexTop3,
                    _                                       => sumTop3,
                };

                result.Add(new CFeature
                {
                    RunIndex  = runIndex,
                    Index     = (idCol >= 0 && idCol < t.Length) ? t[idCol] : r.ToString(Ci),
                    Mass      = mass,
                    RtSeconds = rt * rtScale,
                    Sum       = sum,
                    SumTop3   = sumTop3,
                    ApexTop3  = apexTop3,
                    MinCharge = FeatureAligner.ParseInt(t, cMinCol),
                    MaxCharge = FeatureAligner.ParseInt(t, cMaxCol),
                    Primary   = primary,
                });
            }
            return result;
        }

        // Value at col parsed as double, or 'fallback' if the column is absent/short/unparseable.
        private static double ParseOr(string[] t, int col, double fallback)
            => (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? v : fallback;

        // ====================================================================================================
        //  WRITING THE CONSENSUS TABLE
        // ====================================================================================================

        private static void WriteConsensusTsv(IEnumerable<List<CFeature>> clusters, List<string> labels, string outputPath)
        {
            int runCount = labels.Count;
            var rows = new List<Row>();

            foreach (var cluster in clusters)
            {
                var sum      = new double?[runCount];   // per-run SummedIntensity     (null = no feature in that run)
                var sumTop3  = new double?[runCount];   // per-run SummedIntensityTop3
                var apexTop3 = new double?[runCount];   // per-run ApexIntensityTop3
                var idx      = new string[runCount];    // per-run feature id
                int maxPerRun = 1;

                // In a clean cluster each run contributes exactly one feature. If a run contributes several, they
                // are within tolerance in the SAME run (almost always a deconvolution artifact); keep that run's
                // most-abundant one as its representative and record maxPerRun so such clusters can be flagged.
                foreach (var g in cluster.GroupBy(x => x.RunIndex))
                {
                    if (g.Count() > maxPerRun) maxPerRun = g.Count();
                    var best = g.OrderByDescending(x => x.Primary).First();
                    sum[g.Key]      = best.Sum;
                    sumTop3[g.Key]  = best.SumTop3;
                    apexTop3[g.Key] = best.ApexTop3;
                    idx[g.Key]      = best.Index;
                }

                var withCharge = cluster.Where(x => x.MaxCharge > 0).ToList();
                string chargeRange = withCharge.Count > 0
                    ? $"{withCharge.Min(x => x.MinCharge)}-{withCharge.Max(x => x.MaxCharge)}" : "";

                rows.Add(new Row
                {
                    MassMean    = cluster.Average(x => x.Mass),
                    RtMean      = cluster.Average(x => x.RtSeconds),
                    NumFiles    = sum.Count(q => q.HasValue),
                    MaxPerRun   = maxPerRun,
                    ChargeRange = chargeRange,
                    Sum = sum, SumTop3 = sumTop3, ApexTop3 = apexTop3, Idx = idx
                });
            }

            // Most-shared proteoforms first, then by mass.
            rows.Sort((a, b) =>
            {
                int c = b.NumFiles.CompareTo(a.NumFiles);
                return c != 0 ? c : a.MassMean.CompareTo(b.MassMean);
            });

            var dir = Path.GetDirectoryName(outputPath);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir)) Directory.CreateDirectory(dir);

            using var sw = new StreamWriter(outputPath, append: false);

            // Fixed columns, then a block of each quant measure per run, then the per-run feature ids.
            // max_features_per_run == 1 is clean; > 1 flags a same-run duplicate (suspected deconvolution error).
            var header = new List<string> { "mass_mean", "rt_apex_mean", "n_files", "max_features_per_run", "charge_range" };
            header.AddRange(labels.Select(l => "sum_" + l));
            header.AddRange(labels.Select(l => "sumTop3_" + l));
            header.AddRange(labels.Select(l => "apexTop3_" + l));
            header.AddRange(labels.Select(l => "idx_" + l));
            sw.WriteLine(string.Join("\t", header));

            foreach (var row in rows)
            {
                var cells = new List<string>
                {
                    row.MassMean.ToString("F5", Ci),
                    row.RtMean.ToString("F2", Ci),
                    row.NumFiles.ToString(Ci),
                    row.MaxPerRun.ToString(Ci),
                    row.ChargeRange
                };
                cells.AddRange(row.Sum.Select(q => q.HasValue ? q.Value.ToString("F1", Ci) : ""));
                cells.AddRange(row.SumTop3.Select(q => q.HasValue ? q.Value.ToString("F1", Ci) : ""));
                cells.AddRange(row.ApexTop3.Select(q => q.HasValue ? q.Value.ToString("F1", Ci) : ""));
                cells.AddRange(row.Idx.Select(s => s ?? ""));
                sw.WriteLine(string.Join("\t", cells));
            }
        }

        // ====================================================================================================
        //  HELPERS
        // ====================================================================================================

        // Project an in-memory MassFeature onto the fields the aligner needs. ApexRT (minutes) -> seconds.
        private static CFeature ToCFeature(MassFeature mf, int runIndex, ConsensusQuantMeasure measure)
        {
            double primary = measure switch
            {
                ConsensusQuantMeasure.SummedIntensity   => mf.SummedIntensity,
                ConsensusQuantMeasure.ApexIntensityTop3 => mf.ApexIntensityTop3,
                _                                       => mf.SummedIntensityTop3,
            };
            return new CFeature
            {
                RunIndex  = runIndex,
                Index     = mf.Id.ToString(Ci),
                Mass      = mf.ConsensusMass,
                RtSeconds = mf.ApexRT * 60.0,                       // MassFeature.ApexRT is in minutes
                Sum       = mf.SummedIntensity,
                SumTop3   = mf.SummedIntensityTop3,
                ApexTop3  = mf.ApexIntensityTop3,
                MinCharge = mf.Charges.Count > 0 ? mf.Charges.Min() : 0,
                MaxCharge = mf.Charges.Count > 0 ? mf.Charges.Max() : 0,
                Primary   = primary,
            };
        }

        // A short, unique label per run; duplicates get a ".1", ".2"… suffix (matches FeatureAligner's scheme).
        private static List<string> MakeUniqueLabels(List<string> raw)
        {
            var labels = new List<string>();
            var seen = new Dictionary<string, int>();
            foreach (var r in raw)
            {
                string name = string.IsNullOrWhiteSpace(r) ? "run" : r;
                if (seen.TryGetValue(name, out int count)) { seen[name] = count + 1; labels.Add($"{name}.{count}"); }
                else { seen[name] = 1; labels.Add(name); }
            }
            return labels;
        }

        // ---- record types (private) ------------------------------------------------------------------------

        private sealed class CFeature
        {
            public int RunIndex;      // which run this feature came from (index into the runs/labels list)
            public string Index;      // the feature's id in its source file/run
            public double Mass;       // consensus monoisotopic mass
            public double RtSeconds;  // apex retention time, in SECONDS
            public double Sum;        // SummedIntensity (all charges)
            public double SumTop3;    // SummedIntensityTop3
            public double ApexTop3;   // ApexIntensityTop3
            public int MinCharge;     // lowest charge state observed
            public int MaxCharge;     // highest charge state observed
            public double Primary;    // the chosen quant measure, for anchor order + duplicate resolution
        }

        private sealed class Row
        {
            public double MassMean, RtMean;
            public int NumFiles, MaxPerRun;
            public string ChargeRange;
            public double?[] Sum, SumTop3, ApexTop3;
            public string[] Idx;
        }
    }
}
