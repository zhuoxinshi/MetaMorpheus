using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA.FeatureAlign
{
    /// <summary>
    /// Links each MetaMorpheus identification to the TopFD (or FLASHDeconv) MS1 feature it came from, WITHIN one
    /// run, by matching the identification's precursor to a feature's neutral mass and retention time.
    ///
    /// MetaMorpheus does not write a feature id, so unlike TopPIC (whose <c>Feature ID</c> column is an exact
    /// integer link into the same <c>*_ms1.feature</c> file) there is no direct join. This aligner reconstructs it
    /// the way the TopPIC "combined proteoform to aligned feature" join does: mass within a ppm tolerance, allowing
    /// a small number of C13 monoisotope notches, and retention time inside the feature's own elution window.
    ///
    /// Both sides are already in the SAME units - MetaMorpheus <c>Scan Retention Time</c> and
    /// <see cref="Ms1FeatureRecord"/> retention times are both MINUTES - so no conversion is applied (this is the
    /// seconds-vs-minutes trap that breaks the FeatureAlign-table join, avoided here because the feature file is
    /// read directly).
    ///
    /// One row is written per identification, matched or not, so the output doubles as an ID-level table carrying
    /// the feature's mass, elution window and intensity for every ID that could be placed.
    /// </summary>
    public static class IdToFeatureAligner
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;
        private const double C13 = 1.0033548381;   // 13C-12C mass gap, as used by FeatureAligner / the TopPIC join

        // ---- MetaMorpheus psmtsv column dictionaries (auto-detected, so column order does not matter) ----------
        private static readonly string[] FileNameNames    = { "File Name" };
        private static readonly string[] ScanNames        = { "Scan Number" };
        private static readonly string[] RtNames          = { "Scan Retention Time", "Retention Time" };   // MINUTES
        private static readonly string[] PrecursorMassNms = { "Precursor Mass" };                           // observed, neutral
        private static readonly string[] MonoMassNames    = { "Monoisotopic Mass" };                        // theoretical
        private static readonly string[] ChargeNames      = { "Precursor Charge" };
        private static readonly string[] ScoreNames       = { "Score" };
        private static readonly string[] BaseSeqNames     = { "Base Sequence" };
        private static readonly string[] FullSeqNames     = { "Full Sequence" };
        private static readonly string[] AccNames         = { "Accession", "Protein Accession" };
        private static readonly string[] GeneNames        = { "Gene Name" };
        private static readonly string[] QValNames        = { "QValue" };
        private static readonly string[] DctNames         = { "Decoy/Contaminant/Target" };
        private static readonly string[] DecoyNames       = { "Decoy" };
        private static readonly string[] ContamNames      = { "Contaminant" };

        /// <summary>
        /// Reads the MetaMorpheus psmtsv and the TopFD <c>*_ms1.feature</c>, links every surviving identification to
        /// its best-matching feature, writes the per-identification table, and returns a summary.
        /// </summary>
        /// <param name="metaMorpheusPsmtsv">A MetaMorpheus <c>*_Proteoforms.psmtsv</c> or <c>*_PSMs.psmtsv</c> for
        /// ONE run (a pooled file with a "File Name" column also works, but the feature file must belong to those
        /// same raw scans).</param>
        /// <param name="topFdFeatureFile">The run's <c>*_ms1.feature</c> (TopFD 1.7/1.8 or FLASHDeconv - column
        /// spellings are resolved by <see cref="Ms1FeatureRecord"/>).</param>
        /// <param name="outputPath">tsv to write (one row per identification).</param>
        /// <param name="massTolerance">Mass match tolerance; default 20 ppm.</param>
        /// <param name="retentionTimeToleranceMinutes">Slack added to each end of the feature's elution window
        /// (<c>matchWithinRtWindow = true</c>), or the half-width of the apex window (<c>false</c>). Default 0 -
        /// the feature windows are already wide.</param>
        /// <param name="matchWithinRtWindow">true (default): the ID's RT must fall inside
        /// [feature begin - tol, feature end + tol]. false: |feature apex - ID RT| must be within the tolerance.</param>
        /// <param name="monoisotopicNotches">C13 notches allowed on the mass match; default -2..2 (matches
        /// TopFD's own notch handling). Pass {0} to require an exact monoisotope.</param>
        /// <param name="useTheoreticalMass">false (default): match the observed <c>Precursor Mass</c> to the
        /// feature's observed mass (they are the same ion, so only C13 notches separate them). true: match the
        /// theoretical <c>Monoisotopic Mass</c> instead (needed only if precursor mass is missing).</param>
        /// <param name="requireChargeInFeatureRange">If true, a feature only matches when the ID's precursor charge
        /// lies within the feature's observed charge range.</param>
        /// <param name="qValueThreshold">Drop identifications with QValue above this (default 0.01). Rows with no
        /// QValue column are kept.</param>
        /// <param name="includeDecoys">Keep decoy identifications (default false).</param>
        /// <param name="includeContaminants">Keep contaminant identifications (default false).</param>
        public static AlignmentSummary Align(string metaMorpheusPsmtsv, string topFdFeatureFile, string outputPath,
            Tolerance massTolerance = null, double retentionTimeToleranceMinutes = 0.0, bool matchWithinRtWindow = true,
            IEnumerable<int> monoisotopicNotches = null, bool useTheoreticalMass = false,
            bool requireChargeInFeatureRange = false, double qValueThreshold = 0.01,
            bool includeDecoys = false, bool includeContaminants = false)
        {
            if (string.IsNullOrWhiteSpace(metaMorpheusPsmtsv) || !File.Exists(metaMorpheusPsmtsv))
                throw new MetaMorpheusException("MetaMorpheus psmtsv not found: " + metaMorpheusPsmtsv);
            if (string.IsNullOrWhiteSpace(outputPath))
                throw new ArgumentException("OutputPath is not set.");

            massTolerance ??= new PpmTolerance(20);
            var notches = (monoisotopicNotches ?? new[] { -2, -1, 0, 1, 2 })
                .Distinct().OrderBy(Math.Abs).ThenBy(n => n).ToArray();

            var ids = ReadIdentifications(metaMorpheusPsmtsv, useTheoreticalMass, qValueThreshold,
                includeDecoys, includeContaminants);

            // Features indexed by ascending mass so each ID only scans the mass slice a notch could reach, rather
            // than the whole file - the raw runs carry ~24k features.
            var features = Ms1FeatureRecord.ReadFile(topFdFeatureFile).OrderBy(f => f.Mass).ToList();
            var featureMasses = features.Select(f => f.Mass).ToArray();

            var summary = new AlignmentSummary { TotalIdentifications = ids.Count, FeatureCount = features.Count };
            var notchCounts = new SortedDictionary<int, int>();

            var directory = Path.GetDirectoryName(outputPath);
            if (!string.IsNullOrEmpty(directory)) Directory.CreateDirectory(directory);
            using var writer = new StreamWriter(outputPath, false);
            writer.WriteLine(string.Join("\t",
                "File Name", "Scan Number", "Accession", "Gene", "Base Sequence", "Full Sequence", "QValue",
                "Precursor Charge", "Id Mass", "Scan Retention Time", "Matched",
                "Feature_ID", "Feature_Mass", "Mass_ppm_error", "C13_notch",
                "Feature_RT_begin", "Feature_RT_end", "Feature_RT_apex", "Feature_Intensity",
                "Feature_Charge_min", "Feature_Charge_max"));

            foreach (var id in ids)
            {
                var best = FindBestFeature(id, features, featureMasses, massTolerance, notches,
                    retentionTimeToleranceMinutes, matchWithinRtWindow, requireChargeInFeatureRange);

                if (best.Feature != null)
                {
                    summary.MatchedIdentifications++;
                    notchCounts.TryGetValue(best.Notch, out int c);
                    notchCounts[best.Notch] = c + 1;
                }

                var f = best.Feature;
                writer.WriteLine(string.Join("\t",
                    id.FileName, id.Scan, id.Accession, id.Gene, id.BaseSeq, id.FullSeq,
                    Fmt(id.QValue), id.Charge, Fmt(id.Mass), Fmt(id.Rt), f != null ? "Y" : "N",
                    f?.Id.ToString(Ci) ?? "", f != null ? Fmt(f.Mass) : "",
                    f != null ? best.PpmError.ToString("F2", Ci) : "", f != null ? best.Notch.ToString(Ci) : "",
                    f != null ? Fmt(f.RetentionTimeBegin) : "", f != null ? Fmt(f.RetentionTimeEnd) : "",
                    f != null ? Fmt(f.RetentionTimeApex) : "", f != null ? Fmt(f.Intensity) : "",
                    f != null ? f.ChargeStateMin.ToString(Ci) : "", f != null ? f.ChargeStateMax.ToString(Ci) : ""));
            }

            summary.NotchCounts = notchCounts;
            summary.OutputPath = outputPath;
            return summary;
        }

        // For one identification, the feature that best explains its precursor: the smallest-magnitude C13 notch
        // whose shifted feature mass lands within tolerance and whose elution window covers the ID's RT; ties on
        // notch are broken by smallest ppm error, then by the more intense feature.
        private static (Ms1FeatureRecord Feature, int Notch, double PpmError) FindBestFeature(Identification id,
            List<Ms1FeatureRecord> features, double[] featureMasses, Tolerance massTolerance, int[] notches,
            double rtTol, bool matchWithinRtWindow, bool requireChargeInFeatureRange)
        {
            Ms1FeatureRecord best = null;
            int bestNotch = 0;
            double bestPpm = double.MaxValue;

            // Only a handful of notches, so scan them all and keep the single best candidate under one total order:
            // smallest |notch|, then smallest ppm error, then the more intense feature. No early-exit heuristics.
            foreach (int notch in notches)
            {
                // A feature at mass M read as this proteoform means M ~= idMass - notch*C13, i.e. the ID sits
                // notch C13 units above the feature. Search the feature-mass band that could satisfy the tolerance.
                double centre = id.Mass - notch * C13;
                double window = massTolerance.GetMaximumValue(centre) - centre + 1e-6;
                int lo = LowerBound(featureMasses, centre - window);

                for (int i = lo; i < features.Count && featureMasses[i] <= centre + window; i++)
                {
                    var feature = features[i];
                    if (!massTolerance.Within(id.Mass, feature.Mass + notch * C13))
                        continue;

                    bool rtOk = matchWithinRtWindow
                        ? feature.ContainsRetentionTime(id.Rt, rtTol)
                        : Math.Abs(feature.RetentionTimeApex - id.Rt) <= rtTol;
                    if (!rtOk)
                        continue;

                    if (requireChargeInFeatureRange && id.Charge > 0 && !feature.ContainsCharge(id.Charge))
                        continue;

                    double ppm = Math.Abs((id.Mass - (feature.Mass + notch * C13)) / (feature.Mass + notch * C13)) * 1e6;

                    if (IsBetter(notch, ppm, feature.Intensity, best, bestNotch, bestPpm))
                    {
                        best = feature; bestNotch = notch; bestPpm = ppm;
                    }
                }
            }

            return (best, best == null ? 0 : bestNotch, best == null ? double.NaN : bestPpm);
        }

        // Total order on candidate matches: fewer C13 notches first, then smaller ppm error, then more intense.
        private static bool IsBetter(int notch, double ppm, double intensity,
            Ms1FeatureRecord best, int bestNotch, double bestPpm)
        {
            if (best == null) return true;
            int dNotch = Math.Abs(notch) - Math.Abs(bestNotch);
            if (dNotch != 0) return dNotch < 0;
            if (Math.Abs(ppm - bestPpm) > 1e-9) return ppm < bestPpm;
            return intensity > best.Intensity;
        }

        // First index whose value is >= target, on an ascending array (binary search).
        private static int LowerBound(double[] sorted, double target)
        {
            int lo = 0, hi = sorted.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (sorted[mid] < target) lo = mid + 1; else hi = mid;
            }
            return lo;
        }

        // ---- MetaMorpheus reader -------------------------------------------------------------------------------
        private static List<Identification> ReadIdentifications(string path, bool useTheoreticalMass,
            double qValueThreshold, bool includeDecoys, bool includeContaminants)
        {
            var lines = File.ReadAllLines(path);
            if (lines.Length < 2)
                throw new MetaMorpheusException("MetaMorpheus psmtsv has no rows: " + path);

            var header = lines[0].Split('\t');
            int cFile = Find(header, FileNameNames), cScan = Find(header, ScanNames), cRt = Find(header, RtNames);
            int cPrec = Find(header, PrecursorMassNms), cMono = Find(header, MonoMassNames), cChg = Find(header, ChargeNames);
            int cBase = Find(header, BaseSeqNames), cFull = Find(header, FullSeqNames), cAcc = Find(header, AccNames);
            int cGene = Find(header, GeneNames), cQ = Find(header, QValNames);
            int cDct = Find(header, DctNames), cDecoy = Find(header, DecoyNames), cContam = Find(header, ContamNames);
            int cMass = useTheoreticalMass ? cMono : cPrec;

            if (cAcc < 0 || cRt < 0 || cMass < 0)
                throw new FormatException($"'{Path.GetFileName(path)}' is missing a required column " +
                    $"({(cAcc < 0 ? "Accession " : "")}{(cRt < 0 ? "Scan Retention Time " : "")}" +
                    $"{(cMass < 0 ? (useTheoreticalMass ? "Monoisotopic Mass" : "Precursor Mass") : "")}). " +
                    $"Columns: {string.Join(", ", header)}");

            var ids = new List<Identification>(lines.Length - 1);
            for (int r = 1; r < lines.Length; r++)
            {
                if (lines[r].Length == 0) continue;
                var t = lines[r].Split('\t');

                if (!includeDecoys && IsFlag(t, cDct, "D", cDecoy)) continue;
                if (!includeContaminants && IsFlag(t, cDct, "C", cContam)) continue;

                double q = ParseDouble(t, cQ);
                if (!double.IsNaN(q) && q > qValueThreshold) continue;

                double mass = ParseDouble(t, cMass);
                double rt = ParseDouble(t, cRt);
                if (double.IsNaN(mass) || double.IsNaN(rt)) continue;   // cannot be placed without both

                ids.Add(new Identification
                {
                    FileName = Cell(t, cFile), Scan = Cell(t, cScan),
                    Accession = Cell(t, cAcc), Gene = Cell(t, cGene),
                    BaseSeq = Cell(t, cBase), FullSeq = Cell(t, cFull),
                    Mass = mass, Rt = rt, Charge = ParseInt(t, cChg), QValue = q
                });
            }
            return ids;
        }

        private static bool IsFlag(string[] t, int cDct, string dctLetter, int cBool)
        {
            string dct = Cell(t, cDct);
            if (dct.Length > 0) return string.Equals(dct, dctLetter, StringComparison.OrdinalIgnoreCase);
            string b = Cell(t, cBool);
            return string.Equals(b, "true", StringComparison.OrdinalIgnoreCase) || b == "1";
        }

        private static int Find(string[] header, string[] candidates)
        {
            foreach (var c in candidates)
            {
                int idx = Array.FindIndex(header, h => string.Equals(h.Trim(), c, StringComparison.OrdinalIgnoreCase));
                if (idx >= 0) return idx;
            }
            return -1;
        }

        private static string Cell(string[] t, int col) => (col >= 0 && col < t.Length) ? t[col].Trim() : "";
        private static double ParseDouble(string[] t, int col) =>
            (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? v : double.NaN;
        private static int ParseInt(string[] t, int col) =>
            (col >= 0 && col < t.Length && double.TryParse(t[col], NumberStyles.Any, Ci, out double v)) ? (int)Math.Round(v) : 0;
        private static string Fmt(double v) => double.IsNaN(v) ? "" : v.ToString("R", Ci);

        private sealed class Identification
        {
            public string FileName, Scan, Accession, Gene, BaseSeq, FullSeq;
            public double Mass;    // observed precursor (or theoretical) monoisotopic mass
            public double Rt;      // scan retention time, MINUTES
            public int Charge;
            public double QValue;
        }

        /// <summary>Counts from the most recent <see cref="Align"/> call.</summary>
        public sealed class AlignmentSummary
        {
            public int TotalIdentifications { get; set; }
            public int MatchedIdentifications { get; set; }
            public int FeatureCount { get; set; }
            public string OutputPath { get; set; }
            /// <summary>Matched identifications by the C13 notch that placed them.</summary>
            public SortedDictionary<int, int> NotchCounts { get; set; } = new SortedDictionary<int, int>();
            public double MatchedFraction => TotalIdentifications == 0 ? 0 : (double)MatchedIdentifications / TotalIdentifications;
        }
    }
}
