using EngineLayer.DIA;
using EngineLayer.DIA.FeatureAlign;
using NUnit.Framework;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.DIATests.Manuscript
{
    /// <summary>
    /// Self-test for <see cref="FeatureAligner"/> on THREE synthetic TopFD .feature tables with a known
    /// ground-truth grouping. Not tied to any data on disk, so it runs in the normal suite and documents exactly
    /// what "correct" alignment means here. Design (20 ppm, +/-2 notch, 60 s, ReferenceAnchor):
    ///
    ///   A  present in all 3 runs, identical mass + RT             -> one group, n_files = 3
    ///   B  present in all 3 runs, masses 1 Da apart (monoisotope) -> one group, n_files = 3  (notch test)
    ///   C  present in runs 1 and 2 only                           -> one group, n_files = 2
    ///   D  run 1 only                                             -> singleton, n_files = 1
    ///   G  run 3 only                                             -> singleton, n_files = 1
    ///   H  same mass in runs 1 and 2 but 3 min apart in RT        -> TWO singletons (RT gate), n_files = 1 each
    ///
    /// So: 12 input features -> 7 consensus rows; n_files histogram {3:2, 2:1, 1:4}.
    /// </summary>
    [TestFixture]
    public class FeatureAlignerSelfTest
    {
        private const double C13 = FeatureAligner.NotchStep; // 1.00335... Da

        private const string Header =
            "File_name\tFraction_ID\tFeature_ID\tMass\tIntensity\tMin_time\tMax_time\tMin_scan\tMax_scan\t" +
            "Min_charge\tMax_charge\tApex_time\tApex_scan\tApex_intensity\tRep_charge\tRep_average_mz\t" +
            "Envelope_num\tEC_score";

        // One TopFD row. Only Mass, Intensity, Apex_time (MINUTES), Feature_ID and the charge range are read by
        // the aligner; the rest are filled plausibly. Apex_time is in minutes because TopFD reports minutes and
        // the aligner multiplies that column by 60 to get seconds.
        private static string Row(string file, int id, double mass, double apexMin, double intensity,
                                  int minZ = 5, int maxZ = 10)
        {
            string I(double v) => v.ToString(CultureInfo.InvariantCulture);
            return string.Join("\t",
                file, 0, id, I(mass), I(intensity),
                I(apexMin - 0.2), I(apexMin + 0.2), 100, 120, minZ, maxZ,
                I(apexMin), 110, I(intensity * 0.3), (minZ + maxZ) / 2,
                I(mass / ((minZ + maxZ) / 2)), 50, "0.99");
        }

        [Test]
        public void ThreeRun_KnownGrouping_AlignsAsExpected()
        {
            var dir = Path.Combine(Path.GetTempPath(), "FeatureAlignerSelfTest");
            Directory.CreateDirectory(dir);
            string f1 = Path.Combine(dir, "run1_ms1.feature");
            string f2 = Path.Combine(dir, "run2_ms1.feature");
            string f3 = Path.Combine(dir, "run3_ms1.feature");

            // run1: A, B(8000.0), C, D, H@25.0 min
            File.WriteAllLines(f1, new[] { Header,
                Row("run1", 1, 5000.0,   20.0, 1_000_000),   // A
                Row("run1", 2, 8000.0,   30.0, 2_000_000),   // B anchor (most abundant)
                Row("run1", 3, 12000.0,  40.0,   800_000),   // C
                Row("run1", 4, 15000.0,  50.0,   500_000),   // D  (singleton)
                Row("run1", 5, 20000.0,  25.0,   400_000),   // H run1
            });
            // run2: A, B(+1 Da), C, H@28.0 min (3 min from H run1 -> must NOT group)
            File.WriteAllLines(f2, new[] { Header,
                Row("run2", 1, 5000.0,          20.0, 1_100_000),        // A
                Row("run2", 2, 8000.0 + C13,    30.0, 1_500_000),        // B +1 monoisotope
                Row("run2", 3, 12000.0,         40.0,   900_000),        // C
                Row("run2", 4, 20000.0,         28.0,   380_000),        // H run2 (RT gate -> separate)
            });
            // run3: A, B(-1 Da), G
            File.WriteAllLines(f3, new[] { Header,
                Row("run3", 1, 5000.0,          20.0,   900_000),        // A
                Row("run3", 2, 8000.0 - C13,    30.0, 1_400_000),        // B -1 monoisotope
                Row("run3", 3, 18000.0,         60.0,   300_000),        // G (singleton)
            });

            string outPath = Path.Combine(dir, "consensus.tsv");
            FeatureAligner.Align(
                new[] { f1, f2, f3 }, outPath,
                massTolerance: 20, massUnit: MassToleranceUnit.Ppm, rtTolerance: 60,
                preferredQuantColumn: "Intensity",
                positiveNotch: 2, negativeNotch: 2,
                grouping: GroupingMethod.ReferenceAnchor);

            var lines = File.ReadAllLines(outPath);
            var header = lines[0].Split('\t');
            int nfCol = System.Array.IndexOf(header, "n_files");
            int massCol = System.Array.IndexOf(header, "mass_mean");
            var rows = lines.Skip(1).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
            var nFiles = rows.Select(r => int.Parse(r[nfCol])).ToList();

            TestContext.Out.WriteLine($"consensus rows: {rows.Count}");
            TestContext.Out.WriteLine("n_files histogram: " + string.Join("  ",
                nFiles.GroupBy(x => x).OrderBy(g => g.Key).Select(g => $"{g.Key}:{g.Count()}")));

            // ground truth
            Assert.That(rows.Count, Is.EqualTo(7), "expected 7 consensus features");
            Assert.That(nFiles.Count(n => n == 3), Is.EqualTo(2), "A and B should each span all 3 runs");
            Assert.That(nFiles.Count(n => n == 2), Is.EqualTo(1), "C should span 2 runs");
            Assert.That(nFiles.Count(n => n == 1), Is.EqualTo(4), "D, G, and the two RT-split H's are singletons");

            // the notch group B must be one all-3 group centred near 8000 (mean of 8000, 8001.003, 7998.997)
            bool bGroupedAcross3 = rows.Any(r =>
                int.Parse(r[nfCol]) == 3 &&
                System.Math.Abs(double.Parse(r[massCol], CultureInfo.InvariantCulture) - 8000.0) < 0.01);
            Assert.That(bGroupedAcross3, Is.True,
                "the +/-1 Da monoisotope features must be merged into ONE all-3 group by the notch tolerance");

            // A must be a clean all-3 group at 5000
            var masses = rows.Select(r => double.Parse(r[massCol], CultureInfo.InvariantCulture)).ToList();
            Assert.That(masses.Any(m => System.Math.Abs(m - 5000.0) < 1e-6), Is.True);

            TestContext.Out.WriteLine("PASS: grouping, notch merge, and RT gate all behaved as designed.");
        }

        [Test]
        public void Notch0_SplitsTheMonoisotopeGroup()
        {
            // Same B triplet, but with notch = 0 the +/-1 Da features must NOT merge: 3 separate singletons.
            var dir = Path.Combine(Path.GetTempPath(), "FeatureAlignerSelfTest_notch0");
            Directory.CreateDirectory(dir);
            string f1 = Path.Combine(dir, "r1_ms1.feature");
            string f2 = Path.Combine(dir, "r2_ms1.feature");
            string f3 = Path.Combine(dir, "r3_ms1.feature");
            File.WriteAllLines(f1, new[] { Header, Row("r1", 1, 8000.0,       30.0, 2_000_000) });
            File.WriteAllLines(f2, new[] { Header, Row("r2", 1, 8000.0 + C13, 30.0, 1_500_000) });
            File.WriteAllLines(f3, new[] { Header, Row("r3", 1, 8000.0 - C13, 30.0, 1_400_000) });

            string outPath = Path.Combine(dir, "consensus.tsv");
            FeatureAligner.Align(new[] { f1, f2, f3 }, outPath,
                massTolerance: 20, massUnit: MassToleranceUnit.Ppm, rtTolerance: 60,
                preferredQuantColumn: "Intensity", positiveNotch: 0, negativeNotch: 0,
                grouping: GroupingMethod.ReferenceAnchor);

            var lines = File.ReadAllLines(outPath);
            int nfCol = System.Array.IndexOf(lines[0].Split('\t'), "n_files");
            var nFiles = lines.Skip(1).Where(l => l.Length > 0).Select(l => int.Parse(l.Split('\t')[nfCol])).ToList();
            TestContext.Out.WriteLine("notch=0 n_files: " + string.Join(",", nFiles));
            Assert.That(nFiles.Count, Is.EqualTo(3), "with notch 0 the 1 Da-apart features stay separate");
            Assert.That(nFiles.All(n => n == 1), Is.True, "each is its own singleton at notch 0");
        }

        /// <summary>
        /// Verifies that <see cref="GroupingMethod.OpenMSConsensus"/> reproduces OpenMS TopDownConsensusFeatureGroup
        /// (the ConsensusFeatureGroupDetector exe) on the same 3 FLASHDeconv mix files. OpenMS at 100 ppm / 180 s
        /// produced 151 consensus groups (all present in all 3 runs); the faithful port should match. [Explicit]
        /// (reads real data on disk).
        /// </summary>
        [TestFixture]
        public class OpenMSConsensusPortTest
        {
            [Test, Explicit("checks OpenMSConsensus grouping == OpenMS's 151 groups on the FLASHDeconv mix files")]
            public void OpenMSConsensus_Matches_OpenMS()
            {
                var fd = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\FeatureFiles";
                var files = new[]
                {
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep1.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep2.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep3.tsv",
            }.Select(f => Path.Combine(fd, f)).ToList();
                Assert.That(files.All(File.Exists), "FLASHDeconv FeatureFiles must be present");

                string outPath = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\_compare\ours_openms_port_consensus.tsv";
                FeatureAligner.Align(files, outPath,
                    massTolerance: 100, massUnit: MassToleranceUnit.Ppm, rtTolerance: 180,
                    preferredQuantColumn: "FeatureQuantity",
                    grouping: GroupingMethod.OpenMSConsensus,
                    whenDuplicate: DuplicateResolution.MaxAbundance);

                var lines = File.ReadAllLines(outPath);
                int nf = System.Array.IndexOf(lines[0].Split('\t'), "n_files");
                var nFiles = lines.Skip(1).Where(l => l.Length > 0).Select(l => int.Parse(l.Split('\t')[nf])).ToList();
                int total = nFiles.Count;
                int all3 = nFiles.Count(n => n == 3);

                TestContext.Out.WriteLine($"OpenMSConsensus rows: {total}   all-3: {all3}   (OpenMS produced 151)");
                Assert.That(nFiles.All(n => n == 3), Is.True, "OpenMSConsensus must emit only all-run groups");
                Assert.That(total, Is.EqualTo(151).Within(3),
                    "should match OpenMS's 151 consensus groups (small slack for tie-breaks)");
            }
        }

        /// <summary>
        /// Runs OUR ReferenceAnchor aligner on the same 3 FLASHDeconv mix files that OpenMS
        /// ConsensusFeatureGroupDetector was run on, with OpenMS-matched tolerances (100 ppm, 180 s, notch 0),
        /// so the two consensus tables can be diffed. [Explicit].
        /// </summary>
        [Test, Explicit("head-to-head vs OpenMS ConsensusFeatureGroupDetector on 3 FLASHDeconv mix files")]
        public void OursVsOpenMS_ThreeMixFlashDeconv()
        {
            var fd = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\FeatureFiles";
            var files = new[]
            {
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep1.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep2.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep3.tsv",
            }.Select(f => Path.Combine(fd, f)).ToList();

            string outPath = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\_compare\ours_consensus.tsv";
            FeatureAligner.Align(files, outPath,
                massTolerance: 100, massUnit: MassToleranceUnit.Ppm, rtTolerance: 90,
                preferredQuantColumn: "FeatureQuantity",
                positiveNotch: 0, negativeNotch: 0,
                grouping: GroupingMethod.ReferenceAnchor, whenDuplicate: DuplicateResolution.MaxAbundance);

            var lines = File.ReadAllLines(outPath);
            int nf = System.Array.IndexOf(lines[0].Split('\t'), "n_files");
            var dist = lines.Skip(1).Where(l => l.Length > 0).Select(l => int.Parse(l.Split('\t')[nf]))
                            .GroupBy(x => x).OrderBy(g => g.Key).Select(g => $"{g.Key}:{g.Count()}");
            TestContext.Out.WriteLine($"OURS consensus rows: {lines.Length - 1}");
            TestContext.Out.WriteLine("n_files: " + string.Join("  ", dist));
        }
    }
}
