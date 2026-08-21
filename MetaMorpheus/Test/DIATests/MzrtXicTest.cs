using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// XICs read straight out of TopDIA's mzrt.csv files, for both precursors and fragments.
    ///
    /// The mzrt format stores a per-CYCLE intensity vector of the same length in the MS1 file and in every MS2
    /// window file, so precursor and fragment chromatograms share an axis with no join at all - cycle k is
    /// position k of the scan array. These tests pin that identity down, because everything else depends on it.
    ///
    /// WINDOW-TO-VOLTAGE MAPPING IS NOT THE FILE-NAME ORDER. Each ISD source voltage was given a synthetic DIA
    /// window (+/-373.5 m/z) and TopDIA names MS2 files by the window's LOWER BOUND, which inverts the voltage
    /// order for this run:
    ///     875.5 -> window centre 1249 -> sid=100
    ///     877.5 -> window centre 1251 -> sid=80
    ///     879.5 -> window centre 1253 -> sid=60
    /// </summary>
    [TestFixture]
    public static class MzrtXicTest
    {
        private const string TopDiaDir =
            @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\TopDIA";
        private const string Stem =
            @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_3windows";

        private static string Ms1Mzrt => Path.Combine(TopDiaDir, Stem + "_frac_ms1.mzrt.csv");
        private static string ParentMzml => Path.Combine(TopDiaDir, Stem + ".mzML");

        /// <summary>Fragment window files keyed by the ISD source voltage they actually contain.</summary>
        private static Dictionary<int, string> FragmentMzrtByVoltage => new()
        {
            { 100, Path.Combine(TopDiaDir, Stem + "_875.500000_frac_ms2.mzrt.csv") },
            { 80, Path.Combine(TopDiaDir, Stem + "_877.500000_frac_ms2.mzrt.csv") },
            { 60, Path.Combine(TopDiaDir, Stem + "_879.500000_frac_ms2.mzrt.csv") },
        };

        /// <summary>
        /// Reader-level checks on a real file: rows parse, charge states group under a shared feature id, and the
        /// retained XIC stretch lines up with the specLow/specHi the file reports.
        /// </summary>
        [Test, Explicit]
        public static void ReadsMzrtAndKeepsTheChromatogram()
        {
            string path = FragmentMzrtByVoltage[100];
            Assert.That(File.Exists(path), $"missing {path}");

            var records = TopDiaMzrtFile.ReadFile(path, out int cycleCount);
            TestContext.Progress.WriteLine($"{Path.GetFileName(path)}: {records.Count} rows, " +
                $"{records.Select(r => r.Id).Distinct().Count()} features, {cycleCount} cycles");

            Assert.That(records, Is.Not.Empty);
            Assert.That(cycleCount, Is.GreaterThan(0));

            // the retained stretch must sit inside the reported elution window
            int checkedRows = 0;
            foreach (var record in records.Where(r => r.Intensities.Length > 0).Take(500))
            {
                Assert.That(record.FirstCycle, Is.GreaterThanOrEqualTo(0));
                Assert.That(record.EndCycleExclusive, Is.LessThanOrEqualTo(cycleCount));
                Assert.That(record.Intensities[0], Is.GreaterThan(0), "leading zeros should be trimmed");
                Assert.That(record.Intensities[^1], Is.GreaterThan(0), "trailing zeros should be trimmed");
                checkedRows++;
            }
            TestContext.Progress.WriteLine($"checked {checkedRows} rows: XIC trimmed to its non-zero span, in range");

            var lengths = records.Select(r => r.Intensities.Length).OrderBy(n => n).ToArray();
            TestContext.Progress.WriteLine($"XIC span (cycles): median {lengths[lengths.Length / 2]}, " +
                $"max {lengths[^1]}; non-zero points median " +
                $"{records.Select(r => r.PointCount).OrderBy(n => n).ToArray()[records.Count / 2]}");
        }

        /// <summary>
        /// Full construction over the real run: precursors from the MS1 mzrt, fragments from one mzrt per ISD
        /// channel. Also verifies the cycle axis really is the scan axis, via the apex retention time check.
        /// </summary>
        [Test, Explicit]
        public static void BuildsPrecursorAndFragmentXicsFromMzrt()
        {
            Assert.That(File.Exists(ParentMzml), $"missing {ParentMzml}");

            var dataFile = MsDataFileReader.GetDataFile(ParentMzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);
            TestContext.Progress.WriteLine($"{allScans.Length} scans -> {ms1Scans.Length} MS1 + channels " +
                string.Join(", ", isdMap.OrderBy(k => k.Key).Select(k => $"sid={k.Key} ({k.Value.Count})")));

            // ---- precursors ----
            var ms1Constructor = new MzrtXicConstructor(new PpmTolerance(20), 2, 1, 3, Ms1Mzrt,
                minMass: 4000, maxMass: 70000, minScore: 0);
            var precursorXics = ms1Constructor.GetAllXics(ms1Scans);
            TestContext.Progress.WriteLine($"MS1: {ms1Constructor.LastRecordCount} rows -> " +
                $"{ms1Constructor.LastFeatureCount} features -> {precursorXics.Count} XICs " +
                $"({ms1Constructor.LastBelowMinPeaksCount} below minNumberOfPeaks); " +
                $"RT discrepancy median {ms1Constructor.LastMedianRetentionTimeDiscrepancy:F4} / max {ms1Constructor.LastMaxRetentionTimeDiscrepancy:F3} min");
            Assert.That(precursorXics.Count, Is.GreaterThan(0));

            // the cycle axis IS the scan axis - if it were not, apex RTs would be far from the scans they index
            Assert.That(ms1Constructor.LastMedianRetentionTimeDiscrepancy, Is.LessThan(0.01),
                "elution-start retention times disagree with the scans the cycle indices point at");

            // ---- fragments, one channel at a time ----
            var ms2Constructor = MzrtXicConstructor.BySourceVoltage(new PpmTolerance(20), 2, 1, 3,
                FragmentMzrtByVoltage, maxMass: 70000, minScore: 0);

            var fragmentXics = new List<ExtractedIonChromatogram>();
            foreach (var channel in isdMap.OrderBy(k => k.Key))
            {
                var channelXics = ms2Constructor.GetAllXics(channel.Value.ToArray());
                TestContext.Progress.WriteLine($"MS2 sid={channel.Key}: picked voltage {ms2Constructor.LastSourceVoltage}, " +
                    $"{ms2Constructor.LastRecordCount} rows -> {ms2Constructor.LastFeatureCount} features -> " +
                    $"{channelXics.Count} XICs; RT discrepancy median " +
                    $"{ms2Constructor.LastMedianRetentionTimeDiscrepancy:F4} / max {ms2Constructor.LastMaxRetentionTimeDiscrepancy:F3} min");
                Assert.That(ms2Constructor.LastSourceVoltage, Is.EqualTo(channel.Key));
                Assert.That(ms2Constructor.LastMedianRetentionTimeDiscrepancy, Is.LessThan(0.01));
                fragmentXics.AddRange(channelXics);
            }
            TestContext.Progress.WriteLine($"fragment XICs pooled: {fragmentXics.Count}");
            Assert.That(fragmentXics.Count, Is.GreaterThan(0));

            var lengths = precursorXics.Select(x => x.Peaks.Count).OrderBy(n => n).ToArray();
            TestContext.Progress.WriteLine($"precursor XIC length: median {lengths[lengths.Length / 2]}, max {lengths[^1]}");
        }

        /// <summary>
        /// A scan array that does not cover the same cycles as the mzrt file must fail loudly - silently building
        /// chromatograms on a shifted time axis would be far worse than an exception.
        /// </summary>
        [Test, Explicit]
        public static void MismatchedCycleCountThrows()
        {
            var dataFile = MsDataFileReader.GetDataFile(ParentMzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);

            var constructor = new MzrtXicConstructor(new PpmTolerance(20), 2, 1, 3, Ms1Mzrt);
            var truncated = ms1Scans.Take(ms1Scans.Length - 5).ToArray();

            var ex = Assert.Throws<MetaMorpheusException>(() => constructor.GetAllXics(truncated));
            TestContext.Progress.WriteLine("rejected as expected: " + ex.Message);
            Assert.That(ex.Message, Does.Contain("cycles"));
        }
    }
}
