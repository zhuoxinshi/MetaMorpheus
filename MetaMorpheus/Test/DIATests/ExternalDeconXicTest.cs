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
    /// Tests for reading EXTERNAL deconvolution results (TopFD msAlign) and constructing neutral-mass XICs
    /// from them with the ordinary MetaMorpheus XIC machinery - i.e. external decon in, normal XICs out.
    /// </summary>
    [TestFixture]
    public static class ExternalDeconXicTest
    {
        /// <summary>
        /// Self-contained: writes a small TopFD-style msAlign (4 columns, RT in seconds) covering 6 scans that
        /// contain two co-eluting masses, indexes it against matching fake scans, and checks that the two masses
        /// come back as two XICs of the expected length, mass and intensity. Nothing on E:\ is touched.
        /// </summary>
        [Test]
        public static void ExternalDeconXicsFromWrittenMsAlign()
        {
            const double massA = 8500.1234;
            const double massB = 12750.5678;
            int scanCount = 6;

            // TopFD writes RETENTION_TIME in SECONDS; the reader converts to minutes. Fake scans are built on
            // the same time base so the scan-number join and the RT axis agree.
            double RtMinutes(int oneBasedScanNumber) => 10.0 + (oneBasedScanNumber - 1) * 0.01;

            string msAlignPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ExternalDeconXicTest_ms1.msalign");
            using (var writer = new StreamWriter(msAlignPath))
            {
                writer.WriteLine("#TopFD 1.7.4");
                writer.WriteLine("####################### Parameters ######################");
                writer.WriteLine("#Spectral data type:         \tCentroid");
                writer.WriteLine("####################### Parameters ######################");
                writer.WriteLine();
                for (int scanNumber = 1; scanNumber <= scanCount; scanNumber++)
                {
                    writer.WriteLine("BEGIN IONS");
                    writer.WriteLine("FRACTION_ID=0");
                    writer.WriteLine("SPECTRUM_ID=" + (scanNumber - 1));
                    writer.WriteLine("TITLE=Scan_" + scanNumber);
                    writer.WriteLine("SCANS=" + scanNumber);
                    writer.WriteLine("RETENTION_TIME=" + (RtMinutes(scanNumber) * 60.0).ToString("F2"));
                    writer.WriteLine("LEVEL=1");
                    // 4-column TopFD peak lines: mass, intensity, charge, ECScore
                    writer.WriteLine(massA.ToString("F5") + "\t" + (1000.0 * scanNumber).ToString("F2") + "\t10\t0.85");
                    writer.WriteLine(massB.ToString("F5") + "\t" + (500.0 * scanNumber).ToString("F2") + "\t15\t0.72");
                    writer.WriteLine("END IONS");
                    writer.WriteLine();
                }
            }

            // The raw scans themselves carry no deconvoluted content - all masses come from the msAlign.
            var scans = new MsDataScan[scanCount];
            for (int i = 0; i < scanCount; i++)
            {
                scans[i] = new MsDataScan(new MzSpectrum(new double[] { 700 }, new double[] { 1 }, false),
                    oneBasedScanNumber: i + 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                    retentionTime: RtMinutes(i + 1), scanWindowRange: new MzRange(200, 2000), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1, injectionTime: 1.0, noiseData: null,
                    nativeId: "scan=" + (i + 1));
            }

            var results = ExternalDeconResults.ReadMsAlign(msAlignPath);
            Assert.That(results.ScanCount, Is.EqualTo(scanCount));
            Assert.That(results.PeakCount, Is.EqualTo(scanCount * 2));
            Assert.That(results.PeaksByOneBasedScanNumber[1][0].MonoisotopicMass, Is.EqualTo(massA).Within(1e-4));
            Assert.That(results.PeaksByOneBasedScanNumber[1][0].Charge, Is.EqualTo(10));
            // msAlign RT is in seconds, the reader must have converted it to minutes
            Assert.That(results.RetentionTimeByOneBasedScanNumber[1], Is.EqualTo(RtMinutes(1)).Within(1e-3));

            var constructor = new ExternalDeconXicConstructor(new PpmTolerance(20), maxMissedScansAllowed: 1,
                maxPeakHalfWidth: 5, minNumberOfPeaks: 3, results, maxMass: 30000);
            var xics = constructor.GetAllXics(scans);

            Assert.That(constructor.LastMatchedScanCount, Is.EqualTo(scanCount));
            Assert.That(constructor.LastIndexedPeakCount, Is.EqualTo(scanCount * 2));
            Assert.That(xics.Count, Is.EqualTo(2), "expected one XIC per co-eluting mass");

            var xicA = xics.OrderBy(x => x.ApexPeak.M).Last();
            var xicB = xics.OrderBy(x => x.ApexPeak.M).First();
            Assert.That(xicB.ApexPeak.M, Is.EqualTo(massA).Within(0.01));
            Assert.That(xicA.ApexPeak.M, Is.EqualTo(massB).Within(0.01));

            foreach (var xic in xics)
            {
                Assert.That(xic.Peaks.Count, Is.EqualTo(scanCount));
                Assert.That(xic.Peaks.Select(p => p.ZeroBasedScanIndex), Is.EqualTo(Enumerable.Range(0, scanCount)));
                // intensity rises with scan number in the written file, so the apex is the last scan
                Assert.That(xic.ApexPeak.ZeroBasedScanIndex, Is.EqualTo(scanCount - 1));
            }
            // intensities are carried through from the msAlign, not recomputed from the (empty) raw spectra
            Assert.That(xicB.ApexPeak.Intensity, Is.EqualTo(1000.0 * scanCount).Within(0.1));
            Assert.That(xicA.ApexPeak.Intensity, Is.EqualTo(500.0 * scanCount).Within(0.1));

            File.Delete(msAlignPath);
        }

        /// <summary>
        /// Mass and charge filters drop peaks before indexing, and an isolation window filters on the m/z the
        /// mass was observed at.
        /// </summary>
        [Test]
        public static void ExternalDeconFiltersApply()
        {
            string msAlignPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ExternalDeconFilter_ms1.msalign");
            using (var writer = new StreamWriter(msAlignPath))
            {
                for (int scanNumber = 1; scanNumber <= 4; scanNumber++)
                {
                    writer.WriteLine("BEGIN IONS");
                    writer.WriteLine("SCANS=" + scanNumber);
                    writer.WriteLine("RETENTION_TIME=" + (60.0 * scanNumber));
                    writer.WriteLine("LEVEL=1");
                    writer.WriteLine("1500.00000\t1000.00\t2\t0.9");   // below a 3000 Da minMass
                    writer.WriteLine("9000.00000\t2000.00\t10\t0.9");  // kept; m/z ~ 901
                    writer.WriteLine("END IONS");
                    writer.WriteLine();
                }
            }

            var scans = Enumerable.Range(1, 4).Select(i => new MsDataScan(
                new MzSpectrum(new double[] { 700 }, new double[] { 1 }, false), i, 1, true, Polarity.Positive,
                i * 1.0, new MzRange(200, 2000), "f", MZAnalyzerType.Orbitrap, 1, 1.0, null, "scan=" + i)).ToArray();

            var results = ExternalDeconResults.ReadMsAlign(msAlignPath);

            var filtered = new ExternalDeconXicConstructor(new PpmTolerance(20), 1, 100, 3, results, minMass: 3000);
            var xics = filtered.GetAllXics(scans);
            Assert.That(filtered.LastIndexedPeakCount, Is.EqualTo(4), "only the 9000 Da mass should survive minMass");
            Assert.That(xics.Count, Is.EqualTo(1));

            // isolation window around m/z 901 keeps the 9 kDa 10+ species; a window elsewhere keeps nothing
            var unfiltered = new ExternalDeconXicConstructor(new PpmTolerance(20), 1, 100, 3, results);
            unfiltered.GetAllXics(scans, new MzRange(890, 910));
            Assert.That(unfiltered.LastIndexedPeakCount, Is.EqualTo(4));

            var outsideWindow = new ExternalDeconXicConstructor(new PpmTolerance(20), 1, 100, 3, results);
            Assert.Throws<EngineLayer.MetaMorpheusException>(() => outsideWindow.GetAllXics(scans, new MzRange(1200, 1300)));

            File.Delete(msAlignPath);
        }

        /// <summary>
        /// A scan array whose numbering does not line up with the msAlign must fail loudly under scan-number
        /// matching, and succeed under retention-time matching.
        /// </summary>
        [Test]
        public static void ExternalDeconRetentionTimeMatching()
        {
            string msAlignPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ExternalDeconRt_ms1.msalign");
            using (var writer = new StreamWriter(msAlignPath))
            {
                for (int scanNumber = 1; scanNumber <= 5; scanNumber++)
                {
                    writer.WriteLine("BEGIN IONS");
                    writer.WriteLine("SCANS=" + scanNumber);
                    writer.WriteLine("RETENTION_TIME=" + (300.0 + 6.0 * scanNumber).ToString("F2")); // seconds
                    writer.WriteLine("LEVEL=1");
                    writer.WriteLine("7000.00000\t1000.00\t8\t0.9");
                    writer.WriteLine("END IONS");
                    writer.WriteLine();
                }
            }

            // same retention times, but scan numbers offset by 1000 (as if the deconvoluted file was renumbered)
            var scans = Enumerable.Range(1, 5).Select(i => new MsDataScan(
                new MzSpectrum(new double[] { 700 }, new double[] { 1 }, false), 1000 + i, 1, true, Polarity.Positive,
                (300.0 + 6.0 * i) / 60.0, new MzRange(200, 2000), "f", MZAnalyzerType.Orbitrap, 1, 1.0, null,
                "scan=" + (1000 + i))).ToArray();

            var results = ExternalDeconResults.ReadMsAlign(msAlignPath);

            var byScanNumber = new ExternalDeconXicConstructor(new PpmTolerance(20), 1, 100, 3, results);
            Assert.Throws<EngineLayer.MetaMorpheusException>(() => byScanNumber.GetAllXics(scans));

            var byRetentionTime = new ExternalDeconXicConstructor(new PpmTolerance(20), 1, 100, 3, results,
                scanMatchMode: ExternalDeconScanMatchMode.RetentionTime, scanMatchRetentionTimeTolerance: 0.01);
            var xics = byRetentionTime.GetAllXics(scans);
            Assert.That(byRetentionTime.LastMatchedScanCount, Is.EqualTo(5));
            Assert.That(xics.Count, Is.EqualTo(1));
            Assert.That(xics[0].Peaks.Count, Is.EqualTo(5));

            File.Delete(msAlignPath);
        }

        /// <summary>
        /// A single combined ms2.msalign holding every ISD voltage is split back into per-voltage results by scan
        /// number, both from an explicit scan-to-voltage map and from the repeating acquisition cycle. Also covers
        /// the TopFD 1.8 retention time unit: 1.8 writes MINUTES where earlier versions wrote seconds.
        /// </summary>
        [Test]
        public static void CombinedMsAlignSplitsByVoltage()
        {
            // 3 cycles of MS1, isd60, isd80, isd100 -> scans 1..12; the MS1 scans are not in the ms2 file
            var voltageByScanNumber = new Dictionary<int, int>();
            string msAlignPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "CombinedIsd_ms2.msalign");
            using (var writer = new StreamWriter(msAlignPath))
            {
                writer.WriteLine("#TopFD 1.8.1");
                writer.WriteLine("####################### Parameters ######################");
                writer.WriteLine("#Spectral data type:         \tCentroid");
                writer.WriteLine("####################### Parameters ######################");
                writer.WriteLine();
                for (int scanNumber = 1; scanNumber <= 12; scanNumber++)
                {
                    int cyclePosition = (scanNumber - 1) % 4;
                    if (cyclePosition == 0) { voltageByScanNumber[scanNumber] = 15; continue; } // MS1, not in ms2 file
                    int voltage = cyclePosition == 1 ? 60 : cyclePosition == 2 ? 80 : 100;
                    voltageByScanNumber[scanNumber] = voltage;

                    writer.WriteLine("BEGIN IONS");
                    writer.WriteLine("TITLE=Scan_" + scanNumber);
                    writer.WriteLine("SCANS=" + scanNumber);
                    // TopFD 1.8 writes RETENTION_TIME in MINUTES
                    writer.WriteLine("RETENTION_TIME=" + (0.01 * scanNumber).ToString("F3"));
                    writer.WriteLine("LEVEL=2");
                    // one mass per voltage, so which voltage a split bucket holds is unambiguous
                    writer.WriteLine((1000.0 + voltage).ToString("F5") + "\t1000.00\t5\t0.9");
                    writer.WriteLine("END IONS");
                    writer.WriteLine();
                }
            }

            var combined = ExternalDeconResults.ReadMsAlign(msAlignPath);
            Assert.That(combined.ScanCount, Is.EqualTo(9), "3 cycles x 3 fragment channels");
            // the header says TopFD 1.8, so the minutes in the file must survive as minutes
            Assert.That(combined.RetentionTimeByOneBasedScanNumber[2], Is.EqualTo(0.02).Within(1e-6));

            foreach (var split in new[] { combined.SplitByVoltage(voltageByScanNumber),
                                          combined.SplitByIsdCycle(new[] { 15, 60, 80, 100 }) })
            {
                Assert.That(split.Keys.OrderBy(v => v), Is.EqualTo(new[] { 60, 80, 100 }));
                foreach (var voltage in new[] { 60, 80, 100 })
                {
                    Assert.That(split[voltage].ScanCount, Is.EqualTo(3));
                    // every mass in a bucket must be the one written for that voltage
                    Assert.That(split[voltage].PeaksByOneBasedScanNumber.SelectMany(kvp => kvp.Value)
                        .Select(p => p.MonoisotopicMass), Has.All.EqualTo(1000.0 + voltage).Within(1e-4));
                }
            }

            // explicitly forcing the seconds interpretation reproduces mzLib's divide-by-60
            var asSeconds = ExternalDeconResults.ReadMsAlign(new[] { msAlignPath }, MsAlignRetentionTimeUnit.Seconds);
            Assert.That(asSeconds.RetentionTimeByOneBasedScanNumber[2], Is.EqualTo(0.02 / 60.0).Within(1e-9));

            File.Delete(msAlignPath);
        }

        
        /// <summary>
        /// End-to-end on real data (not CI) for the COMBINED msAlign layout: TopFD 1.8 was run on the whole ISD
        /// file rather than on per-channel splits, so all three source voltages sit in one ms2.msalign and are
        /// told apart only by scan number. MS1 XICs come from the ms1.msalign grouped onto the ms1.feature file.
        /// </summary>
        [Test, Explicit]
        public static void CombinedMsAlignIsdXicsFromRealTopFd()
        {
            string dataDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string stem = @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected";

            string mzmlFile = Path.Combine(dataDir, stem + ".mzML");
            string ms1MsAlign = Path.Combine(dataDir, stem + "_ms1.msalign");
            string ms1Feature = Path.Combine(dataDir, stem + "_ms1.feature");
            string combinedMs2MsAlign = Path.Combine(dataDir, stem + "_ms2.msalign");

            var dataFile = MsDataFileReader.GetDataFile(mzmlFile);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // ---- the combined ms2.msalign, split back into channels by scan number ----
            var combined = ExternalDeconResults.ReadMsAlign(combinedMs2MsAlign);
            var voltageByScanNumber = ExternalDeconResults.GetIsdVoltageByScanNumber(allScans);
            var byVoltage = combined.SplitByVoltage(voltageByScanNumber);

            // the scan-filter split and the pure cycle split (MS1, 60, 80, 100 repeating) must agree
            var byCycle = combined.SplitByIsdCycle(new[] { 15, 60, 80, 100 });

            // retention times must line up with the raw file's, which they only do if the TopFD 1.8
            // minutes-vs-seconds convention was handled
            int probeScan = isdMap.First().Value.First().OneBasedScanNumber;
            double rawRt = allScans.First(s => s.OneBasedScanNumber == probeScan).RetentionTime;
            double msAlignRt = combined.RetentionTimeByOneBasedScanNumber[probeScan];

            // ---- MS1: trace from ms1.msalign, group onto ms1.feature ----
            var ms1Constructor = new FeatureGroupedXicConstructor(new PpmToleranceWithNotch(20, 2, 2), 2, 1, 3,
                ms1MsAlign, ms1Feature, featureMassTolerance: new PpmTolerance(20), minMass: 3000, maxMass: 50000);
            var featureXics = ms1Constructor.GetAllXics(ms1Scans);

            // ---- MS2: one constructor, fed the split combined file ----
            var ms2Constructor = ExternalDeconXicConstructor.BySourceVoltage(new PpmToleranceWithNotch(20, 2, 2), 2, 1, 3,
                byVoltage, combined, maxMass: 50000);
            var fragmentXics = new List<ExtractedIonChromatogram>();
            foreach (var channel in isdMap)
            {
                var channelXics = ms2Constructor.GetAllXics(channel.Value.ToArray());
                fragmentXics.AddRange(channelXics);
            }

            // ---- precursor-fragment grouping, then pseudo-MS2 scans, then the msAlign ----
            // Same three steps ISDEngine.GetPseudoMs2Scans performs internally, done inline here so the grouping
            // thresholds can be changed and an msAlign inspected without running the whole search pipeline.
            var grouper = new XicGroupingEngine(apexRTTolerance: 0.2f, overlapThreshold: 0, correlationThreshold: 0.5,
                maxThreadsForGrouping: Environment.ProcessorCount - 2);
            var pfGroups = grouper.PrecursorFragmentGrouping(featureXics, fragmentXics).ToList();

            var commonParams = new CommonParameters();
            var pseudoScans = new List<Ms2ScanWithSpecificMass>(pfGroups.Count);
            for (int i = 0; i < pfGroups.Count; i++)
            {
                pfGroups[i].PFgroupIndex = i + 1;
                pseudoScans.Add(pfGroups[i].GetPseudoMs2ScanFromPfGroup(
                    PseudoMs2ConstructionType.Mass, commonParams, dataFile.FilePath));
            }

            var fragmentCounts = pseudoScans.Select(s => s.ExperimentalFragments.Length).OrderBy(n => n).ToArray();

            string outDir = Path.Combine(dataDir, "PseudoMs2");
            Directory.CreateDirectory(outDir);
            string msAlignOut = Path.Combine(outDir, stem + "_pseudo_ms2.msalign");
            ISDEngine.WriteMsAlignFile(msAlignOut, pseudoScans);
        }
    }
}
