//using EngineLayer;
//using MassSpectrometry;
//using NUnit.Framework;
//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;
//using Chemistry;
//using MzLibUtil;
//using UsefulProteomicsDatabases;
//using TaskLayer;


//namespace Test.TestDIA
//{
//    public class PeakExtractionTests
//    {
//        // Test cases to be used
//        public static IEnumerable<PeakExtractionTestCase> GetTestCases()
//        {
//            Loaders.LoadElements();

//            // Test Case 1: Fake data with 3 scans, very simple. 
//            var scan1 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 10, 20, 30, 40 }, false),
//                1, 1, true, Polarity.Positive, 0.1, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
//            var scan2 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 100, 200, 300, 400 }, false),
//                3, 1, true, Polarity.Positive, 0.2, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
//            var scan3 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 1000, 2000, 3000, 4000 }, false),
//                5, 1, true, Polarity.Positive, 0.3, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
//            var scans = new MsDataScan[] { scan1, scan2, scan3 };

//            var testCase1 = new PeakExtractionTestCase
//            {
//                TestCaseName = "Fake Data, Simple, 3 scans",
//                FakeScans = scans,
//                NumScansPerCycle = 2,
//                ExpectedPeaksCount = 12,
//                ExpectedIntensity = new List<double>[] { new List<double> { 10, 20, 30, 40 }, new List<double> { 100, 200, 300, 400 },
//                    new List<double> { 1000, 2000, 3000, 4000 } },
//                ExpectedMz = new List<double>[] { new List<double> { 1, 2, 3, 4 }, new List<double> { 1, 2, 3, 4 }, new List<double> { 1, 2, 3, 4 } },
//                ExpectedRt = new List<double>[] { new List<double> { 0.1, 0.1, 0.1, 0.1 }, new List<double> { 0.2, 0.2, 0.2, 0.2 },
//                    new List<double> { 0.3, 0.3, 0.3, 0.3 } },
//                ExpectedZeroBasedIndex = new List<int>[] { new List<int> { 0, 0, 0, 0 }, new List<int> { 1, 1, 1, 1 }, new List<int> { 2, 2, 2, 2 } },
//                PeakTableListCount = 4,
//                PeakTableValidArrays = new List<int> { 100, 200, 300, 400 }
//            };
//            yield return testCase1;

//            //Test Case 2: snip bu-ISD data
//            var testCase2 = new PeakExtractionTestCase();
//            var path = @"E:\ISD Project\ISD_bu\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
//            var myFileManagers = new MyFileManager(true);
//            var dataFile = myFileManagers.LoadFile(path, new CommonParameters());
//            var realScans = dataFile.GetAllScansList().Where(s => s.RetentionTime <= 1).ToArray();
//            testCase2.FakeScans = realScans;
//            testCase2.TestCaseName = "snip bu-ISD";
//            testCase2.NumScansPerCycle = 2;
//            testCase2.ExpectedPeaksCount = realScans.Sum(s => s.MassSpectrum.Size);
//            testCase2.ExpectedMz = new List<double>[realScans.Length];
//            testCase2.ExpectedIntensity = new List<double>[realScans.Length];
//            testCase2.ExpectedRt = new List<double>[realScans.Length];
//            testCase2.ExpectedZeroBasedIndex = new List<int>[realScans.Length];
//            for (int i = 0; i < realScans.Length; i++)
//            {
//                var zeroIndex = (realScans[i].OneBasedScanNumber - 1) / testCase2.NumScansPerCycle;
//                testCase2.ExpectedMz[i] = realScans[i].MassSpectrum.XArray.ToList();
//                testCase2.ExpectedIntensity[i] = realScans[i].MassSpectrum.YArray.ToList();
//                testCase2.ExpectedRt[i] = Enumerable.Repeat(realScans[i].RetentionTime, realScans[i].MassSpectrum.Size).ToList();
//                testCase2.ExpectedZeroBasedIndex[i] = Enumerable.Repeat(zeroIndex, realScans[i].MassSpectrum.Size).ToList();
//            }
//            var allMzs = realScans.SelectMany(s => s.MassSpectrum.XArray).ToList();
//            var entries = allMzs.Select(mz => (int)Math.Round(mz * 100, 0)).Distinct();
//            testCase2.PeakTableListCount = entries.Count();
//            testCase2.PeakTableValidArrays = entries.OrderBy(p => p).ToList();
//            yield return testCase2;
//        }

//        // Structure of test case
//        public class PeakExtractionTestCase
//        {
//            public string TestCaseName { get; set; }
//            public MsDataScan[] FakeScans { get; set; }
//            public int NumScansPerCycle { get; set; }
//            public int ExpectedPeaksCount { get; set; }

//            public List<double>[] ExpectedIntensity { get; set; }
//            public List<double>[] ExpectedMz { get; set; }
//            public List<double>[] ExpectedRt { get; set; }
//            public List<int>[] ExpectedZeroBasedIndex { get; set; }
//            public int PeakTableListCount { get; set; }
//            public List<int> PeakTableValidArrays { get; set; }
//            public override string ToString() => TestCaseName;
//        }


//        [Test, TestCaseSource(nameof(GetTestCases))]
//        public static void TestGetAllPeaks(PeakExtractionTestCase testCase)
//        {
//            var allPeaksByScan = Peak.GetAllPeaksByScan(testCase.FakeScans, testCase.NumScansPerCycle);

//            var allPeakCount = allPeaksByScan.Where(v => v != null).Sum(p => p.Count);
//            Assert.That(allPeakCount, Is.EqualTo(testCase.ExpectedPeaksCount));

//            var validPeaksByScan = allPeaksByScan.Where(v => v != null).ToList();
//            for (int scanIndex = 0; scanIndex < validPeaksByScan.Count; scanIndex++)
//            {
//                List<Peak> scanSpecificPeaks = validPeaksByScan[scanIndex];

//                for (int peakIndex = 0; peakIndex < scanSpecificPeaks.Count; peakIndex++)
//                {
//                    Peak peak = scanSpecificPeaks[peakIndex];
//                    Assert.That(peak.Intensity, Is.EqualTo(testCase.ExpectedIntensity[scanIndex][peakIndex]));
//                    Assert.That(peak.Mz, Is.EqualTo(testCase.ExpectedMz[scanIndex][peakIndex]));
//                    Assert.That(peak.RetentionTime, Is.EqualTo(testCase.ExpectedRt[scanIndex][peakIndex]));
//                    Assert.That(peak.ZeroBasedScanIndex, Is.EqualTo(testCase.ExpectedZeroBasedIndex[scanIndex][peakIndex]));
//                }
//            }
//        }


//        [Test]
//        [TestCaseSource(nameof(GetTestCases))]
//        public static void TestPeakTable(PeakExtractionTestCase testCase)
//        {
//            //Building ground truth; need to past the GetAllPeaksByScan test first
//            var allPeaksByScan = Peak.GetAllPeaksByScan(testCase.FakeScans);
//            var allPeaks = allPeaksByScan.Where(v => v != null).SelectMany(v => v).ToList();
//            var mzGroups = allPeaks.GroupBy(p => Math.Round(p.Mz, 2)).OrderBy(g => g.Key).ToList();

//            //test peakTable
//            var peakTable = Peak.GetPeakTable(allPeaks, 100);
//            Assert.That(peakTable.Where(v => v != null).Count(), Is.EqualTo(testCase.PeakTableListCount));

//            var tableEntries = peakTable.Where(v => v != null).ToList();
//            for (int i = 0; i < tableEntries.Count; i++)
//            {
//                for (int j = 0; j < tableEntries[i].Count; j++)
//                {
//                    Assert.That(tableEntries[i][j], Is.SameAs(mzGroups[i].ElementAt(j)));
//                }
//            }
//            for (int i = 0; i < peakTable.Length; i++)
//            {
//                if (testCase.PeakTableValidArrays.Contains(i))
//                {
//                    Assert.That(peakTable[i], Is.Not.Null);
//                }
//                else
//                {
//                    Assert.That(peakTable[i], Is.Null);
//                }
//            }
//        }
//    }
//}

