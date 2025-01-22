using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MzLibUtil;
using MassSpectrometry;


namespace Test.TestDIA
{
    public class PeakExtractionTests
    {
        public static IEnumerable<PeakExtractionTestCase> TestCases { get; set; }

        public class PeakExtractionTestCase
        {
            public MsDataScan[] FakeScans { get; set; }
            public int NumScansPerCycle { get; set; }
            public int ExpectedPeaksCount { get; set; }

            public List<double>[] ExpectedIntensity { get; set; }
            public List<double>[] ExpectedMz { get; set; }
            public List<double>[] ExpectedRt { get; set; }
            public List<int>[] ExpectedZeroBasedIndex { get; set; }
            public int PeakTableListCount { get; set; }
            public List<double> PeakTableValidArrays { get; set; }
        }

        [OneTimeSetUp]
        public static void OneTimeSetUp()
        {
            List<PeakExtractionTestCase> testCases = new List<PeakExtractionTestCase>();
            MsDataScan[] fakeScans = new MsDataScan[3];

            // TODO: Make fake scans
            var scan1 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 10, 20, 30, 40 }, false),
                1, 1, true, Polarity.Positive, 0.1, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan2 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 100, 200, 300, 400 }, false),
                3, 1, true, Polarity.Positive, 0.2, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan3 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 1000, 2000, 3000, 4000 }, false),
                5, 1, true, Polarity.Positive, 0.3, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scans = new MsDataScan[] { scan1, scan2, scan3 };

            // TODO: Turn fake scans into test cases
            var testCase1 = new PeakExtractionTestCase
            {
                FakeScans = scans,
                NumScansPerCycle = 2,
                ExpectedPeaksCount = 12,
                ExpectedIntensity = new List<double>[] { new List<double> { 10, 20, 30, 40 }, new List<double> { 100, 200, 300, 400 }, 
                    new List<double> { 1000, 2000, 3000, 4000 } },
                ExpectedMz = new List<double>[] { new List<double> { 1, 2, 3, 4 }, new List<double> { 1, 2, 3, 4 }, new List<double> { 1, 2, 3, 4} },
                ExpectedRt = new List<double>[] { new List<double> { 0.1, 0.1, 0.1, 0.1 }, new List<double> { 0.2, 0.2, 0.2, 0.2 }, 
                    new List<double> { 0.3, 0.3, 0.3, 0.3 } },
                ExpectedZeroBasedIndex = new List<int>[] { new List<int> { 0, 0, 0, 0 }, new List<int> { 1, 1, 1, 1 }, new List<int> { 2, 2, 2, 2 } },
                PeakTableListCount = 4,
                PeakTableValidArrays = new List<double> { 100, 200, 300, 400 }
            };

            // TODO: Add more test cases from either fake or real data. 
            testCases.Add(testCase1);
            TestCases = testCases;
        }


        [Test]
        //[TestCaseSource(nameof(TestCases))]
        public static void TestGetAllPeaks(PeakExtractionTestCase testCase)
        {
            var allPeaksByScan = Peak.GetAllPeaksByScan(testCase.FakeScans, testCase.NumScansPerCycle);

            var allPeakCount = allPeaksByScan.Where(v => v != null).Sum(p => p.Count);
            Assert.That(allPeakCount, Is.EqualTo(testCase.ExpectedPeaksCount));

            var validPeaksByScan = allPeaksByScan.Where(v => v != null).ToList();
            for (int scanIndex = 0; scanIndex < validPeaksByScan.Count; scanIndex++)
            {
                List<Peak> scanSpecificPeaks = validPeaksByScan[scanIndex];

                for (int singleScanIndex = 0; singleScanIndex < scanSpecificPeaks.Count; singleScanIndex++)
                {
                    Peak peak = scanSpecificPeaks[singleScanIndex];
                    Assert.That(peak.Intensity, Is.EqualTo(testCase.ExpectedIntensity[scanIndex][singleScanIndex]));
                    Assert.That(peak.Mz, Is.EqualTo(testCase.ExpectedMz[scanIndex][singleScanIndex]));
                    Assert.That(peak.RetentionTime, Is.EqualTo(testCase.ExpectedRt[scanIndex][singleScanIndex]));
                    Assert.That(peak.ZeroBasedScanIndex, Is.EqualTo(testCase.ExpectedZeroBasedIndex[scanIndex][singleScanIndex]));
                }
            }
        }

        [Test]
        public static void TestCasesAll()
        {
            foreach(var testCase in TestCases)
            {
                TestGetAllPeaks(testCase);
                TestPeakTable(testCase);
            }
        }

        [Test]
        public static void TestPeakTable(PeakExtractionTestCase testCase)
        {
            var allPeaksByScan = Peak.GetAllPeaksByScan(testCase.FakeScans);
            var allPeaks = allPeaksByScan.Where(v => v != null).SelectMany(v => v).ToList();
            var allPeaksEntries = allPeaksByScan.Where(v => v != null).ToArray();
            var peakTable = Peak.GetPeakTable(allPeaks, 100);
            Assert.That(peakTable.Where(v => v != null).Count(), Is.EqualTo(testCase.PeakTableListCount));

            var tableEntries = peakTable.Where(v => v != null).ToList();
            for (int i = 0; i < tableEntries.Count; i++)
            {
                for (int j = 0; j < tableEntries[i].Count; j++)
                {
                    Assert.That(tableEntries[i][j], Is.SameAs(allPeaksEntries[j][i]));
                } 
            }
            for (int i = 0; i < peakTable.Length; i++)
            {
                if (testCase.PeakTableValidArrays.Contains(i))
                {
                    Assert.That(peakTable[i], Is.Not.Null);
                }
                else
                {
                    Assert.That(peakTable[i], Is.Null);
                }
            }
        }
    }
}

