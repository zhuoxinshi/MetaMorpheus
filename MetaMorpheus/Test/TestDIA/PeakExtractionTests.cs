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
            public int ExpectedPeaksCount { get; set; }

            public List<double>[] ExpectedIntensity { get; set; }
            public List<double>[] ExpectedMz { get; set; }
            public List<double>[] ExpectedRt { get; set; }
        }


        public static void OneTimeSetUp()
        {
            List<PeakExtractionTestCase> testCases = new List<PeakExtractionTestCase>();
            MsDataScan[] fakeScans = new MsDataScan[3];

            // TODO: Make fake scans
            var spectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 10, 20, 30 }, false);
            var scan1 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 10, 20, 30 }, false), 
                1, 1, true, Polarity.Positive, 0.1, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan2 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 100, 200, 300 }, false),
                3, 1, true, Polarity.Positive, 0.2, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan3 = new MsDataScan(new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 1000, 2000, 3000 }, false),
                5, 1, true, Polarity.Positive, 0.3, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scans = new MsDataScan[] { scan1, scan2, scan3 };

            // TODO: Turn fake scans into test cases
            var testCase1 = new PeakExtractionTestCase
            {
                FakeScans = scans,
                ExpectedPeaksCount = 9,
                ExpectedIntensity = new List<double>[] { new List<double> { 10, 20, 30 }, new List<double> { 100, 200, 300 }, new List<double> { 1000, 2000, 3000 } },
                ExpectedMz = new List<double>[] { new List<double> { 1, 2, 3 }, new List<double> { 1, 2, 3 }, new List<double> { 1, 2, 3 } },
                ExpectedRt = new List<double>[] { new List<double> { 0.1, 0.1, 0.1 }, new List<double> { 0.2, 0.2, 0.2 }, new List<double> { 0.3, 0.3, 0.3 } }
            };

            // TODO: Add more test cases from either fake or real data. 
            testCases.Add(testCase1);
            TestCases = testCases;
        }


        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public static void TestGetAllPeaks(PeakExtractionTestCase testCase)
        {
            var allPeaks = Peak.GetAllPeaksByScan(testCase.FakeScans);

            var allPeakCount = allPeaks.Sum(p => p.Count);
            Assert.That(allPeakCount, Is.EqualTo(testCase.ExpectedPeaksCount));

            for (int scanIndex = 0; scanIndex < allPeaks.Length; scanIndex++)
            {
                List<Peak> scanSpecificPeaks = allPeaks[scanIndex];
                List<double> intensities = testCase.ExpectedIntensity[scanIndex];
                List<double> mzs = testCase.ExpectedMz[scanIndex];
                List<double> rts = testCase.ExpectedRt[scanIndex];

                for (int singleScanIndex = 0; singleScanIndex < scanSpecificPeaks.Count; singleScanIndex++)
                {
                    Peak peak = scanSpecificPeaks[singleScanIndex];
                    Assert.That(peak.Intensity, Is.EqualTo(intensities[singleScanIndex]));
                    Assert.That(peak.Mz, Is.EqualTo(mzs[singleScanIndex]));
                    Assert.That(peak.RetentionTime, Is.EqualTo(rts[singleScanIndex]));
                }
            }
        }
    }
}
