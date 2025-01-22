using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Test.TestDIA.PeakExtractionTests;
using EngineLayer;
using Microsoft.VisualStudio.TestPlatform.ObjectModel;
using TaskLayer;
using Plotly.NET.CSharp;

namespace Test.TestDIA
{
    public class PeakCurveTests
    {
        public static IEnumerable<PeakCurveTestCase> TestCases { get; set; }

        public class PeakCurveTestCase
        {
            public MsDataScan[] FakeScans { get; set; }
            public int NumScansPerCycle { get; set; }
            public int ExpectedPeakCurvesCount { get; set; }
            public List<double> ExpectedAverageIntensity { get; set; }
            public List<double> ExpectedAverageMz { get; set; }
            public List<double> ExpectedApexRt { get; set; }
            public List<double> ExpectedApexCycle { get; set; }
        }

        [OneTimeSetUp]
        public static void OneTimeSetUp()
        {
            List<PeakCurveTestCase> testCases = new List<PeakCurveTestCase>();
            MsDataScan[] fakeScans = new MsDataScan[3];

            // TODO: Make fake scans
            var scan1 = new MsDataScan(new MzSpectrum(new double[] { 400, 500, 1000 }, new double[] { 10, 20, 30 }, false),
                1, 1, true, Polarity.Positive, 0.1, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan2 = new MsDataScan(new MzSpectrum(new double[] { 400.001, 500.001, 1000.001 }, new double[] { 1001, 2001, 3000 }, false),
                3, 1, true, Polarity.Positive, 0.2, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan3 = new MsDataScan(new MzSpectrum(new double[] { 400.002, 500.002, 1000.002 }, new double[] { 100, 200, 300 }, false),
                5, 1, true, Polarity.Positive, 0.3, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan4 = new MsDataScan(new MzSpectrum(new double[] { 400.005, 500.0025, 1000.003 }, new double[] { 1000, 2000, 3001 }, false),
                7, 1, true, Polarity.Positive, 0.4, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scan5 = new MsDataScan(new MzSpectrum(new double[] { 400.0065, 500.0015, 1000.005 }, new double[] { 100, 200, 300 }, false),
                9, 1, true, Polarity.Positive, 0.5, null, "", MZAnalyzerType.Orbitrap, 1, null, null, null);
            var scans = new MsDataScan[] { scan1, scan2, scan3, scan4, scan5};

            var averageMz3 = (400 * 10 + 400.001 * 1001 + 400.002 * 100) / (10 + 1001 + 100);
            var averageMz4 = (400.005 * 1000 + 400.0065 * 100) / (1000 + 100);
            var averageMz2 = (500 * 20 + 500.001 * 2001 + 500.002 * 200 + 500.0025 * 2000 + 500.0015 * 200) / (20 + 2001 + 200 + 2000 + 200);
            var averageMz1 = (1000 * 30 + 1000.001 * 3000 + 1000.002 * 300 + 1000.003 * 3001 + 1000.005 * 300) / (30 + 3000 + 300 + 3001 + 300);


            var testCase1 = new PeakCurveTestCase
            {
                FakeScans = scans,
                NumScansPerCycle = 2,
                ExpectedPeakCurvesCount = 4,
                ExpectedAverageMz = new List<double> {averageMz1, averageMz2, averageMz3, averageMz4},
                ExpectedApexRt = new List<double> { 0.4, 0.2, 0.2, 0.4 },
                ExpectedApexCycle = new List<double> { 3, 1, 1, 3 },
            };

            testCases.Add(testCase1);
            TestCases = testCases;
        }

        [Test]
        //[TestCaseSource(nameof(TestCases))]
        public static void TestPeakFindingMz(PeakCurveTestCase testCase)
        {
            var diaParam = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20));
            var allPeakCurves = ISDEngine_static.GetAllPeakCurves_Peak(testCase.FakeScans, diaParam, new PpmTolerance(5), 2, 
                out List<Peak>[] allPeaksByScan, false);
            Assert.That(allPeakCurves.Count == testCase.ExpectedPeakCurvesCount);

            for (int i = 0; i < allPeakCurves.Count; i++)
            {
                Assert.That(allPeakCurves[i].AveragedMz, Is.EqualTo(testCase.ExpectedAverageMz[i]).Within(0.0000001));
                Assert.That(allPeakCurves[i].ApexRT, Is.EqualTo(testCase.ExpectedApexRt[i]));
                Assert.That(allPeakCurves[i].ApexScanCycle, Is.EqualTo(testCase.ExpectedApexCycle[i]));
            }
        }

        [Test]
        public static void TestPeakFinding()
        {
            foreach (var testCase1 in TestCases)
            {
                TestPeakFindingMz(testCase1);
            }
        }

        [Test]
        public static void SplineTest()
        {
            var path = @"E:\ISD Project\ISD_240606\06-09-24_mix_sample10_5uL_ISD.mzML";
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path,  new CommonParameters());
            var allMs1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToArray();
            var allMs2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var ms1Scans = allMs1Scans.Where(s => s.RetentionTime >=31.96 && s.RetentionTime <= 35).ToArray();
            var ms2Scans = allMs2Scans.Where(s => s.RetentionTime >= 32 && s.RetentionTime <= 35).ToArray();

            var ms1PeakCurves = ISDEngine_static.GetAllPeakCurves_Peak(ms1Scans, new DIAparameters(new PpmTolerance(5), new PpmTolerance(20)), new PpmTolerance(5), 2,
                out List<Peak>[] allPeaksByScan, false);
            var pc1 = ms1PeakCurves.Where(p => Math.Abs(p.AveragedMz - 857.47) < 0.01).First();
            var pc2 = ms1PeakCurves.Where(p => Math.Abs(p.AveragedMz - 845.86) < 0.01).First();

            var ms2PeakCurves = ISDEngine_static.GetAllPeakCurves_Peak(ms2Scans, new DIAparameters(new PpmTolerance(5), new PpmTolerance(20)), new PpmTolerance(5), 2,
                               out List<Peak>[] allPeaksByScan2, false);
            var frag1 = ms2PeakCurves.Where(p => Math.Abs(p.AveragedMz - 1428.27) < 0.01).First();
            var frag2 = ms2PeakCurves.Where(p => Math.Abs(p.AveragedMz - 1427.94) < 0.01).First();
            frag1.GetScanCycleCubicSplineSavgolSmoothedXYData(7, 0.05);
            //frag1.VisualizeGeneral("cycle").Show();

            var rtIndexMap = ISDEngine_static.GetRtIndexMap(ms1Scans);
            var rtMap = ISDEngine_static.GetRtMap(ms1Scans, ms2Scans);

            //pc1.GetScanCycleCubicSplineSavgolSmoothedXYData(7, 0.05);
            //pc2.GetScanCycleCubicSplineSavgolSmoothedXYData(7, 0.05);
            //var pfCorr1 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            //pc1.GetRawXYData();
            //frag1.GetRawXYData();
            //var rawCorr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            //pc1.GetCubicSplineXYData(0.05);
            //frag1.GetMs1SpaceCubicSplineXYData(rtMap, 0.05);
            //var ms1SpaceCorr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            //pc1.GetSavgolSmoothedXYData(rtIndexMap, 7);
            //frag1.GetSavgolSmoothedXYData(rtIndexMap, 7);
            //var sgCorr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            pc1.GetCubicSplineSavgolSmoothedXYData(7, 0.05);
            frag1.GetMs1SpaceCubicSplineSavgolSmoothedXYData(rtMap, 7, 0.05);
            var preX = pc1.XYData.Select(x => x.Item1).ToArray();
            var fragX = frag1.XYData.Select(x => x.Item1).ToArray();
            var corrMs1cubic = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            pc1.GetRawXYData();
            frag1.GetRawXYData();
            var rawPfCorr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            pc1.GetCubicSplineXYData(0.05);
            frag1.GetMs1SpaceCubicSplineXYData(rtMap, 0.05);
            var ms1SpaceCorr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            pc1.GetCubicSplineXYData(0.05);
            frag1.GetCubicSplineXYData(0.05);
            var cubicCorr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);

            pc1.GetScanCycleCubicSplineSavgolSmoothedXYData(7, 0.05);
            frag1.GetScanCycleCubicSplineSavgolSmoothedXYData(7, 0.05);
            var corr1 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, frag1);
            pc1.GetRawXYData();
            pc2.GetRawXYData();
            var corr2 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, pc2);
            pc1.GetSavgolSmoothedXYData(rtIndexMap, 11);
            pc2.GetSavgolSmoothedXYData(rtIndexMap, 11);
            var corr3 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, pc2);
        }

        //[Test]
        //public static Dictionary<SplineType, double> TestAllSplineTypesCorrelation(PeakCurve pc1, PeakCurve pc2, List<SplineType> splineTypes)
        //{
        //    var correlationResults = new Dictionary<SplineType, double>();

        //    foreach (var splineType in splineTypes)
        //    {
        //        switch (splineType)
        //        {
        //            case SplineType.NoSpline:
        //                    pc.GetRawXYData();
        //                break;
        //            case SplineType.CubicSpline:
        //                    pc.GetCubicSplineXYData(diaParam.SplineRtInterval);
        //                break;
        //            case SplineType.ScanCycleCubicSpline:
        //                    pc.GetScanCycleCubicSplineXYData(diaParam.ScanCycleSplineTimeInterval);
        //                break;
        //            case SplineType.CubicSplineSavgolSmoothed:
        //                foreach (var pc in allPeakCurves)
        //                    pc.GetCubicSplineSavgolSmoothedXYData(diaParam.SGfilterWindowSize, diaParam.ScanCycleSplineTimeInterval);
        //                break;
        //            case SplineType.ScanCycleCubicSplineSavgolSmoothed:
        //                foreach (var pc in allPeakCurves)
        //                    pc.GetScanCycleCubicSplineSavgolSmoothedXYData(diaParam.SGfilterWindowSize, diaParam.ScanCycleSplineTimeInterval);
        //                break;
        //            case SplineType.SavgolSmoothedCubicSpline:
        //                foreach (var pc in allPeakCurves)
        //                    pc.GetSavgolSmoothedCubicSplineXYData(rtIndexMap, diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
        //                break;
        //            case SplineType.Ms1SpaceCubicSpline:
        //                foreach (var pc in allPeakCurves)
        //                    pc.GetMs1SpaceCubicSplineXYData(rtMap, diaParam.SplineRtInterval);
        //                break;
        //            case SplineType.Ms1SpaceCubicSplineSavgolSmoothed:
        //                foreach (var pc in allPeakCurves)
        //                    pc.GetMs1SpaceCubicSplineSavgolSmoothedXYData(rtMap, diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
        //                break;
        //        }

        //        var correlation = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, pc2);
        //        correlationResults[splineType] = correlation;
        //    }

        //    VisualizeCorrelationResults(correlationResults);
        //}

        private static void VisualizeCorrelationResults(Dictionary<string, double> correlationResults)
        {
            var figure = new Plotly.NET.Figure();
            figure.AddBar(new Plotly.NET.Bar
            {
                X = correlationResults.Keys.ToArray(),
                Y = correlationResults.Values.ToArray(),
                Name = "Correlation Values"
            });
            figure.WithTitle("Correlation Values for Different Spline Types");
            figure.WithXTitle("Spline Type");
            figure.WithYTitle("Correlation Value");
            figure.Show();
        }
}
