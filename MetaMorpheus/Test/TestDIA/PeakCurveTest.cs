using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using Microsoft.VisualStudio.TestPlatform.ObjectModel;
using TaskLayer;
using ScottPlot;
using MathNet.Numerics.Statistics;
using ScottPlot.Plottables;
using UsefulProteomicsDatabases;
using EngineLayer.DIA.CWT;


namespace Test.TestDIA
{
    [TestFixture]
    public class PeakCurveTests
    {
        static PeakCurveTests()
        {
            Loaders.LoadElements();
        }

        public class AllPeakCurveTestCase
        {
            public MsDataScan[] FakeScans { get; set; }
            public int NumScansPerCycle { get; set; }
            public int ExpectedPeakCurvesCount { get; set; }
            public List<double> ExpectedAverageMz { get; set; }
        }

        public class SinglePeakCurveTestCase
        {
            public PeakCurve PeakCurve { get; set; }
            public double ExpectedAveragedMz { get; set; }
            public double ExpectedAveragedIntensity { get; set; }
            public double ExpectedApexRT { get; set; }
            public int ExpectedApexCycle { get; set; }
            public double ExpectedStartRT { get; set; }
            public double ExpectedEndRT { get; set; }
            public int ExpectedStartCycle { get; set; }
            public int ExpectedEndCycle { get; set; }

            public SinglePeakCurveTestCase(PeakCurve peakCurve, double expectedMz, double expectedApexRt, double expectedStartRt,  double expectedEndRt, int expectedApexCycle, int expectedStartCycle, int expectedEndCycle)
            {
                PeakCurve = peakCurve;
                ExpectedAveragedMz = expectedMz;
                ExpectedApexRT = expectedApexRt;
                ExpectedStartRT = expectedStartRt;
                ExpectedEndRT = expectedEndRt;
                ExpectedApexCycle = expectedApexCycle;
                ExpectedStartCycle = expectedStartCycle;
                ExpectedEndCycle = expectedEndCycle;
            }
        }

        //Get all test cases for extracting all XICs from a set of scans
        public static IEnumerable<AllPeakCurveTestCase> GetAllPeakCurveTestCases()
        {
            var testCases = new List<AllPeakCurveTestCase>();
            AllPeakCurveTestCases = testCases; 

            // Test Case 1: Fake data with 5 scans
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
            var scans = new MsDataScan[] { scan1, scan2, scan3, scan4, scan5 };

            //Theoretical averaged Mz values for each PeakCurve
            var averageMz3 = (400 * 10 + 400.001 * 1001 + 400.002 * 100) / (10 + 1001 + 100);
            var averageMz4 = (400.005 * 1000 + 400.0065 * 100) / (1000 + 100);
            var averageMz2 = (500 * 20 + 500.001 * 2001 + 500.002 * 200 + 500.0025 * 2000 + 500.0015 * 200) / (20 + 2001 + 200 + 2000 + 200);
            var averageMz1 = (1000 * 30 + 1000.001 * 3000 + 1000.002 * 300 + 1000.003 * 3001 + 1000.005 * 300) / (30 + 3000 + 300 + 3001 + 300);

            var testCase1 = new AllPeakCurveTestCase
            {
                FakeScans = scans,
                NumScansPerCycle = 2,
                ExpectedPeakCurvesCount = 4,
                ExpectedAverageMz = new List<double> { averageMz1, averageMz2, averageMz3, averageMz4 }
            };
            testCases.Add(testCase1);

            yield return testCase1;
        }

        public static List<AllPeakCurveTestCase> AllPeakCurveTestCases { get; set; }

        //Test the method to get all peak curves from a set of scans
        [Test]
        [TestCaseSource(nameof(GetAllPeakCurveTestCases))]
        public static void TestGetAllPeakCurves(AllPeakCurveTestCase testCase)
        {
            var diaParam = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20));
            var allPeakCurves = ISDEngine_static.GetAllPeakCurves_Peak(testCase.FakeScans, diaParam, new PpmTolerance(5), 2,
                out List<Peak>[] allPeaksByScan, false);
            Assert.That(allPeakCurves.Count == testCase.ExpectedPeakCurvesCount);
            Assert.That(allPeakCurves.Sum(pc => pc.Peaks.Count) == testCase.FakeScans.Sum(s => s.MassSpectrum.Size));

            for (int i = 0; i < allPeakCurves.Count; i++)
            {
                Assert.That(allPeakCurves[i].AveragedMz, Is.EqualTo(testCase.ExpectedAverageMz[i]).Within(0.0000001));
            }
        }

        //Test single PeakCurve properties with fake data
        [Test]
        public static void TestSinglePeakCurveFakeData()
        {
            var testCases = new List<SinglePeakCurveTestCase>();

            //Test Case 1
            var scans = AllPeakCurveTestCases[0].FakeScans;
            var allPeakCurves = ISDEngine_static.GetAllPeakCurves_Peak(scans, new DIAparameters(new PpmTolerance(5), new PpmTolerance(20)), new PpmTolerance(5), 2,
                               out List<Peak>[] allPeaksByScan, false);
            var fakeCurve = allPeakCurves[0];
            double expectedMz1 = (1000 * 30 + 1000.001 * 3000 + 1000.002 * 300 + 1000.003 * 3001 + 1000.005 * 300) / (30 + 3000 + 300 + 3001 + 300);
            double expectedAveragedIntensity = (30 * 30 + 3000 * 3000 + 300 * 300 + 3001 * 3001 + 300 * 300) / (30 + 3000 + 300 + 3001 + 300);
            var testCase1 = new SinglePeakCurveTestCase(fakeCurve, expectedMz1, 0.4, 0.1, 0.5, 3, 0, 4);
            testCase1.ExpectedAveragedIntensity = expectedAveragedIntensity;
            testCases.Add(testCase1);

            foreach(var testCase in testCases)
            {
                Assert.That(testCase.PeakCurve.AveragedMz, Is.EqualTo(testCase.ExpectedAveragedMz));
                Assert.That(testCase.PeakCurve.AveragedIntensity, Is.EqualTo(testCase.ExpectedAveragedIntensity));
                Assert.That(testCase.PeakCurve.ApexRT, Is.EqualTo(testCase.ExpectedApexRT));
                Assert.That(testCase.PeakCurve.ApexCycle, Is.EqualTo(testCase.ExpectedApexCycle));
                Assert.That(testCase.PeakCurve.StartRT, Is.EqualTo(testCase.ExpectedStartRT));
                Assert.That(testCase.PeakCurve.EndRT, Is.EqualTo(testCase.ExpectedEndRT));
                Assert.That(testCase.PeakCurve.StartCycle, Is.EqualTo(testCase.ExpectedStartCycle));
                Assert.That(testCase.PeakCurve.EndCycle, Is.EqualTo(testCase.ExpectedEndCycle));
            }
        }

        //Test single PeakCurve properties with real data
        [Test]
        public static void TestSinglePeakCurveRealData()
        {
            var testCases = new List<SinglePeakCurveTestCase>();

            var path = @"E:\ISD Project\ISD_bu\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            var myFileManagers = new MyFileManager(true);
            var commonParam = new CommonParameters(trimMsMsPeaks: false);
            var dataFile = myFileManagers.LoadFile(path, commonParam);

            //Test Case 1
            var realScans1 = dataFile.GetAllScansList().Where(s => s.RetentionTime >= 32 && s.RetentionTime <= 39).ToArray();
            var diaParam1 = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), numScanPerCycle: 2, maxNumMissedScan: 1);
            var allPeakCurves1 = ISDEngine_static.GetAllPeakCurves_Peak(realScans1, diaParam1, new PpmTolerance(5), maxRTRange: 0.5, out List<Peak>[] allPeaksByScan2, false);
            var testPC1 = allPeakCurves1.Where(p => Math.Abs(p.AveragedMz - 838.05) < 0.01).First();
            var testCase1 = new SinglePeakCurveTestCase(testPC1, 838.05, 33.91, 33.41, 34.41, 2562, 2523, 2601);
            testCases.Add(testCase1);

            //Test Case 2
            var testPC2 = allPeakCurves1.Where(p => Math.Abs(p.AveragedMz - 569.28) < 0.01).First();
            var testCase2 = new SinglePeakCurveTestCase(testPC2, 569.28, 33.04, 32.75, 33.5, 2495, 2470, 2529);

            //Test Case 3
            var testPC3 = allPeakCurves1.Where(p => Math.Abs(p.AveragedMz - 673.85) < 0.01).First();
            var testCase3 = new SinglePeakCurveTestCase(testPC3, 673.85, 33.48, 33.1, 33.9, 2529, 2495, 2560);

            //loop through all the testCases
            foreach (var testCase in testCases)
            {
                Assert.That(testCase.PeakCurve.ApexRT, Is.EqualTo(testCase.ExpectedApexRT).Within(0.01));
                Assert.That(testCase.PeakCurve.ApexCycle, Is.EqualTo(testCase.ExpectedApexCycle));
                Assert.That(testCase.PeakCurve.StartRT, Is.EqualTo(testCase.ExpectedStartRT).Within(0.25));
                Assert.That(testCase.PeakCurve.EndRT, Is.EqualTo(testCase.ExpectedEndRT).Within(0.25));
                Assert.That(testCase.PeakCurve.StartCycle, Is.EqualTo(testCase.ExpectedStartCycle).Within(12));
                Assert.That(testCase.PeakCurve.EndCycle, Is.EqualTo(testCase.ExpectedEndCycle).Within(12));
            }
        }

        //Test the dictionary used to convert ms2 retention time to the corresponding ms1 retention time
        [Test]
        public static void TestRtMap()
        {
            var path = @"E:\DIA\Data\DIA-bottom-up_241218\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, new CommonParameters());
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var ms1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToArray();
            var rtMap = ISDEngine_static.GetRtMap(ms1Scans, ms2Scans);
            Assert.That(rtMap.Count == ms2Scans.Length);
            foreach(var scan in ms2Scans)
            {
                var ms1Scan = ms1Scans.First(s => s.OneBasedScanNumber == scan.OneBasedPrecursorScanNumber);
                Assert.That(rtMap[scan.RetentionTime], Is.EqualTo(ms1Scan.RetentionTime));
            }
        }

        [Test]
        public static void TestFindPrecursors()
        {
            var path = @"E:\ISD Project\ISD_bu\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            var myFileManagers = new MyFileManager(true);
            var commonParam = new CommonParameters(trimMsMsPeaks: false);
            var dataFile = myFileManagers.LoadFile(path, commonParam);
            var ms1scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToArray();

            var allEnvelopes = new List<IsotopicEnvelope>();
            foreach(var scan in ms1scans)
            {
                var envelopes = Deconvoluter.Deconvolute(scan, commonParam.PrecursorDeconvolutionParameters);
                allEnvelopes.AddRange(envelopes);
            }
        }

        //unfinished
        [Test]
        public static void SplineTest()
        {
            var path = @"E:\ISD Project\ISD_240606\06-09-24_mix_sample10_5uL_ISD.mzML";
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, new CommonParameters());
            var allMs1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToArray();
            var allMs2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var ms1Scans = allMs1Scans.Where(s => s.RetentionTime >= 31.96 && s.RetentionTime <= 35).ToArray();
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
            var corr3 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, pc2);

            var splineTypes = new List<SplineType> { SplineType.NoSpline, SplineType.CubicSpline, SplineType.CubicSplineSavgolSmoothed,
                           SplineType.ScanCycleCubicSplineSavgolSmoothed, SplineType.Ms1SpaceCubicSpline, SplineType.Ms1SpaceCubicSplineSavgolSmoothed };
            var correlations = GetAllCorr(pc1, new List<PeakCurve> { frag1, frag2 }, splineTypes, 0.05, 7, rtMap);
            var plot = CreateDensityPlot(correlations, 0.1);
            plot.SavePng(@"E:\ISD Project\ISD_240606\06-09-24_mix_sample10_5uL_ISD.png", 1600, 1200);
        }

        //unfinished
        [Test]
        public static Dictionary<SplineType, double> TestAllSplineTypesCorrelation(PeakCurve pc1, PeakCurve pc2, List<SplineType> splineTypes,
            double splineRtInterval, int smoothWindowSize, Dictionary<double, double> rtMap)
        {
            var correlationResults = new Dictionary<SplineType, double>();

            foreach (var splineType in splineTypes)
            {
                switch (splineType)
                {
                    case SplineType.NoSpline:
                        pc1.GetRawXYData();
                        pc2.GetRawXYData();
                        break;
                    case SplineType.CubicSpline:
                        pc1.GetCubicSplineXYData(splineRtInterval);
                        pc2.GetCubicSplineXYData(splineRtInterval);
                        break;
                    case SplineType.CubicSplineSavgolSmoothed:
                        pc1.GetCubicSplineSavgolSmoothedXYData(smoothWindowSize, splineRtInterval);
                        pc2.GetCubicSplineSavgolSmoothedXYData(smoothWindowSize, splineRtInterval);
                        break;
                    case SplineType.ScanCycleCubicSplineSavgolSmoothed:
                        pc1.GetScanCycleCubicSplineSavgolSmoothedXYData(smoothWindowSize, splineRtInterval);
                        break;
                    case SplineType.Ms1SpaceCubicSpline:
                        pc1.GetCubicSplineXYData(splineRtInterval);
                        pc2.GetMs1SpaceCubicSplineXYData(rtMap, splineRtInterval);
                        break;
                    case SplineType.Ms1SpaceCubicSplineSavgolSmoothed:
                        pc1.GetCubicSplineSavgolSmoothedXYData(smoothWindowSize, splineRtInterval);
                        pc2.GetMs1SpaceCubicSplineSavgolSmoothedXYData(rtMap, smoothWindowSize, splineRtInterval);
                        break;
                }
                var correlation = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(pc1, pc2);
                correlationResults[splineType] = correlation;
            }
            return correlationResults;
        }

        public static Dictionary<SplineType, List<double>> GetAllCorr(PeakCurve precursorPC, List<PeakCurve> fragmentPCs, List<SplineType> splineTypes,
                       double splineRtInterval, int smoothWindowSize, Dictionary<double, double> rtMap)
        {
            var allCorrelations = new Dictionary<SplineType, List<double>>();
            foreach (var fragmentPC in fragmentPCs)
            {
                var correlationResults = TestAllSplineTypesCorrelation(precursorPC, fragmentPC, splineTypes, splineRtInterval, smoothWindowSize, rtMap);
                foreach (var splineType in correlationResults.Keys)
                {
                    if (!allCorrelations.ContainsKey(splineType))
                    {
                        allCorrelations[splineType] = new List<double>();
                    }
                    allCorrelations[splineType].Add(correlationResults[splineType]);
                }
            }
            return allCorrelations;
        }

        public static ScottPlot.Plot CreateDensityPlot(Dictionary<SplineType, List<double>> data, double bandWidth)
        {
            // Create a new plot
            var plt = new ScottPlot.Plot();

            foreach (var kvp in data)
            {
                string category = kvp.Key.ToString();
                double[] values = kvp.Value.ToArray();

                // Determine the range for x-axis
                double min = Statistics.Minimum(values) - 3 * bandWidth;
                double max = Statistics.Maximum(values) + 3 * bandWidth;
                int points = 100;
                double step = (max - min) / (points - 1);
                double[] x = new double[points];
                double[] y = new double[points];

                for (int i = 0; i < points; i++)
                {
                    x[i] = min + i * step;
                    y[i] = KernelDensity.EstimateGaussian(x[i], bandWidth, values);
                }

                // Add the density curve to the plot
                var scatter = plt.Add.Scatter(x, y); // Create the scatter plot
                scatter.Label = category;
            }
            // Customize the plot
            plt.Title("Density Plot of Correlation Values by Enum Type");
            plt.XLabel("Correlation Values");
            plt.YLabel("Density");
            plt.Legend.IsVisible = true;

            return plt;
        }

        //[Test]
        //[TestCaseSource(nameof(GetSinglePeakCurveTestCasesRealData))]
        //public static void TestPython(SinglePeakCurveTestCase testCase)
        //{
        //    var pc = testCase.PeakCurve;
        //    pc.GetCubicSplineSavgolSmoothedXYData(11, 0.05);
        //    var y = pc.XYData.Select(i => i.Item2).ToArray();
        //    double[] widths = Enumerable.Range(1, 30).Select(p => (double)p).ToArray();
        //    var results = Scipy_signal.FindPeaks_cwt(y, widths);
        //}
    }
}
