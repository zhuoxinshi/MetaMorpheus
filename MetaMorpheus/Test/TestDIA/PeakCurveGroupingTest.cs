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
using TaskLayer;
using Microsoft.VisualStudio.TestPlatform.ObjectModel;
using Org.BouncyCastle.Crypto.Encodings;
using Plotly.NET.CSharp;

namespace Test.TestDIA
{
    public class PeakCurveGroupingTest
    {
        public class PeakCurveGroupingTestCase
        {
            public MsDataScan[] Scans { get; set; }
            public double PrecursorMz { get; set; }
            public int PrecursorCharge { get; set; }
            public double ApexRt { get; set; }
            public List<double> ExpectedFragMz { get; set; }
            public List<double> FalseFragMz { get; set; }

            public PeakCurveGroupingTestCase(MsDataScan[] scans, double precursorMz, int precursorCharge, double apexRt, List<double> expectedFragMz, List<double> falseFragMz = null)
            {
                Scans = scans;
                PrecursorMz = precursorMz;
                PrecursorCharge = precursorCharge;
                ApexRt = apexRt;
                ExpectedFragMz = expectedFragMz;
                FalseFragMz = falseFragMz;
            }
        }

        public static List<PeakCurveGroupingTestCase> GetGroupingTestCases()
        {
            var testCases = new List<PeakCurveGroupingTestCase>();

            //Test case 1: bu-ISD, simple test case, single peptide with minimal coelution
            var filePath = @"E:\ISD Project\ISD_bu\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            var myFileManager = new MyFileManager(true);
            var dataFile1 = myFileManager.LoadFile(filePath, new CommonParameters(trimMsMsPeaks: false));
            var scans1 = dataFile1.GetAllScansList().Where(s => s.OneBasedScanNumber >= 7429 && s.OneBasedScanNumber <= 7550).ToArray();
            var expectedFragMz1 = new List<double> { 429.19, 557.25, 686.29, 815.33, 928.42, 1056.48, 1157.52, 1254.58 };
            var testCase1 = new PeakCurveGroupingTestCase (scans1, 858.40, 2, 49.34, expectedFragMz1);
            testCases.Add(testCase1);

            //Test case 2: bu-ISD, 2 peptides coelution
            var scans2 = dataFile1.GetAllScansList().Where(s => s.OneBasedScanNumber >= 4985 && s.OneBasedScanNumber <= 5370).ToArray();
            var expectedFragMz2 = new List<double> {409.20815, 537.30301, 666.34560, 779.42945, 850.46663, 979.50952, 1094.53650, 1209.56405, 1308.63174, 1437.67639, 
                427.21826, 542.24542, 657.27251, 786.31488, 857.35187 }.Select(x => Math.Round(x, 2)).ToList();
            var testCase2 = new PeakCurveGroupingTestCase(scans2, 545.93, 3, 34.39, expectedFragMz2);
            testCases.Add(testCase2);

            //Test case 3: bu-ISD, 2 peptides coelution
            var expectedFragMz3 = new List<double> { 401.28758, 500.35621, 629.39894, 757.45758, 814.48175, 476.27303, 1008.56041, 1079.59769, 1192.68291, 1307.71100, 1378.74866,
                415.18289, 978.46606, 1107.50964}.Select(x => Math.Round(x, 2)).ToList();
            var testCase3 = new PeakCurveGroupingTestCase(scans2, 536.29, 3, 33.91, expectedFragMz3);
            testCases.Add(testCase3);

            return testCases;
        }

        [Test]
        public static void TestGroupingBottomUpRetentionTime()
        {
            var testCases = GetGroupingTestCases();
            
            foreach(var testCase in testCases)
            {
                var ms1Scans = testCase.Scans.Where(s => s.MsnOrder == 1).ToArray();
                var ms2Scans = testCase.Scans.Where(s => s.MsnOrder == 2).ToArray();
                var diaParam = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5), 1, numScanPerCycle: 2, correlationCutOff: 0.9, apexRtTolerance: 0.2,
                    overlapRatioCutOff: 0.3, pfGroupingType: PFGroupingType.RetentionTime);
                var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, new CommonParameters(), diaParam, XICType.DeconHighestPeak, new PpmTolerance(5), 0.5, out List<Peak>[] peaks);
                var allMs2PeakCurves = ISDEngine_static.GetAllPeakCurves(ms2Scans, new CommonParameters(), diaParam, XICType.Peak, new PpmTolerance(5), 0.5, out List<Peak>[] peaks2);
                var precursorPC = allMs1PeakCurves.Where(p => Math.Abs(p.AveragedMz - testCase.PrecursorMz) < 0.01 && Math.Abs(p.ApexRT - testCase.ApexRt) < 0.01
                                    && p.Charge == testCase.PrecursorCharge).First();
                precursorPC.GetRawXYData();
                foreach (var frag in allMs2PeakCurves)
                {
                    frag.GetRawXYData();
                }
                var pfGroup = ISDEngine_static.PFgrouping(precursorPC, allMs2PeakCurves, diaParam);
                pfGroup.VisualizeNormalized().Show();
                var matchedFragPCs = pfGroup.PFpairs.Select(p => p.FragmentPeakCurve).ToList();
                foreach (var frag in testCase.ExpectedFragMz)
                {
                    Assert.That(matchedFragPCs.Any(p => Math.Abs(p.AveragedMz - frag) < 0.01), Is.True);
                }
            }
        }

        //Test Pearson's correlation calculations with fake PeakCurves
        [Test]
        public static void TestCorrelationCalculation()
        {
            var precursorPeaks = new List<Peak>();
            for (int i = 0; i < 8; i++)
            {
                var peak = new Peak(10, i + 1, i + 1, ZeroBasedScanNumber: i);
                precursorPeaks.Add(peak);
            }
            var precursorPC = new PeakCurve(precursorPeaks);

            var fragmentPeaks1 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(2, i + 1.1, i + 1, ZeroBasedScanNumber: i);
                fragmentPeaks1.Add(peak);
            }
            var fragmentPC1 = new PeakCurve(fragmentPeaks1);

            var fragmentPeaks2 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(5, i + 1.2, (i + 1) * 1000, ZeroBasedScanNumber: i);
                fragmentPeaks2.Add(peak);
            }
            var fragmentPC2 = new PeakCurve(fragmentPeaks2);

            //Test correlation calculations
            //Test using raw data with scan cycles
            precursorPC.GetRawXYData();
            fragmentPC1.GetRawXYData();
            var corr1 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorPC, fragmentPC1);
            Assert.That(corr1 == 1, Is.True);

            //Test using cubic spline
            precursorPC.GetScanCycleCubicSplineXYData(0.05f);
            precursorPC.GetCubicSplineXYData(0.005);
            fragmentPC1.GetCubicSplineXYData(0.005);
            var corr2 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorPC, fragmentPC1);
            Assert.That(corr2, Is.EqualTo(1).Within(0.01));

            //Test using cubic spline
            precursorPC.GetBSplineXYData(0.005, 2);
            fragmentPC2.GetBSplineXYData(0.005, 2);
            var corr4 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorPC, fragmentPC2);
            Assert.That(corr4, Is.EqualTo(1).Within(0.01));

            //Test converting ms2 retention times to ms1 retention times
            var rtMap = new Dictionary<double, double>();
            for (int i = 0; i < 6; i++)
            {
                rtMap[i + 1.1] = 1 + i;
                rtMap[i + 1.2] = 1 + i;
            }
            precursorPC.GetCubicSplineXYData(0.005);
            fragmentPC1.GetMs1SpaceCubicSplineXYData(rtMap, 0.005);
            var corr3 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorPC, fragmentPC1);
            Assert.That(corr3, Is.EqualTo(1));

            fragmentPC2.GetMs1SpaceBSplineXYData(0.005, 2, rtMap);
            var corr5 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorPC, fragmentPC2);
        }

        [Test]
        public static void TestPseudoMs2Construction()
        {
            //make some fake PeakCurves
            var precursorPeaks = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(10, i + 1, i + 1, ZeroBasedScanNumber: i);
                precursorPeaks.Add(peak);
            }
            var precursorPC = new PeakCurve(precursorPeaks);
            precursorPC.MonoisotopicMass = 50;
            precursorPC.Charge = 10;

            var fragmentPeaks1 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(2, i + 1.1, i + 1, ZeroBasedScanNumber: i);
                fragmentPeaks1.Add(peak);
            }
            var fragmentPC1 = new PeakCurve(fragmentPeaks1);
            fragmentPC1.MonoisotopicMass = 10;
            fragmentPC1.Charge = 5;

            var fragmentPeaks2 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(5, i + 1.2, (i + 1) * 1000, ZeroBasedScanNumber: i);
                fragmentPeaks2.Add(peak);
            }
            var fragmentPC2 = new PeakCurve(fragmentPeaks2);
            fragmentPC2.MonoisotopicMass = 10;
            fragmentPC2.Charge = 2;

            //make a fake pfGroup with the fake PeakCurves
            var pfPair1 = new PrecursorFragmentPair(precursorPC, fragmentPC1, 0.9);
            var pfPair2 = new PrecursorFragmentPair(precursorPC, fragmentPC2, 0.9);
            var pfGroup = new PrecursorFragmentsGroup(precursorPC);
            pfGroup.PFpairs.AddRange(new List<PrecursorFragmentPair> { pfPair1, pfPair2 });

            //Test two different ways of making new ms2withspecificmass
            var ms2WithPre_mzPeak = ISDEngine_static.ConstructNewMs2Scans(pfGroup, new CommonParameters(), PseudoMs2ConstructionType.mzPeak, "dataFile");
            Assert.That(ms2WithPre_mzPeak.PrecursorMass == precursorPC.MonoisotopicMass);
            Assert.That(ms2WithPre_mzPeak.TheScan.MassSpectrum.XArray.Length == 2);
            Assert.That(ms2WithPre_mzPeak.TheScan.MassSpectrum.XArray[0] == 2);
            Assert.That(ms2WithPre_mzPeak.TheScan.MassSpectrum.XArray[1] == 5);
            Assert.That(ms2WithPre_mzPeak.TheScan.MassSpectrum.SumOfAllY == fragmentPC1.AveragedIntensity + fragmentPC2.AveragedIntensity);

            //This method is used when DeconHighestPeak is used to find all Ms2 XICs
            var ms2WithPre_mass = ISDEngine_static.ConstructNewMs2Scans(pfGroup, new CommonParameters(), PseudoMs2ConstructionType.neutralMass, "dataFile");
            Assert.That(ms2WithPre_mass.PrecursorMass == precursorPC.MonoisotopicMass);
            Assert.That(ms2WithPre_mass.ExperimentalFragments.Length == 2);
            Assert.That(ms2WithPre_mass.ExperimentalFragments[0].MonoisotopicMass == 10);
            Assert.That(ms2WithPre_mass.ExperimentalFragments[1].MonoisotopicMass == 10);
        }
    }
}
