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

namespace Test.TestDIA
{
    public class PeakCurveGroupingTest
    {

        public class PeakCurveGroupingTestCase
        {
            public MsDataScan[] Scans { get; set; }
            public PeakCurve precursorXIC { get; set; }
            public List<double> ExpectedFragMz { get; set; }
        }

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
            precursorPC.GetScanCycleCubicSplineXYData(0.05);

            var fragmentPeaks1 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(2, i + 1.1, i + 1, ZeroBasedScanNumber: i);
                fragmentPeaks1.Add(peak);
            }
            var fragmentPC1 = new PeakCurve(fragmentPeaks1);
            fragmentPC1.GetScanCycleCubicSplineXYData(0.05);

            var fragmentPeaks2 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(5, i + 1.2, (i + 1) * 1000, ZeroBasedScanNumber: i);
                fragmentPeaks2.Add(peak);
            }
            var fragmentPC2 = new PeakCurve(fragmentPeaks2);
            fragmentPC2.GetScanCycleCubicSplineXYData(0.05);

            //Test correlation calculations
            var corr1 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorPC, fragmentPC2);
            precursorPC.GetScanCycleCubicSplineXYData(0.05f);
            fragmentPC1.GetScanCycleCubicSplineXYData(0.05f);
            fragmentPC2.GetScanCycleCubicSplineXYData(0.05f);
            var corr2 = PrecursorFragmentPair.CalculateCorr_scanCycleSpline_preCalculated(precursorPC, fragmentPC2);
            var corr3 = PrecursorFragmentPair.CalculateCorr_spline(precursorPC, fragmentPC2, "cubic", 0.005);
        }
    }
}
