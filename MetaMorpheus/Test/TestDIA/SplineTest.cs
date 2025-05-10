using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using EngineLayer;
using NUnit.Framework;
using System.IO;
using TaskLayer;
using EngineLayer.DIA;
using MzLibUtil;
using Org.BouncyCastle.Asn1.Cmp;
using Nett;
using Plotly.NET.CSharp;

namespace Test.TestDIA
{
    public class SplineTest
    {
        [Test]
        public void Visualize()
        {
            var path = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YB_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected.mzML";
            var myFileManager = new MyFileManager(false);
            var msDataFile = myFileManager.LoadFile(path, new CommonParameters());
            var scans = msDataFile.GetAllScansList();
            var ms1Scans = scans.Where(x => x.MsnOrder == 1).ToArray();
            var ms2Scans = scans.Where(x => x.MsnOrder == 2).ToArray();
            var isd100 = ms2Scans.Where(s => s.ScanFilter.Contains("100")).ToArray();

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.3, correlationCutOff: 0.25, apexRtTolerance: 0.3,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, combineFragments: true);

            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.MassCurve,
                new PpmTolerance(20), 0.5, out List<Peak>[] peaksByScan);
            var allISD100PeakCurves = ISDEngine_static.GetAllPeakCurves(isd100, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.MassCurve,
                               new PpmTolerance(20), 1, out List<Peak>[] peaksByScan2);
            var targetMass = peaksByScan[3489].Where(p => Math.Abs(p.HighestPeakMz - 1122.67) < 0.01 && p.Charge == 10).FirstOrDefault();
            //var targetMass100 = peaksByScan2[2728].Where(p => Math.Abs(p.HighestPeakMz - 1324.01) < 0.01 && p.Charge == 7).FirstOrDefault();
            var targetPC1 = targetMass.PeakCurve;
            //var targetPC2 = targetMass100.PeakCurve;
            targetPC1.GetRawXYData();
            //targetPC2.GetRawXYData();
            //var corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(targetPC1, targetPC2);
            targetPC1.GetGaussianFitXYData();
            targetPC1.VisualizeGeneral("rt").Show();
            //targetPC2.GetGaussianFitXYData();
            //targetPC2.VisualizeGeneral("rt").Show();
            //var corr2 = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(targetPC1, targetPC2);

        }
    }
}
