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
            var averagedPath = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YB_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_averaged_labelCorrected.mzML";
            
            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig); 
            var myFileManager = new MyFileManager(false);
            var msDataFile = myFileManager.LoadFile(path, task.CommonParameters);
            var averagedFile = myFileManager.LoadFile(averagedPath, task.CommonParameters);
            var scans = msDataFile.GetAllScansList();
            var ms1Scans = scans.Where(x => x.MsnOrder == 1).ToArray();
            var ms2Scans = scans.Where(x => x.MsnOrder == 2).ToArray();
            var isd100 = ms2Scans.Where(s => s.ScanFilter.Contains("100")).ToArray();
            var averagedIsd100Scans = averagedFile.GetAllScansList().Where(s => s.ScanFilter.Contains("100")).ToArray();

            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(500), new PpmTolerance(500),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.3, correlationCutOff: 0.25, apexRtTolerance: 0.3,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, combineFragments: true);

            //var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.MassCurve,
            //    new PpmTolerance(20), 0.5, out List<Peak>[] peaksByScan);
            var allISD100PeakCurves = ISDEngine_static.GetAllPeakCurves(isd100, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.MassCurve,
                               new PpmTolerance(20), 0.5, out List<Peak>[] peaksByScan2).Where(pc => pc.Peaks.Count >= 5).ToList();
            foreach(var pc in allISD100PeakCurves)
            {
                pc.ScanCycleSpline(0.05);
            }
            var allAveragedIsd100PeakCurves = ISDEngine_static.GetAllPeakCurves(averagedIsd100Scans, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.MassCurve,
                               new PpmTolerance(20), 0.5, out List<Peak>[] peaksByScan3);
            //var targetMass = peaksByScan[3089].Where(p => Math.Abs(p.HighestPeakMz - 967.40) < 0.01 && p.Charge == 10).FirstOrDefault();
            //var targetMass2 = peaksByScan[3089].Where(p => Math.Abs(p.HighestPeakMz - 1000.3) < 0.01 && p.Charge == 9).FirstOrDefault();
            var targetMass100 = peaksByScan2[3588].Where(p => Math.Abs(p.HighestPeakMz - 831.39) < 0.01 && p.Charge == 2).FirstOrDefault();
            //var targetPC1 = targetMass.PeakCurve;
            //var targetPC2 = targetMass2.PeakCurve;
            var targetPC3 = targetMass100.PeakCurve;
            targetPC3.VisualizeGeneral("rt").Show();
            var averaged = peaksByScan3[3588].Where(p => Math.Abs(p.HighestPeakMz - 831.39) < 0.01 && p.Charge == 2).FirstOrDefault();
            targetPC3.ScanCycleSpline(0.05);
            targetPC3.VisualizeGeneral("cycle").Show();
            targetPC3.GetExtendedCycleCubicSplineXYData(0.05);
            targetPC3.VisualizeGeneral("cycle").Show();
            averaged.PeakCurve.VisualizeGeneral("rt").Show();
            averaged.PeakCurve.ScanCycleSpline(0.05);
            averaged.PeakCurve.VisualizeGeneral("cycle").Show();

            var scan3588 = peaksByScan2[3588].OrderBy(p => p.HighestPeakMz).ToList();
            var averagedScan3588 = peaksByScan3[3588].OrderBy(p => p.HighestPeakMz).ToList();
            var mass2 = scan3588.Where(p => Math.Abs(p.HighestPeakMz - 882.45) < 0.01 && p.Charge == 3).FirstOrDefault();
            mass2.PeakCurve.ScanCycleSpline(0.05);
            mass2.PeakCurve.VisualizeGeneral("cycle").Show();
            mass2.PeakCurve.GetUmpireBSplineData(150, 2);
            mass2.PeakCurve.VisualizeGeneral("cycle").Show();
            var averagedMass2 = averagedScan3588.Where(p => Math.Abs(p.HighestPeakMz - 882.92) < 0.01 && p.Charge == 3).FirstOrDefault();
            averagedMass2.PeakCurve.VisualizeGeneral("rt").Show();
            int stop3 = 0;
            var mass3 = scan3588.Where(p => Math.Abs(p.HighestPeakMz - 925.15) < 0.01 && p.Charge == 3).FirstOrDefault();
            mass3.PeakCurve.VisualizeGeneral("rt").Show();
            var averagedMass3 = averagedScan3588.Where(p => Math.Abs(p.HighestPeakMz - 925.15) < 0.01 && p.Charge == 3).FirstOrDefault();
            averagedMass3.PeakCurve.VisualizeGeneral("rt").Show();
        }

        [Test]
        public static void SharedXICCalculationTest()
        {
            var peakList = new List<(double, double)> { (1, 0.2), (2, 0.6), (3, 0.2) };
            var area = PrecursorFragmentPair.CalculateNormalizedArea(peakList);
            Assert.That(area, Is.EqualTo(1));
        }

        [Test]
        public static void CycleSpline()
        {
            var filePath1 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected.mzML";
            var filePath2 = @"E:\ISD Project\ISD_250428\0504YB_ISD_rep_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected-calib-averaged.mzML";

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(500), new PpmTolerance(500),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.3, correlationCutOff: 0.25, apexRtTolerance: 0.3,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, combineFragments: true);

            var myFileManager = new MyFileManager(false);
            var msDataFile = myFileManager.LoadFile(filePath1, task.CommonParameters);
            var averagedFile = myFileManager.LoadFile(filePath2, task.CommonParameters);
            var scans = msDataFile.GetAllScansList();
            var ms1Scans = scans.Where(x => x.MsnOrder == 1).ToArray();
            var allMs1PCs = ISDEngine_static.GetAllPeakCurves(ms1Scans, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.MassCurve,
                               new PpmTolerance(20), 0.5, out List<Peak>[] peaksByScan).Where(pc => pc.Peaks.Count >= 5).ToList();
            foreach(var pc in allMs1PCs)
            {
                pc.ScanCycleSpline(0.05);
            }
            var ms2Scans = scans.Where(x => x.MsnOrder == 2).ToArray();
            var isd100 = ms2Scans.Where(s => s.ScanFilter.Contains("100")).ToArray();
            var averagedIsd100Scans = averagedFile.GetAllScansList().Where(s => s.ScanFilter.Contains("100")).ToArray();
        }

        [Test]
        public static void VisualizeFake()
        {
            var peakList = new List<Peak>();
            double[] intensityMultipliers = { 1, 3, 1};

        }
    }
}
