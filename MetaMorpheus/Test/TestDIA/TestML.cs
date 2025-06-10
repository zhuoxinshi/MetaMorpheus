using Nett;
using System.IO.Compression;
using Chemistry;
using MzLibUtil;
using Accord.Math;
using Accord.Statistics.Models.Regression.Fitting;
using Accord.Statistics.Models.Regression;
using Accord.Statistics;
using Microsoft.ML;
using TaskLayer;
using MassSpectrometry;
using NUnit.Framework;
using EngineLayer;
using EngineLayer.DIA;
using System;
using System.Collections.Generic;
using System.Linq;
using Plotly.NET;

namespace Test.TestDIA
{
    public class TestML
    {
        [Test]
        public static void TestML1()
        {
            var psmFilePath = @"E:\ISD Project\ISD_250428\0429YD_ISD&DDA_SearchVariableO\Task1-SearchTask\Individual File Results\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected_PSMs.psmtsv";
            var dataFilePath = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";
            var myDataFileManager = new MyFileManager(true);

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);
            var diaParam = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.2, correlationCutOff: 0.25, apexRtTolerance: 0.2,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 3000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.NoSpline, ms2SplineType: SplineType.NoSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2);
            task.CommonParameters.DIAparameters = diaParam;

            var dataFile = myDataFileManager.LoadFile(dataFilePath, task.CommonParameters);
            var ms1scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToArray();
            var ms2scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var allMs1MassCurves = ISDEngine_static.GetAllPeakCurves(ms1scans, task.CommonParameters, task.CommonParameters.DIAparameters,
                                              diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance, diaParam.MaxRTRangeMS1, out List<Peak>[] ms1PeaksByScan, diaParam.CutMs1Peaks, null,
                    diaParam.MinMS1Mass, diaParam.MinMS1Charge, diaParam.Ms1NumPeaksThreshold);
            var isdScanVoltageMap = ISDEngine_static.ConstructMs2Groups(ms2scans);
            var allMs2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            var allMs2PeaksByScan = new Dictionary<double, List<Peak>[]>();
            foreach (var ms2Group in isdScanVoltageMap)
            {
                allMs2PeakCurves[ms2Group.Key] = ISDEngine_static.GetAllPeakCurves(ms2Group.Value.ToArray(), task.CommonParameters, task.CommonParameters.DIAparameters,
                    diaParam.Ms2XICType, diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out List<Peak>[] peaksByScan2, diaParam.CutMs2Peaks, null, 
                    diaParam.MinMS2Mass, diaParam.MinMS2Charge, diaParam.Ms2NumPeaksThreshold);
                allMs2PeaksByScan[ms2Group.Key] = peaksByScan2;
            }

            var psms = PsmTsvReader.ReadTsv(psmFilePath, out List<string> warnings).ToList();
            var ms1MassTable = DeconvolutedMass.GetMassTable(ms1PeaksByScan.Where(v => v != null).SelectMany(p => p).ToList(), 1);
            var ms2MassTables = new Dictionary<double, List<DeconvolutedMass>[]>();
            foreach(var group in allMs2PeaksByScan)
            {
                ms2MassTables[group.Key] = DeconvolutedMass.GetMassTable(group.Value.Where(v => v != null).SelectMany(p => p).ToList(), 1);
            }

            var pfGroupMetricsList = new List<PFgroupMetrics>();
            var pfPairsList = new List<PFpairMetrics>();
            var groupMetricsFilePath = @"E:\ISD Project\TestIsdDataAnalysis\TestML_pfGroupMetrics_YD_proteoform.tsv";
            var pairMetricsFilePath = @"E:\ISD Project\TestIsdDataAnalysis\TestML_pfPairMetrics_YD_proteoform.tsv";
            foreach (var psm in psms)
            {
                var pfGroup = PrecursorFragmentsGroup.FindPfGroupBasedOnPsm(psm, ms1MassTable, ms2MassTables, task.CommonParameters.DIAparameters);
                if (pfGroup != null)
                {
                    var pfGroupMetrics = new PFgroupMetrics(pfGroup, psm);
                    pfGroupMetricsList.Add(pfGroupMetrics);
                    foreach(var pair in pfGroup.PFpairs)
                    {
                        var pfPair = new PFpairMetrics(pair, pfGroup, psm);
                        pfPairsList.Add(pfPair);
                    }
                }
            }
            //var pfGroupFile = new DIAPFgroupsMetricsFile { FullFileName = Path.GetFileNameWithoutExtension(dataFilePath), Results = pfGroupMetricsList };
            //pfGroupFile.WriteResults(groupMetricsFilePath);
            //var pfPairFile = new PFpairMetricFile {  Results = pfPairsList };
            //pfPairFile.WriteResults(pairMetricsFilePath);

            double[][] featureMatrix = new double[pfPairsList.Count][];
            for (int i = 0; i < pfPairsList.Count; i++)
            {
                featureMatrix[i] = new double[]
                {
                    pfPairsList[i].FragmentIntensity,
                    pfPairsList[i].ApexRtDelta,
                    pfPairsList[i].Correlation,
                    pfPairsList[i].Overlap
                };
            }

            int[] labels = new int[pfPairsList.Count];
            for (int i = 0; i < pfPairsList.Count; i++)
            {
                labels[i] = pfPairsList[i].TagetDecoy == "T" ? 1 : 0; 
            }

            //var teacher = new IterativeReweightedLeastSquares<LogisticRegression>();
            //var model = teacher.Learn(featureMatrix, labels);

            //var lra = new LogisticRegressionAnalysis();
            //var regression = lra.Learn(featureMatrix, labels);

            var indicesTarget = Enumerable.Range(0, labels.Length).Where(i => labels[i] == 1).ToList();
            var indicesDecoy = Enumerable.Range(0, labels.Length).Where(i => labels[i] == 0).ToList();
            indicesTarget.Shuffle();
            indicesDecoy.Shuffle();
            int nTargetTrain = (int)(0.8 * indicesTarget.Count);
            int nDecoyTrain = (int)(0.8 * indicesDecoy.Count);
            var trainIndices = indicesTarget.Take(nTargetTrain).Concat(indicesDecoy.Take(nDecoyTrain)).ToArray();
            var testIndices = indicesTarget.Skip(nTargetTrain).Concat(indicesDecoy.Skip(nDecoyTrain)).ToArray();

            double[][] trainX = trainIndices.Select(i => featureMatrix[i]).ToArray();
            int[] trainY = trainIndices.Select(i => labels[i]).ToArray();
            double[][] testX = testIndices.Select(i => featureMatrix[i]).ToArray();
            int[] testY = testIndices.Select(i => labels[i]).ToArray();

            var teacher = new IterativeReweightedLeastSquares<LogisticRegression>();
            var model = teacher.Learn(trainX, trainY);
            double[] testScores = testX.Select(x => model.Probability(x)).ToArray();

            //// Create ROC curve from scores and labels
            //var roc = new ReceiverOperatingCharacteristic(testScores, testY);

            //// Compute AUC (Area Under Curve)
            //double auc = roc.Area;
            //double[] thresholds = roc.Thresholds;
            //double[] tpr = roc.TruePositiveRates;
            //double[] fpr = roc.FalsePositiveRates;

            //ML.Net
            var dataList = new List<PFPair>();
            foreach(var pfPair in pfPairsList)
            {
                var data = new PFPair
                {
                    Correlation = (float)pfPair.Correlation,
                    ApexRtDelta = (float)pfPair.ApexRtDelta,
                    Overlap = (float)pfPair.Overlap,
                    FragmentIntensity = (float)pfPair.FragmentIntensity,
                    Label = pfPair.TagetDecoy == "T" ? true : false
                };
                dataList.Add(data);
            }
            TrainModel.Train(dataList);
        }

        public static double ComputeAUC(double[] scores, int[] labels)
        {
            var pairs = scores.Zip(labels, (s, l) => new { Score = s, Label = l })
                              .OrderByDescending(p => p.Score)
                              .ToList();

            double tp = 0, fp = 0;
            double prevTp = 0, prevFp = 0;
            double auc = 0.0;

            double pos = labels.Count(l => l == 1);
            double neg = labels.Count(l => l == 0);

            foreach (var p in pairs)
            {
                if (p.Label == 1)
                    tp++;
                else
                {
                    fp++;
                    auc += tp;  // trapezoidal area contribution
                }
            }

            auc /= (pos * neg); // normalize
            return auc;
        }

        [Test]
        public static void TestLearnFromPfPairMetricsFile()
        {
            var metricsFilePath = @"E:\ISD Project\TestSearch\random\B_rep1_ML\search\PFgrouping\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected_PFpairMetrics.tsv";
            var metricsFile = new PFpairMetricFile { FilePath = metricsFilePath };
            metricsFile.LoadResults();
            var results = metricsFile.Results.Where(pf => pf.TagetDecoy == "T" && pf.MatchedIonType != "Internal" && pf.PsmScore >= 8).ToList();
            var pfPairs = TrainModel.GetPFPairsFromPFMetricsFile(metricsFile);
            TrainModel.Train(pfPairs);
        }
    }
}
