using Chemistry;
using Microsoft.ML;
using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using System.Linq.Expressions;
using CsvHelper.Configuration;
using CsvHelper;
using Readers;
using System.Globalization;
using System.IO;
using Microsoft.ML.Trainers.LightGbm;
using BayesianEstimation;
using MathNet.Numerics.Distributions;
using Microsoft.ML.Data;
using System.Data;
using Proteomics.AminoAcidPolymer;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.DIA
{
    public class PFPairFeature
    {
        public float Correlation { get; set; }   // e.g., XIC correlation
        public float ApexRtDelta { get; set; }   // e.g., RT diff
        public float Overlap { get; set; }      // e.g., XIC overlap
        public float FragmentIntensity { get; set; } // e.g., fragment intensity
        public float NormalizedIntensityRank { get; set; } // e.g., normalized intensity rank
        public float PsmScore { get; set; } 
        public float SharedXIC { get; set; } 
        public bool Label { get; set; }       
        public float Weight { get; set; }

        public PFPairFeature() { }

        public PFPairFeature(PeakCurve precursor, PeakCurve fragment)
        {
            Correlation = (float)PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursor, fragment);
            ApexRtDelta = (float)Math.Abs(precursor.ApexRT - fragment.ApexRT);
            Overlap = (float)PrecursorFragmentPair.CalculateOverlapAreaRatio(precursor, fragment);
            FragmentIntensity = (float)fragment.ApexIntensity;
            SharedXIC = (float)PrecursorFragmentPair.CalculateSharedXIC(precursor, fragment);
        }

        public PFPairFeature(PFpairMetrics pfPairMetrics)
        {
            Correlation = (float)pfPairMetrics.Correlation;
            ApexRtDelta = (float)pfPairMetrics.ApexRtDelta;
            Overlap = (float)pfPairMetrics.Overlap;
            FragmentIntensity = (float)pfPairMetrics.FragmentIntensity;
            SharedXIC = (float)pfPairMetrics.SharedXIC;
            NormalizedIntensityRank = (float)pfPairMetrics.NormalizedIntensityRank;
            if (pfPairMetrics.MatchedIonType == "NA") Label = false;
            if (pfPairMetrics.MatchedIonType == "Terminal") Label = true;
            if (pfPairMetrics.MatchedIonType == "NeutralLoss") Label = true;
        }

        public PFPairFeature(PrecursorFragmentPair pfPair)
        {
            Correlation = (float)pfPair.Correlation;
            ApexRtDelta = (float)Math.Abs(pfPair.PrecursorPeakCurve.ApexRT - pfPair.FragmentPeakCurve.ApexRT);
            Overlap = (float)pfPair.Overlap;
            FragmentIntensity = (float)pfPair.FragmentPeakCurve.ApexIntensity;
            SharedXIC = (float)pfPair.SharedXIC;
            NormalizedIntensityRank = (float)pfPair.NormalizedIntensityRank;
        }

        public static IEnumerable<PFPairFeature> ConvertPfPairMetricsToPfPairFeatures(IEnumerable<PFpairMetrics> pfPairMetricsList)
        {
            foreach (var pfPairMetrics in pfPairMetricsList)
            {
                yield return new PFPairFeature(pfPairMetrics);
            }
        }
    }

    public class MyInput
    {
        [VectorType(5)] // Because you have 5 features
        public float[] features { get; set; }

        public MyInput(float[] features)
        {
            this.features = features;
        }

        public MyInput(PeakCurve precursor, PeakCurve fragment)
        {
            features = new float[5];
            features[0] = (float)PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursor, fragment);
            features[1] = (float)Math.Abs(precursor.ApexRT - fragment.ApexRT);
            features[2] = (float)PrecursorFragmentPair.CalculateSharedXIC(precursor, fragment);
            features[3] = (float)fragment.ApexIntensity;
        }

        public MyInput(PrecursorFragmentPair pfPair)
        {
            features = new float[5];
            features[0] = (float)pfPair.Correlation;
            features[1] = (float)Math.Abs(pfPair.PrecursorPeakCurve.ApexRT - pfPair.FragmentPeakCurve.ApexRT);
            features[2] = (float)pfPair.SharedXIC;
            features[3] = (float)pfPair.FragmentPeakCurve.ApexIntensity;
            features[4] = (float)pfPair.NormalizedIntensityRank;
        }
    }

    public class PFpairPrediction
    {
        [ColumnName("PredictedLabel")]
        public bool Prediction;

        public float Probability;

        public float Score;
    }

    public class PFPairFeatureGbm
    {
        public bool Label { get; set; }
        [VectorType(4)]
        public float[] Features { get; set; }

        public PFPairFeatureGbm() { }
        public PFPairFeatureGbm(PFPairFeature pfPairFeature, List<string> featureList)
        {
            Features = new float[4];
            for (int i = 0; i < featureList.Count; i++)
            {
                Features[i] = (float)pfPairFeature.GetType().GetProperty(featureList[i])?.GetValue(pfPairFeature);
            }
        }

        public static IEnumerable<PFPairFeatureGbm> ConvertFeatureSamplesToGbm(IEnumerable<PFPairFeature> pfPairFeatures, List<string> featureList)
        {
            foreach(var pfPairFeature in pfPairFeatures)
            {
                yield return new PFPairFeatureGbm(pfPairFeature, featureList);
            }
        }
    }

    public class PFpairFeatureFile : ResultFile<PFPairFeature>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public PFpairFeatureFile() : base() { }
        public PFpairFeatureFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<PFPairFeature>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<PFPairFeature>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }

    public class TrainModel
    {
        public static ITransformer TrainGradientBoostedTree(IEnumerable<PFpairMetrics> pairMetrics, List<string> featureList, string modelPath = null)
        {
            var samples = PFPairFeature.ConvertPfPairMetricsToPfPairFeatures(pairMetrics);
            var filteredSamples = UnderSample(samples); 
            var gbmSamples = PFPairFeatureGbm.ConvertFeatureSamplesToGbm(filteredSamples, featureList);

            var mlContext = new MLContext();
            IDataView data = mlContext.Data.LoadFromEnumerable(gbmSamples);

            var pipeline = mlContext.BinaryClassification.Trainers.LightGbm(
                new LightGbmBinaryTrainer.Options
                {
                    //LabelColumnName = "Label",
                    //FeatureColumnName = "Features",
                    //NumberOfLeaves = 64,
                    //MinimumExampleCountPerLeaf = 10,
                    //LearningRate = 0.1,
                    //NumberOfIterations = 100,
                    //Seed = 42
                });
            var split = mlContext.Data.TrainTestSplit(data, testFraction: 0.2);
            var trainData = split.TrainSet;
            var testData = split.TestSet;
            var model = pipeline.Fit(data);
            var predictions = model.Transform(testData);
            var metrics = mlContext.BinaryClassification.Evaluate(predictions);
            if (modelPath!=null)
            {
                mlContext.Model.Save(model, data.Schema, modelPath);
            }
            return model;
        }

        public static ITransformer TrainLogisticRegression(IEnumerable<PFpairMetrics> pairMetrics, List<string> featureList, string modelPath = null)
        {
            var samples = PFPairFeature.ConvertPfPairMetricsToPfPairFeatures(pairMetrics);
            var filteredSamples = UnderSample(samples);

            var mlContext = new MLContext();
            IDataView data = mlContext.Data.LoadFromEnumerable(filteredSamples);

            var pipeline = mlContext.Transforms.Concatenate("Features", featureList.ToArray())
                .Append(mlContext.BinaryClassification.Trainers.SdcaLogisticRegression(labelColumnName: "Label", featureColumnName: "Features")); 
            var split = mlContext.Data.TrainTestSplit(data, testFraction: 0.2);
            var trainData = split.TrainSet;
            var testData = split.TestSet;
            var model = pipeline.Fit(trainData);
            var predictions = model.Transform(testData);
            var metrics = mlContext.BinaryClassification.Evaluate(predictions);

            if (modelPath != null)
            {
                mlContext.Model.Save(model, data.Schema, modelPath);
            }
            return model;
        }

        public static ITransformer TrainFastTree(IEnumerable<PFpairMetrics> pairMetrics, List<string> featureList, string modelPath = null)
        {
            var samples = PFPairFeature.ConvertPfPairMetricsToPfPairFeatures(pairMetrics);
            var filteredSamples = UnderSample(samples);

            var mlContext = new MLContext();
            IDataView data = mlContext.Data.LoadFromEnumerable(filteredSamples);

            var pipeline = mlContext.Transforms.Concatenate("Features", featureList.ToArray())
                .Append(mlContext.BinaryClassification.Trainers.FastTree(
                        labelColumnName: "Label",
                        featureColumnName: "Features",
                        numberOfTrees: 100,  // Example: 100 trees
                        learningRate: 0.2f // Example: learning rate of 0.2
                    )); 
            var split = mlContext.Data.TrainTestSplit(data, testFraction: 0.2);
            var trainData = split.TrainSet;
            var testData = split.TestSet;
            var model = pipeline.Fit(trainData);
            var predictions = model.Transform(testData);
            var metrics = mlContext.BinaryClassification.Evaluate(predictions);

            if (modelPath != null)
            {
                mlContext.Model.Save(model, data.Schema, modelPath);
            }
            return model;
        }

        public static IEnumerable<PFPairFeature> RandomSample(IEnumerable<PFPairFeature> pfPairs, int targetCount)
        {
            if (pfPairs.Count() <= targetCount)
            {
                return pfPairs;
            }
            var random = new Random(42);
            var sampledPairs = pfPairs.OrderBy(x => random.Next()).Take(targetCount).ToList();
            return sampledPairs;
        }

        public static IEnumerable<PFPairFeature> UnderSample(IEnumerable<PFPairFeature> allPairFeatures, int targetCount = 0)
        {
            var positives = allPairFeatures.Where(p => p.Label == true);
            var negatives = allPairFeatures.Where(p => p.Label == false);
            int posNum = positives.Count();
            int negNum = negatives.Count();
            if (positives.Count() < negatives.Count())
            {
                if (targetCount == 0) 
                    negatives = RandomSample(negatives, positives.Count());
                else
                    negatives = RandomSample(negatives, targetCount);
            }
            else
            {
                if (targetCount == 0)
                    positives = RandomSample(positives, negatives.Count());
                else
                    positives = RandomSample(positives, targetCount);
            }
            var filteredPairs = positives.Concat(negatives).ToList();
            return filteredPairs;
        }

        public static PFPairFeature GetPFPairFeatureFromPfPairMetrics(PFpairMetrics pfPairMetrics)
        {
            var pfPair = new PFPairFeature
            {
                Correlation = (float)pfPairMetrics.Correlation,
                ApexRtDelta = (float)pfPairMetrics.ApexRtDelta,
                Overlap = (float)pfPairMetrics.Overlap,
                FragmentIntensity = (float)pfPairMetrics.FragmentIntensity,
                PsmScore = (float)pfPairMetrics.PsmScore,
                SharedXIC = (float)pfPairMetrics.SharedXIC,
                Label = pfPairMetrics.MatchedIonType == "NA" ? false : true
            };
            return pfPair;
        }


        public static List<PrecursorFragmentsGroup> GetAllPfGroups_ML(List<PeakCurve> ms1curves, List<PeakCurve> ms2curves, DIAparameters diaParam)
        {
            var allGroups = new List<PrecursorFragmentsGroup>();
            var pseudoGroups = new List<PrecursorFragmentsGroup>();
            foreach (var precursor in ms1curves)
            {
                var pseudoGroup = PrecursorFragmentsGroup.JustPair(precursor, ms2curves, diaParam);
                pseudoGroups.Add(pseudoGroup);
            }

            return null;
        }

        public static List<PFPairFeature> GetPFPairsFromPFMetricsExcludingInternal(IEnumerable<PFpairMetrics> pfPairMetricsList)
        {
            var list = new List<PFPairFeature>();
            foreach (var pfPairMetrics in pfPairMetricsList)
            {
                if (pfPairMetrics.MatchedIonType == "Internal" || pfPairMetrics.MatchedIonType == "NeutralLoss") continue;
                var pfPair = GetPFPairFeatureFromPfPairMetrics(pfPairMetrics);
                list.Add(pfPair);
            }
            return list;
        }

        public static (float, float) ComputeWeights(List<PFPairFeature> pfPairs)
        {
            float w0 = pfPairs.Count / (2 * pfPairs.Where(p => p.Label == false).Count());
            float w1 = pfPairs.Count / (2 * pfPairs.Where(p => p.Label == true).Count());
            return (w0, w1);
        }

        public static void AssignWeightsToPFPairs(IEnumerable<PFPairFeature> pfPairs, (float, float) weights)
        {
            foreach (var pfPair in pfPairs)
            {
                pfPair.Weight = pfPair.Label == true ? weights.Item2 : weights.Item1;
            }
        }

        //this method assumes that Xic construction is already done
        public static PrecursorFragmentsGroup FindPfGroupBasedOnPsm(PsmFromTsv psmTsv, List<DeconvolutedMass>[] ms1PeakTable, Dictionary<double, List<DeconvolutedMass>[]> allMs2MassTables,
            DIAparameters diaParam)
        {
            var cycle = (psmTsv.Ms2ScanNumber - 1) % 4;
            List<DeconvolutedMass>[] ms2PeakTable = null;
            switch (cycle)
            {
                case 1:
                    ms2PeakTable = allMs2MassTables[60];
                    break;
                case 2:
                    ms2PeakTable = allMs2MassTables[80];
                    break;
                case 3:
                    ms2PeakTable = allMs2MassTables[100];
                    break;
            }
            var precursorMass = MassCurve.GetMassFromScan(psmTsv.PrecursorMass, psmTsv.PrecursorCharge, ms1PeakTable, (psmTsv.PrecursorScanNum - 1) / 4,
                diaParam.Ms1PeakFindingTolerance, diaParam.PeakSearchBinSize);
            if (precursorMass == null || precursorMass.PeakCurve == null)
            {
                return null;
            }
            var pfGroup = new PrecursorFragmentsGroup(precursorMass.PeakCurve);
            precursorMass.PeakCurve.GetRawXYData();

            foreach (var matchedIon in psmTsv.MatchedIons)
            {
                var fragmentMass = MassCurve.GetMassFromScan(matchedIon.Mz.ToMass(matchedIon.Charge), matchedIon.Charge, ms2PeakTable, (psmTsv.Ms2ScanNumber - 1) / 4,
                                       diaParam.Ms2PeakFindingTolerance, diaParam.PeakSearchBinSize);
                if (fragmentMass.PeakCurve != null)//fragmentMass == null ||
                {
                    double overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    fragmentMass.PeakCurve.GetRawXYData();
                    double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    var sharedXIC = PrecursorFragmentPair.CalculateSharedXIC(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    var pfPair = new PrecursorFragmentPair(precursorMass.PeakCurve, fragmentMass.PeakCurve, overlap, corr, sharedXIC, psmTsv.DecoyContamTarget);
                    pfGroup.PFpairs.Add(pfPair);
                }
            }
            return pfGroup;
        }



        public static List<PFpairMetrics> GetPFPairMetricsFromPsm(PsmFromTsv psmTsv, List<Peak>[] ms1PeaksByScan, List<Peak>[] ms2PeaksByScan, DIAparameters diaParam)
        {
            var pfPairMetricsList = new List<PFpairMetrics>();
            var precursorMass = SearchPeakByScan(psmTsv.PrecursorMass, psmTsv.PrecursorCharge, ms1PeaksByScan[psmTsv.PrecursorScanNum], diaParam.Ms1PeakFindingTolerance);
            if (precursorMass == null || precursorMass.PeakCurve == null)
            {
                return null;
            }
            var pfGroup = new PrecursorFragmentsGroup(precursorMass.PeakCurve);
            ISDEngine_static.PeakCurveSpline(precursorMass.PeakCurve, diaParam.Ms1SplineType, diaParam);

            var visited = new HashSet<Peak>();
            foreach (var matchedIon in psmTsv.MatchedIons.Where(i => i.IsInternalFragment == false))
            {
                var fragmentMass = SearchPeakByScan(matchedIon.Mz.ToMass(matchedIon.Charge), matchedIon.Charge, ms2PeaksByScan[psmTsv.Ms2ScanNumber],
                                       diaParam.Ms2PeakFindingTolerance);
                if (fragmentMass != null && fragmentMass.PeakCurve != null)
                {
                    double overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    fragmentMass.PeakCurve.GetRawXYData();
                    double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    double sharedXIC = PrecursorFragmentPair.CalculateSharedXIC(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    var pfPair = new PrecursorFragmentPair(precursorMass.PeakCurve, fragmentMass.PeakCurve, overlap, corr, sharedXIC, psmTsv.DecoyContamTarget);
                    var pfPairMetrics = new PFpairMetrics(pfPair);
                    pfPairMetrics.MatchedIonType = "Terminal";
                    visited.Add(fragmentMass);
                    pfPairMetricsList.Add(pfPairMetrics);
                }
            }

            var umatchedMasses = ms2PeaksByScan[psmTsv.Ms2ScanNumber].Where(p => !visited.Contains(p)).ToList();
            foreach(var mass in umatchedMasses)
            {
                double overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursorMass.PeakCurve, mass.PeakCurve);
                mass.PeakCurve.GetRawXYData();
                double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorMass.PeakCurve, mass.PeakCurve);
                double sharedXIC = PrecursorFragmentPair.CalculateSharedXIC(precursorMass.PeakCurve, mass.PeakCurve);
                var pfPair = new PrecursorFragmentPair(precursorMass.PeakCurve, mass.PeakCurve, overlap, corr, sharedXIC, psmTsv.DecoyContamTarget);
                var pfPairMetrics = new PFpairMetrics(pfPair);
                pfPairMetrics.MatchedIonType = "NA";
                pfPairMetricsList.Add(pfPairMetrics);
            }
            return pfPairMetricsList;
        }

        public static Peak SearchPeakByScan(double targetMass, int targetCharge, List<Peak> peaks, Tolerance tolerance)
        {
            var allPossiblePeaks = peaks.Where(p => p.Charge == targetCharge).ToArray();
            if (allPossiblePeaks == null || allPossiblePeaks.Length == 0)
            {
                return null; // No peaks found for the specified charge
            }
            int index = Array.BinarySearch(allPossiblePeaks.Select(p => p.MonoisotopicMass).ToArray(), targetMass);
            if (index < 0)
            {
                index = ~index; 
            }
            if (index == 0)
            {
                if (tolerance.Within(targetMass, allPossiblePeaks[0].MonoisotopicMass))
                    return allPossiblePeaks[0];
                else
                    return null;
            }
            if (index == allPossiblePeaks.Length)
            {
                if (tolerance.Within(targetMass, allPossiblePeaks[allPossiblePeaks.Length - 1].MonoisotopicMass))
                    return allPossiblePeaks[allPossiblePeaks.Length - 1];
                else
                    return null;
            }

            Peak bestMass = null;
            if (tolerance.Within(targetMass, allPossiblePeaks[index - 1].MonoisotopicMass))
            {
                bestMass = allPossiblePeaks[index - 1];
            }
            if (tolerance.Within(targetMass, allPossiblePeaks[index].MonoisotopicMass)
                && (bestMass == null || Math.Abs(bestMass.MonoisotopicMass - targetMass) > Math.Abs(allPossiblePeaks[index].MonoisotopicMass - targetMass)))
            {
                bestMass = allPossiblePeaks[index];
            }
            return bestMass;
        }
    }
    
}
