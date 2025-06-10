using Chemistry;
using Microsoft.ML;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer.DIA;
using System.Reflection.Emit;

namespace EngineLayer.DIA
{
    public class PFPair
    {
        public float Correlation { get; set; }   // e.g., XIC correlation
        public float ApexRtDelta { get; set; }   // e.g., RT diff
        public float Overlap { get; set; }      // e.g., XIC overlap
        public float FragmentIntensity { get; set; } // e.g., fragment intensity ratio
        public bool Label { get; set; }       // true = target, false = decoy
    }

    public class PredictionResult
    {
        public bool PredictedLabel { get; set; }
        public float Probability { get; set; }
        public float Score { get; set; }
    }

    public class TrainModel
    {
        public static void Train(IEnumerable<PFPair> yourDataList)
        {
            var mlContext = new MLContext();

            IDataView data = mlContext.Data.LoadFromEnumerable(yourDataList);

            // Split into train/test
            var split = mlContext.Data.TrainTestSplit(data, testFraction: 0.2);
            var trainData = split.TrainSet;
            var testData = split.TestSet;

            var pipeline = mlContext.Transforms.Concatenate("Features", nameof(PFPair.Correlation), nameof(PFPair.ApexRtDelta),
                nameof(PFPair.Overlap), nameof(PFPair.FragmentIntensity))
                .Append(mlContext.BinaryClassification.Trainers.SdcaLogisticRegression(labelColumnName: "Label", featureColumnName: "Features"));
            var model = pipeline.Fit(trainData);
            var predictions = model.Transform(testData);
            var metrics = mlContext.BinaryClassification.Evaluate(predictions);

            var predictor = mlContext.Model.CreatePredictionEngine<PrecursorFragmentPair, PredictionResult>(model);
            Console.WriteLine($"AUC: {metrics.AreaUnderRocCurve:F4}");
            Console.WriteLine($"F1 Score: {metrics.F1Score:F4}");
            Console.WriteLine($"Accuracy: {metrics.Accuracy:P2}");
        }

        public static PFPair GetPFPairsFromPfPairMetrics(PFpairMetrics pfPairMetrics)
        {
            var pfPair = new PFPair
            {
                Correlation = (float)pfPairMetrics.Correlation,
                ApexRtDelta = (float)pfPairMetrics.ApexRtDelta,
                Overlap = (float)pfPairMetrics.Overlap,
                FragmentIntensity = (float)pfPairMetrics.FragmentIntensity,
                Label = pfPairMetrics.MatchedIonType == "Terminal" ? true : false
            };
            return pfPair;
        }

        public static List<PFPair> GetPFPairsFromPFMetricsFile(PFpairMetricFile pfMetricsFile)
        {
            var list = new List<PFPair>();
            foreach (var pfPairMetrics in pfMetricsFile.Results)
            {
                if (pfPairMetrics.MatchedIonType == "Internal") continue;
                var pfPair = GetPFPairsFromPfPairMetrics(pfPairMetrics);
                list.Add(pfPair);
            }
            return list;
        }
    }
    
}
