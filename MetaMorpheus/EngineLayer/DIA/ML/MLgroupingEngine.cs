using BayesianEstimation;
using MassSpectrometry;
using Microsoft.ML;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class MLgroupingEngine : PfGroupingEngine
    {
        public ITransformer Model { get; set; }
        public double ApexRTTolerance { get; set; } 
        public double ProbabilityThreshold { get; set; }

        public MLgroupingEngine(ITransformer model, double scoreThreshold, double apexRtTolerance)
        {
            Model = model;
            ProbabilityThreshold = scoreThreshold;
            ApexRTTolerance = apexRtTolerance;
        }

        public override List<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, IEnumerable<ExtractedIonChromatogram> fragments)
        {
            var pfGroups = new ConcurrentBag<PrecursorFragmentsGroup>();
            var apexSortedFragmentXics = XicGroupingEngine.BuildApexSortedXics(fragments);

            var mlContext = new MLContext();
            var predictionEngine = mlContext.Model.CreatePredictionEngine<PfPairTrainingSample, PFpairPrediction>(Model);
            Parallel.ForEach(Partitioner.Create(0, precursors.Count), new ParallelOptions { MaxDegreeOfParallelism = 15 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursors[i];
                        var fragmentsInRange = XicGroupingEngine.GetXicsInRange(apexSortedFragmentXics, precursor.ApexRT, ApexRTTolerance);
                        var group = MLgroupingForOnePFgroup(precursor, fragmentsInRange, predictionEngine);
                        if (group != null)
                            pfGroups.Add(group);
                    }
                });
            return pfGroups.ToList();
        }

        public PrecursorFragmentsGroup MLgroupingForOnePFgroup(ExtractedIonChromatogram precursor, IEnumerable<ExtractedIonChromatogram> fragments, PredictionEngine<PfPairTrainingSample, PFpairPrediction> predictionEngine)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragment in fragments)
            {
                var pfPairSample = new PfPairTrainingSample(precursor, fragment);
                var prediction = predictionEngine.Predict(pfPairSample);
                if (prediction.Probability >= ProbabilityThreshold)
                {
                    var pfPair = new PrecursorFragmentPair(precursor, fragment);
                    pfPairs.Add(pfPair);
                }
            }
            if (pfPairs.Count > 0)
            {
                var pfGroup = new PrecursorFragmentsGroup(precursor, pfPairs);
                return pfGroup;
            }
            return null;
        }
    }
}
