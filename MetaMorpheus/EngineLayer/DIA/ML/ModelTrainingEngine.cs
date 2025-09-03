using Microsoft.ML;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using static System.Runtime.InteropServices.JavaScript.JSType;
using MassSpectrometry;

namespace EngineLayer.DIA
{
    public class ModelTrainingEngine
    {
        public int TrainingDataSize { get; set; } 
        public MLbasedDIAparameters MlDIAparams { get; set; }
        public SpectralMatch[] Psms { get; set; }
        public PrecursorFragmentsGroup[] PfGroups { get; set; }
        public Ms2ScanWithSpecificMass[] Ms2Scans { get; set; }

        public ModelTrainingEngine(MLbasedDIAparameters mlDIAparams, SpectralMatch[] psms = null, PrecursorFragmentsGroup[] pfGroups = null, Ms2ScanWithSpecificMass[] ms2Scans = null, int trainingDataSize = 0)
        {
            MlDIAparams = mlDIAparams;
            Psms = psms;
            PfGroups = pfGroups;
            Ms2Scans = ms2Scans;
            TrainingDataSize = trainingDataSize;
        }

        public ITransformer TrainModel()
        {
            var mlContext = new MLContext();
            ITransformer model = null;
            if (MlDIAparams.ExistingModelPath != null)
            {
                model = mlContext.Model.Load(MlDIAparams.ExistingModelPath, out var modelInputSchema);
                return model;
            }

            List<PfPairTrainingSample> trainingSamples = null;
            if (MlDIAparams.ExistingSampleFilePath == null)
            {
                trainingSamples = GetTrainingSamples();
                string newFilePath = @"E:\ISD Project\TestDataForML\SampleFileFolder\umpireXic";
                var sampleFile = new PfPairTrainingSampleFile(newFilePath);
                //sampleFile.WriteResults(newFilePath);
            }
            else
            {
                var sampleFile = new PfPairTrainingSampleFile(MlDIAparams.ExistingSampleFilePath);
                sampleFile.LoadResults();
                trainingSamples = sampleFile.Results.ToList();
            }
            IDataView data = mlContext.Data.LoadFromEnumerable(trainingSamples);
            var split = mlContext.Data.TrainTestSplit(data, testFraction: 0.2);
            var trainData = split.TrainSet;
            var testData = split.TestSet;

            switch (MlDIAparams.ModelType)
            {
                case ModelType.LogisticRegression:
                    var pipeline = mlContext.Transforms.Concatenate("Features", MlDIAparams.Features.ToArray()).Append(mlContext.BinaryClassification.Trainers.SdcaLogisticRegression(labelColumnName: "Label", featureColumnName: "Features"));
                    model = pipeline.Fit(trainData);
                    break;
                case ModelType.FastTree:
                    var pipeline2 = mlContext.Transforms.Concatenate("Features", MlDIAparams.Features.ToArray()).Append(mlContext.BinaryClassification.Trainers.FastTree(
                       labelColumnName: "Label",
                       featureColumnName: "Features",
                       numberOfTrees: 100,  // Example: 100 trees
                       learningRate: 0.2f // Example: learning rate of 0.2
                   ));
                    model = pipeline2.Fit(trainData);
                    break;
            }

            //validation
            var predictions = model.Transform(testData);
            var metrics = mlContext.BinaryClassification.Evaluate(predictions);

            return model;
        }
        
        public List<PfPairTrainingSample> GetTrainingSamples()
        {
            var filteredPsms = Psms.Where(p => p.Score >= MlDIAparams.PsmScoreCutOff).ToArray();
            var peptides = filteredPsms.GroupBy(p => p.FullSequence).ToList();
            var allSamples = new List<PfPairTrainingSample>();
            foreach(var psm in filteredPsms)
            {
                var correspondingPfGroup = PfGroups.FirstOrDefault(g => g.PFgroupIndex == psm.ScanNumber);
                var correspondingMs2Scan = Ms2Scans.FirstOrDefault(s => s.OneBasedScanNumber == psm.ScanNumber);
                var samplesFromThisPsm = GetTrainingSamplesFromPfGroup(psm, correspondingPfGroup, correspondingMs2Scan, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                allSamples.AddRange(samplesFromThisPsm);
            }
            var filteredSamples = BalanceTrainingData(allSamples, TrainingDataSize);
            return filteredSamples.ToList();
        }

        public static List<PfPairTrainingSample> GetTrainingSamplesFromPfGroup(SpectralMatch psm, PrecursorFragmentsGroup pfGroup, Ms2ScanWithSpecificMass ms2WithPrecursor, Tolerance tol)
        {
            var samples = new List<PfPairTrainingSample>();
            var sortedPfPairs = pfGroup.PFpairs.OrderByDescending(pf => pf.FragmentXic.AveragedMassOrMz).ToArray();
            var sortedPfPairMzs = sortedPfPairs.Select(pf => (double)pf.FragmentXic.ApexPeak.M).ToArray();

            var positiveIndices = new List<int>();
            foreach (var ion in psm.MatchedFragmentIons)
            {
                double minMass = tol.GetMinimumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                double maxMass = tol.GetMaximumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                var masses = ms2WithPrecursor.GetClosestExperimentalIsotopicEnvelopeList(minMass, maxMass);
                var targetMass = masses.FirstOrDefault(m => m.Charge == ion.Charge);
                if (targetMass != null)
                {
                    foreach (var peak in targetMass.Peaks)
                    {
                        var closestIndex = FindClosestIndexOfPfPair(sortedPfPairMzs, peak.mz);
                        positiveIndices.Add(closestIndex);
                    }
                }
            }

            for (int i = 0; i < sortedPfPairs.Length; i++)
            {
                var newSample = new PfPairTrainingSample(sortedPfPairs[i]);
                if (positiveIndices.Contains(i))
                {
                    newSample.Label = true;
                }
                else
                {
                    newSample.Label = false;
                }
                samples.Add(newSample);
            }
            return samples;
        }

        public static List<PfPairTrainingSample> GetTrainingSamplesFromDirectSearch(SpectralMatch psm, Ms2ScanWithSpecificMass ms2WithPrecursor, Dictionary<IIndexedPeak, ExtractedIonChromatogram> ms1Xics, Dictionary<IIndexedPeak, ExtractedIonChromatogram> ms2Xics)
        {
            return null;
        }

        public static IEnumerable<PfPairTrainingSample> BalanceTrainingData(IEnumerable<PfPairTrainingSample> allPairFeatures, int targetCount = 0)
        {
            var positives = allPairFeatures.Where(p => p.Label == true);
            var negatives = allPairFeatures.Where(p => p.Label == false);
            //debug
            int posCount = positives.Count();
            int negCount = negatives.Count();

            if (targetCount != 0)
            {
                positives = RandomSample(positives, targetCount);
                negatives = RandomSample(negatives, targetCount);
            }
            else
            {
                if (positives.Count() == negatives.Count()) return positives.Concat(negatives).ToList();
                if (positives.Count() < negatives.Count())
                {
                    negatives = RandomSample(negatives, positives.Count());
                }
                else
                {
                    positives = RandomSample(positives, negatives.Count());
                }
            }
            return positives.Concat(negatives).ToList();
        }

        public static IEnumerable<PfPairTrainingSample> RandomSample(IEnumerable<PfPairTrainingSample> pfPairs, int targetCount)
        {
            if (pfPairs.Count() <= targetCount)
            {
                return pfPairs;
            }
            var random = new Random();
            var sampledPairs = pfPairs.OrderBy(x => random.Next()).Take(targetCount).ToList();
            return sampledPairs;
        }
        private static int FindClosestIndexOfPfPair(double[] sortedPfPairMzs, double targetMz)
        {
            int closestIndex = -1;
            var index = Array.BinarySearch(sortedPfPairMzs, targetMz);
            if (index >= 0)
            {
                return index;
            }
            else
            {
                index = ~index;
                if (index == 0)
                {
                    closestIndex = 0;
                }
                else if (index >= sortedPfPairMzs.Length)
                {
                    closestIndex = sortedPfPairMzs.Length - 1;
                }
                else
                {
                    if (Math.Abs(sortedPfPairMzs[index] - targetMz) < Math.Abs(sortedPfPairMzs[index - 1] - targetMz))
                    {
                        closestIndex = index;
                    }
                    else
                    {
                        closestIndex = index - 1;
                    }
                }
            }
            return closestIndex;
        }
    }
}
