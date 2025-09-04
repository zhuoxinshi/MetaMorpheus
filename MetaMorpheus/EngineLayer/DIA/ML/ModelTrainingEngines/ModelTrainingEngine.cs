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
using FlashLFQ;
using System.Net.Http.Headers;
using EngineLayer.ClassicSearch;
using UsefulProteomicsDatabases;
using Omics.Modifications;
using Proteomics;

namespace EngineLayer.DIA
{
    public abstract class ModelTrainingEngine
    {
        public MLbasedDIAparameters MlDIAparams { get; set; }
        public CommonParameters CommonParameters { get; set; } 
        public SpectralMatch[] Psms { get; set; }
        public Ms2ScanWithSpecificMass[] PseudoSearchMs2Scans { get; set; }

        public ModelTrainingEngine(MLbasedDIAparameters mlDIAparams, CommonParameters commonParameters)
        {
            MlDIAparams = mlDIAparams;
            CommonParameters = commonParameters;
        }

        public abstract IEnumerable<PfPairTrainingSample> GetTrainingSamples();

        public abstract void GeneratePseudoSearchScans();

        public ITransformer TrainModel()
        {
            IEnumerable<PfPairTrainingSample> trainingSamples = null;  
            if (MlDIAparams.ExistingSampleFilePath != null)
            {
                var sampleFile = new PfPairTrainingSampleFile(MlDIAparams.ExistingSampleFilePath);
                sampleFile.LoadResults();
                trainingSamples = sampleFile.Results.ToList();
            }
            else
            {
                trainingSamples = GetTrainingSamples();
            }

            var model = ModelTraining(trainingSamples);
            return model;
        }
        
        public ITransformer ModelTraining(IEnumerable<PfPairTrainingSample> trainingSamples)
        {
            var mlContext = new MLContext();

            //balance training data
            var balancedSamples = BalanceTrainingData(trainingSamples, MlDIAparams.TargetSampleCount);

            IDataView data = mlContext.Data.LoadFromEnumerable(balancedSamples);
            var split = mlContext.Data.TrainTestSplit(data, testFraction: 0.2);
            var trainData = split.TrainSet;
            var testData = split.TestSet;

            ITransformer model = null;
            switch (MlDIAparams.ModelType)
            {
                case ModelType.LogisticRegression:
                    var pipeline = mlContext.Transforms.Concatenate("Features", MlDIAparams.Features.ToArray()).Append(mlContext.BinaryClassification.Trainers.SdcaLogisticRegression(labelColumnName: "Label", featureColumnName: "Features"));
                    model = pipeline.Fit(trainData);
                    break;
                case ModelType.FastTree:
                    var pipeline2 = mlContext.Transforms.Concatenate("Features", MlDIAparams.Features.ToArray()).Append(mlContext.BinaryClassification.Trainers.FastTree(
                                              labelColumnName: "Label",
                                                                    featureColumnName: "Features"));
                    model = pipeline2.Fit(trainData);
                    break;
                default:
                    throw new NotImplementedException($"Model type {MlDIAparams.ModelType} not implemented.");
            }
            return model;
        }

        public SpectralMatch[] RunClassicSearch(Ms2ScanWithSpecificMass[] pseudoSearchScans)
        {
            SpectralMatch[] fileSpecificPsms = new SpectralMatch[pseudoSearchScans.Length];
            var fixedModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var variableModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            var decoyType = MlDIAparams.SearchDecoys ? DecoyType.Reverse : DecoyType.None;
            var proteins = LoadProteinDb(MlDIAparams.ProteinDb, true, decoyType, null, false, out var um, out int count, CommonParameters);
            var newClassicSearchEngine = new ClassicSearchEngine(fileSpecificPsms, pseudoSearchScans, variableModifications, fixedModifications, null, null, null, proteins, new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, CommonParameters.PrecursorMassTolerance), CommonParameters, null, null, null, true);
            var result = newClassicSearchEngine.Run();
            return fileSpecificPsms.Where(p => p != null).ToArray();
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

        protected static List<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, out Dictionary<string, Modification> um,
            out int emptyEntriesCount, CommonParameters commonParameters)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                um = null;
                proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, out dbErrors,
                    ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, commonParameters.MaxThreadsToUsePerFile, addTruncations: commonParameters.AddTruncations);
            }
            else
            {
                List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
                proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, commonParameters.MaxHeterozygousVariants, commonParameters.MinVariantDepth, addTruncations: commonParameters.AddTruncations);
            }

            emptyEntriesCount = proteinList.Count(p => p.BaseSequence.Length == 0);
            return proteinList.Where(p => p.BaseSequence.Length > 0).ToList();
        }
    }
}
