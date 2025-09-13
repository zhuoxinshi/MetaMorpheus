using pepXML.Generated;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class MLbasedDIAparameters : DIAparameters
    {
        public PseudoSearchScanType PseudoSearchType { get; set; }
        public string ProteinDb { get; set; }
        public bool SearchDecoys { get; set; } 
        public ModelType ModelType { get; set; }
        public List<string> Features { get; set; }
        public double PsmScoreCutOff { get; set; } 
        public string ExistingModelPath { get; set; } 
        public string ExistingSampleFilePath { get; set; }
        public string OutputFolder { get; set; }
        public double TestFraction { get; set; }
        public double ApexRtTolerance { get; set; } 
        public double PredictionScoreThreshold { get; set; } 
        public int TargetSampleCount { get; set; } 
        public bool WriteModel { get; set; } 
        public bool WriteTrainingSamples { get; set; }
        public MLbasedDIAparameters(PseudoSearchScanType pseudoSearchScanType, string proteinDb, bool searchDecoys, ModelType trainingModel, List<string> features, double psmScoreCutOff, string existingModelPath, string existingSampleFilePath, string outFolder, double testFraction, double apexRtTolerance, double predictionScoreThreshold, AnalysisType analysisType, XicConstructor ms1XicConstructor, XicConstructor ms2XicConstructor, PfGroupingEngine pfGroupingEngine, PseudoMs2ConstructionType pseudoMs2ConstructionType, bool combineFragments = false, bool writeModel = false, bool writeTrainingSamples = false) : base(analysisType, ms1XicConstructor, ms2XicConstructor, pfGroupingEngine, pseudoMs2ConstructionType, combineFragments)
        {
            PseudoSearchType = pseudoSearchScanType;
            ProteinDb = proteinDb;
            SearchDecoys = searchDecoys;
            ModelType = trainingModel;
            Features = features;
            PsmScoreCutOff = psmScoreCutOff;
            ExistingModelPath = existingModelPath;
            ExistingSampleFilePath = existingSampleFilePath;
            ApexRtTolerance = apexRtTolerance;
            TestFraction = testFraction;
            PredictionScoreThreshold = predictionScoreThreshold;
            OutputFolder = outFolder;
            WriteModel = writeModel;
            WriteTrainingSamples = writeTrainingSamples;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("DIAparameters:");
            sb.AppendLine($"AnalysisType: {AnalysisType}");
            sb.AppendLine($"{Ms1XicConstructor.ToString()}");
            sb.AppendLine($"{Ms2XicConstructor.ToString()}");
            sb.AppendLine($"PseudoSearchType: {PseudoSearchType}");
            sb.AppendLine($"SearchDecoys: {SearchDecoys}");
            sb.AppendLine($"ModelType: {ModelType}");
            sb.AppendLine($"Features: {string.Join(", ", Features)}");
            sb.AppendLine($"PsmScoreCutOff: {PsmScoreCutOff}");
            if (!string.IsNullOrEmpty(ExistingModelPath)) sb.AppendLine($"ExistingModelPath: {ExistingModelPath}");
            if (!string.IsNullOrEmpty(ExistingSampleFilePath)) sb.AppendLine($"ExistingSampleFilePath: {ExistingSampleFilePath}");
            sb.AppendLine($"ApexRtTolerance: {ApexRtTolerance}");
            sb.AppendLine($"TestFraction: {TestFraction}");
            sb.AppendLine($"PredictionScoreThreshold: {PredictionScoreThreshold}");
            sb.AppendLine($"PseudoMs2ConstructionType: {PseudoMs2ConstructionType}");
            sb.AppendLine($"CombineFragments: {CombineFragments}");
            return sb.ToString();
        }
    }
}
