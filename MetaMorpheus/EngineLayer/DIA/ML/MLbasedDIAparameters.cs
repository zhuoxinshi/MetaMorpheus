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
        public double PredictionScoreThreshold { get; set; } 
        public int TargetSampleCount { get; set; } 
        public MLbasedDIAparameters(PseudoSearchScanType pseudoSearchScanType, string proteinDb, bool searchDecoys, ModelType trainingModel, List<string> features, double psmScoreCutOff, string existingModelPath, string existingSampleFilePath, double testFraction, double predictionScoreThreshold, AnalysisType analysisType, XicConstructor ms1XicConstructor, XicConstructor ms2XicConstructor, PfGroupingEngine pfGroupingEngine, PseudoMs2ConstructionType pseudoMs2ConstructionType, bool combineFragments = false) : base(analysisType, ms1XicConstructor, ms2XicConstructor, pfGroupingEngine, pseudoMs2ConstructionType, combineFragments)
        {
            PseudoSearchType = pseudoSearchScanType;
            ProteinDb = proteinDb;
            SearchDecoys = searchDecoys;
            ModelType = trainingModel;
            Features = features;
            PsmScoreCutOff = psmScoreCutOff;
            ExistingModelPath = existingModelPath;
            ExistingSampleFilePath = existingSampleFilePath;
            TestFraction = testFraction;
            PredictionScoreThreshold = predictionScoreThreshold;
        }
    }
}
