using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DIAparameters
    {
        public AnalysisType AnalysisType { get; set; }
        public XicConstructor Ms1XicConstructor { get; set; }
        public XicConstructor Ms2XicConstructor { get; set; } 
        public PfGroupingEngine PfGroupingEngine { get; set; }
        public PseudoMs2ConstructionType PseudoMs2ConstructionType { get; set; }
        public bool CombineFragments { get; set; } 
        public Dictionary<string, Ms2ScanWithSpecificMass[]> PseudoScans { get; set; }
        public bool WritePseudoScans { get; set; } 
        public string DbPath { get; set; }
        public string PsmPath { get; set; }

        public DIAparameters(AnalysisType analysisType, XicConstructor ms1XicConstructor, XicConstructor ms2XicConstructor, PfGroupingEngine pfGroupingEngine, PseudoMs2ConstructionType pseudoMs2ConstructionType, bool combineFragments = false, bool writePseudoScans = false)
        {
            AnalysisType = analysisType;
            Ms1XicConstructor = ms1XicConstructor;
            Ms2XicConstructor = ms2XicConstructor;
            PfGroupingEngine = pfGroupingEngine;
            PseudoMs2ConstructionType = pseudoMs2ConstructionType;
            CombineFragments = combineFragments;
            WritePseudoScans = writePseudoScans;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("DIAparameters:");
            sb.AppendLine($"AnalysisType: {AnalysisType}");
            sb.AppendLine($"{Ms1XicConstructor.ToString()}");
            sb.AppendLine($"{Ms2XicConstructor.ToString()}");
            sb.AppendLine($"{Ms1XicConstructor.XicSplineEngine?.ToString()}");
            sb.AppendLine($"{Ms2XicConstructor.XicSplineEngine?.ToString()}");
            if (PfGroupingEngine != null) sb.AppendLine($"{PfGroupingEngine.ToString()}");
            sb.AppendLine($"PseudoMs2ConstructionType: {PseudoMs2ConstructionType}");
            sb.AppendLine($"CombineFragments: {CombineFragments}");
            return sb.ToString();
        }
    }
}
