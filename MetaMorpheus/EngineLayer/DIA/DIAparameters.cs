using Easy.Common.Extensions;
using EngineLayer.DIA;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DIAparameters
    {
        public Tolerance Ms1PeakFindingTolerance { get; set; }
        public Tolerance Ms2PeakFindingTolerance { get; set; }
        public int MaxNumMissedScan {  get; set; }
        public int PeakSearchBinSize {  get; set; }
        public double OverlapRatioCutOff {  get; set; }
        public double CorrelationCutOff {  get; set; }
        public double ApexRtTolerance { get; set; }
        public int ApexCycleTolerance { get; set; }
        public int FragmentRankCutOff { get; set; }
        public int PrecursorRankCutOff { get; set; }
        public double HighCorrThreshold { get; set; }
        public int NumHighCorrFragments { get; set; }
        public double MaxRTRangeMS1 { get; set; }
        public double MaxRTRangeMS2 { get; set; }
        public double PrecursorIntensityCutOff { get; set; }
        public double MinRTRangeForCWT { get; set; }
        public CwtParameters CwtParameters => new CwtParameters(2f, 150, 0.1f, 0.3f);
        public bool SplitMS2Peak { get; set; } 
        public bool SplitMS1Peak { get; set; }  
        public float SplineTimeInterval { get; set; }
        public double SplineRtInterval { get; set; }
        public double ScanCycleSplineTimeInterval { get; set; }
        public double MinMS1Mass { get; set; }
        public double MinMS2Mass { get; set; }
        public double MaxMass { get; set; }
        public int MinMS1Charge { get; set; }
        public int MinMS2Charge { get; set; }
        public string Type { get; set; }
        public Tolerance PrecursorMassTolerance => new PpmTolerance(200);
        public bool CutMs1Peaks { get; set; }
        public bool CutMs2Peaks { get; set; }
        public bool AverageMs2Scans { get; set; }
        public XICType Ms1XICType {  get; set; } 
        public XICType Ms2XICType { get; set; }
        public PFGroupingType PFGroupingType { get; set; }
        public PseudoMs2ConstructionType PseudoMs2ConstructionType { get; set; }
        public AnalysisType AnalysisType { get; set; }
        public bool CombineFragments { get; set; }
        public CorrelationType CorrelationType { get; set; }
        public int NumScansPerCycle { get; set; }
        public int SGfilterWindowSize { get; set; }
        public SplineType Ms1SplineType { get; set; }
        public SplineType Ms2SplineType { get; set; }
        public bool TrimMs2Peaks { get; set; }
        public double TrimMs2MinSNR { get; set; }

        public int NoPointsPerMin { get; set; }
        public int NumPeaksThreshold { get; set; }

        public Dictionary<string, List<PrecursorFragmentsGroup>> PFgroupsDictionary { get; set; }
        public Dictionary<string, List<PeakCurve>> PeakCurveDictionary { get; set; }

        public DIAparameters(Tolerance ms1PeakFindingTolerance, Tolerance ms2PeakFindingTolerance, int maxNumMissedScan = 2, int binSize = 100, 
            double overlapRatioCutOff = 0.3, double correlationCutOff = 0.5, double apexRtTolerance = 0.1, int fragmentRankCutOff = 5000, int precursorRankCutOff = 1000
            , double maxRTrangeMS1 = 0.5, double maxRTrangeMS2 = 2, double highCorrThreshold = 0.5, int numHighCorrFragments = 0, double precursorIntensityCutOff = 10000, double minRTRangeForCWT = 0.1,
            bool splitMS2Peak = false, bool splitMS1Peak = false, double splineRtInterval = 0.005, double minMS1Mass = 0, double minMS2Mass = 0, double maxMass = 99999, string type = "DIA", int apexCycleTolerance = 2, 
            double scanCycleSplineInterval = 0.025, int minMS1Charge = 1, int minMS2Charge = 1, bool averageMs2Scans = false, XICType ms1XICType = XICType.DeconHighestPeak, 
            XICType ms2XICType = XICType.Peak, PFGroupingType pfGroupingType = PFGroupingType.ScanCycle, PseudoMs2ConstructionType pseudoMs2Type = PseudoMs2ConstructionType.mzPeak, 
            AnalysisType analysisType = AnalysisType.DIAEngine, bool combineFragments = false, CorrelationType correlationType = CorrelationType.CubicSpline_scanCycle,
            bool cutMs1Peaks = false, bool cutMs2Peaks = false, int sgFilterWindowSize = 5, SplineType ms1SplineType = SplineType.NoSpline, SplineType ms2SplineType = SplineType.NoSpline, 
            float splineTimeInterval = 0.05f, int numScanPerCycle = 0, bool trimMs2Peaks = false, double trimMs2MinSNR = 0.01, int noPointsPerMin = 150, int numPeaksThreshold = 4)
        {
            Ms1PeakFindingTolerance = ms1PeakFindingTolerance;
            Ms2PeakFindingTolerance = ms2PeakFindingTolerance;
            MaxNumMissedScan = maxNumMissedScan;
            PeakSearchBinSize = binSize;
            OverlapRatioCutOff = overlapRatioCutOff;
            CorrelationCutOff = correlationCutOff;
            ApexRtTolerance = apexRtTolerance;
            FragmentRankCutOff = fragmentRankCutOff;
            PrecursorRankCutOff = precursorRankCutOff;
            NumHighCorrFragments = numHighCorrFragments;
            HighCorrThreshold = highCorrThreshold;
            MaxRTRangeMS1 = maxRTrangeMS1;
            MaxRTRangeMS2 = maxRTrangeMS2;
            PrecursorIntensityCutOff = precursorIntensityCutOff;
            MinRTRangeForCWT = minRTRangeForCWT;
            SplitMS2Peak = splitMS2Peak;
            SplitMS1Peak = splitMS1Peak;
            SplineTimeInterval = splineTimeInterval;
            MinMS1Mass = minMS1Mass;
            MinMS2Mass = minMS2Mass;
            MaxMass = maxMass;
            Type = type;
            ApexCycleTolerance = apexCycleTolerance;
            ScanCycleSplineTimeInterval = scanCycleSplineInterval;
            CutMs1Peaks = cutMs1Peaks;
            CutMs2Peaks = cutMs2Peaks;
            MinMS1Charge = minMS1Charge;
            MinMS2Charge = minMS2Charge;
            AverageMs2Scans = averageMs2Scans;
            Ms1XICType = ms1XICType;
            Ms2XICType = ms2XICType;
            PFGroupingType = pfGroupingType;
            PseudoMs2ConstructionType = pseudoMs2Type;
            AnalysisType = analysisType;
            CombineFragments = combineFragments;
            CorrelationType = correlationType;
            SGfilterWindowSize = sgFilterWindowSize;
            Ms1SplineType = ms1SplineType;
            Ms2SplineType = ms2SplineType;
            SplineRtInterval = splineRtInterval;
            NumScansPerCycle = numScanPerCycle;
            TrimMs2Peaks = trimMs2Peaks;
            TrimMs2MinSNR = trimMs2MinSNR;
            NoPointsPerMin = noPointsPerMin;
            NumPeaksThreshold = numPeaksThreshold;
        }

        public StringBuilder WriteDIASettings()
        {
            var settings = new StringBuilder();
            settings.AppendLine("DIA Settings:");
            settings.AppendLine("AnalysisType: " + AnalysisType);
            settings.AppendLine("Ms1XICType: " + Ms1XICType);
            settings.AppendLine("Ms2XICType: " + Ms2XICType);
            settings.AppendLine("PFGroupingType: " + PFGroupingType);
            settings.AppendLine("PseudoMs2ConstructionType: " + PseudoMs2ConstructionType);
            settings.AppendLine("CombineFragments: " + CombineFragments);
            settings.AppendLine("CorrelationType: " + CorrelationType);
            settings.AppendLine("Ms1SplineType: " + Ms1SplineType);
            settings.AppendLine("Ms2SplineType: " + Ms2SplineType);
            settings.Append("\n");
            settings.AppendLine("Ms1PeakFindingTolerance: " + Ms1PeakFindingTolerance.ToString());
            settings.AppendLine("Ms2PeakFindingTolerance: " + Ms2PeakFindingTolerance.ToString());
            settings.AppendLine("OverlapRatioCutOff: " + OverlapRatioCutOff);
            settings.AppendLine("CorrelationCutOff: " + CorrelationCutOff);
            settings.AppendLine("ApexRtTolerance: " + ApexRtTolerance);
            settings.AppendLine("ApexCycleTolerance: " + ApexCycleTolerance);
            settings.AppendLine("FragmentRankCutOff: " + FragmentRankCutOff);
            settings.AppendLine("PrecursorRankCutOff: " + PrecursorRankCutOff);
            settings.AppendLine("MaxRTRangeMS1: " + MaxRTRangeMS1);
            settings.AppendLine("MaxRTRangeMS2: " + MaxRTRangeMS2);
            settings.AppendLine("NumScansPerCycle: " + NumScansPerCycle);
            settings.Append("\n");
            settings.AppendLine("TrimMss2Peaks: " + TrimMs2Peaks);
            settings.AppendLine("TrimMs2MinSNR: " + TrimMs2MinSNR);
            settings.AppendLine("PrecursorIntensityCutOff: " + PrecursorIntensityCutOff);
            settings.AppendLine("MinRTRangeForCWT: " + MinRTRangeForCWT);
            settings.AppendLine("SplitMS2Peak: " + SplitMS2Peak);
            settings.AppendLine("SplitMS1Peak: " + SplitMS1Peak);
            settings.AppendLine("SplineRtInterval: " + SplineRtInterval);
            settings.AppendLine("MinMass: " + MinMS1Mass);
            settings.AppendLine("MaxMass: " + MaxMass);
            settings.AppendLine("NumHighCorrFragments: " + NumHighCorrFragments);
            settings.AppendLine("HighCorrThreshold: " + HighCorrThreshold);
            settings.AppendLine("ScanCycleSplineTimeInterval: " + ScanCycleSplineTimeInterval);
            settings.AppendLine("CutMs1Peaks: " + CutMs1Peaks);
            settings.AppendLine("CutMs2Peaks: " + CutMs2Peaks);
            settings.AppendLine("MinCharge: " + MinMS1Charge);
            settings.AppendLine("AverageMs2Scans: " + AverageMs2Scans);
            settings.AppendLine("MaxNumMissedScan: " + MaxNumMissedScan);
            settings.AppendLine("PeakSearchBinSize: " + PeakSearchBinSize);

            return settings;
        }

    }
}
