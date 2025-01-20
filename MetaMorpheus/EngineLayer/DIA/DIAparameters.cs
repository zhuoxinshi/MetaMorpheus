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
        public double SplineTimeIntervalDouble { get; set; }
        public double ScanCycleSplineTimeInterval { get; set; }
        public double MinMass { get; set; }
        public double MaxMass { get; set; }
        public int MinCharge { get; set; }
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

        public DIAparameters(Tolerance ms1PeakFindingTolerance, Tolerance ms2PeakFindingTolerance, int maxNumMissedScan = 2, int binSize = 100, 
            double overlapRatioCutOff = 0.3, double correlationCutOff = 0.5, double apexRtTolerance = 0.1, int fragmentRankCutOff = 5000, int precursorRankCutOff = 1000
            , double maxRTrangeMS1 = 0.5, double maxRTrangeMS2 = 2, double highCorrThreshold = 0.5, int numHighCorrFragments = 0, double precursorIntensityCutOff = 10000, double minRTRangeForCWT = 0.1,
            bool splitMS2Peak = false, bool splitMS1Peak = false, float splineTimeInterval = 0.05f, double minMass = 0, double maxMass = 99999, string type = "DIA", int apexCycleTolerance = 2, 
            double scanCycleSplineInterval = 0.025, int minCharge = 1, bool averageMs2Scans = false, XICType ms1XICType = XICType.DeconHighestPeak, 
            XICType ms2XICType = XICType.Peak, PFGroupingType pfGroupingType = PFGroupingType.ScanCycle, PseudoMs2ConstructionType pseudoMs2Type = PseudoMs2ConstructionType.mzPeak, 
            AnalysisType analysisType = AnalysisType.DIAEngine, bool combineFragments = false, CorrelationType correlationType = CorrelationType.CubicSpline_scanCycle,
            bool cutMs1Peaks = false, bool cutMs2Peaks = false)
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
            MinMass = minMass;
            MaxMass = maxMass;
            Type = type;
            ApexCycleTolerance = apexCycleTolerance;
            ScanCycleSplineTimeInterval = scanCycleSplineInterval;
            CutMs1Peaks = cutMs1Peaks;
            CutMs2Peaks = cutMs2Peaks;
            MinCharge = minCharge;
            AverageMs2Scans = averageMs2Scans;
            Ms1XICType = ms1XICType;
            Ms2XICType = ms2XICType;
            PFGroupingType = pfGroupingType;
            PseudoMs2ConstructionType = pseudoMs2Type;
            AnalysisType = analysisType;
            CombineFragments = combineFragments;
            CorrelationType = correlationType;
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
            settings.Append("\n");
            settings.AppendLine("PrecursorIntensityCutOff: " + PrecursorIntensityCutOff);
            settings.AppendLine("MinRTRangeForCWT: " + MinRTRangeForCWT);
            settings.AppendLine("SplitMS2Peak: " + SplitMS2Peak);
            settings.AppendLine("SplitMS1Peak: " + SplitMS1Peak);
            settings.AppendLine("SplineTimeInterval: " + SplineTimeInterval);
            settings.AppendLine("MinMass: " + MinMass);
            settings.AppendLine("MaxMass: " + MaxMass);
            settings.AppendLine("NumHighCorrFragments: " + NumHighCorrFragments);
            settings.AppendLine("HighCorrThreshold: " + HighCorrThreshold);
            settings.AppendLine("ScanCycleSplineTimeInterval: " + ScanCycleSplineTimeInterval);
            settings.AppendLine("CutMs1Peaks: " + CutMs1Peaks);
            settings.AppendLine("CutMs2Peaks: " + CutMs2Peaks);
            settings.AppendLine("MinCharge: " + MinCharge);
            settings.AppendLine("AverageMs2Scans: " + AverageMs2Scans);
            settings.AppendLine("MaxNumMissedScan: " + MaxNumMissedScan);
            settings.AppendLine("PeakSearchBinSize: " + PeakSearchBinSize);

            return settings;
        }

    }
}
