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
        public double ScanCycleSplineTimeInterval { get; set; }
        public double MinMass { get; set; }
        public double MaxMass { get; set; }
        public double MinCharge { get; set; }
        public string Type { get; set; }
        public bool TopDown { get; set; }
        public Tolerance PrecursorMassTolerance => new PpmTolerance(200);
        public bool CutPeaks { get; set; }

        public DIAparameters(Tolerance ms1PeakFindingTolerance, Tolerance ms2PeakFindingTolerance,int maxNumMissedScan, int binSize, 
            double overlapRatioCutOff, double correlationCutOff, double apexRtTolerance, int fragmentRankCutOff = 5000, int precursorRankCutOff = 1000
            , double maxRTrangeMS1 = 0.5, double maxRTrangeMS2 = 2, double highCorrThreshold = 0.5, int numHighCorrFragments = 0, double precursorIntensityCutOff = 10000, double minRTRangeForCWT = 0.1,
            bool splitMS2Peak = true, bool splitMS1Peak = false, float splineTimeInterval = 0.05f, double minMass = 0, double maxMass = 99999, string type = "DIA", 
            bool topdown = false, int apexCycleTolerance = 2, double scanCycleSplineInterval = 0.025, bool cutPeaks = false, double minCharge = 1)
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
            TopDown = topdown;
            ApexCycleTolerance = apexCycleTolerance;
            ScanCycleSplineTimeInterval = scanCycleSplineInterval;
            CutPeaks = cutPeaks;
            MinCharge = minCharge;
        }

    }
}
