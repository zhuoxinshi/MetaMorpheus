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
        public int FragmentRankCutOff { get; set; }
        public int PrecursorRankCutOff { get; set; }
        public double HighCorrThreshold { get; set; }
        public int NumHighCorrFragments { get; set; }
        public double MaxRTRange { get; set; }
        public double PrecursorIntensityCutOff { get; set; }
        public double MinRTRangeForCWT { get; set; }
        public CwtParameters CwtParameters => new CwtParameters(2f, 150, 0.1f, 0.3f);

        public DIAparameters(Tolerance ms1PeakFindingTolerance, Tolerance ms2PeakFindingTolerance,int maxNumMissedScan, int binSize, 
            double overlapRatioCutOff, double correlationCutOff, double apexRtTolerance, int fragmentRankCutOff = 5000, int precursorRankCutOff = 1000
            , double maxRTrange = 2, double highCorrThreshold = 0.5, int numHighCorrFragments = 0, double precursorIntensityCutOff = 10000, double minRTRangeForCWT = 0.1)
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
            MaxRTRange = maxRTrange;
            PrecursorIntensityCutOff = precursorIntensityCutOff;
            MinRTRangeForCWT = minRTRangeForCWT;
        }

    }
}
