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

        public DIAparameters(Tolerance ms1PeakFindingTolerance, Tolerance ms2PeakFindingTolerance,int maxNumMissedScan, int binSize, 
            double overlapRatioCutOff, double correlationCutOff, double apexRtTolerance, int fragmentRankCutOff = 5000, int precursorRankCutOff = 1000)
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
        }

    }
}
