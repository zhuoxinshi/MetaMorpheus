using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class CwtParameters
    {
        public float MaxCurveRTRange { get; set; }
        public int NoPeakPerMin { get; set; }

        public float MaxRTDiff { get; set; }
        public float SymThreshold { get; set; }
        public float SNRThreshold { get; set; }
        public float MinPeakWidth { get; set; }
        public float SplineTimeInterval { get; set; }

        public CwtParameters(float maxCurveRTRange, int noPeakPerMin, float minRTRange, float symThreshold, float sNRThreshold = 1.5f, float minPeakWidth = 0.1f, float splineTimeInterval = 0.005f)
        {
            MaxCurveRTRange = maxCurveRTRange;
            NoPeakPerMin = noPeakPerMin;
            MaxRTDiff = minRTRange;
            SymThreshold = symThreshold;
            SNRThreshold = sNRThreshold;
            MinPeakWidth = minPeakWidth;
            SplineTimeInterval = splineTimeInterval;
        }

    }
}
