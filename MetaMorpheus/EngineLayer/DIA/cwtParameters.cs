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

        public float MinRTRange { get; set; }
        public float SymThreshold { get; set; }

        public CwtParameters(float maxCurveRTRange, int noPeakPerMin, float minRTRange, float symThreshold)
        {
            MaxCurveRTRange = maxCurveRTRange;
            NoPeakPerMin = noPeakPerMin;
            MinRTRange = minRTRange;
            SymThreshold = symThreshold;
        }

    }
}
