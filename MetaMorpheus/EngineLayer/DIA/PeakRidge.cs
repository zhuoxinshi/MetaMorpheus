using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PeakRidge
    {
        public float RT { get; set; }
        public float Intensity { get; set; }
        public int lowScale { get; set; }
        public int ContinuousLevel { get; set; }

        public PeakRidge(float RT, float Intensity, int lowScale, int ContinuousLevel = 0)
        {
            this.RT = RT;
            this.Intensity = Intensity;
            this.lowScale = lowScale;
            this.ContinuousLevel = ContinuousLevel;
        }
    }
}
