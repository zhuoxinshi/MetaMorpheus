using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class Feature
    {
        public List<PeakCurve> PeakCurves;
        public PeakCurve PeakCurve;

        public Feature(List<PeakCurve> peakCurves)
        {
            PeakCurves = peakCurves;
        }

        public void GetPeakCurve()
        {
            var peakCurve = new PeakCurve();
            var startCycle = PeakCurves.Min(p => p.StartCycle);
            var endCycle = PeakCurves.Max(p => p.EndCycle);
            for (int i = startCycle; i <= endCycle; i++)
            {
                foreach(var peak in PeakCurves)
                {

                }
            }
        }
    }
}
