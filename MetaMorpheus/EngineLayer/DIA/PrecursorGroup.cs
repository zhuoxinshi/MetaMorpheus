using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PrecursorGroup
    {
        public List<PeakCurve> PeakCurves;
        public PeakCurve CombinedPeakCurve;

        public PrecursorGroup()
        {
            PeakCurves = new List<PeakCurve>();
        }

        public PrecursorGroup(List<PeakCurve> peakCurves)
        {
            PeakCurves = peakCurves;
        }

        public void CombinePeakCurves()
        {
            CombinedPeakCurve = new PeakCurve();
            var start = PeakCurves.Select(p => p.StartCycle).Min();
            var end = PeakCurves.Select(p => p.EndCycle).Max();
            for (int i = start; i <= end; i++)
            {
                var peaksAtCycle = new List<Peak>();
                foreach (var pc in PeakCurves)
                {
                    int index = Array.BinarySearch(pc.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray(), i);
                    if (index >= 0)
                    {
                        peaksAtCycle.Add(pc.Peaks[index]);
                    }
                }
                var combinedPeak = new Peak(rt: peaksAtCycle.First().RetentionTime, ZeroBasedScanNumber: i, intensity: peaksAtCycle.Sum(p => p.Intensity));
                CombinedPeakCurve.Peaks.Add(combinedPeak);
            }
        }
    }

    
}
