using FlashLFQ;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class EnvelopeCurve
    {
        public EnvelopeCurve()
        {

        }

        public int MassIndex {  get; set; }
        public List<Peak>[] Peaks { get; set; }
        public double MonoisotopicMass {  get; set; }
        public int Charge { get; set; }

        public PeakCurve GetEnvelopePeakCurve()
        {
            var fakePeaks = new List<Peak>();
            foreach(var envelopePeaks in Peaks)
            {
                var mz = envelopePeaks.OrderByDescending(p => p.Intensity).First().Mz;
                var rt = envelopePeaks.First().RetentionTime;
                var intensity = envelopePeaks.Sum(p => p.Intensity);
                var fakePeak = new Peak(mz, rt, intensity, 1, envelopePeaks.First().ScanNumber, envelopePeaks.First().ZeroBasedScanIndex);
                fakePeaks.Add(fakePeak);
            }
            var newPC = new PeakCurve(fakePeaks, 1, null, MonoisotopicMass, Charge = Charge);
            return newPC;
        }
    }
}
