using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class IndexedChargeEnvelope : IIndexedPeak
    {
        public float Intensity { get; set; }
        public float RetentionTime { get; set; }
        public int ZeroBasedScanIndex { get; set; }
        public float M { get; set; }
        public int MinCharge { get; set; }
        public int MaxCharge { get; set; }
        public int MsLevel { get; set; }

        public IndexedChargeEnvelope(double mass, double intensity, double retentionTime, int minCharge, int maxCharge, int zeroBasedScanIndex)
        {
            M = (float)mass;
            Intensity = (float)intensity;
            RetentionTime = (float)retentionTime;
            ZeroBasedScanIndex = zeroBasedScanIndex;
        }
    }
}
