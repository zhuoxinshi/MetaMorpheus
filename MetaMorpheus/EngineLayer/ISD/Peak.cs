using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class Peak
    {
        public Peak(double mz, double rt, double intensity, int msLevel, int scanNumber = 0, int index = 0, PeakCurve curve = null)
        {
            Mz = mz;
            Intensity = intensity;
            RT = rt;
            ScanNumber = scanNumber;
            Index = index;
            XIC = curve;
            MsLevel = msLevel;
        }

        public double Mz { get; set; }
        public double Intensity { get; set; }
        public double RT { get; set; }
        public int ScanNumber { get; set; }
        public int Index { get; set; }
        public PeakCurve XIC { get; set; }   
        public int MsLevel { get; set; }
    }
}
