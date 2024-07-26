using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.ISD
{
    public class Precursor
    {
        public double MonoPeakMz { get; set; }
        public int Charge {  get; set; }
        public double RT {  get; set; }
        public double HighestPeakMz { get; set; }
        public double HighestPeakIntensity {  get; set; }
        public double MonoisotopicMass { get; set; }
        public XIC PrecursorXIC { get; set; }

        public int ScanNumber {  get; set; }

        public Precursor(double monoPeakMz, int charge, double rt, double highestPeakMz, double highestPeakIntensity, double monoisotopicMass, int scanNumber)
        {
            MonoPeakMz = monoPeakMz;
            Charge = charge;
            RT = rt;
            HighestPeakMz = highestPeakMz;
            HighestPeakIntensity = highestPeakIntensity;
            MonoisotopicMass = monoisotopicMass;
            ScanNumber = scanNumber;
        }
    }
}
