using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DeconvolutedMass
    {
        public int Charge { get; set; }
        public double RetentionTime { get; set; }
        public double HighestPeakMz { get; set; }
        public double HighestPeakIntensity { get; set; }
        public double MonoisotopicMass { get; set; }

        public int ScanNumber { get; set; }
        public int ZeroBasedScanIndex {  get; set; }
        public IsotopicEnvelope Envelope { get; set; }
        public int MassIndex => Envelope.MassIndex;
        public int Index;

        public DeconvolutedMass(IsotopicEnvelope envelope, int charge, double rt, double highestPeakMz, double highestPeakIntensity, double monoisotopicMass, int scanNumber, int zeroBasedScanNum)
        {
            Envelope = envelope;
            Charge = charge;
            RetentionTime = rt;
            HighestPeakMz = highestPeakMz;
            HighestPeakIntensity = highestPeakIntensity;
            MonoisotopicMass = monoisotopicMass;
            ScanNumber = scanNumber;
            ZeroBasedScanIndex = zeroBasedScanNum;
        }

        //public static List<Precursor> FindPrecursor(Precursor precursor, List<Precursor> allPrecursors, int maxMissedScans, double massTolerance)
        //{
        //    var precursorList = new List<Precursor>(); 
        //    precursorList.Add(precursor);
        //    int numScans = allPrecursors.Max(p => p.ZeroBasedScanIndex) + 1;

            
        //}
    }
}
