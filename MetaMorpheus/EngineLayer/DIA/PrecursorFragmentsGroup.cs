using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PrecursorFragmentsGroup
    {
        public PrecursorFragmentsGroup(PeakCurve precursorPeakCurve)
        {
            PrecursorPeakCurve = precursorPeakCurve;
            PFpairs = new List<PrecursorFragmentPair>();
        }

        public PeakCurve PrecursorPeakCurve { get; set; }
        public List<PrecursorFragmentPair> PFpairs { get; set; }
        public int Index => PrecursorPeakCurve.Index;
        public int NumHighCorrFragments { get; set; }
        public double MonoPeakMz { get; set; }

        //TODO: finish this
        public static PrecursorFragmentsGroup GroupPF(PeakCurve prePeakCurve, List<PeakCurve> fragPeakCurves)
        {
            return new PrecursorFragmentsGroup(prePeakCurve);
        }

        public void GetNumberOfHighCorrFragments(DIAparameters diaParam)
        {
            NumHighCorrFragments = PFpairs.Where(p => p.Correlation >= diaParam.HighCorrThreshold).Count();
        }
    }
}
