using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PrecursorFragmentPair
    {
        public PrecursorFragmentPair(PeakCurve pre, PeakCurve frag, double corr)
        {
            PrecursorPeakCurve = pre;
            FragmentPeakCurve = frag;
            Correlation = corr;
        }

        public PeakCurve PrecursorPeakCurve { get; set; }
        public PeakCurve FragmentPeakCurve { get; set; }    
        public double Correlation { get; set; }
        public int FragmentRank { get; set; }
        public int PrecursorRank { get; set; }

    }
}
