using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA.ML
{
    public class MLgroupingEngine : PfGroupingEngine
    {
        public override List<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, IEnumerable<ExtractedIonChromatogram> fragments)
        {
            throw new NotImplementedException();
        }
    }
}
