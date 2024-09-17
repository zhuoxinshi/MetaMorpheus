
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DeconvolutionResult
    {
        public IsotopicEnvelope PrecursorEnvelope { get; set; }
        public MsDataScan Ms2Scan { get; set; }
    }
}
