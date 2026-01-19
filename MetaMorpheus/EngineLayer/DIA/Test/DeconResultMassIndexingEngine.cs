using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace EngineLayer.DIA
{
    public class DeconResultMassIndexingEngine : MassIndexingEngine
    {
        public DeconResultMassIndexingEngine(IEnumerable<IndexedMass> indexedMasses)
        {
            int maxMass = (int)indexedMasses.Max(im => im.IsotopicEnvelope.MonoisotopicMass);
            IndexedPeaks = new List<IndexedMass>[maxMass + 1];

        }
    }
}
