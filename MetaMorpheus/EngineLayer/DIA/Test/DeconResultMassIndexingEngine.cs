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
        public DeconResultMassIndexingEngine()
        {
            
        }

        public bool IndexPeaks(MsDataScan[] scanArray, IEnumerable<IndexedMass> indexedMasses)
        {
            if (indexedMasses == null)
                return false;

            //store metadata for each scan
            ScanInfoArray = new ScanInfo[scanArray.Length];
            for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
            {
                var scan = scanArray[scanIndex];
                ScanInfoArray[scanIndex] = new ScanInfo(scan.OneBasedScanNumber, scanIndex, scan.RetentionTime, scan.MsnOrder);
            }

            MaxMass = indexedMasses.Max(p => (int)Math.Round(p.M, 0));
            IndexedPeaks = new List<IndexedMass>[MaxMass + 1];
            bool anyPeaksIndexed = false;
            foreach (var indexedMass in indexedMasses)
            {
                // Calculate the mass bin index for the envelope
                int roundedMass = (int)Math.Round(indexedMass.M * BinsPerDalton, 0);

                // Skip if the mass bin index is out of range
                if (roundedMass >= IndexedPeaks.Length)
                    continue;

                // Initialize the list for this mass bin if it doesn't exist
                IndexedPeaks[roundedMass] ??= new List<IndexedMass>();

                // Add the indexed mass to the appropriate mass bin
                IndexedPeaks[roundedMass].Add(indexedMass);
                anyPeaksIndexed = true;
            }
            return anyPeaksIndexed;
        }

    }
}
