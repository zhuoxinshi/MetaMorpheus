using Easy.Common.Extensions;
using MzLibUtil;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Nett.TomlObjectFactory;
using static Plotly.NET.StyleParam.Range;

namespace MassSpectrometry
{
    public class ChargeEnvelopeIndexingEngine : IndexingEngine<IndexedChargeEnvelope>
    {
        protected override int BinsPerDalton => 1;
        public int MaxMass { get; set; }
        public ChargeEnvelopeIndexingEngine()
        {
        }

        public bool IndexPeaks(MsDataScan[] scanArray, IEnumerable<IndexedChargeEnvelope> envelopes)
        {
            if (envelopes == null)
                return false;

            //store metadata for each scan
            ScanInfoArray = new ScanInfo[scanArray.Length];
            for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
            {
                var scan = scanArray[scanIndex];
                ScanInfoArray[scanIndex] = new ScanInfo(scan.OneBasedScanNumber, scanIndex, scan.RetentionTime, scan.MsnOrder);
            }

            MaxMass = envelopes.Max(p => (int)Math.Round(p.M, 0));
            IndexedPeaks = new List<IndexedChargeEnvelope>[MaxMass + 1];
            bool anyPeaksIndexed = false;
            foreach (var envelope in envelopes)
            {
                // Calculate the mass bin index for the envelope
                int roundedMass = (int)Math.Round(envelope.M * BinsPerDalton, 0);

                // Skip if the mass bin index is out of range
                if (roundedMass >= IndexedPeaks.Length)
                    continue;

                // Initialize the list for this mass bin if it doesn't exist
                IndexedPeaks[roundedMass] ??= new List<IndexedChargeEnvelope>();

                // Add the indexed mass to the appropriate mass bin
                IndexedPeaks[roundedMass].Add(envelope);
                anyPeaksIndexed = true;
            }
            return anyPeaksIndexed;
        }

        public static ChargeEnvelopeIndexingEngine? GetChargeEnvelopeIndexingEngine(IEnumerable<MsDataScan> scans, string flashDeconvResultPath)
        {
            var envelopes = new List<IndexedChargeEnvelope>();
            var flashResultFile = new FlashDeconvFile(flashDeconvResultPath);
            foreach (var result in flashResultFile.Results)
            {
                int zeroBasedScanIndex = (result.ZeroBasedScanNumber - 1)/4;
                var indexedChargeEnvelope = new IndexedChargeEnvelope(mass: result.MonoisotopicMass, intensity: result.SumIntensity, retentionTime: result.RetentionTime / 60.0, result.MinCharge, result.MaxCharge, zeroBasedScanIndex);
                envelopes.Add(indexedChargeEnvelope);
            }

            var engine = new ChargeEnvelopeIndexingEngine();
            if (engine.IndexPeaks(scans.ToArray(), envelopes))
                return engine;
            return null;
        }
    }
}
