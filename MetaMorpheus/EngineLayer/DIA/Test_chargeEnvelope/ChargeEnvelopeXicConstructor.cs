using MassSpectrometry;
using MzLibUtil;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Plotly.NET.StyleParam.Range;

namespace EngineLayer.DIA
{
    public class ChargeEnvelopeXicConstructor : XicConstructor
    {
        private readonly string FlashDeconvResultPath;
        public ChargeEnvelopeXicConstructor(string flashDeconvResultPath, Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, XicSpline? xicSpline = null) : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            FlashDeconvResultPath = flashDeconvResultPath;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, out Dictionary<IIndexedPeak, ExtractedIonChromatogram> matchedPeaks, out object indexingEngine, MzRange isolationRange = null)
        {
            var envelopes = new List<IndexedChargeEnvelope>();
            var flashResultFile = new FlashDeconvFile(FlashDeconvResultPath);
            foreach (var result in flashResultFile.Results)
            {
                int zeroBasedScanIndex = (result.ZeroBasedScanNumber - 1) / 4;
                var indexedChargeEnvelope = new IndexedChargeEnvelope(mass: result.MonoisotopicMass, intensity: result.SumIntensity, retentionTime: result.RetentionTime / 60.0, result.MinCharge, result.MaxCharge, zeroBasedScanIndex);
                envelopes.Add(indexedChargeEnvelope);
            }

            var chargeEnvelopeIndexingEngine = new ChargeEnvelopeIndexingEngine();
            indexingEngine = chargeEnvelopeIndexingEngine;
            if (chargeEnvelopeIndexingEngine.IndexPeaks(scans, envelopes))
            {
                return chargeEnvelopeIndexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks, out matchedPeaks);
            }
            else
            {
                throw new MetaMorpheusException("XIC construction failed.");
            }
        }
    }
}
