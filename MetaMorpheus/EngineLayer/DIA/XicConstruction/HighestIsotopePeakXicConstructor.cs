using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MzLibUtil;
using FlashLFQ;
using IsotopicEnvelope = MassSpectrometry.IsotopicEnvelope;

namespace EngineLayer.DIA.XicConstruction
{
    public class HighestIsotopePeakXicConstructor : XicConstructor
    {
        public DeconvolutionParameters DeconParameters { get; set; }

        public HighestIsotopePeakXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, DeconvolutionParameters deconParameters, XicSpline? xicSpline = null)
            : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            DeconParameters = deconParameters;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            //find all XICs in mz space
            var mzPeakIndexingEngine = PeakIndexingEngine.InitializeIndexingEngine(scans);
            var allMzXics = mzPeakIndexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks, out var matchedPeaks);
            var foundXics = new HashSet<ExtractedIonChromatogram>();

            //deconvolute everything and order the envelopes by descending intensity of the highest isotope peak
            var deconvolutedMasses = new List<(IsotopicEnvelope envelope, double rt, int scanIndex)>();
            for (int i = 0; i < scans.Length; i++)
            {
                var envelopes = Deconvoluter.Deconvolute(scans[i], DeconParameters, isolationRange);
                deconvolutedMasses.AddRange(envelopes.Select(e => (e, scans[i].RetentionTime, i)));
            }
            deconvolutedMasses.Sort((a, b) => b.envelope.Peaks.Max(p => p.intensity).CompareTo(a.envelope.Peaks.Max(p => p.intensity)));

            //go down the list of envelopes to find the XIC of highest isotope peak from matchedPeaks output from mzPeakIndexingEngine
            foreach (var deconMass in deconvolutedMasses)
            {
                var highestPeak = deconMass.envelope.Peaks.MaxBy(p => p.intensity);
                var indexedPeak = mzPeakIndexingEngine.GetIndexedPeak(highestPeak.mz, deconMass.scanIndex, PeakFindingTolerance);
                if (indexedPeak != null && matchedPeaks.ContainsKey(indexedPeak))
                {
                    var foundXic = matchedPeaks[indexedPeak];
                    if (foundXic != null && !foundXics.Contains(foundXic))
                    {
                        //create an indexedMass which contains mass and charge info to replace the mz peak 
                        var indexedMass = new IndexedMass(deconMass.envelope, deconMass.rt, deconMass.scanIndex, 1);
                        int index = foundXic.Peaks.IndexOf(indexedPeak);
                        if (index >= 0)
                        {
                            foundXic.Peaks[index] = indexedMass;
                            foundXic.SetXicInfo();
                        }
                        foundXics.Add(foundXic);
                    }
                }
            }
            return foundXics.ToList();
        }
    }
}
