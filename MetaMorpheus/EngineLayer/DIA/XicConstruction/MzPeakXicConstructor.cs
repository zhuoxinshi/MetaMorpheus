using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer.DIA
{
    public class MzPeakXicConstructor : XicConstructor
    {
        public MzPeakXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, XicSpline? xicSpline = null)
            : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, out Dictionary<IIndexedPeak, ExtractedIonChromatogram> matchedPeaks, out object indexingEngine, MzRange isolationRange = null)
        {
            var mzPeakIndexingEngine = PeakIndexingEngine.InitializeIndexingEngine(scans);
            indexingEngine = mzPeakIndexingEngine;
            return mzPeakIndexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks, out matchedPeaks);
        }

    }
}
