using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// A <see cref="MassIndexingEngine"/> that is filled from EXTERNAL deconvolution results instead of by
    /// deconvoluting each scan itself.
    ///
    /// mzLib's <see cref="MassIndexingEngine.IndexPeaks(MsDataScan[], DeconvolutionParameters, MzRange, double, int)"/>
    /// hard-wires the call to <see cref="Deconvoluter"/>. This subclass swaps out only that step: it fills the
    /// same mass-binned <c>IndexedPeaks</c> jagged array and the same <c>ScanInfoArray</c> from masses that were
    /// already deconvoluted by TopFD (or FLASHDeconv / IsoDec), so every peak-tracing method on the base class -
    /// <c>GetXicByScanIndex</c>, <c>GetAllXics</c> - behaves exactly as it does for MetaMorpheus's own
    /// deconvolution. Nothing about how XICs are traced changes; only where the masses came from.
    ///
    /// Peaks are indexed by their 0-based position in the scan array passed to <see cref="IndexPeaks"/>, and are
    /// appended in ascending scan order, which is the ordering the base class's binary search over each mass bin
    /// requires.
    /// </summary>
    public class ExternalMassIndexingEngine : MassIndexingEngine
    {
        /// <summary>Number of scans in the passed array that had a matching entry in the external results.
        /// Zero means the join failed - most often because the deconvoluted file was renumbered from scan 1
        /// and needs <see cref="ExternalDeconScanMatchMode.RetentionTime"/> matching instead.</summary>
        public int MatchedScanCount { get; private set; }

        /// <summary>Number of external peaks actually indexed after the mass, charge and isolation filters.</summary>
        public int IndexedPeakCount { get; private set; }

        /// <summary>Number of peaks dropped as duplicates of another peak in the same scan (same charge, mass
        /// within tolerance). TopFD 1.8 emits repeated envelope entries; see <see cref="CollapseDuplicates"/>.</summary>
        public int DuplicatePeakCount { get; private set; }

        /// <summary>
        /// Indexes external deconvolution results against <paramref name="scanArray"/>, joining each scan to its
        /// entry in <paramref name="externalResults"/>.
        /// </summary>
        /// <param name="scanArray">Scans from the raw file, in the order that defines the XIC time axis.</param>
        /// <param name="externalResults">Deconvolution results read from msAlign file(s).</param>
        /// <param name="matchMode">How to join a scan to its external entry.</param>
        /// <param name="retentionTimeTolerance">RT half-window (minutes) when matching on retention time.</param>
        /// <param name="mzRange">Optional isolation window; peaks whose m/z falls outside it are dropped.</param>
        /// <param name="minMass">Minimum monoisotopic mass to keep.</param>
        /// <param name="minCharge">Minimum charge state to keep.</param>
        /// <returns>True if at least one peak was indexed.</returns>
        /// <param name="duplicateMassTolerance">When set, peaks in the same scan that share a charge and agree in
        /// mass within this tolerance are collapsed to the most intense one. Pass null to index every peak.</param>
        public bool IndexPeaks(MsDataScan[] scanArray, ExternalDeconResults externalResults,
            ExternalDeconScanMatchMode matchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double retentionTimeTolerance = 0.01, MzRange mzRange = null, double minMass = 0, int minCharge = 1,
            Tolerance duplicateMassTolerance = null)
        {
            if (scanArray == null || scanArray.Length == 0)
                return false;
            if (externalResults == null)
                throw new MetaMorpheusException("No external deconvolution results were supplied.");

            // The index is a jagged array with one bucket per DALTON of neutral mass (MassIndexingEngine overrides
            // BinsPerDalton to 1), each holding the peaks whose mass rounds into that bucket. The tracer later
            // binary-searches these buckets by scan index, which only works because we append in ascending scan
            // order below. MaxMass is the array length, so anything at or above it cannot be indexed at all.
            IndexedPeaks = new List<IndexedMass>[MaxMass];
            ScanInfoArray = new ScanInfo[scanArray.Length];
            MatchedScanCount = 0;
            IndexedPeakCount = 0;
            DuplicatePeakCount = 0;

            // walk the scans in the order given: position in THIS array becomes the peak's ZeroBasedScanIndex,
            // which is the x-axis every XIC built from these scans will share
            for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
            {
                var scan = scanArray[scanIndex];
                if (scan == null)
                    continue;

                // ScanInfo is recorded for EVERY scan, matched or not, so that gaps in the external results
                // are treated as missed scans by the tracer rather than collapsing the time axis.
                ScanInfoArray[scanIndex] = new ScanInfo(scan.OneBasedScanNumber, scanIndex, scan.RetentionTime, scan.MsnOrder);

                // the join: this scan's deconvoluted masses, looked up by its one-based scan number (or by
                // retention time). Null means the external software wrote no entry for it - a skipped scan, or a
                // scan belonging to an ISD channel this results object does not cover.
                var externalPeaks = externalResults.GetPeaksForScan(scan, matchMode, retentionTimeTolerance);
                if (externalPeaks == null)
                    continue;
                MatchedScanCount++;

                foreach (var externalPeak in CollapseDuplicates(externalPeaks, duplicateMassTolerance))
                {
                    // drop peaks outside the mass/charge window, and outside the isolation window when one is
                    // given - external results carry no m/z, so it is reconstructed from mass and charge
                    if (externalPeak.MonoisotopicMass < minMass || externalPeak.Charge < minCharge)
                        continue;
                    if (mzRange != null && !mzRange.Contains(externalPeak.Mz))
                        continue;

                    int roundedMass = (int)Math.Round(externalPeak.MonoisotopicMass * BinsPerDalton, 0);
                    if (roundedMass < 0 || roundedMass >= IndexedPeaks.Length)
                        continue;

                    // one IndexedMass per msAlign peak line - charge states are NOT merged here, exactly as
                    // Deconvoluter yields one envelope per charge, so the tracer sees the same shape either way
                    IndexedPeaks[roundedMass] ??= new List<IndexedMass>();
                    IndexedPeaks[roundedMass].Add(new IndexedMass(
                        externalPeak.ToIsotopicEnvelope(), scan.RetentionTime, scanIndex, scan.MsnOrder));
                    IndexedPeakCount++;
                }
            }

            if (MatchedScanCount == 0)
                throw new MetaMorpheusException(
                    "None of the " + scanArray.Length + " scans matched an entry in the external deconvolution results (" +
                    string.Join(", ", externalResults.SourceFilePaths) + "). " +
                    (matchMode == ExternalDeconScanMatchMode.OneBasedScanNumber
                        ? "Scan numbers do not line up; if the deconvoluted file was renumbered from 1, match on retention time instead."
                        : "No entry fell within the retention time tolerance of " + retentionTimeTolerance + " minutes."));

            return IndexedPeakCount > 0;
        }

        /// <summary>
        /// Collapses peaks of one scan that describe the SAME species - same charge, mass within
        /// <paramref name="tolerance"/> - down to the most intense one.
        ///
        /// TopFD 1.8 repeats envelope entries: one real ISD ms1.msalign holds 3088 such repeats, most of them
        /// byte-identical mass and charge. Beyond being redundant, those duplicates break mzLib's shared tracer:
        /// <c>IndexingEngine.GetAllXics</c> only checks its SEED peak against the already-matched set, while the
        /// peak <c>GetBestPeakFromBins</c> hands back for that seed's scan can be an equal-mass twin that a
        /// previous XIC already claimed - which is then added to the matched dictionary a second time and throws
        /// "An item with the same key has already been added". Removing the duplicate on the way in keeps tracing
        /// deterministic and sidesteps that.
        ///
        /// Intensities are NOT summed: a repeated entry is one species reported twice, so summing would inflate it.
        /// </summary>
        private List<ExternalDeconPeak> CollapseDuplicates(List<ExternalDeconPeak> peaks, Tolerance tolerance)
        {
            if (tolerance == null || peaks.Count < 2)
                return peaks;

            var collapsed = new List<ExternalDeconPeak>(peaks.Count);
            foreach (var chargeGroup in peaks.GroupBy(p => p.Charge))
            {
                ExternalDeconPeak best = null;
                foreach (var peak in chargeGroup.OrderBy(p => p.MonoisotopicMass))
                {
                    if (best != null && tolerance.Within(peak.MonoisotopicMass, best.MonoisotopicMass))
                    {
                        // same species reported twice: keep whichever entry is more intense
                        DuplicatePeakCount++;
                        if (peak.Intensity > best.Intensity)
                        {
                            collapsed[collapsed.Count - 1] = peak;
                            best = peak;
                        }
                        continue;
                    }

                    collapsed.Add(peak);
                    best = peak;
                }
            }

            return collapsed;
        }
    }
}
