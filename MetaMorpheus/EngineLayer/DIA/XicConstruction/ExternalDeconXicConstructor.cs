using MassSpectrometry;
using MzLibUtil;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace EngineLayer.DIA
{
    /// <summary>
    /// A <see cref="XicConstructor"/> that builds neutral-mass XICs from EXTERNAL deconvolution results
    /// (TopFD / FLASHDeconv / IsoDec msAlign files) rather than from MetaMorpheus's own per-scan deconvolution.
    ///
    /// It is the drop-in counterpart of <see cref="NeutralMassXicConstructor"/>: same mass-binned indexing, same
    /// charge-aware peak tracing across scans, same <see cref="ExtractedIonChromatogram"/> output consumed
    /// unchanged by <see cref="XicGroupingEngine"/> and pseudo-MS2 generation. The ONLY difference is the source
    /// of the per-scan neutral masses - <see cref="ExternalDeconResults"/> read from msAlign, instead of
    /// <see cref="Deconvoluter"/> run on the raw spectra. Set it as the MS1 and/or MS2 XIC constructor on
    /// <see cref="DIAparameters"/> to run the ISD/DIA workflow on TopFD deconvolution.
    ///
    /// The join to the raw scans is by one-based scan number, so a single instance can be loaded with every
    /// msAlign of a run (the MS1 msAlign plus one per ISD source-voltage channel) and handed to both the
    /// precursor and the fragment channel; only the scans actually passed to <see cref="GetAllXics"/> are
    /// indexed on any given call.
    /// </summary>
    public class ExternalDeconXicConstructor : XicConstructor
    {
        /// <summary>The deconvolution results read from msAlign file(s). Used whenever
        /// <see cref="ResultsBySourceVoltage"/> is not set.</summary>
        public ExternalDeconResults ExternalResults { get; set; }

        /// <summary>
        /// Optional: external results keyed by ISD SOURCE VOLTAGE (the <c>sid=NN</c> in the scan filter), so each
        /// fragment channel is traced strictly against its own deconvolution file. When set, the voltage is read
        /// off the scans handed to <see cref="GetAllXics"/> - which <see cref="ISDEngine"/> already partitions one
        /// channel at a time - and the matching results are used; <see cref="ExternalResults"/> is ignored.
        ///
        /// Matching on scan number alone already separates the channels whenever the channel splits preserved the
        /// parent scan numbering (channel k of cycle n is raw scan 4n+k, so the channels never collide). Keying on
        /// voltage makes the separation explicit and holds even if a split renumbered its scans from 1, which
        /// would otherwise let one channel silently pull another channel's masses.
        /// </summary>
        public Dictionary<int, ExternalDeconResults> ResultsBySourceVoltage { get; set; }

        /// <summary>ISD source voltage of the channel used by the most recent <see cref="GetAllXics"/> call, or 0
        /// when voltage keying is not in use.</summary>
        public int LastSourceVoltage { get; private set; }

        /// <summary>How a raw scan is joined to its entry in the external results.</summary>
        public ExternalDeconScanMatchMode ScanMatchMode { get; set; }

        /// <summary>RT half-window in MINUTES, used only when <see cref="ScanMatchMode"/> is
        /// <see cref="ExternalDeconScanMatchMode.RetentionTime"/>.</summary>
        public double ScanMatchRetentionTimeTolerance { get; set; }

        /// <summary>Minimum monoisotopic mass to keep. Use e.g. 3000 on the precursor channel to keep only
        /// intact proteoforms; leave 0 on the fragment channel.</summary>
        public double MinMass { get; set; }

        /// <summary>Minimum charge state to keep.</summary>
        public int MinCharge { get; set; }

        /// <summary>Upper bound of the mass index, in Da. Masses at or above this are dropped. Raise it above the
        /// mzLib default of 30000 when the external software was run with a higher maximum mass (TopFD's default
        /// is 70000).</summary>
        public int MaxMass { get; set; }

        /// <summary>Within one scan, peaks sharing a charge whose masses agree within this tolerance are collapsed
        /// to the most intense - they are one species the external software reported more than once. Set to null
        /// to index every peak as written. Deliberately tighter than <see cref="XicConstructor.PeakFindingTolerance"/>
        /// and never a notched tolerance, so only true repeats are merged.</summary>
        public Tolerance DuplicateMassTolerance { get; set; } = new PpmTolerance(5);

        /// <summary>Diagnostics from the most recent <see cref="GetAllXics"/> call: how many of the passed scans
        /// had a matching external entry, and how many external peaks survived the filters and were indexed.</summary>
        public int LastMatchedScanCount { get; private set; }
        public int LastIndexedPeakCount { get; private set; }

        /// <summary>Peaks dropped by the most recent call as repeats of another peak in the same scan.</summary>
        public int LastDuplicatePeakCount { get; private set; }

        public ExternalDeconXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, ExternalDeconResults externalResults,
            double minMass = 0, int minCharge = 1, int maxMass = 30000,
            ExternalDeconScanMatchMode scanMatchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double scanMatchRetentionTimeTolerance = 0.01, XicSpline xicSpline = null)
            : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            ExternalResults = externalResults;
            MinMass = minMass;
            MinCharge = minCharge;
            MaxMass = maxMass;
            ScanMatchMode = scanMatchMode;
            ScanMatchRetentionTimeTolerance = scanMatchRetentionTimeTolerance;
        }

        /// <summary>
        /// Convenience overload that reads the msAlign file(s) itself.
        /// </summary>
        public ExternalDeconXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, IEnumerable<string> msAlignFilePaths,
            double minMass = 0, int minCharge = 1, int maxMass = 30000,
            ExternalDeconScanMatchMode scanMatchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double scanMatchRetentionTimeTolerance = 0.01, XicSpline xicSpline = null)
            : this(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks,
                ExternalDeconResults.ReadMsAlign(msAlignFilePaths), minMass, minCharge, maxMass,
                scanMatchMode, scanMatchRetentionTimeTolerance, xicSpline)
        {
        }

        /// <summary>
        /// Constructs one <see cref="ExternalDeconXicConstructor"/> that keeps each ISD source voltage on its own
        /// deconvolution file. <paramref name="msAlignFilePathsByVoltage"/> maps the voltage (60, 80, 100 ...) to
        /// that channel's msAlign file(s).
        /// </summary>
        public static ExternalDeconXicConstructor BySourceVoltage(Tolerance peakFindingTolerance,
            int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks,
            Dictionary<int, List<string>> msAlignFilePathsByVoltage, double minMass = 0, int minCharge = 1,
            int maxMass = 30000, ExternalDeconScanMatchMode scanMatchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double scanMatchRetentionTimeTolerance = 0.01, XicSpline xicSpline = null)
        {
            var resultsByVoltage = msAlignFilePathsByVoltage.ToDictionary(
                kvp => kvp.Key, kvp => ExternalDeconResults.ReadMsAlign(kvp.Value));

            // ExternalResults is still populated (with everything) so the class works even if a scan array
            // arrives without a recognisable sid= filter
            return BySourceVoltage(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks,
                resultsByVoltage, ExternalDeconResults.ReadMsAlign(msAlignFilePathsByVoltage.SelectMany(kvp => kvp.Value)),
                minMass, minCharge, maxMass, scanMatchMode, scanMatchRetentionTimeTolerance, xicSpline);
        }

        /// <summary>
        /// Constructs one <see cref="ExternalDeconXicConstructor"/> from results ALREADY split by ISD source
        /// voltage. This is the path for a single combined <c>*_ms2.msalign</c> holding every voltage: read it
        /// once, split it on scan number with <see cref="ExternalDeconResults.SplitByVoltage"/> (or
        /// <see cref="ExternalDeconResults.SplitByIsdCycle"/>), and hand the result here.
        /// </summary>
        /// <param name="fallbackResults">Used when a scan array arrives with no recognisable sid= filter.
        /// Defaults to the union of the per-voltage results.</param>
        public static ExternalDeconXicConstructor BySourceVoltage(Tolerance peakFindingTolerance,
            int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks,
            Dictionary<int, ExternalDeconResults> resultsByVoltage, ExternalDeconResults fallbackResults = null,
            double minMass = 0, int minCharge = 1, int maxMass = 30000,
            ExternalDeconScanMatchMode scanMatchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double scanMatchRetentionTimeTolerance = 0.01, XicSpline xicSpline = null)
        {
            if (resultsByVoltage == null || resultsByVoltage.Count == 0)
                throw new MetaMorpheusException("No per-voltage external deconvolution results were supplied.");

            var constructor = new ExternalDeconXicConstructor(peakFindingTolerance, maxMissedScansAllowed,
                maxPeakHalfWidth, minNumberOfPeaks,
                fallbackResults ?? resultsByVoltage.OrderBy(kvp => kvp.Key).First().Value,
                minMass, minCharge, maxMass, scanMatchMode, scanMatchRetentionTimeTolerance, xicSpline);
            constructor.ResultsBySourceVoltage = resultsByVoltage;
            return constructor;
        }

        /// <summary>
        /// Convenience path for the combined-msAlign layout: read <paramref name="combinedMsAlignFilePath"/> once
        /// and split it into per-voltage results using the scan filters of <paramref name="allScans"/> (the raw
        /// file's scans). Every voltage present in the scans is keyed, MS1 included, so the caller can use the
        /// same object for whichever channels it traces.
        /// </summary>
        public static ExternalDeconXicConstructor FromCombinedMsAlign(Tolerance peakFindingTolerance,
            int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks,
            string combinedMsAlignFilePath, IEnumerable<MsDataScan> allScans, double minMass = 0, int minCharge = 1,
            int maxMass = 30000, ExternalDeconScanMatchMode scanMatchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double scanMatchRetentionTimeTolerance = 0.01, XicSpline xicSpline = null)
        {
            var combined = ExternalDeconResults.ReadMsAlign(combinedMsAlignFilePath);
            var resultsByVoltage = combined.SplitByVoltage(ExternalDeconResults.GetIsdVoltageByScanNumber(allScans));
            return BySourceVoltage(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks,
                resultsByVoltage, combined, minMass, minCharge, maxMass, scanMatchMode,
                scanMatchRetentionTimeTolerance, xicSpline);
        }

        /// <summary>
        /// Picks the external results to trace <paramref name="scans"/> against: the channel-specific results when
        /// <see cref="ResultsBySourceVoltage"/> is set, otherwise <see cref="ExternalResults"/>.
        /// </summary>
        protected ExternalDeconResults GetResultsForScans(MsDataScan[] scans)
        {
            LastSourceVoltage = 0;
            if (ResultsBySourceVoltage == null || ResultsBySourceVoltage.Count == 0)
                return ExternalResults;

            var voltages = scans
                .Where(s => s != null && !string.IsNullOrEmpty(s.ScanFilter))
                .Select(s => ExternalDeconResults.SourceVoltagePattern.Match(s.ScanFilter))
                .Where(m => m.Success)
                .Select(m => int.Parse(m.Groups[1].Value))
                .Distinct()
                .ToList();

            if (voltages.Count == 0)
                return ExternalResults; // no sid= in the filter (e.g. the MS1 channel): fall back to everything
            if (voltages.Count > 1)
                throw new MetaMorpheusException(
                    "XIC construction was handed scans from more than one ISD source voltage (" +
                    string.Join(", ", voltages.OrderBy(v => v)) + "). Each voltage must be traced separately.");

            LastSourceVoltage = voltages[0];
            if (!ResultsBySourceVoltage.TryGetValue(LastSourceVoltage, out var results))
                throw new MetaMorpheusException(
                    "No external deconvolution results were supplied for ISD source voltage " + LastSourceVoltage +
                    ". Voltages available: " + string.Join(", ", ResultsBySourceVoltage.Keys.OrderBy(v => v)) + ".");

            return results;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            // 1. pick which deconvolution results these scans should be traced against. With voltage keying this
            //    reads the sid= off the scans and returns just that ISD channel's file; otherwise it is simply
            //    ExternalResults. Nothing is read from disk here - the msAlign was parsed once, at construction.
            var results = GetResultsForScans(scans);

            // 2. fill the mass index from those results instead of by deconvoluting the spectra. This is the ONLY
            //    step that differs from NeutralMassXicConstructor: same jagged mass-binned array, same ScanInfo,
            //    but each scan's neutral masses come from the file rather than from Deconvoluter. Peaks are
            //    filtered by mass/charge/isolation window and same-scan duplicates collapsed on the way in.
            var indexingEngine = new ExternalMassIndexingEngine { MaxMass = MaxMass };
            bool indexed = indexingEngine.IndexPeaks(scans, results, ScanMatchMode,
                ScanMatchRetentionTimeTolerance, isolationRange, MinMass, MinCharge, DuplicateMassTolerance);

            LastMatchedScanCount = indexingEngine.MatchedScanCount;
            LastIndexedPeakCount = indexingEngine.IndexedPeakCount;
            LastDuplicatePeakCount = indexingEngine.DuplicatePeakCount;

            // IndexPeaks throws if NO scan matched (a broken join); reaching here with false means the join worked
            // but every peak was filtered out, which is a parameter problem rather than a wiring one
            if (!indexed)
                throw new MetaMorpheusException(
                    "XIC construction from external deconvolution results failed: " + indexingEngine.MatchedScanCount +
                    " of " + scans.Length + " scans matched, but no peaks passed the mass/charge/isolation filters.");

            // 3. hand off to mzLib's ordinary tracer, untouched: it sorts every indexed peak by intensity, seeds an
            //    XIC from the most intense unclaimed one, and walks backwards then forwards through scan indices
            //    collecting the nearest-mass peak OF THE SAME CHARGE within PeakFindingTolerance, stopping after
            //    MaxMissedScansAllowed consecutive empty scans or once MaxPeakHalfWidth minutes from the seed.
            //    Traces shorter than MinNumberOfPeaks are discarded; peaks in a kept XIC are claimed exclusively.
            return indexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks);
        }
    }
}
