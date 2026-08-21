using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// Traces neutral-mass XICs from an external <c>*_ms1.msalign</c> (per-scan deconvoluted masses) and then
    /// GROUPS those XICs onto the LC-MS features of the matching <c>*_ms1.feature</c> file, emitting one
    /// <see cref="ExtractedIonChromatogram"/> per feature.
    ///
    /// The two TopFD outputs describe the same run at different levels: the msAlign holds a mass per charge state
    /// per scan, the feature file holds the assembled LC-MS features (a mass, an RT window and a charge range).
    /// Tracing gives the elution profile the feature file does not carry; grouping onto the feature file gives the
    /// cross-charge assembly the msAlign does not carry. The result is one XIC per TopFD feature, with the
    /// per-charge traces of that feature summed scan by scan.
    ///
    /// This is the external-decon analogue of <see cref="ConsensusMassXicConstructor"/> with
    /// <c>AggregateCharges = true</c>: same one-XIC-per-feature output consumed unchanged by
    /// <see cref="XicGroupingEngine"/>, except that TopFD's feature file defines the features instead of mzLib's
    /// <c>MassFeatureBuilder</c>.
    /// </summary>
    public class FeatureGroupedXicConstructor : ExternalDeconXicConstructor
    {
        /// <summary>The LC-MS features the traced XICs are grouped onto, already restricted to
        /// [<see cref="ExternalDeconXicConstructor.MinMass"/>, <see cref="ExternalDeconXicConstructor.MaxMass"/>).
        /// Traced XICs are filtered to the same mass range before they are ever indexed, so features outside it
        /// could never receive one and would only inflate the denominator when judging how many features were hit.</summary>
        public List<Ms1FeatureRecord> Features { get; set; }

        /// <summary>How many features the feature file held before the mass filter.</summary>
        public int FeatureCountBeforeMassFilter { get; private set; }

        /// <summary>Tolerance for matching a traced XIC's mass to a feature's mass. Falls back to
        /// <see cref="XicConstructor.PeakFindingTolerance"/> when not set.</summary>
        public Tolerance FeatureMassTolerance { get; set; }

        /// <summary>Padding in MINUTES added to each end of a feature's retention time window when deciding
        /// whether a traced XIC belongs to it. 0 means the XIC apex must fall strictly inside the window.</summary>
        public double FeatureRetentionTimeTolerance { get; set; }

        /// <summary>If true (default), a traced XIC only joins a feature whose charge range covers the XIC's
        /// charge state. TopFD reports the charge range it assembled the feature from, so this stops an
        /// unrelated co-eluting species of similar mass being folded in.</summary>
        public bool RequireChargeWithinFeatureRange { get; set; } = true;

        /// <summary>If true, traced XICs that match no feature are passed through as their own XICs instead of
        /// being dropped. Default false - the point of this constructor is to work at TopFD feature level.</summary>
        public bool KeepUnassignedXics { get; set; }

        /// <summary>Diagnostics from the most recent <see cref="GetAllXics"/> call.</summary>
        public int LastTracedXicCount { get; private set; }
        public int LastAssignedXicCount { get; private set; }
        public int LastFeaturesWithXicsCount { get; private set; }

        /// <summary>Why traced XICs failed to land on a feature, from the most recent call: no feature within
        /// <see cref="FeatureMassTolerance"/> at all, versus a mass match that was then rejected on retention time
        /// or charge. Separating the two says whether to loosen the mass tolerance or the RT/charge constraints.</summary>
        public int LastUnassignedNoMassMatchCount { get; private set; }
        public int LastUnassignedRejectedOnRtOrChargeCount { get; private set; }

        /// <summary>Features that cannot receive an XIC no matter the data, because their elution spans fewer than
        /// <see cref="XicConstructor.MinNumberOfPeaks"/> of the scans handed to <see cref="GetAllXics"/>. A large
        /// value means the run's features are too short for the current minNumberOfPeaks, not that matching failed.</summary>
        public int LastFeaturesTooShortToTraceCount { get; private set; }

        /// <summary>The feature each emitted XIC came from, for the most recent call. Null-valued for
        /// pass-through XICs when <see cref="KeepUnassignedXics"/> is set.</summary>
        public Dictionary<ExtractedIonChromatogram, Ms1FeatureRecord> LastFeatureByXic { get; private set; }
            = new Dictionary<ExtractedIonChromatogram, Ms1FeatureRecord>();

        public FeatureGroupedXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, ExternalDeconResults externalResults,
            List<Ms1FeatureRecord> features, Tolerance featureMassTolerance = null,
            double featureRetentionTimeTolerance = 0, double minMass = 0, int minCharge = 1, int maxMass = 30000,
            ExternalDeconScanMatchMode scanMatchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double scanMatchRetentionTimeTolerance = 0.01, XicSpline xicSpline = null)
            : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, externalResults,
                minMass, minCharge, maxMass, scanMatchMode, scanMatchRetentionTimeTolerance, xicSpline)
        {
            FeatureMassTolerance = featureMassTolerance ?? peakFindingTolerance;
            FeatureRetentionTimeTolerance = featureRetentionTimeTolerance;

            // Apply the same mass window the traced XICs get. Without this, features outside [MinMass, MaxMass)
            // sit in the list unmatchable forever and make the "features hit" fraction meaningless.
            FeatureCountBeforeMassFilter = features?.Count ?? 0;
            Features = features?.Where(f => f.Mass >= MinMass && f.Mass < MaxMass).ToList()
                       ?? new List<Ms1FeatureRecord>();
        }

        /// <summary>
        /// Convenience overload reading both files itself: the msAlign the XICs are traced from and the feature
        /// file they are grouped onto.
        /// </summary>
        public FeatureGroupedXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, string msAlignFilePath, string featureFilePath,
            Tolerance featureMassTolerance = null, double featureRetentionTimeTolerance = 0, double minMass = 0,
            int minCharge = 1, int maxMass = 30000,
            ExternalDeconScanMatchMode scanMatchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double scanMatchRetentionTimeTolerance = 0.01, XicSpline xicSpline = null)
            : this(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks,
                ExternalDeconResults.ReadMsAlign(msAlignFilePath), Ms1FeatureRecord.ReadFile(featureFilePath),
                featureMassTolerance, featureRetentionTimeTolerance, minMass, minCharge, maxMass, scanMatchMode,
                scanMatchRetentionTimeTolerance, xicSpline)
        {
        }

        /// <summary>The feature file, as a scalar, so <see cref="DIAparameters.GetSettingsString"/> records how
        /// many features a given search was grouped onto.</summary>
        public int FeatureCount => Features?.Count ?? 0;

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            // 1. trace XICs from the msAlign exactly as ExternalDeconXicConstructor does. These come back ONE PER
            //    CHARGE STATE, because mzLib's tracer never crosses charges - collapsing them is step 3's job.
            var tracedXics = base.GetAllXics(scans, isolationRange);
            LastTracedXicCount = tracedXics.Count;
            LastFeatureByXic = new Dictionary<ExtractedIonChromatogram, Ms1FeatureRecord>();

            if (Features == null || Features.Count == 0)
                throw new MetaMorpheusException("No LC-MS features were supplied to group the traced XICs onto.");

            // 2. assign each traced XIC to at most one feature: mass within tolerance, apex RT inside the
            //    feature's elution window, charge covered by the feature's charge range. Features are sorted by
            //    mass so the candidates for a given XIC can be found by binary search rather than a full scan.
            var featuresByMass = Features.OrderBy(f => f.Mass).ToList();
            var featureMasses = featuresByMass.Select(f => f.Mass).ToArray();
            var xicsByFeature = new Dictionary<Ms1FeatureRecord, List<ExtractedIonChromatogram>>();
            var unassigned = new List<ExtractedIonChromatogram>();

            LastUnassignedNoMassMatchCount = 0;
            LastUnassignedRejectedOnRtOrChargeCount = 0;

            foreach (var xic in tracedXics)
            {
                // each XIC lands on at most ONE feature, but a feature collects many XICs - typically one per
                // charge state of the same proteoform, which is exactly what step 3 wants to sum
                var best = FindBestFeature(xic, featuresByMass, featureMasses, out bool anyMassMatch);
                if (best == null)
                {
                    // split the misses: "nothing near this mass" argues for a wider mass tolerance (or a notch,
                    // since TopFD flips the monoisotope between scans), whereas "right mass, wrong place" argues
                    // for padding the RT window or relaxing RequireChargeWithinFeatureRange
                    if (anyMassMatch) LastUnassignedRejectedOnRtOrChargeCount++;
                    else LastUnassignedNoMassMatchCount++;
                    unassigned.Add(xic);
                    continue;
                }

                if (!xicsByFeature.TryGetValue(best, out var list))
                {
                    list = new List<ExtractedIonChromatogram>();
                    xicsByFeature[best] = list;
                }
                list.Add(xic);
            }

            LastAssignedXicCount = xicsByFeature.Sum(kvp => kvp.Value.Count);
            LastFeaturesWithXicsCount = xicsByFeature.Count;

            // How many features are simply too brief to yield MinNumberOfPeaks points on THIS scan axis? Counted
            // from the scans actually passed in, so it reflects the real acquisition rate rather than an assumed one.
            var scanRetentionTimes = scans.Where(s => s != null).Select(s => s.RetentionTime).OrderBy(t => t).ToArray();
            LastFeaturesTooShortToTraceCount = Features.Count(f =>
                CountScansInWindow(scanRetentionTimes, f.RetentionTimeBegin - FeatureRetentionTimeTolerance,
                    f.RetentionTimeEnd + FeatureRetentionTimeTolerance) < MinNumberOfPeaks);

            // 3. one XIC per feature: sum the assigned per-charge traces scan by scan, at the feature's mass.
            //    This is the cross-charge collapse - TopFD's feature file decides WHICH traces belong together,
            //    where ConsensusMassXicConstructor would instead work that out with MassFeatureBuilder.
            var featureXics = new List<ExtractedIonChromatogram>();
            foreach (var kvp in xicsByFeature)
            {
                // MergeXics returns null when the union of the assigned traces still has too few time points
                var merged = MergeXics(kvp.Key, kvp.Value);
                if (merged == null)
                    continue;
                featureXics.Add(merged);
                LastFeatureByXic[merged] = kvp.Key;
            }

            if (KeepUnassignedXics)
            {
                foreach (var xic in unassigned)
                {
                    featureXics.Add(xic);
                    LastFeatureByXic[xic] = null;
                }
            }

            return featureXics;
        }

        /// <summary>
        /// The feature an XIC belongs to: among features whose mass is within tolerance, whose RT window contains
        /// the XIC apex and whose charge range covers the XIC charge, the one closest in mass; ties broken by
        /// apex-RT proximity. Null when nothing qualifies.
        /// </summary>
        /// <param name="anyMassMatch">Set true when at least one feature was within mass tolerance, whatever the
        /// outcome on retention time and charge. Distinguishes "nothing near this mass" from "right mass, wrong
        /// place".</param>
        private Ms1FeatureRecord FindBestFeature(ExtractedIonChromatogram xic, List<Ms1FeatureRecord> featuresByMass,
            double[] featureMasses, out bool anyMassMatch)
        {
            anyMassMatch = false;
            double xicMass = xic.ApexPeak.M;
            double apexRetentionTime = xic.ApexRT;
            int xicCharge = xic.ApexPeak is IndexedMass indexedMass ? indexedMass.Charge : 0;

            double minMass = FeatureMassTolerance.GetMinimumValue(xicMass);
            double maxMass = FeatureMassTolerance.GetMaximumValue(xicMass);

            int start = System.Array.BinarySearch(featureMasses, minMass);
            if (start < 0) start = ~start;

            Ms1FeatureRecord best = null;
            double bestMassError = double.MaxValue;
            double bestRetentionTimeError = double.MaxValue;

            for (int i = start; i < featuresByMass.Count && featureMasses[i] <= maxMass; i++)
            {
                var feature = featuresByMass[i];
                anyMassMatch = true;
                if (!feature.ContainsRetentionTime(apexRetentionTime, FeatureRetentionTimeTolerance))
                    continue;
                if (RequireChargeWithinFeatureRange && xicCharge > 0 && !feature.ContainsCharge(xicCharge))
                    continue;

                double massError = System.Math.Abs(feature.Mass - xicMass);
                double retentionTimeError = System.Math.Abs(feature.RetentionTimeApex - apexRetentionTime);
                if (massError < bestMassError || (massError == bestMassError && retentionTimeError < bestRetentionTimeError))
                {
                    best = feature;
                    bestMassError = massError;
                    bestRetentionTimeError = retentionTimeError;
                }
            }

            return best;
        }

        /// <summary>
        /// Number of scan retention times falling in [begin, end]. <paramref name="sortedRetentionTimes"/> must be
        /// ascending.
        /// </summary>
        private static int CountScansInWindow(double[] sortedRetentionTimes, double begin, double end)
        {
            int low = System.Array.BinarySearch(sortedRetentionTimes, begin);
            if (low < 0) low = ~low;
            int high = System.Array.BinarySearch(sortedRetentionTimes, end);
            if (high < 0) high = ~high; else high++;
            return System.Math.Max(0, high - low);
        }

        /// <summary>
        /// Collapses the XICs assigned to one feature into a single XIC at the feature's mass: one time point per
        /// scan the feature was seen in, intensities summed over the charge states present in that scan, and the
        /// most intense charge kept so downstream m/z conversion is representative. Mirrors how
        /// <see cref="ConsensusMassXicConstructor"/> collapses a cross-charge MassFeature.
        /// </summary>
        private ExtractedIonChromatogram MergeXics(Ms1FeatureRecord feature, List<ExtractedIonChromatogram> xics)
        {
            var peaks = new List<IIndexedPeak>();
            foreach (var scanGroup in xics.SelectMany(x => x.Peaks).GroupBy(p => p.ZeroBasedScanIndex))
            {
                double intensity = scanGroup.Sum(p => (double)p.Intensity);
                var dominant = scanGroup.OrderByDescending(p => p.Intensity).First();
                int charge = dominant is IndexedMass indexedMass && indexedMass.Charge != 0 ? indexedMass.Charge : 1;
                var envelope = new IsotopicEnvelope(
                    new List<(double mz, double intensity)> { (feature.Mass.ToMz(charge), intensity) },
                    feature.Mass, charge, intensity, 0);
                peaks.Add(new IndexedMass(envelope, dominant.RetentionTime, scanGroup.Key, 1));
            }

            if (peaks.Count < MinNumberOfPeaks)
                return null;

            return new ExtractedIonChromatogram(peaks);
        }
    }
}
