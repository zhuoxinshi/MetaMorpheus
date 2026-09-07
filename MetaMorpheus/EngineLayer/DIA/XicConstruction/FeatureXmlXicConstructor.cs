using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// How the charge states of one TopFD feature are turned into XICs by <see cref="FeatureXmlXicConstructor"/>.
    /// A feature is one neutral mass observed at several charges, each with its own envelope list, so there is
    /// more than one defensible way to collapse that into a single elution profile.
    /// </summary>
    public enum FeatureChargeMode
    {
        /// <summary>One XIC per feature, intensities SUMMED across charge states within each scan. The default,
        /// and the analogue of what the msAlign path carries.</summary>
        SumAllCharges,

        /// <summary>One XIC per feature, built from the SINGLE most abundant charge state only - the charge whose
        /// in-range envelopes carry the most intensity. Scans in which only weaker charges were observed are not
        /// part of the XIC, so the profile is one charge's measured elution rather than a blend of several.</summary>
        MostAbundantCharge,

        /// <summary>One XIC per charge state, charges left separate.</summary>
        SeparateCharges
    }

    /// <summary>
    /// A <see cref="XicConstructor"/> that reads feature XICs STRAIGHT OUT of TopFD's <c>*_feature.xml</c>,
    /// with no peak finding of any kind.
    ///
    /// The other external-decon constructors rebuild elution profiles by tracing per-scan masses across scans
    /// (<see cref="ExternalDeconXicConstructor"/>) and optionally grouping the traces onto features
    /// (<see cref="FeatureGroupedXicConstructor"/>). TopFD already did both: <c>*_feature.xml</c> stores, for
    /// every proteoform feature, a per-charge <c>envelope_list</c> giving scan, retention time and intensity for
    /// each point of the elution. This class simply converts those lists into
    /// <see cref="ExtractedIonChromatogram"/> objects, so nothing is lost to tracing tolerance, missed-scan
    /// limits, monoisotope flips or peak claiming.
    ///
    /// Envelopes are joined to the scans passed to <see cref="GetAllXics"/> by one-based scan number, exactly as
    /// the msAlign path joins its entries, so an ISD channel array picks up only its own scans.
    /// </summary>
    public class FeatureXmlXicConstructor : XicConstructor
    {
        /// <summary>The features read from <c>*_feature.xml</c>, already restricted to [MinMass, MaxMass). Used
        /// whenever <see cref="FeaturesBySourceVoltage"/> is not set.</summary>
        public List<TopFdFeature> Features { get; set; }

        /// <summary>
        /// Optional: features keyed by ISD SOURCE VOLTAGE (the <c>sid=NN</c> in the scan filter), for when each
        /// fragment channel was deconvoluted separately and so has its own <c>*_feature.xml</c>. When set, the
        /// voltage is read off the scans handed to <see cref="GetAllXics"/> - which <see cref="ISDEngine"/>
        /// already partitions one channel at a time - and that channel's features are used.
        ///
        /// This is what lets a single <see cref="FeatureXmlXicConstructor"/> serve as
        /// <see cref="DIAparameters.Ms2XicConstructor"/> even though the fragment channels live in different
        /// files: ISDEngine has only one MS2 constructor slot, so the per-channel dispatch has to happen in here.
        /// </summary>
        public Dictionary<int, List<TopFdFeature>> FeaturesBySourceVoltage { get; set; }

        /// <summary>ISD source voltage used by the most recent <see cref="GetAllXics"/> call, or 0 when voltage
        /// keying is not in use.</summary>
        public int LastSourceVoltage { get; private set; }

        /// <summary>How many features the file held before the mass, charge and score filters.</summary>
        public int FeatureCountBeforeFilters { get; private set; }

        /// <summary>Minimum monoisotopic mass to keep. Use e.g. 3000 on the precursor channel.</summary>
        public double MinMass { get; set; }

        /// <summary>Maximum monoisotopic mass to keep.</summary>
        public double MaxMass { get; set; }

        /// <summary>Minimum number of charge states a feature must have been observed at.</summary>
        public int MinChargeCount { get; set; }

        /// <summary>Minimum TopFD ECScore. TopFD's own default cutoff when writing the file is only 0.1, so this
        /// is the main lever for dropping low-confidence features.</summary>
        public double MinEcScore { get; set; }

        /// <summary>How the charge states of a feature become XICs; see <see cref="FeatureChargeMode"/>. Defaults
        /// to <see cref="FeatureChargeMode.SumAllCharges"/>.</summary>
        public FeatureChargeMode ChargeMode { get; set; } = FeatureChargeMode.SumAllCharges;

        /// <summary>
        /// Back-compatible view of <see cref="ChargeMode"/>: true for either one-XIC-per-feature mode
        /// (<see cref="FeatureChargeMode.SumAllCharges"/> and <see cref="FeatureChargeMode.MostAbundantCharge"/>),
        /// false for <see cref="FeatureChargeMode.SeparateCharges"/>. Because a bool cannot name three modes,
        /// SETTING it can only choose between summing and separating - assign <see cref="ChargeMode"/> directly to
        /// ask for <see cref="FeatureChargeMode.MostAbundantCharge"/>.
        /// </summary>
        public bool AggregateCharges
        {
            get => ChargeMode != FeatureChargeMode.SeparateCharges;
            set => ChargeMode = value ? FeatureChargeMode.SumAllCharges : FeatureChargeMode.SeparateCharges;
        }

        /// <summary>If true, use each envelope's tallest isotopic peak (<see cref="TopFdEnvelope.MaxIntensity"/>)
        /// as the XIC intensity instead of the summed envelope intensity. Summed is the default because it is the
        /// analogue of what the msAlign path carries.</summary>
        public bool UseMaxIntensityInsteadOfSum { get; set; }

        /// <summary>Diagnostics from the most recent <see cref="GetAllXics"/> call.</summary>
        public int LastFeatureCount { get; private set; }
        public int LastXicCount { get; private set; }
        /// <summary>Features whose envelopes fell entirely outside the scan array handed in - for ISD this is
        /// normal on a fragment channel, since the features were detected on the MS1 channel.</summary>
        public int LastFeaturesWithNoScansInRangeCount { get; private set; }
        /// <summary>Features dropped because too few of their envelopes landed on scans in this array.</summary>
        public int LastFeaturesBelowMinPeaksCount { get; private set; }

        /// <summary>The feature each emitted XIC came from, for the most recent call.</summary>
        public Dictionary<ExtractedIonChromatogram, TopFdFeature> LastFeatureByXic { get; private set; }
            = new Dictionary<ExtractedIonChromatogram, TopFdFeature>();

        /// <summary>Under <see cref="FeatureChargeMode.MostAbundantCharge"/>, the charge state chosen to represent
        /// each feature in the most recent call. Empty in the other modes.</summary>
        public Dictionary<TopFdFeature, int> LastMostAbundantChargeByFeature { get; private set; }
            = new Dictionary<TopFdFeature, int>();

        /// <summary>How many of the charges in <see cref="LastMostAbundantChargeByFeature"/> matched TopFD's own
        /// <c>rep_charge</c> for the feature. TopFD picks its representative charge over the whole run, while this
        /// class picks it over the scans handed in, so the two need not agree on an ISD channel.</summary>
        public int LastMostAbundantChargeMatchingTopFdCount { get; private set; }

        /// <param name="peakFindingTolerance">Not used for peak finding - there is none - but kept on the base
        /// class so this constructor is interchangeable with the others and its settings are recorded.</param>
        /// <param name="aggregateCharges">Legacy two-way switch, used only when <paramref name="chargeMode"/> is
        /// null: true sums the charge states, false keeps them separate.</param>
        /// <param name="chargeMode">Names the charge handling outright and overrides
        /// <paramref name="aggregateCharges"/>. Pass <see cref="FeatureChargeMode.MostAbundantCharge"/> to build
        /// each feature's XIC from its strongest charge state alone.</param>
        public FeatureXmlXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, List<TopFdFeature> features, double minMass = 0,
            double maxMass = double.MaxValue, int minChargeCount = 1, double minEcScore = 0,
            bool aggregateCharges = true, XicSpline xicSpline = null, FeatureChargeMode? chargeMode = null)
            : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            MinMass = minMass;
            MaxMass = maxMass;
            MinChargeCount = minChargeCount;
            MinEcScore = minEcScore;
            ChargeMode = chargeMode ??
                (aggregateCharges ? FeatureChargeMode.SumAllCharges : FeatureChargeMode.SeparateCharges);

            FeatureCountBeforeFilters = features?.Count ?? 0;
            Features = (features ?? new List<TopFdFeature>())
                .Where(f => f.Mass >= MinMass && f.Mass < MaxMass
                            && f.ChargeCount >= MinChargeCount && f.EcScore >= MinEcScore)
                .ToList();
        }

        /// <summary>Convenience overload that reads the feature xml itself.</summary>
        public FeatureXmlXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, string featureXmlFilePath, double minMass = 0,
            double maxMass = double.MaxValue, int minChargeCount = 1, double minEcScore = 0,
            bool aggregateCharges = true, XicSpline xicSpline = null, FeatureChargeMode? chargeMode = null)
            : this(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks,
                TopFdFeatureXmlFile.ReadFile(featureXmlFilePath), minMass, maxMass, minChargeCount, minEcScore,
                aggregateCharges, xicSpline, chargeMode)
        {
        }

        public int FeatureCount => Features?.Count ?? 0;

        /// <summary>
        /// Builds one constructor that serves every ISD fragment channel, each from its own feature xml.
        /// <paramref name="featureXmlPathsByVoltage"/> maps the source voltage (60, 80, 100 ...) to that
        /// channel's file. Hand the result to <see cref="DIAparameters.Ms2XicConstructor"/> and
        /// <see cref="ISDEngine"/> drives it unchanged.
        /// </summary>
        public static FeatureXmlXicConstructor BySourceVoltage(Tolerance peakFindingTolerance,
            int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks,
            Dictionary<int, string> featureXmlPathsByVoltage, double minMass = 0, double maxMass = double.MaxValue,
            int minChargeCount = 1, double minEcScore = 0, bool aggregateCharges = true, XicSpline xicSpline = null,
            FeatureChargeMode? chargeMode = null)
        {
            if (featureXmlPathsByVoltage == null || featureXmlPathsByVoltage.Count == 0)
                throw new MetaMorpheusException("No per-voltage feature xml files were supplied.");

            var byVoltage = new Dictionary<int, List<TopFdFeature>>();
            foreach (var kvp in featureXmlPathsByVoltage)
            {
                // filter each channel the same way the single-file constructor filters its own list
                byVoltage[kvp.Key] = TopFdFeatureXmlFile.ReadFile(kvp.Value)
                    .Where(f => f.Mass >= minMass && f.Mass < maxMass
                                && f.ChargeCount >= minChargeCount && f.EcScore >= minEcScore)
                    .ToList();
            }

            // Features is seeded with the lowest voltage so the object is usable even if a scan array arrives
            // without a recognisable sid= filter
            var constructor = new FeatureXmlXicConstructor(peakFindingTolerance, maxMissedScansAllowed,
                maxPeakHalfWidth, minNumberOfPeaks, byVoltage.OrderBy(k => k.Key).First().Value,
                minMass, maxMass, minChargeCount, minEcScore, aggregateCharges, xicSpline, chargeMode);
            constructor.FeaturesBySourceVoltage = byVoltage;
            constructor.FeatureCountBeforeFilters = byVoltage.Sum(kvp => kvp.Value.Count);
            return constructor;
        }

        /// <summary>
        /// Picks the feature list for <paramref name="scans"/>: the channel-specific one when
        /// <see cref="FeaturesBySourceVoltage"/> is set, otherwise <see cref="Features"/>.
        /// </summary>
        private List<TopFdFeature> GetFeaturesForScans(MsDataScan[] scans)
        {
            LastSourceVoltage = 0;
            if (FeaturesBySourceVoltage == null || FeaturesBySourceVoltage.Count == 0)
                return Features;

            var voltages = scans
                .Where(s => s != null && !string.IsNullOrEmpty(s.ScanFilter))
                .Select(s => ExternalDeconResults.SourceVoltagePattern.Match(s.ScanFilter))
                .Where(m => m.Success)
                .Select(m => int.Parse(m.Groups[1].Value))
                .Distinct()
                .ToList();

            if (voltages.Count == 0)
                return Features; // no sid= in the filter (e.g. an MS1 channel): fall back
            if (voltages.Count > 1)
                throw new MetaMorpheusException(
                    "XIC construction was handed scans from more than one ISD source voltage (" +
                    string.Join(", ", voltages.OrderBy(v => v)) + "). Each voltage must be built separately.");

            LastSourceVoltage = voltages[0];
            if (!FeaturesBySourceVoltage.TryGetValue(LastSourceVoltage, out var features))
                throw new MetaMorpheusException(
                    "No feature xml was supplied for ISD source voltage " + LastSourceVoltage +
                    ". Voltages available: " + string.Join(", ", FeaturesBySourceVoltage.Keys.OrderBy(v => v)) + ".");

            return features;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            // when voltage keyed, swap in this channel's features before anything else reads Features
            Features = GetFeaturesForScans(scans);

            if (Features == null || Features.Count == 0)
                throw new MetaMorpheusException("No TopFD features were supplied to build XICs from.");

            // 1. build the scan-number -> position lookup. TopFD's envelopes name RAW one-based scan numbers,
            //    while an XIC peak is addressed by its position in the array handed in here. For an ISD channel
            //    that array holds only every 4th raw scan, so this map is also what confines a feature to the
            //    channel it was detected on: envelopes naming a scan absent from the map are simply skipped.
            var scanIndexByScanNumber = new Dictionary<int, int>(scans.Length);
            var retentionTimeByScanIndex = new double[scans.Length];
            for (int i = 0; i < scans.Length; i++)
            {
                if (scans[i] == null) continue;
                scanIndexByScanNumber[scans[i].OneBasedScanNumber] = i;
                retentionTimeByScanIndex[i] = scans[i].RetentionTime;
            }

            LastFeatureCount = Features.Count;
            LastFeaturesWithNoScansInRangeCount = 0;
            LastFeaturesBelowMinPeaksCount = 0;
            LastFeatureByXic = new Dictionary<ExtractedIonChromatogram, TopFdFeature>();
            LastMostAbundantChargeByFeature = new Dictionary<TopFdFeature, int>();
            LastMostAbundantChargeMatchingTopFdCount = 0;

            // 2. one pass over the features - there is no seeding, no walking and no tolerance anywhere below,
            //    because TopFD already decided which envelopes form which feature. This loop only reshapes its
            //    answer into the ExtractedIonChromatogram the rest of the pipeline consumes.
            var xics = new List<ExtractedIonChromatogram>();
            foreach (var feature in Features)
            {
                // features carry no m/z of their own, so an isolation window is tested against the per-charge
                // mono m/z: keep the feature if ANY of its charge states was observed inside the window
                if (isolationRange != null && !feature.SingleChargeFeatures.Any(c => isolationRange.Contains(c.MonoisotopicMz)))
                    continue;

                // The mode decides how the feature's charge states are partitioned into XICs:
                //   SumAllCharges      - one group holding every charge -> one XIC per feature, summed per scan
                //   MostAbundantCharge - one group holding only the strongest charge -> one XIC per feature, but
                //                        carrying that charge's own measured profile rather than a blend
                //   SeparateCharges    - one group per charge -> one XIC per charge state, which is what a
                //                        fragment channel wants if there is no cross-charge proteoform to rebuild
                IEnumerable<TopFdSingleChargeFeature>[] chargeGroups;
                switch (ChargeMode)
                {
                    case FeatureChargeMode.SeparateCharges:
                        chargeGroups = feature.SingleChargeFeatures.Select(c => new[] { c }.AsEnumerable()).ToArray();
                        break;

                    case FeatureChargeMode.MostAbundantCharge:
                        var dominantCharge = SelectMostAbundantCharge(feature, scanIndexByScanNumber);
                        // no charge has an envelope on these scans: leave chargeGroups empty so the feature falls
                        // through to the "no scans in range" tally below, exactly as in the other modes
                        chargeGroups = dominantCharge == null
                            ? new IEnumerable<TopFdSingleChargeFeature>[0]
                            : new[] { new[] { dominantCharge }.AsEnumerable() };
                        if (dominantCharge != null)
                        {
                            LastMostAbundantChargeByFeature[feature] = dominantCharge.Charge;
                            if (dominantCharge.Charge == feature.RepresentativeCharge)
                                LastMostAbundantChargeMatchingTopFdCount++;
                        }
                        break;

                    default:
                        chargeGroups = new[] { feature.SingleChargeFeatures.AsEnumerable() };
                        break;
                }

                bool anyScanInRange = false;
                foreach (var group in chargeGroups)
                {
                    var peaks = BuildPeaks(feature, group, scanIndexByScanNumber, retentionTimeByScanIndex,
                        out bool groupHadScansInRange);
                    anyScanInRange |= groupHadScansInRange;

                    if (peaks == null || peaks.Count < MinNumberOfPeaks)
                    {
                        // only count it as "too short" if the feature actually overlapped these scans - a feature
                        // from another ISD channel contributing nothing is expected, not a failure
                        if (groupHadScansInRange) LastFeaturesBelowMinPeaksCount++;
                        continue;
                    }

                    var xic = new ExtractedIonChromatogram(peaks);
                    if (xic.EndRT - xic.StartRT > 2) continue;
                    xics.Add(xic);
                    LastFeatureByXic[xic] = feature;
                }

                if (!anyScanInRange)
                    LastFeaturesWithNoScansInRangeCount++;
            }

            LastXicCount = xics.Count;
            return xics;
        }

        /// <summary>
        /// Picks the charge state that carries the most intensity ON THE SCANS HANDED IN, which is what
        /// <see cref="FeatureChargeMode.MostAbundantCharge"/> builds its XIC from.
        ///
        /// The ranking is deliberately computed from the in-range envelopes rather than read off the feature xml's
        /// per-charge <c>intensity</c> field or TopFD's <c>rep_charge</c>: both of those describe the charge over
        /// the WHOLE run, and an ISD fragment channel sees only every fourth scan, so the strongest charge there
        /// need not be the strongest overall. It also uses the same intensity measure the XIC will carry
        /// (<see cref="UseMaxIntensityInsteadOfSum"/>), so the charge that wins is the charge that would have
        /// contributed most to the summed XIC.
        ///
        /// Returns null when no charge state has an envelope on any of these scans - i.e. the feature does not
        /// elute within the array at all.
        /// </summary>
        private TopFdSingleChargeFeature SelectMostAbundantCharge(TopFdFeature feature,
            Dictionary<int, int> scanIndexByScanNumber)
        {
            TopFdSingleChargeFeature best = null;
            double bestIntensity = double.NegativeInfinity;

            foreach (var chargeFeature in feature.SingleChargeFeatures)
            {
                double total = 0;
                bool anyInRange = false;
                foreach (var envelope in chargeFeature.Envelopes)
                {
                    if (!scanIndexByScanNumber.ContainsKey(envelope.OneBasedScanNumber))
                        continue;
                    anyInRange = true;
                    total += UseMaxIntensityInsteadOfSum ? envelope.MaxIntensity : envelope.IntensitySum;
                }

                // a charge with no envelope here cannot represent the feature on these scans; one whose envelopes
                // happen to sum to zero still can, which is why eligibility is anyInRange and not total > 0
                if (!anyInRange)
                    continue;

                // strictly greater, and the xml lists charges in ascending order, so ties resolve to the LOWER
                // charge - the choice has to be reproducible run to run
                if (total > bestIntensity)
                {
                    bestIntensity = total;
                    best = chargeFeature;
                }
            }

            return best;
        }

        /// <summary>
        /// Turns a feature's envelopes into XIC peaks at the feature's mass. Envelopes of the same scan (i.e.
        /// different charge states, when aggregating) are summed, and the most intense charge in that scan is kept
        /// so the m/z conversion downstream is representative.
        /// </summary>
        private List<IIndexedPeak> BuildPeaks(TopFdFeature feature, IEnumerable<TopFdSingleChargeFeature> chargeFeatures,
            Dictionary<int, int> scanIndexByScanNumber, double[] retentionTimeByScanIndex,
            out bool hadScansInRange)
        {
            hadScansInRange = false;

            // Accumulate one entry per scan. Intensity is the running total over charge states in that scan;
            // Charge/ChargeIntensity track which charge contributed the most there, so the winner can be kept as
            // the representative charge without a second pass.
            var byScan = new Dictionary<int, (double Intensity, int Charge, double ChargeIntensity)>();
            foreach (var chargeFeature in chargeFeatures)
            {
                foreach (var envelope in chargeFeature.Envelopes)
                {
                    if (!scanIndexByScanNumber.TryGetValue(envelope.OneBasedScanNumber, out int scanIndex))
                        continue; // envelope belongs to a scan outside this array (e.g. another ISD channel)
                    hadScansInRange = true;

                    // IntensitySum is the whole isotopic envelope in this scan (the msAlign path's intensity is
                    // the same quantity); MaxIntensity is just its tallest isotope peak, ~6-7x smaller
                    double intensity = UseMaxIntensityInsteadOfSum ? envelope.MaxIntensity : envelope.IntensitySum;
                    if (byScan.TryGetValue(scanIndex, out var existing))
                    {
                        // always add to the total; swap the representative charge only if this one is stronger
                        byScan[scanIndex] = intensity > existing.ChargeIntensity
                            ? (existing.Intensity + intensity, chargeFeature.Charge, intensity)
                            : (existing.Intensity + intensity, existing.Charge, existing.ChargeIntensity);
                    }
                    else
                    {
                        byScan[scanIndex] = (intensity, chargeFeature.Charge, intensity);
                    }
                }
            }

            if (byScan.Count == 0)
                return null; // the feature elutes entirely outside the scans handed in

            // emit in scan order - ExtractedIonChromatogram expects its peaks ordered along the time axis
            var peaks = new List<IIndexedPeak>(byScan.Count);
            foreach (var kvp in byScan.OrderBy(k => k.Key))
            {
                int charge = kvp.Value.Charge != 0 ? kvp.Value.Charge : 1;
                // the raw file's retention time, so the XIC shares an axis with everything else built from
                // these scans rather than mixing in the feature xml's own clock
                double retentionTime = retentionTimeByScanIndex[kvp.Key];

                // Every point sits at the FEATURE's mass, not at the per-scan deconvoluted mass. That is
                // deliberate: it is what makes the XIC represent one proteoform, and it sidesteps the
                // monoisotope flipping that breaks mass-based tracing on these files. The synthetic single-peak
                // envelope exists only so IndexedMass can carry mass, charge and intensity together.
                var envelope = new IsotopicEnvelope(
                    new List<(double mz, double intensity)> { (feature.Mass.ToMz(charge), kvp.Value.Intensity) },
                    feature.Mass, charge, kvp.Value.Intensity, 0);
                peaks.Add(new IndexedMass(envelope, retentionTime, kvp.Key, 1));
            }

            return peaks;
        }
    }
}
