using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace EngineLayer.DIA
{
    /// <summary>
    /// An <see cref="XicConstructor"/> that reads XICs directly from TopFD feature files (<c>*_feature.xml</c>)
    /// instead of tracing peaks or masses across scans. Each TopFD feature becomes one XIC at the feature's
    /// monoisotopic mass. For every scan in which the feature was observed, the envelope intensities of all its
    /// charge states are summed into one XIC point.
    ///
    /// Envelopes are joined to the scans handed to <see cref="GetAllXics"/> by one-based scan number, so the
    /// feature file must use the scan numbers of the spectra file those scans come from. Envelopes on scans that
    /// are not in the array are ignored.
    ///
    /// For ISD data, where TopFD is run separately on the MS1 channel and on each source-voltage channel, build
    /// the fragment constructor with <see cref="BySourceVoltage"/>: it reads the source voltage (<c>sid=</c>) from
    /// the filter of the scans it is handed and uses that channel's features.
    /// </summary>
    public class FeatureXmlXicConstructor : XicConstructor
    {
        private static readonly Regex SourceVoltagePattern = new Regex(@"sid=(\d+)");

        /// <summary>Features used when the constructor is not keyed by source voltage; null otherwise.</summary>
        public List<TopFdFeature> Features { get; private set; }

        /// <summary>Features of each ISD source-voltage channel, keyed by voltage; null when not keyed by voltage.</summary>
        public Dictionary<int, List<TopFdFeature>> FeaturesBySourceVoltage { get; private set; }

        /// <summary>Minimum monoisotopic mass of a feature (inclusive).</summary>
        public double MinMass { get; }

        /// <summary>Maximum monoisotopic mass of a feature (exclusive).</summary>
        public double MaxMass { get; }

        /// <summary>Minimum TopFD ECScore of a feature.</summary>
        public double MinEcScore { get; }

        /// <summary>XICs spanning more than this many minutes, first to last point, are discarded.</summary>
        public double MaxXicWidth { get; }

        /// <summary>The feature each XIC returned by the most recent <see cref="GetAllXics"/> call came from.</summary>
        public Dictionary<ExtractedIonChromatogram, TopFdFeature> LastFeatureByXic { get; private set; }
            = new Dictionary<ExtractedIonChromatogram, TopFdFeature>();

        /// <param name="features">TopFD features, e.g. from <see cref="TopFdFeatureXmlFile.ReadFile"/>.</param>
        /// <param name="minNumberOfPeaks">Minimum number of XIC points, i.e. scans of the array handed to
        /// <see cref="GetAllXics"/> in which the feature was observed.</param>
        /// <param name="minMass">Minimum monoisotopic mass of a feature (inclusive).</param>
        /// <param name="maxMass">Maximum monoisotopic mass of a feature (exclusive).</param>
        /// <param name="minEcScore">Minimum TopFD ECScore of a feature.</param>
        /// <param name="maxXicWidth">XICs spanning more than this many minutes are discarded.</param>
        /// <param name="xicSpline">Optional spline applied by <see cref="XicConstructor.GetAllXicsWithXicSpline"/>.</param>
        public FeatureXmlXicConstructor(List<TopFdFeature> features, int minNumberOfPeaks, double minMass = 0,
            double maxMass = double.MaxValue, double minEcScore = 0, double maxXicWidth = double.MaxValue,
            XicSpline xicSpline = null)
            : this(minNumberOfPeaks, minMass, maxMass, minEcScore, maxXicWidth, xicSpline)
        {
            Features = ApplyFeatureFilters(features);
        }

        /// <summary>Reads the features from a TopFD <c>*_feature.xml</c>.</summary>
        public FeatureXmlXicConstructor(string featureXmlFilePath, int minNumberOfPeaks, double minMass = 0,
            double maxMass = double.MaxValue, double minEcScore = 0, double maxXicWidth = double.MaxValue,
            XicSpline xicSpline = null)
            : this(TopFdFeatureXmlFile.ReadFile(featureXmlFilePath), minNumberOfPeaks, minMass, maxMass,
                minEcScore, maxXicWidth, xicSpline)
        {
        }

        // The base-class peak-finding settings (tolerance, missed scans, peak half width) do not apply: TopFD has
        // already assigned every envelope to its feature, so no peaks are traced here.
        private FeatureXmlXicConstructor(int minNumberOfPeaks, double minMass, double maxMass, double minEcScore,
            double maxXicWidth, XicSpline xicSpline)
            : base(null, 0, 0, minNumberOfPeaks, xicSpline)
        {
            MinMass = minMass;
            MaxMass = maxMass;
            MinEcScore = minEcScore;
            MaxXicWidth = maxXicWidth;
        }

        /// <summary>
        /// Builds one constructor for all ISD fragment channels, each read from its own TopFD feature xml.
        /// <paramref name="featureXmlPathsByVoltage"/> maps the source voltage (the <c>sid=</c> value in the scan
        /// filter, e.g. 60, 80, 100) to that channel's file. Every scan array handed to <see cref="GetAllXics"/>
        /// must come from a single voltage, which is how <see cref="ISDEngine"/> calls it.
        /// </summary>
        public static FeatureXmlXicConstructor BySourceVoltage(Dictionary<int, string> featureXmlPathsByVoltage,
            int minNumberOfPeaks, double minMass = 0, double maxMass = double.MaxValue, double minEcScore = 0,
            double maxXicWidth = double.MaxValue, XicSpline xicSpline = null)
        {
            if (featureXmlPathsByVoltage == null || featureXmlPathsByVoltage.Count == 0)
                throw new MetaMorpheusException("No feature xml files were supplied for the ISD source voltages.");

            var constructor = new FeatureXmlXicConstructor(minNumberOfPeaks, minMass, maxMass, minEcScore, maxXicWidth, xicSpline);
            constructor.FeaturesBySourceVoltage = featureXmlPathsByVoltage.ToDictionary(
                kvp => kvp.Key,
                kvp => constructor.ApplyFeatureFilters(TopFdFeatureXmlFile.ReadFile(kvp.Value)));
            return constructor;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            var features = GetFeaturesForScans(scans);

            // TopFD envelopes name one-based scan numbers; an XIC point is addressed by its position in this array
            var scanIndexByScanNumber = new Dictionary<int, int>(scans.Length);
            for (int i = 0; i < scans.Length; i++)
            {
                if (scans[i] != null)
                    scanIndexByScanNumber[scans[i].OneBasedScanNumber] = i;
            }

            var xics = new List<ExtractedIonChromatogram>();
            var featureByXic = new Dictionary<ExtractedIonChromatogram, TopFdFeature>();
            foreach (var feature in features)
            {
                // a feature has no m/z of its own; keep it if any of its charge states lies in the isolation window
                if (isolationRange != null && !feature.SingleChargeFeatures.Any(c => isolationRange.Contains(c.MonoisotopicMz)))
                    continue;

                var peaks = BuildPeaks(feature, scans, scanIndexByScanNumber);
                if (peaks.Count == 0 || peaks.Count < MinNumberOfPeaks)
                    continue;

                var xic = new ExtractedIonChromatogram(peaks);
                if (xic.EndRT - xic.StartRT > MaxXicWidth)
                    continue;

                xics.Add(xic);
                featureByXic[xic] = feature;
            }

            LastFeatureByXic = featureByXic;
            return xics;
        }

        private List<TopFdFeature> ApplyFeatureFilters(IEnumerable<TopFdFeature> features)
        {
            return (features ?? Enumerable.Empty<TopFdFeature>())
                .Where(f => f.Mass >= MinMass && f.Mass < MaxMass && f.EcScore >= MinEcScore)
                .ToList();
        }

        /// <summary>
        /// The feature list for <paramref name="scans"/>: <see cref="Features"/>, or, when keyed by source voltage,
        /// the list of the one voltage found in the scan filters.
        /// </summary>
        private List<TopFdFeature> GetFeaturesForScans(MsDataScan[] scans)
        {
            if (FeaturesBySourceVoltage == null)
                return Features;

            var voltages = scans
                .Where(s => s != null && !string.IsNullOrEmpty(s.ScanFilter))
                .Select(s => SourceVoltagePattern.Match(s.ScanFilter))
                .Where(m => m.Success)
                .Select(m => int.Parse(m.Groups[1].Value))
                .Distinct()
                .OrderBy(v => v)
                .ToList();

            if (voltages.Count != 1)
                throw new MetaMorpheusException("Expected scans from exactly one ISD source voltage (sid= in the scan filter) but found "
                    + (voltages.Count == 0 ? "none" : string.Join(", ", voltages)) + ".");

            if (!FeaturesBySourceVoltage.TryGetValue(voltages[0], out var features))
                throw new MetaMorpheusException("No feature xml was supplied for ISD source voltage " + voltages[0] + ".");

            return features;
        }

        /// <summary>
        /// One XIC point per scan in which the feature was observed. The envelope intensities of all charge states
        /// in that scan are summed, and the charge state with the most intense envelope in that scan is kept as the
        /// point's charge. Every point sits at the feature's monoisotopic mass, so the XIC follows one feature
        /// rather than per-scan deconvoluted masses.
        /// </summary>
        private static List<IIndexedPeak> BuildPeaks(TopFdFeature feature, MsDataScan[] scans, Dictionary<int, int> scanIndexByScanNumber)
        {
            var byScanIndex = new Dictionary<int, (double Intensity, int Charge, double ChargeIntensity)>();
            foreach (var chargeFeature in feature.SingleChargeFeatures)
            {
                foreach (var envelope in chargeFeature.Envelopes)
                {
                    if (!scanIndexByScanNumber.TryGetValue(envelope.OneBasedScanNumber, out int scanIndex))
                        continue;

                    double intensity = envelope.IntensitySum;
                    if (byScanIndex.TryGetValue(scanIndex, out var existing))
                    {
                        // add to the total; replace the representative charge only if this envelope is more intense
                        byScanIndex[scanIndex] = intensity > existing.ChargeIntensity
                            ? (existing.Intensity + intensity, chargeFeature.Charge, intensity)
                            : (existing.Intensity + intensity, existing.Charge, existing.ChargeIntensity);
                    }
                    else
                    {
                        byScanIndex[scanIndex] = (intensity, chargeFeature.Charge, intensity);
                    }
                }
            }

            var peaks = new List<IIndexedPeak>(byScanIndex.Count);
            foreach (var kvp in byScanIndex.OrderBy(k => k.Key))
            {
                int charge = kvp.Value.Charge != 0 ? kvp.Value.Charge : 1;
                // a one-peak envelope, so that IndexedMass carries the feature's mass, charge and intensity
                var isotopicEnvelope = new IsotopicEnvelope(
                    new List<(double mz, double intensity)> { (feature.Mass.ToMz(charge), kvp.Value.Intensity) },
                    feature.Mass, charge, kvp.Value.Intensity, 0);
                peaks.Add(new IndexedMass(isotopicEnvelope, scans[kvp.Key].RetentionTime, kvp.Key, 1));
            }
            return peaks;
        }
    }
}
