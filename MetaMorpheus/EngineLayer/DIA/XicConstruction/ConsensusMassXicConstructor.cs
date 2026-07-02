using Chemistry;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using MzLibUtil;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    /// <summary>
    /// A <see cref="XicConstructor"/> that builds XICs from the mzLib CONSENSUS mass tracer
    /// (<see cref="MassTraceBuilder"/> → <see cref="TraceCorrector"/> → <see cref="MassFeatureBuilder"/>) —
    /// the same "consensus feature tracing" the project's deconscan pipeline uses — instead of
    /// MetaMorpheus's per-scan neutral-mass indexing (<see cref="NeutralMassXicConstructor"/>).
    ///
    /// Pipeline per call: deconvolute every scan → build charge-locked mass traces → off-by-one
    /// correct each trace → stitch per-charge traces into cross-charge <see cref="MassFeature"/>s.
    /// Each MassFeature becomes ONE <see cref="ExtractedIonChromatogram"/> whose per-scan
    /// <see cref="IndexedMass"/> peaks carry the feature's consensus mass and summed intensity, so the
    /// existing precursor–fragment grouping (<see cref="XicGroupingEngine"/>) and pseudo-MS2 generation
    /// (<see cref="PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup"/>) consume them unchanged.
    ///
    /// Peaks are indexed by their 0-based position within the scan array passed to <see cref="GetAllXics"/>.
    /// For a strict ISD cycle each voltage channel has one scan per cycle, so cycle k has the same scan
    /// index in every channel and precursor↔fragment XICs correlate on a common axis.
    /// </summary>
    public class ConsensusMassXicConstructor : XicConstructor
    {
        public DeconvolutionParameters DeconParameters { get; set; }
        /// <summary>Da tolerance for linking envelopes into a charge-locked trace across scans.</summary>
        public double TraceToleranceDa { get; set; }
        /// <summary>Max consecutive missed scans before a trace is closed.</summary>
        public int MaxGap { get; set; }
        /// <summary>ppm tolerance for stitching per-charge traces into a cross-charge feature.</summary>
        public double FeatureMassPpm { get; set; }
        /// <summary>Minimum consensus (neutral) mass to keep a feature. Use e.g. 3000 on the precursor
        /// channel to keep only intact proteoforms and drop small noise; leave 0 on the fragment channel.</summary>
        public double MinMass { get; set; }
        /// <summary>Minimum number of distinct charge states a feature must be observed at (charge
        /// multiplicity as a confidence signal). Use e.g. 3 on the precursor channel; 1 on the fragment channel.</summary>
        public int MinChargeCount { get; set; }
        /// <summary>If true (default), cross-charge-aggregate corrected traces into <see cref="MassFeature"/>s
        /// (one XIC per feature). If false, emit one XIC per charge-locked <see cref="CorrectedTrace"/> — i.e.
        /// mass tracing only, NO charge aggregation. Use false on the fragment channel to keep each fragment
        /// charge state as its own trace.</summary>
        public bool AggregateCharges { get; set; } = true;

        public ConsensusMassXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, DeconvolutionParameters deconParameters,
            double traceToleranceDa = 0.02, int maxGap = 1, double featureMassPpm = 15.0, double minMass = 0,
            int minChargeCount = 1, XicSpline? xicSpline = null)
            : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            DeconParameters = deconParameters;
            TraceToleranceDa = traceToleranceDa;
            MaxGap = maxGap;
            FeatureMassPpm = featureMassPpm;
            MinMass = minMass;
            MinChargeCount = minChargeCount;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            // 1. deconvolute every scan (parallel; order preserved by index -> aligns 1:1 with scans)
            var perScan = new IReadOnlyList<IsotopicEnvelope>[scans.Length];
            Parallel.For(0, scans.Length, i =>
            {
                perScan[i] = Deconvoluter.Deconvolute(scans[i], DeconParameters).ToList();
            });

            // 2. consensus mass tracing: charge-locked traces -> off-by-one correction (-> cross-charge features)
            var traces = MassTraceBuilder.BuildTraces(scans, perScan, TraceToleranceDa, MaxGap);
            var corrected = traces.SelectMany(t => TraceCorrector.Correct(t)).ToList();

            var xics = new List<ExtractedIonChromatogram>();

            // 3b. mass tracing only (no charge aggregation): one XIC per charge-locked corrected trace.
            if (!AggregateCharges)
            {
                foreach (var trace in corrected)
                {
                    if (trace.Envelopes.Count < MinNumberOfPeaks) continue;
                    if (trace.ConsensusMass < MinMass) continue;
                    int z = trace.Charge != 0 ? trace.Charge : 1;
                    double tmass = trace.ConsensusMass;
                    var tpeaks = new List<IIndexedPeak>();
                    foreach (var e in trace.Envelopes)
                    {
                        var env = new IsotopicEnvelope(
                            new List<(double mz, double intensity)> { (tmass.ToMz(z), e.Intensity) },
                            tmass, z, e.Intensity, 0);
                        tpeaks.Add(new IndexedMass(env, e.RT, e.ScanIndex, 1));
                    }
                    if (tpeaks.Count < MinNumberOfPeaks) continue;
                    xics.Add(new ExtractedIonChromatogram(tpeaks));
                }
                return xics;
            }

            var features = MassFeatureBuilder.BuildFeatures(corrected, FeatureMassPpm);

            // 3. each cross-charge MassFeature -> one XIC of per-scan IndexedMass peaks at the consensus mass
            foreach (var feature in features)
            {
                feature.Finalise();
                if (feature.MaxTraceLength < MinNumberOfPeaks) continue;
                if (feature.ConsensusMass < MinMass) continue;
                if (feature.ChargeCount < MinChargeCount) continue;

                double mass = feature.ConsensusMass;
                var peaks = new List<IIndexedPeak>();
                // one time point per scan index the feature was observed in; sum intensity over charges,
                // keep the dominant charge so downstream m/z conversion is representative
                foreach (var scanGroup in feature.Traces
                             .SelectMany(t => t.Envelopes)
                             .GroupBy(e => e.ScanIndex))
                {
                    double intensity = scanGroup.Sum(e => e.Intensity);
                    var dominant = scanGroup.OrderByDescending(e => e.Intensity).First();
                    int charge = dominant.Charge != 0 ? dominant.Charge : 1;
                    double mz = mass.ToMz(charge);
                    var env = new IsotopicEnvelope(
                        new List<(double mz, double intensity)> { (mz, intensity) },
                        mass, charge, intensity, 0);
                    peaks.Add(new IndexedMass(env, dominant.RT, scanGroup.Key, 1));
                }

                if (peaks.Count < MinNumberOfPeaks) continue;
                xics.Add(new ExtractedIonChromatogram(peaks));
            }
            return xics;
        }
    }
}
