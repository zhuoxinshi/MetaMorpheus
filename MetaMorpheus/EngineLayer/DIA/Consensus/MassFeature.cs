using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// A cross-charge-state consensus feature: a group of
    /// <see cref="CorrectedTrace"/> entries (each at one charge state) whose
    /// consensus masses agree within a ppm tolerance and whose RT windows
    /// overlap.
    ///
    /// A real proteoform or peptide produces envelopes at several charge
    /// states (BU peptides commonly +2/+3; TD proteoforms across +10..+15
    /// or wider). The trace builder is charge-locked, so each charge gets
    /// its own trace. <see cref="MassFeatureBuilder"/> stitches those
    /// per-charge traces back together. Charge multiplicity then becomes
    /// a confidence signal: a feature seen at multiple charges is
    /// corroborated by independent charge calculations; a single-charge,
    /// single-envelope feature is more likely to be noise.
    ///
    /// Mutation policy: callers append to <see cref="Traces"/> during
    /// construction, then call <see cref="Finalise"/> exactly once to
    /// derive the aggregate fields. After Finalise, the feature is
    /// treated as read-only by downstream consumers.
    /// </summary>
    public sealed class MassFeature
    {
        public int Id;
        public List<CorrectedTrace> Traces = new();
        public double ConsensusMass;
        public HashSet<int> Charges = new();
        public int ChargeCount => Charges.Count;
        public int MaxTraceLength;
        public double RTStart;
        public double RTEnd;

        /// <summary>
        /// Retention time (minutes) of the feature apex scan — the scan index where intensity summed across ALL
        /// charge traces is highest (the same apex scan used for <see cref="ApexIntensityTop3"/>). This is a true
        /// apex, unlike the RT-range midpoint written as "Apex_time" in the `_ms1.feature` file.
        /// </summary>
        public double ApexRT;

        public double SummedIntensity;

        /// <summary>
        /// Summed intensity of the three most abundant charge-state traces: the top-three
        /// <see cref="CorrectedTrace"/>s ranked by <see cref="CorrectedTrace.TotalIntensity"/>,
        /// each summed over its own scans, then added together. A feature with fewer than three
        /// traces uses all it has. This is the noise-suppressed counterpart of
        /// <see cref="SummedIntensity"/> (which sums every charge, including weak/spurious ones).
        /// </summary>
        public double SummedIntensityTop3;

        /// <summary>
        /// Apex intensity of the three most abundant charge-state traces, read at a single shared
        /// apex scan. The feature apex scan is the scan index where intensity summed across ALL
        /// charge traces is greatest; this value is the sum of the top-three traces' intensities
        /// at that one scan (co-eluting charges measured at the same retention time). A top-three
        /// trace with no envelope at the apex scan contributes 0. Fewer than three traces uses all
        /// available.
        /// </summary>
        public double ApexIntensityTop3;

        /// <summary>
        /// Populate derived fields from the current <see cref="Traces"/>
        /// list. Idempotent; safe to re-run after the trace list changes.
        /// <see cref="ConsensusMass"/> is the intensity-weighted mean of
        /// per-trace consensus masses (heavier traces contribute more).
        /// </summary>
        public void Finalise()
        {
            if (Traces.Count == 0)
                throw new System.InvalidOperationException("Cannot Finalise a MassFeature with no traces.");

            Charges = new HashSet<int>(Traces.Select(t => t.Charge));
            SummedIntensity = Traces.Sum(t => t.TotalIntensity);
            RTStart = Traces.Min(t => t.FirstRT);
            RTEnd = Traces.Max(t => t.LastRT);
            MaxTraceLength = Traces.Max(t => t.Envelopes.Count);

            // Intensity-weighted mean of per-trace consensus masses. Heavier
            // traces (more intense, more confident) contribute more.
            double w = Traces.Sum(t => t.TotalIntensity);
            ConsensusMass = w == 0
                ? Traces.Average(t => t.ConsensusMass)
                : Traces.Sum(t => t.ConsensusMass * t.TotalIntensity) / w;

            // --- Top-3 charge-state quantities ------------------------------------------------
            // Rank the charge-state traces by total (summed-over-scans) intensity and keep the
            // three most abundant. Ties are broken deterministically (charge, then mass) so the
            // metrics are reproducible run-to-run on identical input.
            var top3 = Traces
                .OrderByDescending(t => t.TotalIntensity)
                .ThenBy(t => t.Charge)
                .ThenBy(t => t.ConsensusMass)
                .Take(3)
                .ToList();

            // (1) Summed intensity: each top-3 trace's sum-over-scans, added together.
            SummedIntensityTop3 = top3.Sum(t => t.TotalIntensity);

            // (2) Apex intensity at a shared apex scan. First locate the feature apex scan: the
            // scan index where intensity summed across ALL charge traces is greatest.
            var perScanTotal = new Dictionary<int, double>();
            foreach (var t in Traces)
                foreach (var e in t.Envelopes)
                {
                    perScanTotal.TryGetValue(e.ScanIndex, out double cur);
                    perScanTotal[e.ScanIndex] = cur + e.Intensity;
                }

            int apexScanIndex = perScanTotal
                .OrderByDescending(kv => kv.Value)
                .ThenBy(kv => kv.Key)          // deterministic tie-break on scan index
                .First().Key;

            // Retention time of that apex scan (every trace's envelope at the scan carries the same RT).
            ApexRT = Traces.SelectMany(t => t.Envelopes).First(e => e.ScanIndex == apexScanIndex).RT;

            // Sum the top-3 traces' intensities read at that single apex scan (a trace with no
            // envelope at the apex scan contributes 0).
            ApexIntensityTop3 = top3.Sum(t =>
                t.Envelopes.FirstOrDefault(e => e.ScanIndex == apexScanIndex)?.Intensity ?? 0.0);
        }
    }
}
