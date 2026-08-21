using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// A <see cref="XicConstructor"/> that reads XICs straight out of TopDIA's <c>*.mzrt.csv</c> files.
    ///
    /// Of the four external-decon routes this is the most direct. The msAlign route rebuilds chromatograms by
    /// tracing per-scan masses; the feature.xml route reads TopFD's per-charge envelope lists and maps their scan
    /// numbers onto the passed scan array. TopDIA's mzrt already stores a per-CYCLE intensity vector of identical
    /// length in the MS1 file and every MS2 window file, so the chromatogram needs neither tracing nor a join -
    /// cycle <c>k</c> of the vector IS position <c>k</c> of the scan array, for precursors and fragments alike.
    ///
    /// That direct indexing is also the one assumption this class makes, so it is checked: the XIC vector length
    /// must equal the number of scans handed in, and the retention times are cross-checked against those scans
    /// (see <see cref="LastMaxRetentionTimeDiscrepancy"/>). Both fail loudly rather than silently producing
    /// chromatograms on the wrong time axis.
    ///
    /// One mzrt file covers one DIA window. For ISD, each source voltage was given its own synthetic window, so
    /// use <see cref="BySourceVoltage"/> to map voltage to file. Beware that TopDIA names MS2 files by the window's
    /// LOWER BOUND, which need not order the same way as the voltages - pass the mapping explicitly rather than
    /// inferring it from file names.
    /// </summary>
    public class MzrtXicConstructor : XicConstructor
    {
        /// <summary>Rows for the current file, one per charge state. Used when
        /// <see cref="RecordsBySourceVoltage"/> is not set.</summary>
        public List<TopDiaMzrtRecord> Records { get; set; }

        /// <summary>Length of the XIC vector, i.e. the number of cycles in the run.</summary>
        public int CycleCount { get; set; }

        /// <summary>Optional: records keyed by ISD source voltage, so one constructor can serve every fragment
        /// channel from <see cref="DIAparameters.Ms2XicConstructor"/>. The voltage is read off the sid= in the
        /// scan filter of whichever channel is passed in.</summary>
        public Dictionary<int, (List<TopDiaMzrtRecord> Records, int CycleCount)> RecordsBySourceVoltage { get; set; }

        /// <summary>ISD source voltage used by the most recent call, or 0 when voltage keying is off.</summary>
        public int LastSourceVoltage { get; private set; }

        /// <summary>Minimum / maximum neutral monoisotopic mass to keep.</summary>
        public double MinMass { get; set; }
        public double MaxMass { get; set; }

        /// <summary>Minimum TopDIA feature score. The files carry a large low-score tail, so this is the main
        /// quality lever - on the fragment windows roughly 40-60% of rows score below 0.5.</summary>
        public double MinScore { get; set; }

        /// <summary>Minimum number of charge states a feature must show. Fragment features are usually observed
        /// at one charge, so leave this at 1 for fragment channels.</summary>
        public int MinChargeCount { get; set; }

        /// <summary>If true (default), sum a feature's charge states into ONE XIC per feature id. If false, emit
        /// one XIC per charge state.</summary>
        public bool AggregateCharges { get; set; } = true;

        /// <summary>Diagnostics from the most recent <see cref="GetAllXics"/> call.</summary>
        public int LastRecordCount { get; private set; }
        public int LastFeatureCount { get; private set; }
        public int LastXicCount { get; private set; }
        public int LastBelowMinPeaksCount { get; private set; }

        /// <summary>
        /// Disagreement, in minutes, between each record's reported elution start (<c>rtLo</c>) and the retention
        /// time of the scan its start cycle (<c>specLow</c>) indexes — a direct check that cycle k really is
        /// scans[k].
        ///
        /// Judge the axis on the MEDIAN, not the max: a max over tens of thousands of rows is set by a single
        /// outlier and says nothing about the mapping. This deliberately uses the elution START rather than the
        /// apex, because the apex depends on how TopDIA defines it and can legitimately sit cycles away from the
        /// tallest point of a broad or multi-modal trace.
        /// </summary>
        public double LastMedianRetentionTimeDiscrepancy { get; private set; }
        public double LastMaxRetentionTimeDiscrepancy { get; private set; }

        /// <summary>The record each emitted XIC came from (the most intense charge state when aggregating).</summary>
        public Dictionary<ExtractedIonChromatogram, TopDiaMzrtRecord> LastRecordByXic { get; private set; }
            = new Dictionary<ExtractedIonChromatogram, TopDiaMzrtRecord>();

        /// <param name="peakFindingTolerance">Unused - nothing is matched by mass here - but kept so this class
        /// stays interchangeable with the other constructors and its settings reach the settings dump.</param>
        public MzrtXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth,
            int minNumberOfPeaks, List<TopDiaMzrtRecord> records, int cycleCount, double minMass = 0,
            double maxMass = double.MaxValue, double minScore = 0, int minChargeCount = 1,
            bool aggregateCharges = true, XicSpline xicSpline = null)
            : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            MinMass = minMass;
            MaxMass = maxMass;
            MinScore = minScore;
            MinChargeCount = minChargeCount;
            AggregateCharges = aggregateCharges;
            CycleCount = cycleCount;
            Records = Filter(records);
        }

        /// <summary>Convenience overload that reads the mzrt file itself.</summary>
        public MzrtXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth,
            int minNumberOfPeaks, string mzrtFilePath, double minMass = 0, double maxMass = double.MaxValue,
            double minScore = 0, int minChargeCount = 1, bool aggregateCharges = true, XicSpline xicSpline = null)
            : this(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks,
                ReadWithCycleCount(mzrtFilePath, out int cycles), cycles, minMass, maxMass, minScore,
                minChargeCount, aggregateCharges, xicSpline)
        {
        }

        private static List<TopDiaMzrtRecord> ReadWithCycleCount(string path, out int cycleCount) =>
            TopDiaMzrtFile.ReadFile(path, out cycleCount);

        private List<TopDiaMzrtRecord> Filter(List<TopDiaMzrtRecord> records)
        {
            if (records == null)
                return new List<TopDiaMzrtRecord>();

            var kept = records.Where(r => r.Mass >= MinMass && r.Mass < MaxMass && r.Score >= MinScore).ToList();
            if (MinChargeCount > 1)
            {
                var chargesPerFeature = kept.GroupBy(r => r.Id).ToDictionary(g => g.Key, g => g.Select(r => r.Charge).Distinct().Count());
                kept = kept.Where(r => chargesPerFeature[r.Id] >= MinChargeCount).ToList();
            }
            return kept;
        }

        /// <summary>
        /// Builds one constructor covering every ISD fragment channel. <paramref name="mzrtPathsByVoltage"/> maps
        /// source voltage to that channel's mzrt file - supply it explicitly, because TopDIA names the files by
        /// window lower bound and that ordering can invert relative to voltage.
        /// </summary>
        public static MzrtXicConstructor BySourceVoltage(Tolerance peakFindingTolerance, int maxMissedScansAllowed,
            double maxPeakHalfWidth, int minNumberOfPeaks, Dictionary<int, string> mzrtPathsByVoltage,
            double minMass = 0, double maxMass = double.MaxValue, double minScore = 0, int minChargeCount = 1,
            bool aggregateCharges = true, XicSpline xicSpline = null)
        {
            if (mzrtPathsByVoltage == null || mzrtPathsByVoltage.Count == 0)
                throw new MetaMorpheusException("No per-voltage TopDIA mzrt files were supplied.");

            MzrtXicConstructor constructor = null;
            var byVoltage = new Dictionary<int, (List<TopDiaMzrtRecord>, int)>();
            foreach (var kvp in mzrtPathsByVoltage.OrderBy(k => k.Key))
            {
                var records = TopDiaMzrtFile.ReadFile(kvp.Value, out int cycles);
                constructor ??= new MzrtXicConstructor(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth,
                    minNumberOfPeaks, records, cycles, minMass, maxMass, minScore, minChargeCount, aggregateCharges,
                    xicSpline);
                byVoltage[kvp.Key] = (constructor.Filter(records), cycles);
            }

            constructor.RecordsBySourceVoltage = byVoltage;
            return constructor;
        }

        /// <summary>Picks the records for these scans: the channel's when voltage keyed, otherwise
        /// <see cref="Records"/>.</summary>
        private List<TopDiaMzrtRecord> GetRecordsForScans(MsDataScan[] scans)
        {
            LastSourceVoltage = 0;
            if (RecordsBySourceVoltage == null || RecordsBySourceVoltage.Count == 0)
                return Records;

            var voltages = scans
                .Where(s => s != null && !string.IsNullOrEmpty(s.ScanFilter))
                .Select(s => ExternalDeconResults.SourceVoltagePattern.Match(s.ScanFilter))
                .Where(m => m.Success)
                .Select(m => int.Parse(m.Groups[1].Value))
                .Distinct()
                .ToList();

            if (voltages.Count == 0)
                return Records;
            if (voltages.Count > 1)
                throw new MetaMorpheusException(
                    "XIC construction was handed scans from more than one ISD source voltage (" +
                    string.Join(", ", voltages.OrderBy(v => v)) + "). Each voltage must be built separately.");

            LastSourceVoltage = voltages[0];
            if (!RecordsBySourceVoltage.TryGetValue(LastSourceVoltage, out var entry))
                throw new MetaMorpheusException(
                    "No TopDIA mzrt file was supplied for ISD source voltage " + LastSourceVoltage +
                    ". Voltages available: " + string.Join(", ", RecordsBySourceVoltage.Keys.OrderBy(v => v)) + ".");

            CycleCount = entry.CycleCount;
            return entry.Records;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            Records = GetRecordsForScans(scans);
            if (Records == null || Records.Count == 0)
                throw new MetaMorpheusException("No TopDIA mzrt records were supplied to build XICs from.");

            // The whole design rests on cycle k == scans[k]; if the lengths differ that identity is broken and
            // every chromatogram would land on the wrong time axis, so refuse rather than guess.
            if (CycleCount != scans.Length)
                throw new MetaMorpheusException(
                    "TopDIA XIC vectors hold " + CycleCount + " cycles but " + scans.Length + " scans were passed" +
                    (LastSourceVoltage > 0 ? " for sid=" + LastSourceVoltage : "") +
                    ". The mzrt file and the scan array must cover the same cycles.");

            LastRecordCount = Records.Count;
            LastBelowMinPeaksCount = 0;
            LastMaxRetentionTimeDiscrepancy = 0;
            LastMedianRetentionTimeDiscrepancy = 0;
            _retentionTimeDiscrepancies = new List<double>(Records.Count);
            LastRecordByXic = new Dictionary<ExtractedIonChromatogram, TopDiaMzrtRecord>();

            var groups = AggregateCharges
                ? Records.GroupBy(r => r.Id).Select(g => g.ToList()).ToList()
                : Records.Select(r => new List<TopDiaMzrtRecord> { r }).ToList();
            LastFeatureCount = groups.Count;

            var xics = new List<ExtractedIonChromatogram>(groups.Count);
            foreach (var group in groups)
            {
                if (isolationRange != null && !group.Any(r => isolationRange.Contains(r.MonoisotopicMz)))
                    continue;

                var peaks = BuildPeaks(group, scans);
                if (peaks == null || peaks.Count < MinNumberOfPeaks)
                {
                    LastBelowMinPeaksCount++;
                    continue;
                }

                var xic = new ExtractedIonChromatogram(peaks);
                xics.Add(xic);
                LastRecordByXic[xic] = group.OrderByDescending(r => r.Intensity).First();
            }

            if (_retentionTimeDiscrepancies.Count > 0)
            {
                _retentionTimeDiscrepancies.Sort();
                LastMedianRetentionTimeDiscrepancy = _retentionTimeDiscrepancies[_retentionTimeDiscrepancies.Count / 2];
                LastMaxRetentionTimeDiscrepancy = _retentionTimeDiscrepancies[^1];
            }

            LastXicCount = xics.Count;
            return xics;
        }

        private List<double> _retentionTimeDiscrepancies = new List<double>();

        /// <summary>
        /// Turns a feature's rows into XIC points at the feature mass. Charge states are summed cycle by cycle
        /// when aggregating, keeping the most intense charge of each cycle as the representative so the
        /// downstream m/z conversion is sensible.
        /// </summary>
        private List<IIndexedPeak> BuildPeaks(List<TopDiaMzrtRecord> group, MsDataScan[] scans)
        {
            var byCycle = new Dictionary<int, (double Intensity, int Charge, float ChargeIntensity)>();
            foreach (var record in group)
            {
                for (int i = 0; i < record.Intensities.Length; i++)
                {
                    float intensity = record.Intensities[i];
                    if (intensity <= 0)
                        continue; // interior zeros are real gaps, not points

                    int cycle = record.FirstCycle + i;
                    if (cycle < 0 || cycle >= scans.Length || scans[cycle] == null)
                        continue;

                    if (byCycle.TryGetValue(cycle, out var existing))
                    {
                        byCycle[cycle] = intensity > existing.ChargeIntensity
                            ? (existing.Intensity + intensity, record.Charge, intensity)
                            : (existing.Intensity + intensity, existing.Charge, existing.ChargeIntensity);
                    }
                    else
                    {
                        byCycle[cycle] = (intensity, record.Charge, intensity);
                    }
                }

                // Cross-check the cycle axis against the scans using the elution START: specLow is a cycle index
                // and rtLo is the retention time TopDIA assigns it, so scans[specLow] must report that same time.
                // Deliberately not the apex - argmax of the trace and TopDIA's rtApex legitimately differ by a
                // few cycles on broad peaks, which would make this look broken when it is not.
                if (record.SpectrumLow >= 0 && record.SpectrumLow < scans.Length
                    && scans[record.SpectrumLow] != null && record.RetentionTimeBegin > 0)
                {
                    _retentionTimeDiscrepancies.Add(
                        Math.Abs(scans[record.SpectrumLow].RetentionTime - record.RetentionTimeBegin));
                }
            }

            if (byCycle.Count == 0)
                return null;

            double mass = group.OrderByDescending(r => r.Intensity).First().Mass;
            var peaks = new List<IIndexedPeak>(byCycle.Count);
            foreach (var kvp in byCycle.OrderBy(k => k.Key))
            {
                int charge = kvp.Value.Charge != 0 ? kvp.Value.Charge : 1;
                var envelope = new IsotopicEnvelope(
                    new List<(double mz, double intensity)> { (mass.ToMz(charge), kvp.Value.Intensity) },
                    mass, charge, kvp.Value.Intensity, 0);
                peaks.Add(new IndexedMass(envelope, scans[kvp.Key].RetentionTime, kvp.Key, 1));
            }

            return peaks;
        }
    }
}
