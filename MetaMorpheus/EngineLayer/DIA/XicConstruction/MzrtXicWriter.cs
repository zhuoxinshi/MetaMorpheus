using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer.DIA
{
    /// <summary>
    /// Writes precursor and fragment XICs to a TopDIA-style <c>.mzrt.csv</c> so the co-elution ("alignment") of a
    /// precursor with its fragments can be visualised the way TopDIA does it.
    ///
    /// The idea borrowed from TopDIA (see <see cref="TopDiaMzrtFile"/>): each XIC is written as ONE row whose
    /// <c>XIC</c> column is a DENSE per-cycle intensity vector - one value per MS1 cycle of the run, zeros where the
    /// feature is absent. Because every row is written on the SAME cycle axis (the same length), a precursor and its
    /// fragments align index-for-index with no scan-number join; plotting each row's XIC vector against the cycle (or
    /// its retention time) overlays the traces directly.
    ///
    /// The columns are TopDIA's (<c>ID,Fraction_ID,Envelope_num,Mass,MonoMz,Charge,Intensity,mzLo,mzHi,rtLo,rtHi,</c>
    /// <c>specLow,specHi,rtApex,Score,...,XIC,Envelope</c>) so the file round-trips through <see cref="TopDiaMzrtFile"/>,
    /// plus three extra columns - <c>role</c> (precursor/fragment), <c>correlation</c> (to the precursor) and
    /// <c>fragment_rank</c> - that carry OUR precursor-fragment grouping. <c>ID</c> is the group index, so every row
    /// sharing an id is one precursor-fragment group (analogous to TopDIA nesting charge states under one id).
    /// </summary>
    public static class MzrtXicWriter
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

        /// <summary>
        /// Write the precursor and fragment XICs of <paramref name="groups"/> to <paramref name="outputPath"/>.
        /// <paramref name="maxGroups"/> caps how many groups are written (0 = all). <paramref name="cycleCount"/> sets
        /// the XIC vector length (0 = auto: one past the largest cycle index seen); pass the run's MS1 cycle count to
        /// match TopDIA's full-length vectors exactly. Returns the output path.
        /// </summary>
        public static string Write(IEnumerable<PrecursorFragmentsGroup> groups, string outputPath,
            int maxGroups = 0, int cycleCount = 0)
        {
            if (groups == null) throw new ArgumentNullException(nameof(groups));
            if (string.IsNullOrWhiteSpace(outputPath)) throw new ArgumentException("OutputPath is not set.");

            // First pass: collect the rows and find the cycle axis (dense vectors must all share one length).
            var rows = new List<Row>();
            int written = 0, maxCycle = 0;
            foreach (var g in groups)
            {
                if (g?.PrecursorXic == null) continue;
                if (maxGroups > 0 && written >= maxGroups) break;
                written++;

                AddRow(rows, ref maxCycle, g.PFgroupIndex, "precursor", g.PrecursorXic, null, null);
                foreach (var pf in g.PFpairs ?? Enumerable.Empty<PrecursorFragmentPair>())
                    if (pf?.FragmentXic != null)
                        AddRow(rows, ref maxCycle, g.PFgroupIndex, "fragment", pf.FragmentXic, pf.Correlation, pf.FragmentRank);
            }

            int cycles = cycleCount > 0 ? cycleCount : maxCycle + 1;

            var dir = Path.GetDirectoryName(outputPath);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir)) Directory.CreateDirectory(dir);

            using var sw = new StreamWriter(outputPath, append: false);
            sw.WriteLine("ID,Fraction_ID,Envelope_num,Mass,MonoMz,Charge,Intensity,mzLo,mzHi,rtLo,rtHi," +
                         "specLow,specHi,rtApex,Score,role,correlation,fragment_rank,XIC,Envelope");

            foreach (var r in rows)
                sw.WriteLine(FormatRow(r, cycles));

            return outputPath;
        }

        private static void AddRow(List<Row> rows, ref int maxCycle, int groupId, string role,
            ExtractedIonChromatogram xic, double? correlation, int? fragmentRank)
        {
            if (xic?.Peaks == null || xic.Peaks.Count == 0) return;
            if (xic.EndScanIndex > maxCycle) maxCycle = xic.EndScanIndex;
            rows.Add(new Row { GroupId = groupId, Role = role, Correlation = correlation, FragmentRank = fragmentRank, Xic = xic });
        }

        private static string FormatRow(Row r, int cycles)
        {
            var xic = r.Xic;
            int charge = (xic.ApexPeak as IndexedMass)?.Charge ?? 0;
            double m = xic.ApexPeak?.M ?? xic.AveragedMassOrMz;     // neutral mass (mass XICs) or m/z (mz-peak XICs)
            double monoMz = charge > 0 ? m.ToMz(charge) : m;
            double intensity = xic.Peaks.Sum(p => (double)p.Intensity);

            // dense per-cycle intensity vector: intensity at each peak's cycle index, 0 elsewhere
            var vec = new float[cycles];
            foreach (var p in xic.Peaks)
                if (p.ZeroBasedScanIndex >= 0 && p.ZeroBasedScanIndex < cycles)
                    vec[p.ZeroBasedScanIndex] = p.Intensity;

            var xicSb = new StringBuilder(cycles * 3);
            for (int i = 0; i < cycles; i++)
            {
                if (i > 0) xicSb.Append(';');
                xicSb.Append(vec[i] == 0f ? "0" : vec[i].ToString(Ci));
            }

            // We do not carry the isotope envelope on an XIC, so Envelope is a single apex mz&intensity pair
            // (the reader does not require this column; it is written so the schema is complete).
            string envelope = monoMz.ToString("F4", Ci) + "&" +
                              (xic.ApexPeak != null ? ((double)xic.ApexPeak.Intensity).ToString(Ci) : "0");

            string corr = r.Correlation.HasValue && !double.IsNaN(r.Correlation.Value)
                ? r.Correlation.Value.ToString("F4", Ci) : "";
            string rank = r.FragmentRank?.ToString(Ci) ?? "";

            return string.Join(",",
                r.GroupId.ToString(Ci),                               // ID (= group index; rows share it within a PF group)
                "0",                                                  // Fraction_ID
                xic.Peaks.Count.ToString(Ci),                         // Envelope_num (points in the trace)
                m.ToString("F5", Ci),                                 // Mass
                monoMz.ToString("F5", Ci),                            // MonoMz
                (charge > 0 ? charge : 1).ToString(Ci),               // Charge (1 when the peaks carry no charge)
                intensity.ToString(Ci),                               // Intensity (summed over the trace)
                monoMz.ToString("F5", Ci),                            // mzLo (envelope span not tracked -> mono m/z)
                monoMz.ToString("F5", Ci),                            // mzHi
                xic.StartRT.ToString("F2", Ci),                       // rtLo (minutes)
                xic.EndRT.ToString("F2", Ci),                         // rtHi
                xic.StartScanIndex.ToString(Ci),                      // specLow (first cycle of the elution)
                xic.EndScanIndex.ToString(Ci),                        // specHi (last cycle)
                xic.ApexRT.ToString("F2", Ci),                        // rtApex
                "1.00",                                               // Score (placeholder; our metric is `correlation`)
                r.Role,                                               // role: precursor / fragment
                corr,                                                 // correlation to the precursor (fragments only)
                rank,                                                 // fragment_rank
                xicSb.ToString(),                                     // XIC (dense per-cycle vector)
                envelope);                                            // Envelope
        }

        private sealed class Row
        {
            public int GroupId;
            public string Role;
            public double? Correlation;
            public int? FragmentRank;
            public ExtractedIonChromatogram Xic;
        }
    }
}
