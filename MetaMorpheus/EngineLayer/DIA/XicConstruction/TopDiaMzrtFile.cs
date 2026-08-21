using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// One row of a TopDIA <c>*_frac_ms1.mzrt.csv</c> / <c>*_&lt;window&gt;_frac_ms2.mzrt.csv</c> file: a single
    /// CHARGE STATE of a feature, together with its extracted ion chromatogram.
    ///
    /// <see cref="Id"/> groups charge states into one feature, the same way TopFD's feature xml nests
    /// <c>single_charge_feature</c> inside <c>frac_feature</c> - here the nesting is flattened into repeated rows.
    ///
    /// The <c>XIC</c> column is the reason this format is worth reading: it is a full per-CYCLE intensity vector,
    /// the same length in the MS1 file and every MS2 window file (one entry per MS1 cycle of the run). Precursor
    /// and fragment chromatograms therefore arrive already on a common axis, with no scan-number join needed -
    /// unlike the msAlign and feature.xml routes, which both have to reconstruct that alignment.
    ///
    /// Only the non-zero stretch of the vector is retained (see <see cref="FirstCycle"/> /
    /// <see cref="Intensities"/>); keeping all ~1545 entries for every row would cost hundreds of MB per file for
    /// values that are almost entirely zero.
    /// </summary>
    public class TopDiaMzrtRecord
    {
        /// <summary>Feature id. Rows sharing an id are charge states of the same species.</summary>
        public int Id { get; set; }
        public int FractionId { get; set; }
        public int EnvelopeNumber { get; set; }
        /// <summary>Neutral monoisotopic mass.</summary>
        public double Mass { get; set; }
        /// <summary>Monoisotopic m/z of this charge state.</summary>
        public double MonoisotopicMz { get; set; }
        public int Charge { get; set; }
        /// <summary>Summed intensity TopDIA reports for this charge state.</summary>
        public double Intensity { get; set; }
        public double MzMin { get; set; }
        public double MzMax { get; set; }
        /// <summary>Elution bounds in MINUTES.</summary>
        public double RetentionTimeBegin { get; set; }
        public double RetentionTimeEnd { get; set; }
        public double RetentionTimeApex { get; set; }
        /// <summary>First and last CYCLE index of the elution, as reported in the file. These index the XIC
        /// vector directly.</summary>
        public int SpectrumLow { get; set; }
        public int SpectrumHigh { get; set; }
        /// <summary>TopDIA's feature score, in [0, 1].</summary>
        public double Score { get; set; }

        /// <summary>Cycle index of <see cref="Intensities"/>[0].</summary>
        public int FirstCycle { get; set; }
        /// <summary>The non-zero stretch of the XIC, one entry per cycle starting at <see cref="FirstCycle"/>.
        /// Interior zeros are preserved so the trace keeps its shape.</summary>
        public float[] Intensities { get; set; } = Array.Empty<float>();

        /// <summary>Number of cycles with non-zero signal.</summary>
        public int PointCount => Intensities.Count(v => v > 0);

        /// <summary>Cycle index just past the end of <see cref="Intensities"/>.</summary>
        public int EndCycleExclusive => FirstCycle + Intensities.Length;

        public override string ToString() =>
            "Feature " + Id + " z" + Charge + ": " + Mass.ToString("F4") + " Da, cycles " +
            FirstCycle + "-" + (EndCycleExclusive - 1) + ", RT " + RetentionTimeBegin.ToString("F2") +
            "-" + RetentionTimeEnd.ToString("F2");
    }

    /// <summary>
    /// Reader for TopDIA mzrt.csv files.
    ///
    /// NOTE this is not the same format as mzLib's <c>TopFdMzrt</c>, which describes an older visualization mzrt
    /// (columns <c>color</c> / <c>opacity</c> / <c>promex_score</c>) and carries no chromatogram. The TopDIA
    /// variant read here has <c>specLow</c> / <c>specHi</c> / <c>rtApex</c> / <c>Score</c> / <c>XIC</c> /
    /// <c>Envelope</c> instead, so the two cannot share a reader.
    /// </summary>
    public static class TopDiaMzrtFile
    {
        /// <summary>
        /// Reads every row. <paramref name="cycleCount"/> returns the length of the XIC vector - the number of
        /// cycles in the run - which every row in a file shares and which callers must match their scan array to.
        /// </summary>
        public static List<TopDiaMzrtRecord> ReadFile(string filePath, out int cycleCount)
        {
            if (!File.Exists(filePath))
                throw new MetaMorpheusException("TopDIA mzrt file not found: " + filePath);

            var records = new List<TopDiaMzrtRecord>();
            cycleCount = -1;

            using var reader = new StreamReader(filePath);
            string headerLine = reader.ReadLine();
            if (headerLine == null)
                throw new MetaMorpheusException("TopDIA mzrt file is empty: " + filePath);

            var columns = headerLine.Split(',');
            var index = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
            for (int i = 0; i < columns.Length; i++)
                index[columns[i].Trim()] = i;

            foreach (var required in new[] { "ID", "Mass", "Charge", "specLow", "specHi", "XIC" })
                if (!index.ContainsKey(required))
                    throw new MetaMorpheusException("TopDIA mzrt file " + filePath + " has no '" + required +
                        "' column. Columns present: " + string.Join(", ", columns) + ".");

            string line;
            while ((line = reader.ReadLine()) != null)
            {
                if (string.IsNullOrWhiteSpace(line))
                    continue;

                // no field in this format contains a comma, so a plain split is safe and avoids a CSV parser
                var splits = line.Split(',');
                if (splits.Length < columns.Length)
                    continue;

                var record = new TopDiaMzrtRecord
                {
                    Id = ReadInt(splits, index, "ID"),
                    FractionId = ReadInt(splits, index, "Fraction_ID"),
                    EnvelopeNumber = ReadInt(splits, index, "Envelope_num"),
                    Mass = ReadDouble(splits, index, "Mass"),
                    MonoisotopicMz = ReadDouble(splits, index, "MonoMz"),
                    Charge = ReadInt(splits, index, "Charge"),
                    Intensity = ReadDouble(splits, index, "Intensity"),
                    MzMin = ReadDouble(splits, index, "mzLo"),
                    MzMax = ReadDouble(splits, index, "mzHi"),
                    RetentionTimeBegin = ReadDouble(splits, index, "rtLo"),
                    RetentionTimeEnd = ReadDouble(splits, index, "rtHi"),
                    RetentionTimeApex = ReadDouble(splits, index, "rtApex"),
                    SpectrumLow = ReadInt(splits, index, "specLow"),
                    SpectrumHigh = ReadInt(splits, index, "specHi"),
                    Score = ReadDouble(splits, index, "Score"),
                };

                ParseXic(splits[index["XIC"]], record, ref cycleCount);
                records.Add(record);
            }

            if (records.Count == 0)
                throw new MetaMorpheusException("TopDIA mzrt file has no rows: " + filePath);

            return records;
        }

        /// <summary>Convenience overload for callers that do not need the cycle count.</summary>
        public static List<TopDiaMzrtRecord> ReadFile(string filePath) => ReadFile(filePath, out _);

        /// <summary>
        /// Parses the semicolon-delimited XIC vector, keeping only the stretch between the first and last
        /// non-zero value. Written as two passes over the raw string with no substring allocation - a file holds
        /// ~17k rows x ~1545 values, so splitting each row would churn ~26M strings.
        /// </summary>
        private static void ParseXic(string xic, TopDiaMzrtRecord record, ref int cycleCount)
        {
            if (string.IsNullOrEmpty(xic))
                return;

            var span = xic.AsSpan();

            int count = 1;
            for (int i = 0; i < span.Length; i++)
                if (span[i] == ';') count++;

            if (cycleCount < 0)
                cycleCount = count;
            else if (count != cycleCount)
                throw new MetaMorpheusException(
                    "Inconsistent XIC length in a TopDIA mzrt file: row " + record.Id + " has " + count +
                    " cycles, earlier rows had " + cycleCount + ".");

            var values = new float[count];
            int firstNonZero = -1, lastNonZero = -1, cycle = 0, start = 0;
            for (int i = 0; i <= span.Length; i++)
            {
                if (i != span.Length && span[i] != ';')
                    continue;

                var field = span.Slice(start, i - start);
                if (double.TryParse(field, NumberStyles.Any, CultureInfo.InvariantCulture, out double v) && v > 0)
                {
                    values[cycle] = (float)v;
                    if (firstNonZero < 0) firstNonZero = cycle;
                    lastNonZero = cycle;
                }
                cycle++;
                start = i + 1;
            }

            if (firstNonZero < 0)
            {
                record.FirstCycle = 0;
                record.Intensities = Array.Empty<float>();
                return;
            }

            // keep the non-zero span, interior zeros included so gaps in the elution survive
            int length = lastNonZero - firstNonZero + 1;
            var trimmed = new float[length];
            Array.Copy(values, firstNonZero, trimmed, 0, length);
            record.FirstCycle = firstNonZero;
            record.Intensities = trimmed;
        }

        private static double ReadDouble(string[] splits, Dictionary<string, int> index, string name) =>
            index.TryGetValue(name, out int i) && i < splits.Length &&
            double.TryParse(splits[i], NumberStyles.Any, CultureInfo.InvariantCulture, out double v) ? v : 0;

        private static int ReadInt(string[] splits, Dictionary<string, int> index, string name) =>
            index.TryGetValue(name, out int i) && i < splits.Length &&
            double.TryParse(splits[i], NumberStyles.Any, CultureInfo.InvariantCulture, out double v)
                ? (int)Math.Round(v) : 0;
    }
}
