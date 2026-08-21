using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// One LC-MS feature row from a TopFD / FLASHDeconv <c>*_ms1.feature</c> file: a neutral mass observed over a
    /// retention time window across a range of charge states.
    ///
    /// Read here rather than through mzLib's <c>Ms1Feature</c> because the mass column has been renamed across
    /// TopFD versions - 1.7.x writes <c>Monoisotopic_mass</c>, 1.8.x writes <c>Mass</c> - and mzLib only aliases
    /// the latter, so a 1.7.x file fails to map. This reader resolves every column by a list of accepted header
    /// names, so both schemas (and FLASHDeconv's) load.
    /// </summary>
    public class Ms1FeatureRecord
    {
        public int Id { get; set; }
        public double Mass { get; set; }
        public double Intensity { get; set; }
        public double RetentionTimeBegin { get; set; }
        public double RetentionTimeEnd { get; set; }
        public double RetentionTimeApex { get; set; }
        public double ApexIntensity { get; set; }
        public int ChargeStateMin { get; set; }
        public int ChargeStateMax { get; set; }

        public bool ContainsRetentionTime(double retentionTime, double tolerance = 0) =>
            retentionTime >= RetentionTimeBegin - tolerance && retentionTime <= RetentionTimeEnd + tolerance;

        public bool ContainsCharge(int charge) => charge >= ChargeStateMin && charge <= ChargeStateMax;

        public override string ToString() =>
            "Feature " + Id + ": " + Mass.ToString("F4") + " Da, z " + ChargeStateMin + "-" + ChargeStateMax +
            ", RT " + RetentionTimeBegin.ToString("F2") + "-" + RetentionTimeEnd.ToString("F2");

        // accepted header spellings, most recent first. TopFD 1.8.x / TopFD 1.7.x / FLASHDeconv
        private static readonly string[] IdNames = { "Feature_ID", "ID" };
        private static readonly string[] MassNames = { "Mass", "Monoisotopic_mass" };
        private static readonly string[] IntensityNames = { "Intensity" };
        private static readonly string[] BeginNames = { "Min_time", "Time_begin" };
        private static readonly string[] EndNames = { "Max_time", "Time_end" };
        private static readonly string[] ApexNames = { "Apex_time", "Time_apex" };
        private static readonly string[] ApexIntensityNames = { "Apex_intensity", "Intensity_Apex" };
        private static readonly string[] MinChargeNames = { "Min_charge", "Minimum_charge_state" };
        private static readonly string[] MaxChargeNames = { "Max_charge", "Maximum_charge_state" };

        /// <summary>
        /// Reads every row of a tab-delimited <c>*_ms1.feature</c> file. Retention times are in MINUTES in both
        /// TopFD and FLASHDeconv output, matching <see cref="MassSpectrometry.MsDataScan.RetentionTime"/>.
        /// </summary>
        public static List<Ms1FeatureRecord> ReadFile(string filePath)
        {
            if (!File.Exists(filePath))
                throw new MetaMorpheusException("Feature file not found: " + filePath);

            var lines = File.ReadAllLines(filePath);
            if (lines.Length == 0)
                throw new MetaMorpheusException("Feature file is empty: " + filePath);

            // the BOM some TopFD builds emit would otherwise corrupt the first column name
            var header = lines[0].TrimStart('﻿').Split('\t');
            var columnIndexByName = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
            for (int i = 0; i < header.Length; i++)
                columnIndexByName[header[i].Trim()] = i;

            int massColumn = ResolveColumn(columnIndexByName, MassNames, filePath, required: true);
            int beginColumn = ResolveColumn(columnIndexByName, BeginNames, filePath, required: true);
            int endColumn = ResolveColumn(columnIndexByName, EndNames, filePath, required: true);
            int minChargeColumn = ResolveColumn(columnIndexByName, MinChargeNames, filePath, required: true);
            int maxChargeColumn = ResolveColumn(columnIndexByName, MaxChargeNames, filePath, required: true);
            int idColumn = ResolveColumn(columnIndexByName, IdNames, filePath, required: false);
            int intensityColumn = ResolveColumn(columnIndexByName, IntensityNames, filePath, required: false);
            int apexColumn = ResolveColumn(columnIndexByName, ApexNames, filePath, required: false);
            int apexIntensityColumn = ResolveColumn(columnIndexByName, ApexIntensityNames, filePath, required: false);

            var records = new List<Ms1FeatureRecord>(lines.Length - 1);
            for (int lineIndex = 1; lineIndex < lines.Length; lineIndex++)
            {
                if (string.IsNullOrWhiteSpace(lines[lineIndex]))
                    continue;
                var splits = lines[lineIndex].Split('\t');
                if (splits.Length <= massColumn)
                    continue;

                var record = new Ms1FeatureRecord
                {
                    Id = ReadInt(splits, idColumn, lineIndex - 1),
                    Mass = ReadDouble(splits, massColumn, 0),
                    Intensity = ReadDouble(splits, intensityColumn, 0),
                    RetentionTimeBegin = ReadDouble(splits, beginColumn, 0),
                    RetentionTimeEnd = ReadDouble(splits, endColumn, 0),
                    ChargeStateMin = ReadInt(splits, minChargeColumn, 1),
                    ChargeStateMax = ReadInt(splits, maxChargeColumn, 1),
                };
                record.RetentionTimeApex = ReadDouble(splits, apexColumn,
                    (record.RetentionTimeBegin + record.RetentionTimeEnd) / 2.0);
                record.ApexIntensity = ReadDouble(splits, apexIntensityColumn, record.Intensity);
                records.Add(record);
            }

            return records;
        }

        private static int ResolveColumn(Dictionary<string, int> columnIndexByName, string[] acceptedNames,
            string filePath, bool required)
        {
            foreach (var name in acceptedNames)
            {
                if (columnIndexByName.TryGetValue(name, out int index))
                    return index;
            }

            if (required)
                throw new MetaMorpheusException("Feature file " + filePath + " has none of the expected columns [" +
                    string.Join(", ", acceptedNames) + "]. Columns present: " +
                    string.Join(", ", columnIndexByName.Keys) + ".");
            return -1;
        }

        private static double ReadDouble(string[] splits, int column, double fallback) =>
            column >= 0 && column < splits.Length &&
            double.TryParse(splits[column], NumberStyles.Any, CultureInfo.InvariantCulture, out double value)
                ? value : fallback;

        private static int ReadInt(string[] splits, int column, int fallback)
        {
            if (column < 0 || column >= splits.Length)
                return fallback;
            // some columns are written as floats ("10.000000") even where an integer is meant
            return double.TryParse(splits[column], NumberStyles.Any, CultureInfo.InvariantCulture, out double value)
                ? (int)Math.Round(value) : fallback;
        }
    }
}
