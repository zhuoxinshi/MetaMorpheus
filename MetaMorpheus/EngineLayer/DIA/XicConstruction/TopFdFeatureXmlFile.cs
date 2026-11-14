using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml;

namespace EngineLayer.DIA
{
    /// <summary>
    /// One point of a TopFD feature's elution profile: the feature was observed, at one charge state, in this
    /// scan with this envelope intensity.
    /// </summary>
    public class TopFdEnvelope
    {
        public int OneBasedScanNumber { get; set; }
        /// <summary>Retention time in minutes (the file stores seconds; the reader converts).</summary>
        public double RetentionTime { get; set; }
        /// <summary>Summed intensity of the isotopic envelope in this scan.</summary>
        public double IntensitySum { get; set; }
        /// <summary>Intensity of the tallest isotopic peak of the envelope in this scan.</summary>
        public double MaxIntensity { get; set; }
    }

    /// <summary>
    /// One charge state of a TopFD feature, with its own elution window and per-scan envelope list.
    /// </summary>
    public class TopFdSingleChargeFeature
    {
        public int Charge { get; set; }
        public double RetentionTimeBegin { get; set; }
        public double RetentionTimeEnd { get; set; }
        public int ScanBegin { get; set; }
        public int ScanEnd { get; set; }
        public double Intensity { get; set; }
        public int EnvelopeNumber { get; set; }
        public double MonoisotopicMz { get; set; }
        public double AverageMz { get; set; }
        public double ReferenceMz { get; set; }
        public double MaxIntensity { get; set; }
        public List<TopFdEnvelope> Envelopes { get; } = new List<TopFdEnvelope>();
    }

    /// <summary>
    /// One TopFD feature (a <c>frac_feature</c> element) with its per-charge elution profiles.
    /// </summary>
    public class TopFdFeature
    {
        public int Id { get; set; }
        public int FractionId { get; set; }
        public string FileName { get; set; }
        public double Mass { get; set; }
        public double Intensity { get; set; }
        public double RetentionTimeBegin { get; set; }
        public double RetentionTimeEnd { get; set; }
        public double RetentionTimeApex { get; set; }
        public int ScanBegin { get; set; }
        public int ScanEnd { get; set; }
        public int ApexScan { get; set; }
        public double ApexIntensity { get; set; }
        public int ChargeStateMin { get; set; }
        public int ChargeStateMax { get; set; }
        public int RepresentativeCharge { get; set; }
        public double RepresentativeAverageMz { get; set; }
        public int EnvelopeNumber { get; set; }
        /// <summary>TopFD's ECScore for the feature, in [0, 1].</summary>
        public double EcScore { get; set; }
        public List<TopFdSingleChargeFeature> SingleChargeFeatures { get; } = new List<TopFdSingleChargeFeature>();

        public override string ToString() =>
            "Feature " + Id + ": " + Mass.ToString("F4", CultureInfo.InvariantCulture) + " Da, z " + ChargeStateMin + "-" + ChargeStateMax +
            ", RT " + RetentionTimeBegin.ToString("F2", CultureInfo.InvariantCulture) + "-" + RetentionTimeEnd.ToString("F2", CultureInfo.InvariantCulture) +
            ", " + EnvelopeNumber + " envelopes";
    }

    /// <summary>
    /// Streaming reader for TopFD's <c>*_feature.xml</c>. Unlike the <c>*_ms1.feature</c> table, which reduces
    /// each feature's elution to begin/end/apex values, the xml lists every envelope of every charge state
    /// (scan, retention time, intensity), so it yields feature XICs directly (see
    /// <see cref="FeatureXmlXicConstructor"/>). Retention times are converted from seconds to minutes to match
    /// <see cref="MassSpectrometry.MsDataScan.RetentionTime"/>.
    /// </summary>
    public static class TopFdFeatureXmlFile
    {
        /// <summary>Reads every feature into memory.</summary>
        public static List<TopFdFeature> ReadFile(string filePath) => EnumerateFeatures(filePath).ToList();

        /// <summary>Streams the features one at a time.</summary>
        public static IEnumerable<TopFdFeature> EnumerateFeatures(string filePath)
        {
            if (!File.Exists(filePath))
                throw new MetaMorpheusException("TopFD feature xml not found: " + filePath);

            const double secondsToMinutes = 1.0 / 60.0;
            var settings = new XmlReaderSettings { IgnoreWhitespace = true, IgnoreComments = true };
            using var reader = XmlReader.Create(filePath, settings);

            TopFdFeature feature = null;
            TopFdSingleChargeFeature chargeFeature = null;
            TopFdEnvelope envelope = null;
            string currentElement = null;

            // Values are read from Text nodes. ReadElementContentAsString would consume the end tag and, combined
            // with the loop's own Read(), skip the following element.
            while (reader.Read())
            {
                switch (reader.NodeType)
                {
                    case XmlNodeType.Element:
                        switch (reader.Name)
                        {
                            case "frac_feature":
                                feature = new TopFdFeature();
                                chargeFeature = null;
                                envelope = null;
                                break;
                            case "single_charge_feature":
                                chargeFeature = new TopFdSingleChargeFeature();
                                envelope = null;
                                break;
                            case "envelope":
                                envelope = new TopFdEnvelope();
                                break;
                        }
                        currentElement = reader.Name;
                        break;

                    case XmlNodeType.Text:
                    case XmlNodeType.CDATA:
                        if (feature == null || currentElement == null)
                            break;
                        // the same element names (time_begin, scan_begin, intensity, envelope_num, max_intensity)
                        // appear at several levels, so assign to the innermost open object
                        if (envelope != null)
                            AssignEnvelope(envelope, currentElement, reader.Value, secondsToMinutes);
                        else if (chargeFeature != null)
                            AssignChargeFeature(chargeFeature, currentElement, reader.Value, secondsToMinutes);
                        else
                            AssignFeature(feature, currentElement, reader.Value, secondsToMinutes);
                        break;

                    case XmlNodeType.EndElement:
                        switch (reader.Name)
                        {
                            case "envelope":
                                if (chargeFeature != null && envelope != null) chargeFeature.Envelopes.Add(envelope);
                                envelope = null;
                                break;
                            case "single_charge_feature":
                                if (feature != null && chargeFeature != null) feature.SingleChargeFeatures.Add(chargeFeature);
                                chargeFeature = null;
                                break;
                            case "frac_feature":
                                if (feature != null) yield return feature;
                                feature = null;
                                break;
                        }
                        currentElement = null;
                        break;
                }
            }
        }

        private static void AssignEnvelope(TopFdEnvelope envelope, string name, string value, double timeScale)
        {
            switch (name)
            {
                case "scan": envelope.OneBasedScanNumber = ParseInt(value); break;
                case "rt": envelope.RetentionTime = ParseDouble(value) * timeScale; break;
                case "intensity_sum": envelope.IntensitySum = ParseDouble(value); break;
                case "max_intensity": envelope.MaxIntensity = ParseDouble(value); break;
            }
        }

        private static void AssignChargeFeature(TopFdSingleChargeFeature chargeFeature, string name, string value, double timeScale)
        {
            switch (name)
            {
                case "charge": chargeFeature.Charge = ParseInt(value); break;
                case "time_begin": chargeFeature.RetentionTimeBegin = ParseDouble(value) * timeScale; break;
                case "time_end": chargeFeature.RetentionTimeEnd = ParseDouble(value) * timeScale; break;
                case "scan_begin": chargeFeature.ScanBegin = ParseInt(value); break;
                case "scan_end": chargeFeature.ScanEnd = ParseInt(value); break;
                case "intensity": chargeFeature.Intensity = ParseDouble(value); break;
                case "envelope_num": chargeFeature.EnvelopeNumber = ParseInt(value); break;
                case "mono_mz": chargeFeature.MonoisotopicMz = ParseDouble(value); break;
                case "average_mz": chargeFeature.AverageMz = ParseDouble(value); break;
                case "refer_mz": chargeFeature.ReferenceMz = ParseDouble(value); break;
                case "max_intensity": chargeFeature.MaxIntensity = ParseDouble(value); break;
            }
        }

        private static void AssignFeature(TopFdFeature feature, string name, string value, double timeScale)
        {
            switch (name)
            {
                case "file_name": feature.FileName = value; break;
                case "frac_id": feature.FractionId = ParseInt(value); break;
                case "feat_id": feature.Id = ParseInt(value); break;
                case "mono_mass": feature.Mass = ParseDouble(value); break;
                case "intensity": feature.Intensity = ParseDouble(value); break;
                case "time_begin": feature.RetentionTimeBegin = ParseDouble(value) * timeScale; break;
                case "time_end": feature.RetentionTimeEnd = ParseDouble(value) * timeScale; break;
                case "scan_begin": feature.ScanBegin = ParseInt(value); break;
                case "scan_end": feature.ScanEnd = ParseInt(value); break;
                case "min_charge": feature.ChargeStateMin = ParseInt(value); break;
                case "max_charge": feature.ChargeStateMax = ParseInt(value); break;
                case "apex_time": feature.RetentionTimeApex = ParseDouble(value) * timeScale; break;
                case "apex_scan": feature.ApexScan = ParseInt(value); break;
                case "apex_inte": feature.ApexIntensity = ParseDouble(value); break;
                case "rep_charge": feature.RepresentativeCharge = ParseInt(value); break;
                case "rep_avg_mz": feature.RepresentativeAverageMz = ParseDouble(value); break;
                case "envelope_num": feature.EnvelopeNumber = ParseInt(value); break;
                case "ec_score": feature.EcScore = ParseDouble(value); break;
            }
        }

        private static double ParseDouble(string value) =>
            double.TryParse(value, NumberStyles.Any, CultureInfo.InvariantCulture, out double result) ? result : 0;

        // some integer fields are written as decimals ("626.00"), so parse as double and round
        private static int ParseInt(string value) =>
            double.TryParse(value, NumberStyles.Any, CultureInfo.InvariantCulture, out double result)
                ? (int)Math.Round(result) : 0;
    }
}
