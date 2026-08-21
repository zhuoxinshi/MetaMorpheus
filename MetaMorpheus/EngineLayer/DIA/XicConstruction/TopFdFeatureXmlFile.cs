using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml;

namespace EngineLayer.DIA
{
    /// <summary>
    /// One point of a TopFD proteoform feature's elution profile: the feature was observed, at one charge state,
    /// in this scan with this summed envelope intensity.
    /// </summary>
    public class TopFdEnvelope
    {
        public int OneBasedScanNumber { get; set; }
        /// <summary>Retention time in MINUTES (the file stores seconds; the reader converts).</summary>
        public double RetentionTime { get; set; }
        /// <summary>Summed intensity of the isotopic envelope in this scan.</summary>
        public double IntensitySum { get; set; }
        /// <summary>Intensity of the tallest isotopic peak of the envelope in this scan.</summary>
        public double MaxIntensity { get; set; }
    }

    /// <summary>
    /// One charge state of a TopFD proteoform feature, with its own elution window and envelope list. This is
    /// already an XIC: <see cref="Envelopes"/> is the per-scan intensity series TopFD measured for this charge.
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
    /// One TopFD proteoform feature (a <c>frac_feature</c>) with its per-charge elution profiles.
    ///
    /// This is the same feature the <c>*_ms1.feature</c> tsv summarises, but the XML additionally carries the
    /// ENVELOPE LISTS - scan, retention time and intensity for every point of every charge state - which the tsv
    /// reduces to a single <c>Envelope_num</c> count.
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
        /// <summary>TopFD's ECScore for the feature, in [0, 1]. TopFD's own default cutoff is 0.1, which is
        /// permissive; raise it here to keep only confident features.</summary>
        public double EcScore { get; set; }
        public List<TopFdSingleChargeFeature> SingleChargeFeatures { get; } = new List<TopFdSingleChargeFeature>();

        public int ChargeCount => SingleChargeFeatures.Count;

        /// <summary>Every envelope of every charge state, flattened.</summary>
        public IEnumerable<TopFdEnvelope> AllEnvelopes => SingleChargeFeatures.SelectMany(c => c.Envelopes);

        /// <summary>Distinct scans the feature was observed in, across all charge states.</summary>
        public int DistinctScanCount => AllEnvelopes.Select(e => e.OneBasedScanNumber).Distinct().Count();

        public override string ToString() =>
            "Feature " + Id + ": " + Mass.ToString("F4") + " Da, z " + ChargeStateMin + "-" + ChargeStateMax +
            ", RT " + RetentionTimeBegin.ToString("F2") + "-" + RetentionTimeEnd.ToString("F2") +
            ", " + EnvelopeNumber + " envelopes";
    }

    /// <summary>
    /// Streaming reader for TopFD's <c>*_feature.xml</c>.
    ///
    /// This file is the richest of TopFD's three MS1 outputs. <c>*_ms1.msalign</c> holds per-scan deconvoluted
    /// masses with no feature assignment; <c>*_ms1.feature</c> holds features but collapses their elution to a
    /// begin/end/apex triple; <c>*_feature.xml</c> holds the features AND the full per-charge, per-scan envelope
    /// lists. Reading it therefore yields feature XICs directly, with no peak tracing - see
    /// <see cref="FeatureXmlXicConstructor"/>.
    ///
    /// Parsed with <see cref="XmlReader"/> rather than a DOM because these files run to ~100 MB for a single ISD
    /// run (24k features, 280k envelopes).
    ///
    /// TIME UNITS: the XML stores seconds while the companion tsv stores minutes (the same feature reads
    /// <c>apex_time</c> 821.94 in the XML and <c>Apex_time</c> 13.699 in the tsv). This reader converts to
    /// MINUTES so everything matches <see cref="MassSpectrometry.MsDataScan.RetentionTime"/>.
    /// </summary>
    public static class TopFdFeatureXmlFile
    {
        /// <summary>Reads every feature into memory.</summary>
        public static List<TopFdFeature> ReadFile(string filePath, bool convertSecondsToMinutes = true) =>
            EnumerateFeatures(filePath, convertSecondsToMinutes).ToList();

        /// <summary>
        /// Streams features one at a time, so a caller that only needs a subset never materialises the whole file.
        /// </summary>
        public static IEnumerable<TopFdFeature> EnumerateFeatures(string filePath, bool convertSecondsToMinutes = true)
        {
            if (!File.Exists(filePath))
                throw new MetaMorpheusException("TopFD feature xml not found: " + filePath);

            double timeScale = convertSecondsToMinutes ? 1.0 / 60.0 : 1.0;
            var settings = new XmlReaderSettings { IgnoreWhitespace = true, IgnoreComments = true };

            using var reader = XmlReader.Create(filePath, settings);

            TopFdFeature feature = null;
            TopFdSingleChargeFeature chargeFeature = null;
            TopFdEnvelope envelope = null;
            string currentElement = null;

            // Driven off Text nodes rather than ReadElementContentAsString: that helper consumes the element's
            // end tag and leaves the reader on the FOLLOWING node, so combining it with the loop's own Read()
            // silently skips every other element.
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
                        // the SAME element names appear at more than one level (time_begin, scan_begin,
                        // intensity, envelope_num, max_intensity), so assign to the innermost object open
                        if (envelope != null)
                            AssignEnvelope(envelope, currentElement, reader.Value, timeScale);
                        else if (chargeFeature != null)
                            AssignChargeFeature(chargeFeature, currentElement, reader.Value, timeScale);
                        else
                            AssignFeature(feature, currentElement, reader.Value, timeScale);
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

        // some integer-valued fields are written as floats ("626.00"), so round rather than fail
        private static int ParseInt(string value) =>
            double.TryParse(value, NumberStyles.Any, CultureInfo.InvariantCulture, out double result)
                ? (int)Math.Round(result) : 0;
    }
}
