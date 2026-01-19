using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using System.Globalization;
using System.Text;
using Chemistry;
using MassSpectrometry;
using System.Collections.Generic;

namespace Readers
{
    /// <summary>
    /// A class representing a single entry in FlashDeconv's ms1.tsv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class FlashDeconvRecord : IHasMass
    {
        [Ignore]
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t",
        };

        [Optional]
        [Name("Index")]
        public int Index { get; set; }

        [Name("FileName")]
        public string FileName { get; set; }

        [Name("ScanNum")]
        public int ZeroBasedScanNumber { get; set; }

        [Ignore]
        public int OneBasedScanNumber => ZeroBasedScanNumber + 1;

        /// <summary>
        /// Target Decoy of FlashDeconv
        ///  0 -> Target
        ///  1 -> Isotope Decoy
        ///  2 -> Noise Decoy
        ///  3 -> Charge Decoy
        /// </summary>
        [Optional]
        [Name("Decoy")]
        public int Decoy { get; set; }

        [Name("RetentionTime")]
        public double RetentionTime { get; set; }

        /// <summary>
        /// Number of m/z in the spectrum after deconvolution
        /// </summary>
        [Name("MassCountInSpec")]
        public int MassCountInSpec { get; set; }

        [Name("AverageMass")]
        public double AverageMass { get; set; }

        [Name("MonoisotopicMass")]
        public double MonoisotopicMass { get; set; }

        [Name("SumIntensity")]
        public double SumIntensity { get; set; }

        /// <summary>
        /// Min charge of all charge states identified
        /// </summary>
        [Name("MinCharge")]
        public int MinCharge { get; set; }

        /// <summary>
        /// Max charge of all charge states identified
        /// </summary>
        [Name("MaxCharge")]
        public int MaxCharge { get; set; }

        /// <summary>
        /// The number of isotopic peaks from all detected charge states 
        /// </summary>
        [Name("PeakCount")]
        public int PeakCount { get; set; }

        [Name("IsotopeCosine")]
        public double IsotopeCosine { get; set; }

        [Name("ChargeScore")]
        public double ChargeScore { get; set; }

        /// <summary>
        /// Calculation can be found here
        /// https://www.nature.com/articles/s41467-022-31922-z
        /// </summary>
        [Name("MassSNR")]
        public double MassSNR { get; set; }

        /// <summary>
        /// Calculation can be found here
        /// https://www.nature.com/articles/s41467-022-31922-z
        /// </summary>
        [Name("ChargeSNR")]
        public double ChargeSNR { get; set; }

        /// <summary>
        /// Charge of the entries representative peak
        /// </summary>
        [Name("RepresentativeCharge")]
        public int RepresentativeCharge { get; set; }

        /// <summary>
        /// Start mz of the entries representative peak
        /// </summary>
        [Name("RepresentativeMzStart")]
        public double RepresentativeMzMin { get; set; }

        /// <summary>
        /// End mz of the entries representative peak
        /// </summary>
        [Name("RepresentativeMzEnd")]
        public double RepresentativeMzMax { get; set; }
        [Optional]
        [Name("QScore")]
        public double QScore { get; set; }

        [Optional]
        /// <summary>
        /// IMPORTANT: Decoys will always have their Q value set to 1
        /// </summary>
        [Name("Qvalue")]
        public double Qvalue { get; set; }

        [Optional]
        [Name("QvalueWithIsotopeDecoyOnly")]
        public double QvalueWithIsotopeDecoyOnly { get; set; }
        [Optional]
        [Name("QvalueWithNoiseDecoyOnly")]
        public double QvalueWithNoiseDecoyOnly { get; set; }
        [Optional]
        [Name("QvalueWithChargeDecoyOnly")]
        public double QvalueWithChargeDecoyOnly { get; set; }
    }
}
