using MassSpectrometry;
using Microsoft.ML.Data;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using CsvHelper.Configuration;
using CsvHelper;
using System.Globalization;
using System.IO;

namespace EngineLayer.DIA
{
    public class PfPairTrainingSample
    {
        public float Correlation { get; set; }   // e.g., XIC correlation
        public float ApexRtDelta { get; set; }   // e.g., RT diff
        public float Overlap { get; set; }      // e.g., XIC overlap
        public float FragmentIntensity { get; set; } // e.g., fragment intensity
        public float NormalizedIntensityRank { get; set; } // e.g., normalized intensity rank
        public float PsmScore { get; set; }
        public float SharedXIC { get; set; }
        public bool Label { get; set; }
        public float Weight { get; set; }
        public float FragmentRank { get; set; }
        public float PrecursorRank { get; set; }

        public PfPairTrainingSample(ExtractedIonChromatogram precursor, ExtractedIonChromatogram fragment, bool label = false, SpectralMatch psm = null)
        {
            Correlation = (float)PrecursorFragmentsGroup.CalculateXicCorrelationXYData(precursor, fragment);
            ApexRtDelta = (float)Math.Abs(precursor.ApexRT - fragment.ApexRT);
            Overlap = (float)PrecursorFragmentsGroup.CalculateXicOverlapRatio(precursor, fragment);
            FragmentIntensity = (float)fragment.ApexPeak.Intensity;
            SharedXIC = (float)PrecursorFragmentPair.CalculateSharedXIC(precursor, fragment);
            PsmScore = psm != null ? (float)psm.Score : 0;
            Label = label; // default label is false, should be set to true for positive samples
        }

        public PfPairTrainingSample(PrecursorFragmentPair pfPair)
        {
            Correlation = pfPair.Correlation.HasValue? (float)pfPair.Correlation : (float)PrecursorFragmentsGroup.CalculateXicCorrelationXYData(pfPair.PrecursorXic, pfPair.FragmentXic);
            ApexRtDelta = (float)Math.Abs(pfPair.PrecursorXic.ApexRT - pfPair.FragmentXic.ApexRT);
            Overlap = pfPair.Overlap.HasValue ? (float)pfPair.Overlap : (float)PrecursorFragmentsGroup.CalculateXicOverlapRatio(pfPair.PrecursorXic, pfPair.FragmentXic);
            FragmentIntensity = (float)pfPair.FragmentXic.ApexPeak.Intensity;
            //SharedXIC = (float)PrecursorFragmentPair.CalculateSharedXIC(pfPair.PrecursorXic, pfPair.FragmentXic);
            FragmentRank = (float)pfPair.FragmentRank;
            PrecursorRank = (float)pfPair.PrecursorRank;
        }

        public PfPairTrainingSample() { }
    }

    public class PfPairTrainingSampleFile : ResultFile<PfPairTrainingSample>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public PfPairTrainingSampleFile() : base() { }
        public PfPairTrainingSampleFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<PfPairTrainingSample>().ToList();
        }

        public string FullFileName { get; set; }
        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<PfPairTrainingSampleFile>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }

    }

    public class PFpairPrediction
    {
        [ColumnName("PredictedLabel")]
        public bool Prediction;

        public float Probability;

        public float Score;
    }
}
