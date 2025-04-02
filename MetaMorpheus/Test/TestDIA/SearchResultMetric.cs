using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using System.Globalization;
using EngineLayer;
using Omics.Fragmentation;
using Org.BouncyCastle.Bcpg;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using LinqStatistics;
using System.Windows.Markup;
using NUnit.Framework;
using static iText.StyledXmlParser.Jsoup.Select.Evaluator;
using Omics.SpectrumMatch;
using Easy.Common.Extensions;
using ThermoFisher.CommonCore.Data.FilterEnums;

namespace Test.TestDIA
{
    public class FragmentMetric
    {
        public string Identifier { get; set; }
        public string PsmId { get; set; }
        public string PsmFullSeq { get; set; }
        public string FragmentId { get; set; }
        public int Charge { get; set; }
        public double Mz { get; set; }
        public int ISDEnergy { get; set; }
        public double RoundedMz => Math.Round(Mz, 0);
        public double Intensity { get; set; }
        public double NormalizedIntensity { get; set; }
        public double NeutralMass { get; set; }
        public bool IsTerminalFragment { get; set; }

        public static IEnumerable<FragmentMetric> FromPsm(PsmFromTsv psm, string identifer = "")
        {
            foreach (var matchedFragment in psm.MatchedIons)
            {

                var result = new FragmentMetric
                {
                    Identifier = identifer,
                    Charge = matchedFragment.Charge,
                    Mz = matchedFragment.Mz,
                    Intensity = matchedFragment.Intensity,
                    NormalizedIntensity = matchedFragment.Intensity / psm.TotalIonCurrent.Value,
                    NeutralMass = matchedFragment.NeutralTheoreticalProduct.NeutralMass,
                    IsTerminalFragment = matchedFragment.NeutralTheoreticalProduct.SecondaryProductType == null,
                    FragmentId = matchedFragment.Annotation,
                    PsmId = psm.BaseSeq,
                    PsmFullSeq = psm.FullSequence,
                    ISDEnergy = FindEnergy(psm),
                };

                yield return result;
            }
        }

        public static int FindEnergy(PsmFromTsv psm)
        {
            int remainder = (psm.Ms2ScanNumber - 1) % 4;
            switch (remainder)
            {
                case 1:
                    return 60;
                case 2:
                    return 80;
                case 3:
                    return 100;
            }
            return 0;
        }

        public static IEnumerable<FragmentMetric> FromPsmFile(string path, string identifer = "")
        {
            var psms = PsmTsvReader.ReadTsv(path, out _).Where(psm => psm.PassesConfidenceFilter());
            return psms.SelectMany(p => FromPsm(p, identifer));
        }

        public FragmentMetric() { }
    }

    public class FragmentMetricFile : ResultFile<FragmentMetric>, IResultFile
    {

        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public FragmentMetricFile(string filePath) : base(filePath, Software.Unspecified) { }
        public FragmentMetricFile(string filePath, IEnumerable<FragmentMetric> results) : base(filePath, Software.Unspecified)
        {
            Results = results.ToList();
        }
        public FragmentMetricFile() : base() { }


        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<FragmentMetric>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<FragmentMetric>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }



    public static class PsmExtensions
    {
        public static bool PassesConfidenceFilter(this PsmFromTsv psm) => psm.DecoyContamTarget == "T" && psm.QValue <= 0.01 && psm.QValueNotch <= 0.01;
    }

    // Generate one of these for each PSM after the search
    public class SearchResultMetric
    {
        [Name("Terminal Fragment Count")]
        public double TerminalFragmentCount { get; set; }

        [Name("Internal Fragment Count")]
        public double InternalFragmentCount { get; set; }

        [Name("Median Fragment Intensity")]
        public double MedianFragmentIntensity { get; set; }

        [Name("Matched Intensity Fraction")]
        public double MatchedIntensityFraction { get; set; }

        [Name("Mean Fragment Charge")]
        public double MeanFragmentCharge { get; set; }

        public SearchResultMetric(PsmFromTsv psm)
        {
            TerminalFragmentCount = psm.MatchedIons.Count(p => p.NeutralTheoreticalProduct.SecondaryProductType == null);
            InternalFragmentCount = psm.MatchedIons.Count(p => p.NeutralTheoreticalProduct.SecondaryProductType != null);

            // TODO: set values
            MedianFragmentIntensity = psm.MatchedIons.Select(p => p.Intensity).Median();
            MatchedIntensityFraction = psm.MatchedIons.Sum(p => p.Intensity) / psm.TotalIonCurrent.Value;
            MeanFragmentCharge = psm.MatchedIons.Select(p => p.Charge).Average();
        }

        public SearchResultMetric() { }
    }

    public class SearchResultsMetricsFile : ResultFile<SearchResultMetric>, IResultFile
    {

        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public SearchResultsMetricsFile(string filePath) : base (filePath, Software.Unspecified) { }
        public SearchResultsMetricsFile(string filePath, IEnumerable<SearchResultMetric> results) : base (filePath, Software.Unspecified) 
        { 
            Results = results.ToList();
        }
        public SearchResultsMetricsFile() : base () { }

        public static SearchResultsMetricsFile GetFromPsmFilePath(string psmFromTsvPath)
        {
            List<SearchResultMetric> results = new List<SearchResultMetric>();
            if (File.Exists(psmFromTsvPath))
            {
                var psms = PsmTsvReader.ReadTsv(psmFromTsvPath, out _).Where(psm => psm.PassesConfidenceFilter()).ToArray();

                foreach (var psm in psms)
                {
                    results.Add(new SearchResultMetric(psm));
                }
            }
            else
            {
                throw new FileNotFoundException("File not found", psmFromTsvPath);
            }

            var searchResultsMetricsFile = new SearchResultsMetricsFile()
            {
                Results = results,
            };
            return searchResultsMetricsFile;
        }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<SearchResultMetric>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<SearchResultMetric>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }

    public class Tests
    {
        [Test]
        public static void WriteResults()
        {
            var psmPath = @"E:\DIA\Data\250313_DIA\0317_5pro_ISD\Task1-SearchTask\Individual File Results\03-17-25_CE_5pro_90min_ISD60-80-100_overlap-equal_labelCorrected_PSMs.psmtsv";
            var frags = FragmentMetric.FromPsmFile(psmPath, "overlap");

            //var largePath = @"E:\DIA\Data\250313_DIA\0318_5pro\Task1-SearchTask\Individual File Results\03-18-25_CE_5pro_90min_ISD60-80-100_equal_labelCorrected_PSMs.psmtsv";
            //var largeFrags = FragmentMetric.FromPsmFile(largePath, "equal");

            //var allFrags = frags.Concat(largeFrags).ToList();
            var outPath = @"E:\ISD Project\TestIsdDataAnalysis\overlap_FragmentsMetrics.tsv";
            var fragFile = new FragmentMetricFile(outPath, frags);
            fragFile.WriteResults(outPath);

            //var searchResults = SearchResultsMetricsFile.GetFromPsmFilePath(psmPath);
            //var outPath = @"E:\ISD Project\TestIsdDataAnalysis\SearchResultsMetricsFile_03-17-25_CE_5pro_90min_ISD60-80-100_overlap-equal_labelCorrected_PSMs.tsv";
            //searchResults.WriteResults(outPath);
        }
    }

}
