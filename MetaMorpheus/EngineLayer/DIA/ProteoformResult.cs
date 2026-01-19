using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CsvHelper.Configuration;
using Readers;

namespace EngineLayer.DIA
{
    public class ProteoformResult
    {
        public string Accession { get; set; }
        public double MonoMass { get; set; }
        public string Modifications { get; set; }
        public string BaseSequence { get; set; }
        public string FullSequence { get; set; }
        public int UniqueFragmentCount { get; set; }
        public int UniqueTerminalFragmentCount { get; set; }
        public string MatchedIons { get; set; }

        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t"
        };

        public ProteoformResult() { }
    }

    public class ProteoformResultFile : ResultFile<ProteoformResult>, IResultFile
    {
        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }

        public ProteoformResultFile(string filePath) : base(filePath) { }

        public ProteoformResultFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvHelper.CsvReader(new StreamReader(FilePath), ProteoformResult.CsvConfiguration);
            Results = csv.GetRecords<ProteoformResult>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(FilePath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvHelper.CsvWriter(new StreamWriter(File.Create(outputPath)), ProteoformResult.CsvConfiguration))
            {
                csv.WriteHeader<ProteoformResult>();
                foreach (var result in Results)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }

        public static void WriteProteoformResults(string outputPath, IEnumerable<SpectralMatch> psms)
        {
            var groupByProteoform = psms.GroupBy(p => p.FullSequence);
            var allResults = new List<ProteoformResult>();
            foreach (var group in groupByProteoform)
            {
                var allMatchedIons = group.SelectMany(g => g.MatchedFragmentIons).ToList();
                int allUniqueFragCount = allMatchedIons.Select(i => i.NeutralTheoreticalProduct.MonoisotopicMass).Distinct().Count();
                int uniqueTerminalFragmentCount = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false)
                    .Select(i => i.MonoisotopicMass).Distinct().Count();
                var annotations = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false).OrderBy(i => i.FragmentNumber).Select(i => i.Annotation).Distinct().ToArray();
                var proteoformResult = new ProteoformResult
                {
                    Accession = group.First().Accession,
                    MonoMass = group.First().ScanPrecursorMass,
                    Modifications = group.First().ModsIdentified.ToString(),
                    BaseSequence = group.First().BaseSequence,
                    FullSequence = group.First().FullSequence,
                    UniqueFragmentCount = allUniqueFragCount,
                    UniqueTerminalFragmentCount = uniqueTerminalFragmentCount,
                    MatchedIons = string.Join(",", annotations)
                };
                allResults.Add(proteoformResult);
            }
            var resultFile = new ProteoformResultFile
            {
                Results = allResults.OrderByDescending(r => r.UniqueTerminalFragmentCount).ToList()
            };
            resultFile.WriteResults(outputPath);
        }

        public static void WriteProteoformResults(string outputPath, IEnumerable<PsmFromTsv> psmTsvs)
        {
            var groupByProteoform = psmTsvs.GroupBy(p => p.FullSequence);
            var allResults = new List<ProteoformResult>();
            foreach (var group in groupByProteoform)
            {
                var allMatchedIons = group.SelectMany(g => g.MatchedIons).ToList();
                int allUniqueFragCount = allMatchedIons.Select(i => i.NeutralTheoreticalProduct.MonoisotopicMass).Distinct().Count();
                int uniqueTerminalFragmentCount = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false)
                    .Select(i => i.MonoisotopicMass).Distinct().Count();
                var annotations = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false).OrderBy(i => i.FragmentNumber).Select(i => i.Annotation).Distinct().ToArray();
                var proteoformResult = new ProteoformResult
                {
                    Accession = group.First().Accession,
                    MonoMass = group.First().MonoisotopicMass,
                    BaseSequence = group.First().BaseSequence,
                    FullSequence = group.First().FullSequence,
                    UniqueFragmentCount = allUniqueFragCount,
                    UniqueTerminalFragmentCount = uniqueTerminalFragmentCount,
                    MatchedIons = string.Join(",", annotations)
                };
                allResults.Add(proteoformResult);
            }
            var resultFile = new ProteoformResultFile
            {
                Results = allResults.OrderByDescending(r => r.UniqueTerminalFragmentCount).ToList()
            };
            resultFile.WriteResults(outputPath);
        }
    }

}
