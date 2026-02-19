using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using Easy.Common.Extensions;
using NWaves.Transforms;
using Omics.Fragmentation;
using Omics.Modifications;
using Readers;
using static Microsoft.FSharp.Core.ByRefKinds;

namespace EngineLayer.DIA
{
    public class ProteoformResult
    {
        public string ProrteinName { get; set; }
        public string Accession { get; set; }
        public double MonoMass { get; set; }
        public string Modifications { get; set; }
        public string BaseSequence { get; set; }
        public string FullSequence { get; set; }
        public int UniqueFragmentCount { get; set; }
        public int UniqueTerminalFragmentCount { get; set; }
        public string MatchedTerminalIons { get; set; }
        //public string MatchedInternalIons { get; set; }
        public string ModSiteTerminalIons { get; set; }
        public string ModSiteInternalIons { get; set; }

        [Optional]
        public string Isd60UniqueTerminal { get; set; }
        [Optional]
        public string Isd60FragmentCount { get; set; }
        [Optional]
        public string Isd80UniqueTerminal { get; set; }
        [Optional]
        public string Isd80FragmentCount { get; set; }
        [Optional]
        public string Isd100UniqueTerminal { get; set; }
        [Optional]
        public string Isd100FragmentCount { get; set; }

        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t"
        };

        public ProteoformResult() { }

        public void FindModSiteCoverage()
        {
            
        }
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
                int allUniqueFragCount = allMatchedIons.Select(i => i.NeutralTheoreticalProduct.Annotation).Distinct().Count();
                int uniqueTerminalFragmentCount = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false)
                    .Select(i => i.MonoisotopicMass).Distinct().Count();
                var annotations = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false).OrderBy(i => i.FragmentNumber).Select(i => i.Annotation).Distinct().ToArray();
                var internalIons = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == true).OrderBy(i => i.FragmentNumber).Select(i => i.Annotation).Distinct().ToArray();
                var proteoformResult = new ProteoformResult
                {
                    Accession = group.First().Accession,
                    MonoMass = group.First().ScanPrecursorMass,
                    Modifications = group.First().ModsIdentified.ToString(),
                    BaseSequence = group.First().BaseSequence,
                    FullSequence = group.First().FullSequence,
                    UniqueFragmentCount = allUniqueFragCount,
                    UniqueTerminalFragmentCount = uniqueTerminalFragmentCount,
                    MatchedTerminalIons = string.Join(",", annotations),
                    //MatchedInternalIons = string.Join(",", internalIons)
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
                int allUniqueFragCount = allMatchedIons.Select(i => i.NeutralTheoreticalProduct.Annotation).Distinct().Count();
                var terminalIons = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false).OrderBy(i => i.FragmentNumber);
                var internalIons = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == true && n.SecondaryFragmentNumber - n.FragmentNumber > 10).OrderBy(i => i.FragmentNumber);
                var uniqueTerminalIons = terminalIons.Select(i => i.Annotation).Distinct().ToArray();
                var uniqueInternalIons = internalIons.Select(i => i.Annotation).Distinct().ToArray();

                var proteoformResult = new ProteoformResult
                {
                    ProrteinName = group.First().Name,
                    Accession = group.First().Accession,
                    MonoMass = group.First().MonoisotopicMass,
                    BaseSequence = group.First().BaseSequence,
                    FullSequence = group.First().FullSequence,
                    UniqueFragmentCount = allUniqueFragCount,
                    UniqueTerminalFragmentCount = uniqueTerminalIons.Length,
                    MatchedTerminalIons = string.Join(",", uniqueTerminalIons),
                    //MatchedInternalIons = string.Join(",", uniqueInternalIons)
                };

                var mods = SpectrumMatchFromTsv.ParseModifications(group.First().FullSequence).Where(m => !(m.Value.Count() == 1 && m.Value.First().Contains("Fixed")));
                if (mods.Count() == 0)
                {
                    proteoformResult.Modifications = "NA";
                }
                else
                {
                    proteoformResult.Modifications = string.Join("; ", mods.Select(kvp => $"{kvp.Key}: [{string.Join(", ", kvp.Value)}]"));
                    var outputStringTerminal = new StringBuilder();
                    var outputStringInternal = new StringBuilder();
                    foreach (var kvp in mods.Where(kvp => kvp.Key != 0))
                    {
                        var b_ions = terminalIons.Where(i => i.ProductType == ProductType.b && i.AminoAcidPosition > kvp.Key);
                        var y_ions = terminalIons.Where(i => i.ProductType == ProductType.y && i.AminoAcidPosition < kvp.Key);
                        var terminalIonsWithModSite = b_ions.Concat(y_ions).Select(i => i.Annotation).Distinct().ToArray();
                        outputStringTerminal.AppendLine(kvp.Key + $":({terminalIonsWithModSite.Length}) " + string.Join(",", terminalIonsWithModSite));

                        var internalIonsWithModSite = internalIons.Where(i => i.FragmentNumber <= kvp.Key + 1 && i.SecondaryFragmentNumber >= kvp.Key + 1).Select(i => i.Annotation).Distinct().ToArray();
                        outputStringInternal.AppendLine(kvp.Key + $":({internalIonsWithModSite.Length}) " + string.Join(",", internalIonsWithModSite));
                    }
                    proteoformResult.ModSiteTerminalIons = outputStringTerminal.ToString();
                    proteoformResult.ModSiteInternalIons = outputStringInternal.ToString();
                }

                if (outputPath.Contains("ISD"))
                {
                    var psms60Fragments = group.Where(g => g.Ms2ScanNumber % 4 == 2).SelectMany(g => g.MatchedIons).Select(m => m.NeutralTheoreticalProduct).ToList();
                    var psms80Fragments = group.Where(g => g.Ms2ScanNumber % 4 == 3).SelectMany(g => g.MatchedIons).Select(m => m.NeutralTheoreticalProduct).ToList();
                    var psms100Fragments = group.Where(g => g.Ms2ScanNumber % 4 == 0).SelectMany(g => g.MatchedIons).Select(m => m.NeutralTheoreticalProduct).ToList();

                    var isd60Terminal = psms60Fragments.Where(n => n.IsInternalFragment == false).Select(i => i.Annotation).Distinct().ToArray();
                    var isd80Terminal = psms80Fragments.Where(n => n.IsInternalFragment == false).Select(i => i.Annotation).Distinct().ToArray();
                    var isd100Terminal = psms100Fragments.Where(n => n.IsInternalFragment == false).Select(i => i.Annotation).Distinct().ToArray();

                    var isd60UniqueTerminal = isd60Terminal.Where(i => !isd80Terminal.Contains(i) && !isd100Terminal.Contains(i)).ToList();
                    var isd80UniqueTerminal = isd80Terminal.Where(i => !isd60Terminal.Contains(i) && !isd100Terminal.Contains(i)).ToList();
                    var isd100UniqueTerminal = isd100Terminal.Where(i => !isd60Terminal.Contains(i) && !isd80Terminal.Contains(i)).ToList();

                    proteoformResult.Isd60UniqueTerminal = string.Join(",", isd60UniqueTerminal);
                    proteoformResult.Isd60FragmentCount = isd60Terminal.Length.ToString();
                    proteoformResult.Isd80UniqueTerminal = string.Join(",", isd80UniqueTerminal);
                    proteoformResult.Isd80FragmentCount = isd80Terminal.Length.ToString();
                    proteoformResult.Isd100UniqueTerminal = string.Join(",", isd100UniqueTerminal);
                    proteoformResult.Isd100FragmentCount = isd100Terminal.Length.ToString();
                }
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
