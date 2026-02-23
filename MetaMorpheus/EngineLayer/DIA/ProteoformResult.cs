using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.WebSockets;
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

        public static Dictionary<int, string> ParseModifications(string fullSeq)
        {
            // use a regex to get all modifications
            // "-?": checks if the mod starts with an optional "-", which marks the mod as an C-Terminus for position tracking.
            // "\[": indicates the start of a mod, and the end of the mod is indicated by "]"
            // "(.+?)": captures the content of the mod, which can be anything except for a closing bracket
            // "(?<!\[I+)": negative lookbehind to ensure that the closing bracket match does not correspond to a cation charge state (also defined with brackets).
            // "\]": indicates the end of the mod
            string pattern = @"-?\[(.+?)(?<!\[I+)\]";
            Regex regex = new(pattern);

            Dictionary<int, string> modDict = new();

            MatchCollection matches = regex.Matches(fullSeq);
            int totalCaptureLength = 0;
            foreach (Match match in matches)
            {
                GroupCollection group = match.Groups;
                string val = group[1].Value;
                int startIndex = group[0].Index;

                // int to add is startIndex - current position
                int positionToAddToDict = startIndex - totalCaptureLength;
                if (group[0].Value.StartsWith('-')) //if the mod starts with a "-", it is a C-Terminus mod, and the position should be incremented by 1
                {
                    positionToAddToDict++;
                }

                modDict.Add(positionToAddToDict, val);
                totalCaptureLength += group[0].Length;
            }
            return modDict;
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
            var groupByProteoform = psmTsvs.GroupBy(p => p.FullSequence).Where(g => !g.Key.Contains("|"));
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

                var mods = ProteoformResult.ParseModifications(group.First().FullSequence).Where(m => !m.Value.Contains("Fixed"));
                if (mods.Count() == 0)
                {
                    proteoformResult.Modifications = "NA";
                }
                else
                {
                    proteoformResult.Modifications = string.Join("; ", mods.Select(kvp => $"{kvp.Key}: [{string.Join(", ", kvp.Value)}]"));
                    var modSites = mods.GroupBy(kvp => kvp.Key).Select(g => g.Key).ToArray();
                    var outputStringTerminal = new StringBuilder();
                    var outputStringInternal = new StringBuilder();
                    foreach (var site in modSites)
                    {
                        var b_ions = terminalIons.Where(i => i.ProductType == ProductType.b && i.AminoAcidPosition > site);
                        var y_ions = terminalIons.Where(i => i.ProductType == ProductType.y && i.AminoAcidPosition < site);
                        var terminalIonsWithModSite = b_ions.Concat(y_ions).Select(i => i.Annotation).Distinct().ToArray();
                        outputStringTerminal.AppendLine(site + $":({terminalIonsWithModSite.Length}) " + string.Join(",", terminalIonsWithModSite));

                        var internalIonsWithModSite = internalIons.Where(i => i.FragmentNumber <= site + 1 && i.SecondaryFragmentNumber >= site + 1).Select(i => i.Annotation).Distinct().ToArray();
                        outputStringInternal.AppendLine(site + $":({internalIonsWithModSite.Length}) " + string.Join(",", internalIonsWithModSite));
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


    public class CombinedProteoform
    {
        public string ProrteinName { get; set; }
        public double MonoMass { get; set; }
        public string Modifications { get; set; }
        public string BaseSequence { get; set; }
        public string FullSequence { get; set; }
        public string TerminalIonCoverage { get; set; }
        public string InternalIonCoverage { get; set; }
        public string ModSiteTerminalIonCoverage { get; set; }
        public string ModSiteInternalIonCoverage { get; set; }

        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t"
        };

        public CombinedProteoform() { }
    }

    public class CombinedProteoformFile : ResultFile<CombinedProteoform>, IResultFile
    {
        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }

        public CombinedProteoformFile(string filePath) : base(filePath) { }

        public CombinedProteoformFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvHelper.CsvReader(new StreamReader(FilePath), CombinedProteoform.CsvConfiguration);
            Results = csv.GetRecords<CombinedProteoform>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(FilePath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvHelper.CsvWriter(new StreamWriter(File.Create(outputPath)), CombinedProteoform.CsvConfiguration))
            {
                csv.WriteHeader<CombinedProteoform>();
                foreach (var result in Results)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }

        public static void WriteCombinedProteoformResults(string outputPath, IEnumerable<string> individualResultFiles)
        {
            var allPsmFiles = individualResultFiles.Select(f => new PsmFromTsvFile(f)).ToList();
            var allProteoforms = allPsmFiles.SelectMany(p => p.Results).GroupBy(p => p.FullSequence).Where(g => !g.Key.Contains("|"));
            var dict = allPsmFiles.SelectMany(p => p.Results).GroupBy(p => p.FileNameWithoutExtension).ToDictionary(g => g.Key, g => g.ToList());
            var allResults = new List<CombinedProteoform>();
            foreach (var group in allProteoforms)
            {
                var combinedProteoform = new CombinedProteoform
                {
                    ProrteinName = group.First().Name,
                    MonoMass = group.First().MonoisotopicMass,
                    BaseSequence = group.First().BaseSequence,
                    FullSequence = group.First().FullSequence,
                };

                var terminalIonSb = new StringBuilder();
                var internalIonSb = new StringBuilder();
                foreach (var file in dict.Keys)
                {
                    var allPsms = dict[file].Where(p => p.FullSequence == group.Key).ToList();
                    var fileTerminalIons = allPsms.SelectMany(p => p.MatchedIons).Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false).Select(i => i.Annotation).Distinct().ToArray();
                    var fileInternalIons = allPsms.SelectMany(p => p.MatchedIons).Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == true && n.SecondaryFragmentNumber - n.FragmentNumber > 10).Select(i => i.Annotation).Distinct().ToArray();
                    terminalIonSb.AppendLine(file + ": " + string.Join(",", fileTerminalIons));
                    internalIonSb.AppendLine(file + ": " + string.Join(",", fileInternalIons));
                }

                var mods = ProteoformResult.ParseModifications(group.First().FullSequence).Where(m => !m.Value.Contains("Fixed"));
                if (mods.Count() == 0)
                {
                    combinedProteoform.Modifications = "NA";
                    continue;
                }
                else
                {
                    combinedProteoform.Modifications = string.Join("; ", mods.Select(kvp => $"{kvp.Key}: [{string.Join(", ", kvp.Value)}]"));
                    var modSites = mods.GroupBy(kvp => kvp.Key).Select(g => g.Key).ToArray();
                    var modSiteTerminalSb = new StringBuilder();
                    var modSiteInternalSb = new StringBuilder();
                    foreach (var site in modSites)
                    {
                        modSiteTerminalSb.AppendLine(site.ToString() + " :");
                        modSiteInternalSb.AppendLine(site.ToString() + " :");

                        foreach (var file in dict.Keys)
                        {
                            var allPsms = dict[file].Where(p => p.FullSequence == group.Key).ToList();
                            var allMatchedIons = allPsms.SelectMany(g => g.MatchedIons).ToList();
                            var terminalIons = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == false).OrderBy(i => i.FragmentNumber);
                            var internalIons = allMatchedIons.Select(i => i.NeutralTheoreticalProduct).Where(n => n.IsInternalFragment == true && n.SecondaryFragmentNumber - n.FragmentNumber > 10).OrderBy(i => i.FragmentNumber);
                            var b_ions = terminalIons.Where(i => i.ProductType == ProductType.b && i.AminoAcidPosition > site);
                            var y_ions = terminalIons.Where(i => i.ProductType == ProductType.y && i.AminoAcidPosition < site);
                            var terminalIonsWithModSite = b_ions.Concat(y_ions).Select(i => i.Annotation).Distinct().ToArray();
                            modSiteTerminalSb.AppendLine(file + $" :({terminalIonsWithModSite.Length}) " + string.Join(",", terminalIonsWithModSite));
                            var internalIonsWithModSite = internalIons.Where(i => i.FragmentNumber <= site + 1 && i.SecondaryFragmentNumber >= site + 1).Select(i => i.Annotation).Distinct().ToArray();
                            modSiteInternalSb.AppendLine(file + $" :({internalIonsWithModSite.Length}) " + string.Join(",", internalIonsWithModSite));
                        }
                        modSiteTerminalSb.AppendLine(" ");
                        modSiteInternalSb.AppendLine(" ");
                    }
                    combinedProteoform.ModSiteTerminalIonCoverage = modSiteTerminalSb.ToString();
                    combinedProteoform.ModSiteInternalIonCoverage = modSiteInternalSb.ToString();
                }
                allResults.Add(combinedProteoform);
            }

            var resultFile = new CombinedProteoformFile { Results = allResults.OrderBy(r => r.MonoMass).ToList() };
            resultFile.WriteResults(outputPath);
        }
    }
}
