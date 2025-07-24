using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using System.Globalization;
using EngineLayer;
using Omics.Fragmentation;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using LinqStatistics;
using System.Windows.Markup;
using Omics.SpectrumMatch;
using Easy.Common.Extensions;
using System.Security.Cryptography.X509Certificates;
using Chemistry;
using MzLibUtil;

namespace EngineLayer.DIA
{

    public class PFpairMetrics
    {
        public int PFgroupIndex { get; set; }
        public double PrecursorMass { get; set; }
        public int PrecursorCharge { get; set; }
        public double PrecursorApexRt { get; set; }
        public double PrecursorIntensity { get; set; }
        public double FragmentMass { get; set; }
        public int FragmentCharge { get; set; } 
        public double FragmentIntensity { get; set; }
        [Optional]
        public double FragmentIonMz { get; set; }
        public double Correlation { get; set; }
        public double Overlap { get; set; }
        public double ApexRtDelta { get; set; }
        public string TargetDecoy { get; set; }
        public string MatchedIonType { get; set; }
        public double RoundedApexRtDelta => Math.Round(ApexRtDelta, 2);
        public double RoundedCorrelation => Math.Round(Correlation, 1);
        public double RoundedOverlap => Math.Round(Overlap, 1);
        public int Ms2Group { get; set; }
        public double PsmScore { get; set; }
        [Optional]
        public double RoundedPsmScore => Math.Round(PsmScore, 0);
        public int Ms2ScanNumber { get; set; }
        public double FragmentFractionalIntensity { get; set; }
        [Optional]
        public int NormalizedIntensityRank { get; set; }
        [Optional]
        public double SharedXIC { get; set; }
        [Optional]
        public int PrecursorRank { get; set; }
        [Optional]
        public int FragmentRank { get; set; } 

        public PFpairMetrics(PrecursorFragmentPair pfPair, PrecursorFragmentsGroup pfGroup = null, SpectralMatch psm = null,PsmFromTsv psmFromTsv = null)
        {
            Correlation = pfPair.Correlation;
            Overlap = pfPair.Overlap;
            FragmentMass = pfPair.FragmentPeakCurve.MonoisotopicMass;
            FragmentCharge = pfPair.FragmentPeakCurve.Charge;
            FragmentIntensity = pfPair.FragmentPeakCurve.TotalIntensity;
            ApexRtDelta = Math.Abs(pfPair.FragmentPeakCurve.ApexRT - pfPair.PrecursorPeakCurve.ApexRT);
            PrecursorMass = pfPair.PrecursorPeakCurve.MonoisotopicMass;
            PrecursorCharge = pfPair.PrecursorPeakCurve.Charge;
            PrecursorApexRt = pfPair.PrecursorPeakCurve.ApexRT;
            if (pfPair.FragmentPeakCurve.MonoisotopicMass > 0)
            {
                FragmentIonMz = pfPair.FragmentPeakCurve.MonoisotopicMass.ToMz(pfPair.FragmentPeakCurve.Charge);
            }
            else
            {
                FragmentIonMz = pfPair.FragmentPeakCurve.AveragedMz;
            }
            
            Ms2ScanNumber = pfPair.FragmentPeakCurve.Peaks.First().ScanNumber;
            NormalizedIntensityRank = pfPair.NormalizedIntensityRank;
            SharedXIC = pfPair.SharedXIC;
            PrecursorRank = pfPair.PrecursorRank;
            FragmentRank = pfPair.FragmentRank;

            if (pfGroup != null)
            {
                PFgroupIndex = pfGroup.PFgroupIndex;
                FragmentFractionalIntensity = pfPair.FragmentPeakCurve.ApexIntensity / pfGroup.PFpairs.Sum(pf => pf.FragmentPeakCurve.ApexIntensity);
            }
            if (psm != null)
            {
                SetPsmInfo(psm);
            }

            if (psmFromTsv != null)
            {
                SetPsmTsvInfo(psmFromTsv);
            }
        }

        public void SetMs2Group(int ms2Group)
        {
            Ms2Group = ms2Group;
        }

        public void SetPsmInfo(SpectralMatch psm)
        {
            switch (psm.IsDecoy)
            {
                case true:
                    TargetDecoy = "D";
                    break;
                case false:
                    TargetDecoy = "T";
                    break;
            }
            PsmScore = psm.Score;

            var sortedMatchedIons = psm.MatchedFragmentIons.OrderBy(ion => ion.Mz).ToArray();
            var roundedMz = sortedMatchedIons.Select(i => Math.Round(i.Mz, 1)).ToArray();
            var index = Array.BinarySearch(roundedMz, Math.Round(FragmentIonMz, 1));
            if (index >= 0)
            {
                MatchedIonType = sortedMatchedIons[index].NeutralTheoreticalProduct.SecondaryProductType == null? "Terminal" : "Internal";
                if (FragmentCharge == 0)
                {
                    FragmentCharge = sortedMatchedIons[index].Charge;
                }
            }
            else
            {
                MatchedIonType = "NA";
            }
        }

        public void SetPsmTsvInfo(PsmFromTsv psmFromTsv)
        {
            if (psmFromTsv.DecoyContamTarget == "T")
            {
                TargetDecoy = psmFromTsv.DecoyContamTarget;
            }
            else if (psmFromTsv.DecoyContamTarget == "D")
            {
                TargetDecoy = psmFromTsv.DecoyContamTarget;
            }
            else
            {
                TargetDecoy = "NA";
            }
            PsmScore = psmFromTsv.Score;
        }

        public static void NeutralLossReSearchFromPFMetrics_neutralMass(IEnumerable<PFpairMetrics> pFpairMetricsList, List<double> neutralLosses, Tolerance ppmTolerance)
        {
            var allTerminalFragments = pFpairMetricsList.Where(p => p.MatchedIonType == "Terminal").ToList();
            var unmatchedFragments = pFpairMetricsList.Where(p => p.MatchedIonType == "NA").ToList();
            foreach (var frag in allTerminalFragments)
            {
                foreach (var loss in neutralLosses)
                {
                    double massToSearch = frag.FragmentMass - loss;
                    var bestFrag = PFpairMetrics.FindPfPair(unmatchedFragments, massToSearch, frag.FragmentCharge, ppmTolerance);
                    if (bestFrag != null)
                    {
                        bestFrag.MatchedIonType = "NeutralLoss";
                    }
                }
            }
        }

        public static void IsotopeReSearch(List<PFpairMetrics> pfPairMetrics, DIAparameters diaParam)
        {
            var allMatchedIons = pfPairMetrics.Where(pf => pf.MatchedIonType == "Terminal").ToList();
            var unmatchedIons = pfPairMetrics.Where(pf => pf.MatchedIonType == "NA").ToList();
            foreach(var pair in allMatchedIons)
            {
                var isotopes = SearchIsotopesOfOnePeak(unmatchedIons, pair.FragmentIonMz, pair.FragmentCharge, diaParam.NumIsotopesToSearch, diaParam.Ms2PeakFindingTolerance);
                if (isotopes.Count > 0)
                    foreach(var isotope in isotopes) isotope.MatchedIonType = "Isotope";
            }
        }
        public static List<PFpairMetrics> SearchIsotopesOfOnePeak(IEnumerable<PFpairMetrics> allPairs, double knownMz, int knownCharge, int numIsotopes, Tolerance ppmTolerance)
        {
            var foundPairs = new List<PFpairMetrics>();
            for (int i = 1; i <= numIsotopes; i++)
            {
                double isotopeMz = (knownMz.ToMass(knownCharge) + i * Constants.C13MinusC12).ToMz(knownCharge);
                var pair = FindPfPair(allPairs, isotopeMz, 0, ppmTolerance);
                if (pair != null)
                {
                    foundPairs.Add(pair);
                }
            }
            return foundPairs;
        }
        public static PFpairMetrics FindPfPair(IEnumerable<PFpairMetrics> allPairs, double targetM, int targetCharge, Tolerance ppmTolerance)
        {
            PFpairMetrics frag = null;
            if (targetCharge > 0)
            {
                frag = allPairs.Where(f => f.FragmentCharge == targetCharge && ppmTolerance.Within(targetM, f.FragmentMass))
                                      .OrderBy(f => Math.Abs(f.FragmentMass - targetM)).FirstOrDefault();
            }
            else
            {
                frag = allPairs.Where(f => ppmTolerance.Within(targetM, f.FragmentIonMz)).OrderByDescending(f => f.FragmentIntensity).FirstOrDefault();
            }
            return frag;
        }

        public PFpairMetrics() { }
    }

    public class PFpairMetricFile : ResultFile<PFpairMetrics>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public PFpairMetricFile() : base() { }
        public PFpairMetricFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<PFpairMetrics>().ToList();
        }

        //public static PFpairMetricFile GetPFpairsFromPfGroupAndPsms(List<PrecursorFragmentsGroup> pfGroups, SpectralMatch[] sortedPsms)
        //{
        //    List<PFpairMetrics> results = new List<PFpairMetrics>();
        //    var sortedScanNumberArray = sortedPsms.Select(psm => psm.ScanNumber).ToArray();
        //    foreach (var pfGroup in pfGroups)
        //    {
        //        int index = Array.BinarySearch(sortedScanNumberArray, pfGroup.PFgroupIndex);
        //        if (index >= 0)
        //        {
        //            foreach (var pfPair in pfGroup.PFpairs)
        //            {
        //                var pfPairMetrics = new PFpairMetrics(pfPair, pfGroup);
        //                pfPairMetrics.SetPsmInfo(sortedPsms[index]);
        //                results.Add(pfPairMetrics);
        //            }
        //        }
        //    }
        //    var pfPairMetricFile = new PFpairMetricFile()
        //    {
        //        Results = results
        //    };
        //    return pfPairMetricFile;
        //}

        //only write out PFpairs that have a PSM associated with the pfGroup
        public static PFpairMetricFile GetAllPFpairsFromPfGroupAndPsms(List<PrecursorFragmentsGroup> pfGroups, SpectralMatch[] sortedPsms, DIAparameters diaParam)
        {
            List<PFpairMetrics> results = new List<PFpairMetrics>();
            var sortedScanNumberArray = sortedPsms.Select(psm => psm.ScanNumber).ToArray();
            foreach (var pfGroup in pfGroups)
            {
                var pfPairMetricsForThisGroup = new List<PFpairMetrics>();
                int index = Array.BinarySearch(sortedScanNumberArray, pfGroup.PFgroupIndex);
                if (index >= 0)
                {
                    foreach (var pfPair in pfGroup.PFpairs)
                    {
                        var pfPairMetrics = new PFpairMetrics(pfPair, pfGroup, sortedPsms[index]);
                        results.Add(pfPairMetrics);
                        pfPairMetricsForThisGroup.Add(pfPairMetrics);
                    }
                    //if we want to research isotopes
                    if (diaParam.NumIsotopesToSearch != 0)
                    {
                        PFpairMetrics.IsotopeReSearch(pfPairMetricsForThisGroup, diaParam);
                    }
                    //if we want to research the neutral loss fragments
                    if (diaParam.NeutralLossSearch && index >= 0 && sortedPsms[index].IsDecoy == false)
                    {
                        //var allQualifiedPFpairs = pfPairMetricsForThisGroup.Where(p => p.TargetDecoy == "T" && p.MatchedIonType == "Terminal").ToList();
                        PFpairMetrics.NeutralLossReSearchFromPFMetrics_neutralMass(pfPairMetricsForThisGroup, new List<double> { 18.010564684, 17.0265491 }, diaParam.Ms2PeakFindingTolerance);
                    }
                }
            }
            var pfPairMetricFile = new PFpairMetricFile()
            {
                Results = results
            };
            return pfPairMetricFile;
        }

        public static PFpairMetricFile NeutralLossResearchFromPFpairMetricsFile(PFpairMetricFile pfPairMetricsFile, List<double> neutralLosses)
        {
            var results = pfPairMetricsFile.Results;
            var groupedResults = results.GroupBy(p => p.PFgroupIndex).ToList();
            foreach (var group in groupedResults)
            {
                if (group.First().TargetDecoy == "D")
                {
                    continue;
                }
                var allTerminalFragments = group.Where(p => p.MatchedIonType == "Terminal").ToList();
                var sortedFragments = group.Where(p => p.MatchedIonType == "NA").OrderBy(p => p.FragmentMass).ToList();
                foreach (var frag in allTerminalFragments)
                {
                    foreach (var loss in neutralLosses)
                    {
                        double massToSearch = frag.FragmentMass - loss;
                        var bestFrag = PFpairMetrics.FindPfPair(sortedFragments, massToSearch, frag.FragmentCharge, new PpmTolerance(10));
                        if (bestFrag != null)
                        {
                            bestFrag.MatchedIonType = "NeutralLoss";
                        }
                    }
                }
            }
            var newPfPairMetricsFile = new PFpairMetricFile();
            newPfPairMetricsFile.Results = results;
            return newPfPairMetricsFile;
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<PFpairMetrics>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }

    //PFgroup
    public class PFgroupMetrics
    {
        public double PrecursorMass { get; set; }
        public double PrecursorCharge{ get; set; }
        public double PrecursorMz { get; set; }
        public double PrecursorIntensity { get; set; }
        public double PrecursorApexRt { get; set; }
        public double PrecursorNumOfPoints { get; set; }
        public int NumberOfFragments { get; set; }
        public double MedianFragmentsCorrelation { get; set; }
        public double AverageFragmentsCorrelation { get; set; }
        public double AverageFragmentIntensity { get; set; }
        public int NumberOfHighCorrelationFragments { get; set; }
        public double MedianApexRtDelta { get; set; }
        public double MedianOverlap { get; set; }
        public int PFgroupIndex { get; set; }
        public string TargetDecoy { get; set; }
        public string FullSequence { get; set; }
        public double PsmScore { get; set; }
        public double PsmQValue { get; set; }
        public int Ms2Group { get; set; }

        public PFgroupMetrics(PrecursorFragmentsGroup pfGroup, PsmFromTsv psmTsv = null)
        {
            PrecursorMass = pfGroup.PrecursorPeakCurve.MonoisotopicMass;
            PrecursorCharge = pfGroup.PrecursorPeakCurve.Charge;
            PrecursorMz = pfGroup.PrecursorPeakCurve.Peaks.OrderByDescending(p => p.Intensity).First().HighestPeakMz;
            PrecursorIntensity = pfGroup.PrecursorPeakCurve.ApexIntensity;
            PrecursorApexRt = pfGroup.PrecursorPeakCurve.ApexRT;
            PrecursorNumOfPoints = pfGroup.PrecursorPeakCurve.Peaks.Count;
            NumberOfFragments = pfGroup.PFpairs.Count;
            try
            {
                MedianFragmentsCorrelation = pfGroup.PFpairs.Select(pf => pf.Correlation).Median();
            } catch (Exception e)
            {
                int stop = 1;
            }
            
            AverageFragmentsCorrelation = pfGroup.PFpairs.Select(pf => pf.Correlation).Average();
            AverageFragmentIntensity = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.TotalIntensity).Average();
            NumberOfHighCorrelationFragments = pfGroup.PFpairs.Count(pf => pf.Correlation > 0.9);
            MedianApexRtDelta = pfGroup.PFpairs.Select(pf => Math.Abs(pf.FragmentPeakCurve.ApexRT - pf.PrecursorPeakCurve.ApexRT)).Median();
            MedianOverlap = pfGroup.PFpairs.Select(pf => pf.Overlap).Median();
            PFgroupIndex = pfGroup.PFgroupIndex;
            Ms2Group = pfGroup.PFpairs.First().FragmentPeakCurve.Peaks.First().ScanNumber % 4;
            if (psmTsv != null)
            {
                TargetDecoy = psmTsv.DecoyContamTarget;
                PsmScore = psmTsv.Score;
                PsmQValue = psmTsv.QValue;
                FullSequence = psmTsv.FullSequence;
            }
        }

        public void SetTargetDecoy(SpectralMatch psm)
        {
            switch (psm.IsDecoy)
            {
                case true:
                    TargetDecoy = "D";
                    break;
                case false:
                    TargetDecoy = "T";
                    break;
            }
            FullSequence = psm.FullSequence;
            PsmScore = psm.Score;
            PsmQValue = psm.FdrInfo.QValue;
        }
        public PFgroupMetrics() { }
    }

    public class DIAPFgroupsMetricsFile: ResultFile<PFgroupMetrics>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public DIAPFgroupsMetricsFile() : base() { }
        public DIAPFgroupsMetricsFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<PFgroupMetrics>().ToList();
        }

        public string FullFileName { get; set; }
        public static DIAPFgroupsMetricsFile GetFromPFgroups(List<PrecursorFragmentsGroup> pfGroups, string fullFileName, SpectralMatch[] psms = null)
        {
            List<PFgroupMetrics> results = new List<PFgroupMetrics>();
            foreach (var group in pfGroups)
            {
                var pfGroupMetrics = new PFgroupMetrics(group);
                results.Add(pfGroupMetrics);
            }
            var diaPFgroupsMetricsFile = new DIAPFgroupsMetricsFile()
            {
                Results = results,
                FullFileName = fullFileName
            };
            return diaPFgroupsMetricsFile;
        }

        public void SetTargetDecoyForPFgroupsFromPsms(SpectralMatch[] sortedPsms)
        {
            var sortedPsmScanNumbers = sortedPsms.Select(psm => psm.ScanNumber).ToArray();
            foreach (var pfGroupMetrics in Results)
            {
                int index = Array.BinarySearch(sortedPsmScanNumbers, pfGroupMetrics.PFgroupIndex);
                if (index >= 0)
                {
                    pfGroupMetrics.SetTargetDecoy(sortedPsms[index]);
                }
                else
                {
                    pfGroupMetrics.TargetDecoy = "NA";
                }
            }
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<PFgroupMetrics>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }

}
