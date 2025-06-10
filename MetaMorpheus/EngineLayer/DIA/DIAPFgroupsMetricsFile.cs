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

namespace EngineLayer.DIA
{

    public class PFpairMetrics
    {
        public int PFgroupIndex { get; set; }
        public double PrecursorMass { get; set; }
        public int PrecursorCharge { get; set; }
        public double PrecursorApexRt { get; set; }
        public double PrecursorIntensity { get; set; }
        public double FragmentIntensity { get; set; }
        public double FragmentIonMz { get; set; }
        public double Correlation { get; set; }
        public double Overlap { get; set; }
        public double ApexRtDelta { get; set; }
        public string TagetDecoy { get; set; }
        public string MatchedIonType { get; set; }
        public double RoundedApexRtDelta => Math.Round(ApexRtDelta, 2);
        public double RoundedCorrelation => Math.Round(Correlation, 1);
        public double RoundedOverlap => Math.Round(Overlap, 1);
        public int Ms2Group { get; set; }
        public double PsmScore { get; set; }
        public int Ms2ScanNumber { get; set; }
        public double FragmentFractionalIntensity { get; set; } 

        public PFpairMetrics(PrecursorFragmentPair pfPair, PrecursorFragmentsGroup pfGroup, PsmFromTsv psmFromTsv = null)
        {
            Correlation = pfPair.Correlation;
            Overlap = pfPair.Overlap;
            FragmentIntensity = pfPair.FragmentPeakCurve.TotalIntensity;
            ApexRtDelta = Math.Abs(pfPair.FragmentPeakCurve.ApexRT - pfPair.PrecursorPeakCurve.ApexRT);
            PFgroupIndex = pfGroup.PFgroupIndex;
            PrecursorMass = pfGroup.PrecursorPeakCurve.MonoisotopicMass;
            PrecursorCharge = pfGroup.PrecursorPeakCurve.Charge;
            PrecursorApexRt = pfGroup.PrecursorPeakCurve.ApexRT;
            FragmentIonMz = pfPair.FragmentPeakCurve.MonoisotopicMass.ToMz(pfPair.FragmentPeakCurve.Charge);
            Ms2ScanNumber = pfPair.FragmentPeakCurve.Peaks.First().ScanNumber;
            FragmentFractionalIntensity = pfPair.FragmentPeakCurve.ApexIntensity/ pfGroup.PFpairs.Sum(pf => pf.FragmentPeakCurve.ApexIntensity);

            if (psmFromTsv != null)
            {
                if (psmFromTsv.DecoyContamTarget == "T")
                {
                    TagetDecoy = psmFromTsv.DecoyContamTarget;
                }
                else if (psmFromTsv.DecoyContamTarget == "D")
                {
                    TagetDecoy = psmFromTsv.DecoyContamTarget;
                }
                else
                {
                    TagetDecoy = "NA";
                }

                PsmScore = psmFromTsv.Score;
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
                    TagetDecoy = "D";
                    break;
                case false:
                    TagetDecoy = "T";
                    break;
            }
            PsmScore = psm.Score;

            var sortedMatchedIons = psm.MatchedFragmentIons.OrderBy(ion => ion.Mz).ToArray();
            var roundedMz = sortedMatchedIons.Select(i => Math.Round(i.Mz, 2)).ToArray();
            var index = Array.BinarySearch(roundedMz, Math.Round(FragmentIonMz, 2));
            if (index >= 0)
            {
                MatchedIonType = sortedMatchedIons[index].NeutralTheoreticalProduct.SecondaryProductType == null? "Terminal" : "Internal";
            }
            else
            {
                MatchedIonType = "NA";
            }
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

        public static PFpairMetricFile GetPFpairsFromPsms(List<PrecursorFragmentsGroup> pfGroups, SpectralMatch[] sortedPsms)
        {
            List<PFpairMetrics> results = new List<PFpairMetrics>();
            var sortedScanNumberArray = sortedPsms.Select(psm => psm.ScanNumber).ToArray();
            foreach (var pfGroup in pfGroups)
            {
                int index = Array.BinarySearch(sortedScanNumberArray, pfGroup.PFgroupIndex);
                if (index >= 0)
                {
                    foreach (var pfPair in pfGroup.PFpairs)
                    {
                        var pfPairMetrics = new PFpairMetrics(pfPair, pfGroup);
                        pfPairMetrics.SetPsmInfo(sortedPsms[index]);
                        results.Add(pfPairMetrics);
                    }
                }
            }
            var pfPairMetricFile = new PFpairMetricFile()
            {
                Results = results
            };
            return pfPairMetricFile;
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
