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
using MassSpectrometry;
using static Plotly.NET.StyleParam.Range;

namespace EngineLayer.DIA
{
    public class DIAQuantResult
    {
        public int PFGroupIndex { get; set; }
        public double PrecursorMass { get; set; }
        public double PrecursorCharge { get; set; }
        public double PrecursorHighestMz { get; set; }
        public double PrecursorApexIntensity { get; set; }
        public double PrecursorApexRt { get; set; }
        public double PrecursorNumOfPoints { get; set; }
        public string PrecursorBaseSequence { get; set; }
        public string PrecursorFullSequence { get; set; }
        public double PsmScore { get; set; }
        public double PsmQvalue { get; set; }
        public string PrecursorRTs { get; set; }
        public string PrecursorHighestIsotopeIntensities { get; set; }
        public string PrecursorEnvelopeTotalIntensities { get; set; }
        public double IntegratedMs1Intensity { get; set; }
        public double IntegratedHighestIsoPeakIntensity { get; set; }

        public DIAQuantResult (SpectralMatch psm, PeakCurve precursorPeakCurve)
        {
            PFGroupIndex = psm.ScanNumber;
            PrecursorMass = psm.ScanPrecursorMass;
            PrecursorHighestMz = psm.PrecursorHighestPeakMz;
            PrecursorCharge = psm.ScanPrecursorCharge;
            PrecursorApexIntensity = precursorPeakCurve.ApexIntensity;
            PrecursorApexRt = precursorPeakCurve.ApexRT;
            PrecursorNumOfPoints = precursorPeakCurve.Peaks.Count;
            PrecursorBaseSequence = psm.BaseSequence;
            PrecursorFullSequence = psm.FullSequence;
            PsmScore = psm.Score;
            PsmQvalue = psm.PsmFdrInfo.QValue;
            PrecursorRTs = string.Join(", ", precursorPeakCurve.Peaks.Select(p => p.RetentionTime).ToList());
            PrecursorHighestIsotopeIntensities = string.Join(", ", precursorPeakCurve.Peaks.Select(p => p.HighestPeakIntensity).ToList());
            PrecursorEnvelopeTotalIntensities = string.Join(", ", precursorPeakCurve.Peaks.Select(p => p.TotalIntensity).ToList());
            IntegratedMs1Intensity = precursorPeakCurve.Peaks.Sum(p => p.TotalIntensity);
            IntegratedHighestIsoPeakIntensity = precursorPeakCurve.Peaks.Sum(p => p.HighestPeakIntensity);
        } 
        public DIAQuantResult() { }
    }

    public class DIAQuantFile : ResultFile<DIAQuantResult>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public DIAQuantFile() : base() { }
        public DIAQuantFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<DIAQuantResult>().ToList();
        }

        public string FullFileName { get; set; }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<DIAQuantResult>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public static DIAQuantFile CreateDIAQuantFile(SpectralMatch[] psms, PrecursorFragmentsGroup[] sortedPfGroups)
        {
            var results = new List<DIAQuantResult>();
            var sortedPfGroupIndex = sortedPfGroups.Select(g => g.PFgroupIndex).ToArray();
            foreach (var peptide in psms)
            {
                int index = Array.BinarySearch(sortedPfGroupIndex, peptide.ScanNumber);
                if (index >= 0)
                {
                    var quantResult = new DIAQuantResult(peptide, sortedPfGroups[index].PrecursorPeakCurve);
                    results.Add(quantResult);
                }
            }
            var diaQuantFile = new DIAQuantFile { Results = results };
            return diaQuantFile;
        }

        public static DIAQuantFile DDAQuant(SpectralMatch[] psms, MsDataScan[] ms1Scans, CommonParameters commonParameters)
        {
            var results = new List<DIAQuantResult>();
            var allMassesByScan = DeconvolutedMass.GetAllNeutralMassesByScan(ms1Scans, commonParameters.PrecursorDeconvolutionParameters);
            var allMasses = allMassesByScan.Where(v => v != null).SelectMany(p => p).ToList();
            var massTable = DeconvolutedMass.GetMassTable(allMasses, commonParameters.DIAparameters.PeakSearchBinSize);
            foreach (var peptide in psms)
            {
                var targetMass = allMassesByScan[peptide.PrecursorScanNumber.Value].FirstOrDefault(p => Math.Round(p.MonoisotopicMass, 3) == Math.Round(peptide.ScanPrecursorMass, 3) && p.Charge == peptide.ScanPrecursorCharge);
                if (targetMass == null)
                {
                    continue;
                }
                var newMassCurve = MassCurve.FindMassCurve((DeconvolutedMass)targetMass, massTable, ms1Scans, null, commonParameters.DIAparameters.MaxNumMissedScan,
                                  commonParameters.DIAparameters.Ms1PeakFindingTolerance, commonParameters.DIAparameters.PeakSearchBinSize, commonParameters.DIAparameters.MaxRTRangeMS1);
                var quantResult = new DIAQuantResult(peptide, newMassCurve);
                results.Add(quantResult);
            }
            var diaQuantFile = new DIAQuantFile { Results = results };
            return diaQuantFile;
        }
        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }

}
