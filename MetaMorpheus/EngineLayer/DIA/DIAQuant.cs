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
using Accord;

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
        public double IntegratedMs1TotalIntensity { get; set; }
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
            IntegratedMs1TotalIntensity = precursorPeakCurve.Peaks.Sum(p => p.TotalIntensity);
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

        //Take in psms to write and all pfGroups sorted
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

        //public static DIAQuantFile DDAQuantFromExistingSearch(string dataFilePath, string psmFilePath, CommonParameters commonParameters)
        //{
        //    var allPsms = PsmTsvReader.ReadTsv(psmFilePath, out List<string> warnings);
        //    var filteredPsms = allPsms.Where(psm => psm.DecoyContamTarget == "T" && psm.QValue <= 0.01 && psm.QValueNotch <= 0.01).ToList();
        //    var results = new List<DIAQuantResult>();
        //    var allMassesByScan = DeconvolutedMass.GetAllNeutralMassesByScan(ms1Scans, commonParameters.PrecursorDeconvolutionParameters);
        //    var allMasses = allMassesByScan.Where(v => v != null).SelectMany(p => p).ToList();
        //    var massTable = DeconvolutedMass.GetMassTable(allMasses, commonParameters.DIAparameters.PeakSearchBinSize);
        //    foreach (var peptide in psms)
        //    {
        //        var targetMass = allMassesByScan[peptide.PrecursorScanNumber.Value].FirstOrDefault(p => Math.Round(p.MonoisotopicMass, 3) == Math.Round(peptide.ScanPrecursorMass, 3) && p.Charge == peptide.ScanPrecursorCharge);
        //        if (targetMass == null)
        //        {
        //            continue;
        //        }
        //        var newMassCurve = MassCurve.FindMassCurve((DeconvolutedMass)targetMass, massTable, ms1Scans, null, commonParameters.DIAparameters.MaxNumMissedScan,
        //                          commonParameters.DIAparameters.Ms1PeakFindingTolerance, commonParameters.DIAparameters.PeakSearchBinSize, commonParameters.DIAparameters.MaxRTRangeMS1);
        //        var quantResult = new DIAQuantResult(peptide, newMassCurve);
        //        results.Add(quantResult);
        //    }
        //    var diaQuantFile = new DIAQuantFile { Results = results };
        //    return diaQuantFile;
        //}

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }

    public class DIAProteoformQuant
    {
        public string FileName { get; set; }
        public double PrecursorMass { get; set; }
        public string ChargeStates { get; set; }
        public string PrecursorBaseSequence { get; set; }
        public string PrecursorFullSequence { get; set; }
        public double PsmScore { get; set; }
        public string PrecursorRTs { get; set; }
        public string ElutionTimeWindow { get; set; } 
        public double IntegratedMs1TotalIntensityAllCharges { get; set; }
        public double IntegratedMs1TotalIntensityHighestThreeCharges { get; set; }
        public double IntegratedMs1TotalIntensityHighestCharge { get; set; }
        public double IntegratedHighestIsoPeakIntensityAllCharges { get; set; }
        public double IntegratedHighestIsoPeakIntensityHighestThreeCharges { get; set; }
        public double IntegratedHighestIsoPeakIntensityHighestCharge { get; set; }

        public DIAProteoformQuant(string fileName, List<DIAQuantResult> quantResults)
        {
            FileName = fileName;
            PrecursorMass = quantResults.First().PrecursorMass;
            ChargeStates = string.Join(", ", quantResults.OrderByDescending(q => q.IntegratedMs1TotalIntensity).Select(q => q.PrecursorCharge).ToList());
            PrecursorBaseSequence = quantResults.First().PrecursorBaseSequence;
            PrecursorFullSequence = quantResults.First().PrecursorFullSequence;
            PsmScore = quantResults.Max(q => q.PsmScore);
            PrecursorRTs = string.Join(", ", quantResults.Select(q => q.PrecursorRTs).ToList());
            ElutionTimeWindow = string.Join(PrecursorRTs.Min(t => t),", ", PrecursorRTs.Max(t => t));
            IntegratedMs1TotalIntensityAllCharges = quantResults.Sum(q => q.IntegratedMs1TotalIntensity);
            IntegratedMs1TotalIntensityHighestThreeCharges = quantResults.OrderByDescending(q => q.IntegratedMs1TotalIntensity).Take(3).Sum(q => q.IntegratedMs1TotalIntensity);
            IntegratedMs1TotalIntensityHighestCharge = quantResults.Max(q => q.IntegratedMs1TotalIntensity);
            IntegratedHighestIsoPeakIntensityAllCharges = quantResults.Sum(q => q.IntegratedHighestIsoPeakIntensity);
            IntegratedHighestIsoPeakIntensityHighestThreeCharges = quantResults.OrderByDescending(q => q.IntegratedHighestIsoPeakIntensity).Take(3).Sum(q => q.IntegratedHighestIsoPeakIntensity);
            IntegratedHighestIsoPeakIntensityHighestCharge = quantResults.Max(q => q.IntegratedMs1TotalIntensity);
        }

        public DIAProteoformQuant() { }
    }

    public class DIAProteoformQuantFile : ResultFile<DIAProteoformQuant>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public DIAProteoformQuantFile() : base() { }
        public DIAProteoformQuantFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<DIAProteoformQuant>().ToList();
        }

        public string FullFileName { get; set; }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<DIAProteoformQuant>();
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
