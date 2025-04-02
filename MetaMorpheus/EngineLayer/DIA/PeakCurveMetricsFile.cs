using CsvHelper.Configuration;
using CsvHelper;
using Readers;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace EngineLayer.DIA
{
    public class PeakCurveMetrics
    {
        public double MonoisotopicMass { get; set; }
        public double Charge { get; set; }
        public double Mz { get; set; }
        public int MsLevel { get; set; }
        public double ApexRt { get; set; }
        public double ApexScanNumber { get; set; }
        public double ApexIntensity { get; set; }
        public double StartRt { get; set; }
        public double EndRt { get; set; }
        public double StartIntensity { get; set; }
        public double EndIntensity { get; set; }
        public int NumOfPoints { get; set; }

        public PeakCurveMetrics(PeakCurve pc)
        {
            MonoisotopicMass = pc.MonoisotopicMass;
            Charge = pc.Charge;
            Mz = pc.AveragedMz == 0? pc.MonoisotopicMass.ToMz(pc.Charge) : pc.AveragedMz;
            MsLevel = pc.MsLevel;
            ApexRt = pc.ApexRT;
            ApexScanNumber = pc.Peaks.OrderByDescending(p => p.Intensity).First().ScanNumber;
            ApexIntensity = pc.ApexIntensity;
            StartRt = pc.StartRT;
            EndRt = pc.EndRT;
            StartIntensity = pc.Peaks.First().Intensity;
            EndIntensity = pc.Peaks.Last().Intensity;
            NumOfPoints = pc.Peaks.Count;
        }
        public PeakCurveMetrics() { }
    }
    public class PeakCurveMetricsFile : ResultFile<PeakCurveMetrics>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public PeakCurveMetricsFile() : base() { }
        public PeakCurveMetricsFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<PeakCurveMetrics>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<PeakCurveMetrics>();
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
