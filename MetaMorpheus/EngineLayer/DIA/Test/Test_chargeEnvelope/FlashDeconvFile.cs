using CsvHelper;
using MassSpectrometry;
using System.Collections.Generic;
using System.Data.SQLite;
using System.IO;
using System.Linq;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a FlashDeconv's ms1.tsv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class FlashDeconvFile : ResultFile<FlashDeconvRecord>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.Ms1Tsv_FlashDeconv;
        public override Software Software { get; set; }

        public FlashDeconvFile(string filePath) : base(filePath, Software.FLASHDeconv) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public FlashDeconvFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), FlashDeconvRecord.CsvConfiguration);
            Results = csv.GetRecords<FlashDeconvRecord>().ToList();
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), FlashDeconvRecord.CsvConfiguration);

            csv.WriteHeader<FlashDeconvRecord>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        //public ChargeEnvelopeIndexingEngine? GetChargeEnvelopeIndexingEngine()
        //{
        //    var envelopes = new List<IndexedChargeEnvelope>();
        //    var resultsByScan = Results.GroupBy(r => r.ZeroBasedScanNumber).OrderBy(g => g.Key);

        //    int zeroBasedScanIndex = 0;
        //    foreach (var group in resultsByScan)
        //    {
        //        foreach (var result in group)
        //        {
        //            var indexedChargeEnvelope = new IndexedChargeEnvelope(mass: result.MonoisotopicMass, intensity: result.SumIntensity, retentionTime: result.RetentionTime, result.MinCharge, result.MaxCharge, zeroBasedScanIndex);
        //            envelopes.Add(indexedChargeEnvelope);
        //        }
        //        zeroBasedScanIndex++;
        //    }

        //    int maxMass = Results.Max(r => (int)r.MonoisotopicMass);
        //    var engine = new ChargeEnvelopeIndexingEngine(maxMass);
        //    if (engine.IndexPeaks(envelopes))
        //        return engine;
        //    return null;
        //}
    }
}

