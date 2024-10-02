using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using PsmFromTsv = EngineLayer.PsmFromTsv;
using Chemistry;
using System.IO;
using TaskLayer;

namespace Test
{
    public class TestPrecursorHighestPeakMz
    {
        [Test]
        public static void TestPreHighestPeakMz()
        {
            string file = @"E:\ISD Project\Targeted\fd_incList3.psmtsv";
            var psms = PsmTsvReader.ReadTsv(file, out var warnings);
            var outputList = new List<(string sequence, double mz, double? startRT, double? EndRT)>();
            foreach(var psm in psms)
            {
                var highestPeakMz = psm.PrecursorHighestPeakMz;
                var startRT = psm.RetentionTime - 2.5;
                var endRT = psm.RetentionTime + 2.5;
                var seq = psm.FullSequence;
                if (psm.RetentionTime < 3)
                {
                    startRT = 0;
                }
                var charges = new List<int> { psm.PrecursorCharge - 1, psm.PrecursorCharge, psm.PrecursorCharge + 1 };
                foreach (var charge in charges)
                {
                    outputList.Add((seq, highestPeakMz.ToMass(psm.PrecursorCharge).ToMz(charge), startRT, endRT));
                }
            }

            var outputTsv = @"E:\ISD Project\Targeted\fd_inclusionList3.csv";
            using (var sw = new StreamWriter(File.Create(outputTsv)))
            {
                sw.WriteLine("Compound, Mz, t start (min), t stop (min)");
                foreach (var id in outputList)
                {
                    sw.WriteLine($"{id.sequence},{id.mz},{id.startRT},{id.EndRT}");
                }
            }
        }

        [Test]
        public static void TestISDFileReading()
        {
            var file = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionB_orbiMS1_ISD60-80-100.mzML";
            var myFM = new MyFileManager(false);
            var commonParameters = new CommonParameters();
            var myMSDataFile = myFM.LoadFile(file, commonParameters);
            var ms2scans = myMSDataFile.GetAllScansList().FindAll(p => p.MsnOrder == 2);
        }
    }
}
