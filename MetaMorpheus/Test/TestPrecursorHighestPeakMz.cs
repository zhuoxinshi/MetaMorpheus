using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using PsmFromTsv = EngineLayer.PsmFromTsv;
using Chemistry;

namespace Test
{
    public class TestPrecursorHighestPeakMz
    {
        [Test]
        public static void TestPreHighestPeakMz()
        {
            string file = @"E:\ISD Project\Targeted\unique_isd_fb_new.tsv";
            var psms = PsmTsvReader.ReadTsv(file, out var warnings);
            var outputList = new List<(string sequence, double mz, double? startRT, double? EndRT)>();
            foreach(var psm in psms)
            {
                var highestPeakMz = psm.PrecursorHighestPeakMz;
                var startRT = psm.RetentionTime - 5;
                var endRT = psm.RetentionTime + 5;
                var seq = psm.FullSequence;
                if (psm.RetentionTime < 5)
                {
                    startRT = 0;
                }
                var charges = new List<int> { psm.PrecursorCharge - 1, psm.PrecursorCharge, psm.PrecursorCharge + 1 };
                foreach (var charge in charges)
                {
                    outputList.Add((seq, highestPeakMz.ToMass(psm.PrecursorCharge).ToMz(charge), startRT, endRT));
                }
            }
        }
    }
}
