using EngineLayer.DIA.FeatureAlign;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// 8-run FeatureAlign over TopDIA's own MS1 feature tables (*_3windows_ms1.feature), so that
    /// PRECURSOR_FEATURE_ID in the TopDIA msalign joins exactly to an aligned row - the same construction
    /// used for the XicGrouping split features, letting the two pipelines be compared on equal footing.
    /// HCD family: 2 pure HCD + 6 HCD-YBC mixtures. Tolerances match: 20 ppm, 120 s, +/-2 notch.
    /// </summary>
    [TestFixture]
    public static class TopDiaFeatureAlignTest
    {
        [Test]
        [Explicit("manual runner over local TopDIA output")]
        public static void AlignTopDiaMs1Features()
        {
            string raw = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\TopDIA\Raw";
            string outDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\TopDIA\FeatureAlign";
            Directory.CreateDirectory(outDir);

            var all = Directory.GetFiles(raw, "*_3windows_ms1.feature").ToList();
            Assert.That(all.Count, Is.EqualTo(10));
            var mix = all.Where(f => Path.GetFileName(f).Contains("HCD-YBC")).OrderBy(f => f).ToList();
            var files = all.Where(f => Path.GetFileName(f).Contains("_HCD_2-3_")).OrderBy(f => f).Concat(mix).ToList();
            Assert.That(files.Count, Is.EqualTo(8));

            string outPath = Path.Combine(outDir, "FeatureAlign_TopDIA_HCD_8files_20ppm_120s_notch2.tsv");
            string written = FeatureAligner.Align(files, outPath, 20, MassToleranceUnit.Ppm, 120,
                grouping: GroupingMethod.ReferenceAnchor, positiveNotch: 2, negativeNotch: 2);
            var lines = File.ReadAllLines(written);
            int all8 = lines.Skip(1).Count(l => l.Split('\t')[2] == "8");
            TestContext.Progress.WriteLine($"TopDIA HCD: {lines.Length - 1:N0} rows, {all8:N0} in all 8 -> {written}");
        }
    }
}
