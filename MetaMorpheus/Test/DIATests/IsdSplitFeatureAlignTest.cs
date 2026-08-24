using EngineLayer.DIA.FeatureAlign;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Builds 8-run FeatureAlign consensus tables from the MS1-CHANNEL SPLIT TopFD runs
    /// (SeparateISD\TopFD\*_ms1_ms1.feature).
    ///
    /// Why not reuse the existing FeatureAlign_*_ISD_8files tables: those were aligned over the
    /// WHOLE-PARENT ISD TopFD run (max feature idx 24,278), while the pseudo-MS2 precursors - and hence
    /// every PRECURSOR_FEATURE_ID written into the msalign - come from the MS1-split run (5,430 features
    /// in rep1). Joining the two gives a median mass deviation of ~5.5e5 ppm, i.e. nonsense. Aligning the
    /// split features instead makes feature id -> aligned row an exact integer join, matching how the DDA
    /// analysis works.
    ///
    /// Tolerances match the DDA/ISD precedent: 20 ppm, 120 s, +/-2 monoisotope notches, ReferenceAnchor.
    /// </summary>
    [TestFixture]
    public static class IsdSplitFeatureAlignTest
    {
        [Test]
        [Explicit("manual runner over local TopFD output")]
        public static void AlignSplitMs1Features()
        {
            string topfdDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\SeparateISD\TopFD";
            string outDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\SeparateISD\FeatureAlign";
            Directory.CreateDirectory(outDir);

            const double massTol = 20;
            const double rtTol = 120;   // seconds

            var all = Directory.GetFiles(topfdDir, "*_ms1_ms1.feature").ToList();
            Assert.That(all.Count, Is.EqualTo(10), "expected 10 split MS1 feature tables");

            // the six mixtures appear in BOTH families, exactly as in the DDA design
            var mix = all.Where(f => Path.GetFileName(f).Contains("HCD-YBC")).OrderBy(f => f).ToList();
            Assert.That(mix.Count, Is.EqualTo(6));

            var families = new (string Name, List<string> Files)[]
            {
                ("HCD_ISD_split_8files",
                 all.Where(f => Path.GetFileName(f).Contains("_HCD_2-3_")).OrderBy(f => f).Concat(mix).ToList()),
                ("YBC_ISD_split_8files",
                 all.Where(f => Path.GetFileName(f).Contains("_YBC_1-1_")).OrderBy(f => f).Concat(mix).ToList()),
            };

            foreach (var (name, files) in families)
            {
                Assert.That(files.Count, Is.EqualTo(8), $"{name}: expected 8 files, got {files.Count}");
                string outPath = Path.Combine(outDir, $"FeatureAlign_{name}_20ppm_120s_notch2.tsv");
                string written = FeatureAligner.Align(files, outPath, massTol, MassToleranceUnit.Ppm, rtTol,
                    grouping: GroupingMethod.ReferenceAnchor, positiveNotch: 2, negativeNotch: 2);

                var lines = File.ReadAllLines(written);
                int all8 = lines.Skip(1).Count(l => l.Split('\t')[2] == "8");
                TestContext.Progress.WriteLine($"{name}: {lines.Length - 1:N0} rows, {all8:N0} in all 8 -> {written}");
            }
        }
    }
}
