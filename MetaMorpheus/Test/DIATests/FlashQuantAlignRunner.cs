using EngineLayer.DIA;
using EngineLayer.DIA.FeatureAlign;
using NUnit.Framework;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Runs OUR <see cref="FeatureAligner"/> on the three FLASHQuant feature tables and writes the consensus into
    /// the FlashQuant folder, so it can be diffed against OpenMS ConsensusFeatureGroupDetector's
    /// <c>flashquant_consensus.tsv</c>.
    ///
    /// Tolerances are set to MATCH OpenMS ConsensusFeatureGroupDetector's defaults so the two are comparable:
    ///   mass 100 ppm, RT 180 s, ReferenceAnchor grouping, keep-most-abundant on a duplicate, and notch 0
    ///   (OpenMS has no monoisotope-notch concept). Both tools align on MedianApexRetentionTime.
    ///
    /// Remember the two tools report different scopes: OpenMS emits ONLY groups present in every file, so compare
    /// its rows against this file's rows with n_files == 3.
    ///
    /// [Explicit] so it only runs when you select it. Requires the three *_flashquant.tsv to exist first.
    /// </summary>
    [TestFixture]
    public class FlashQuantAlignRunner
    {
        [Test, Explicit("aligns the 3 FLASHQuant tables with our aligner (OpenMS-matched params) -> FlashQuant folder")]
        public void AlignFlashQuant_MatchOpenMS()
        {
            var dir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashQuant";
            var files = Directory.GetFiles(dir, "HCD-YBC_1-1_rep*_flashquant.tsv").OrderBy(f => f).ToList();
            Assert.That(files.Count, Is.EqualTo(3),
                "expected the 3 HCD-YBC_1-1_rep*_flashquant.tsv in " + dir + " (run FLASHQuant first)");
            foreach (var f in files) TestContext.Out.WriteLine("input: " + Path.GetFileName(f));

            string outPath = Path.Combine(dir, "ours_reference_anchor_consensus.tsv");
            FeatureAligner.Align(
                files, outPath,
                massTolerance: 100, massUnit: MassToleranceUnit.Ppm, rtTolerance: 180,
                preferredQuantColumn: "FeatureGroupQuantity",   // OpenMS quant_method default
                positiveNotch: 0, negativeNotch: 0,             // OpenMS has no notch
                grouping: GroupingMethod.ReferenceAnchor,
                whenDuplicate: DuplicateResolution.MaxAbundance);

            var lines = File.ReadAllLines(outPath);
            int nf = System.Array.IndexOf(lines[0].Split('\t'), "n_files");
            var nFiles = lines.Skip(1).Where(l => l.Length > 0).Select(l => int.Parse(l.Split('\t')[nf])).ToList();
            var dist = nFiles.GroupBy(x => x).OrderBy(g => g.Key).Select(g => $"{g.Key}:{g.Count()}");

            TestContext.Out.WriteLine("wrote " + outPath);
            TestContext.Out.WriteLine($"consensus rows: {lines.Length - 1}   n_files: {string.Join("  ", dist)}");
            TestContext.Out.WriteLine($"present in all 3 (compare to OpenMS): {nFiles.Count(n => n == 3)}");
        }
    }
}
