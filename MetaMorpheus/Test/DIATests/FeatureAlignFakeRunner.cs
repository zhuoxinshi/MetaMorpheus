using EngineLayer.DIA;
using EngineLayer.DIA.FeatureAlign;
using NUnit.Framework;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Aligns the 3 synthetic .feature files in the FeatureAlign_test folder and writes the consensus tsv INTO
    /// that same folder, so the fake-file alignment result is easy to find. [Explicit].
    /// </summary>
    [TestFixture]
    public class FeatureAlignFakeRunner
    {
        [Test, Explicit("aligns the fake files in FeatureAlign_test and writes the consensus there")]
        public void AlignFakeFilesInTestFolder()
        {
            var dir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\TopFD\FeatureAlign_test";
            var files = Directory.GetFiles(dir, "*_ms1.feature").OrderBy(f => f).ToList();
            Assert.That(files.Count, Is.EqualTo(3), "expected run1/run2/run3 _ms1.feature in " + dir);

            string outPath = Path.Combine(dir, "FeatureAlign_fake_20ppm_60s_notch2.tsv");
            FeatureAligner.Align(files, outPath,
                massTolerance: 20, massUnit: MassToleranceUnit.Ppm, rtTolerance: 60,
                preferredQuantColumn: "Intensity",
                positiveNotch: 2, negativeNotch: 2,
                grouping: GroupingMethod.ReferenceAnchor);

            var lines = File.ReadAllLines(outPath);
            int nf = System.Array.IndexOf(lines[0].Split('\t'), "n_files");
            var dist = lines.Skip(1).Where(l => l.Length > 0).Select(l => int.Parse(l.Split('\t')[nf]))
                            .GroupBy(x => x).OrderBy(g => g.Key).Select(g => $"{g.Key}:{g.Count()}");
            TestContext.Out.WriteLine("wrote " + outPath);
            TestContext.Out.WriteLine($"rows: {lines.Length - 1}   n_files: {string.Join("  ", dist)}");
        }
    }
}
