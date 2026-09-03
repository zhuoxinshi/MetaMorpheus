using EngineLayer.DIA;
using EngineLayer.DIA.FeatureAlign;
using NUnit.Framework;
using System.IO;
using System.Linq;

namespace Test.DIATests.Manuscript
{
    /// <summary>
    /// Verifies that <see cref="GroupingMethod.OpenMSConsensus"/> reproduces OpenMS TopDownConsensusFeatureGroup
    /// (the ConsensusFeatureGroupDetector exe) on the same 3 FLASHDeconv mix files. OpenMS at 100 ppm / 180 s
    /// produced 151 consensus groups (all present in all 3 runs); the faithful port should match. [Explicit]
    /// (reads real data on disk).
    /// </summary>
    [TestFixture]
    public class OpenMSConsensusPortTest
    {
        [Test, Explicit("checks OpenMSConsensus grouping == OpenMS's 151 groups on the FLASHDeconv mix files")]
        public void OpenMSConsensus_Matches_OpenMS()
        {
            var fd = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\FeatureFiles";
            var files = new[]
            {
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep1.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep2.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep3.tsv",
            }.Select(f => Path.Combine(fd, f)).ToList();
            Assert.That(files.All(File.Exists), "FLASHDeconv FeatureFiles must be present");

            string outPath = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\_compare\ours_openms_port_consensus.tsv";
            FeatureAligner.Align(files, outPath,
                massTolerance: 100, massUnit: MassToleranceUnit.Ppm, rtTolerance: 180,
                preferredQuantColumn: "FeatureQuantity",
                grouping: GroupingMethod.OpenMSConsensus,
                whenDuplicate: DuplicateResolution.MaxAbundance);

            var lines = File.ReadAllLines(outPath);
            int nf = System.Array.IndexOf(lines[0].Split('\t'), "n_files");
            var nFiles = lines.Skip(1).Where(l => l.Length > 0).Select(l => int.Parse(l.Split('\t')[nf])).ToList();
            int total = nFiles.Count;
            int all3 = nFiles.Count(n => n == 3);

            TestContext.Out.WriteLine($"OpenMSConsensus rows: {total}   all-3: {all3}   (OpenMS produced 151)");
            Assert.That(nFiles.All(n => n == 3), Is.True, "OpenMSConsensus must emit only all-run groups");
            Assert.That(total, Is.EqualTo(151).Within(3),
                "should match OpenMS's 151 consensus groups (small slack for tie-breaks)");
        }
    }
}
