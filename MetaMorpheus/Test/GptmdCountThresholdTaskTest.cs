using EngineLayer;
using EngineLayer.DatabaseLoading;
using Nett;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class GptmdCountThresholdTaskTest
    {
        /// <summary>
        /// Deterministic coverage of the core algorithm: turning a distribution of localized-mod
        /// counts into a threshold, selecting the mods above it, and restricting the G-PTM-D mod list.
        /// </summary>
        [Test]
        public static void TestCountThresholdSelectionLogic()
        {
            var localizedModCounts = new Dictionary<string, int>
            {
                { "Phosphorylation on S", 100 },
                { "Acetylation on K", 50 },
                { "Oxidation on M", 10 },
                { "Methylation on R", 2 },
                { "Carbamyl on K", 1 },
            };

            // median of {1, 2, 10, 50, 100} is 10
            double median = GptmdCountThresholdTask.DetermineCountThreshold(
                localizedModCounts.Values.ToList(), GptmdCountThresholdTask.ThresholdMethod.Median, 50);
            Assert.That(median, Is.EqualTo(10));

            // mean of {1, 2, 10, 50, 100} is 32.6
            double mean = GptmdCountThresholdTask.DetermineCountThreshold(
                localizedModCounts.Values.ToList(), GptmdCountThresholdTask.ThresholdMethod.Mean, 0);
            Assert.That(mean, Is.EqualTo((100 + 50 + 10 + 2 + 1) / 5.0));

            // strictly above the median (10) keeps Phospho (100) and Acetyl (50)
            var selected = GptmdCountThresholdTask.SelectModsAboveThreshold(localizedModCounts, median);
            Assert.That(selected, Is.EquivalentTo(new[] { "Phosphorylation on S", "Acetylation on K" }));

            // the (ModificationType, IdWithMotif) shape is preserved and only the kept mods survive
            var originalGptmdMods = new List<(string, string)>
            {
                ("Common Biological", "Phosphorylation on S"),
                ("Common Biological", "Acetylation on K"),
                ("Common Variable", "Oxidation on M"),
                ("Metal", "Sodium on D"),
            };
            var filtered = GptmdCountThresholdTask.FilterGptmdModList(originalGptmdMods, selected);
            Assert.That(filtered, Is.EquivalentTo(new[]
            {
                ("Common Biological", "Phosphorylation on S"),
                ("Common Biological", "Acetylation on K"),
            }));

            // empty distribution -> threshold 0, nothing selected
            Assert.That(GptmdCountThresholdTask.DetermineCountThreshold(
                new List<int>(), GptmdCountThresholdTask.ThresholdMethod.Median, 50), Is.EqualTo(0));
            Assert.That(GptmdCountThresholdTask.SelectModsAboveThreshold(new Dictionary<string, int>(), 0), Is.Empty);
        }

        /// <summary>
        /// End-to-end run of the two-round workflow: default G-PTM-D + search, threshold the localized-mod
        /// counts, write a reduced G-PTM-D toml, then re-run G-PTM-D + search. The spectra file path is a
        /// variable so this can be pointed at other data.
        /// </summary>
        [Test]
        public static void TestGptmdCountThresholdWorkflow()
        {
            // ---- input file paths (left as variables so this can be repointed at other data) ----
            string spectraFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\ModificationAnalysis\modificationAnalysis.mzML");
            string databaseFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\ModificationAnalysis\modificationAnalysis.fasta");

            // reuse the known-good parameters that already produce localized mods on this test data
            var gptmdLoaded = Toml.ReadFile<GptmdTask>(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\ModificationAnalysis\GPTMDTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            var searchLoaded = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\ModificationAnalysis\SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);

            var task = new GptmdCountThresholdTask
            {
                CommonParameters = gptmdLoaded.CommonParameters,
                GptmdParameters = gptmdLoaded.GptmdParameters,
                SearchParameters = searchLoaded.SearchParameters,
                ThresholdSelectionMethod = GptmdCountThresholdTask.ThresholdMethod.Median,
            };
            int originalModCount = task.GptmdParameters.ListOfModsGptmd.Count;

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestGptmdCountThreshold");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(outputFolder);

            var results = task.Run(outputFolder,
                new List<DbForTask> { new DbForTask(databaseFilePath, false) },
                new List<string> { spectraFilePath });

            // round 1 produced a G-PTM-D database and observed some localized mods
            Assert.That(results.Round1GptmdDatabase, Is.Not.Null);
            Assert.That(File.Exists(results.Round1GptmdDatabase));
            Assert.That(results.LocalizedModCounts.Count, Is.GreaterThan(0));

            // threshold is a real value within the observed count range
            var counts = results.LocalizedModCounts.Values.ToList();
            Assert.That(results.CountThreshold, Is.GreaterThanOrEqualTo(counts.Min()));
            Assert.That(results.CountThreshold, Is.LessThanOrEqualTo(counts.Max()));

            // every selected mod is strictly above the threshold (selection is consistent with the threshold)
            foreach (var modId in results.SelectedModIdsWithMotif)
                Assert.That(results.LocalizedModCounts[modId], Is.GreaterThan(results.CountThreshold));

            // the reduced G-PTM-D mod list is a genuine subset of the original default list
            Assert.That(results.FilteredGptmdModList.Count, Is.LessThanOrEqualTo(originalModCount));
            foreach (var mod in results.FilteredGptmdModList)
                Assert.That(task.GptmdParameters.ListOfModsGptmd, Does.Contain(mod));

            // a new, count-filtered G-PTM-D toml was written and round-trips with the same mod list
            Assert.That(File.Exists(results.FilteredGptmdTomlPath));
            var reloaded = Toml.ReadFile<GptmdTask>(results.FilteredGptmdTomlPath, MetaMorpheusTask.tomlConfig);
            Assert.That(reloaded.GptmdParameters.ListOfModsGptmd.Count, Is.EqualTo(results.FilteredGptmdModList.Count));

            // round 2 re-ran G-PTM-D on the filtered list and searched the resulting database
            Assert.That(results.Round2GptmdDatabase, Is.Not.Null);
            Assert.That(File.Exists(results.Round2GptmdDatabase));
            Assert.That(results.Round2SearchResults, Is.Not.Null);

            Directory.Delete(outputFolder, true);
        }
    }
}
