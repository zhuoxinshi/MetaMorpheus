using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.ModificationAnalysis;
using Nett;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    /// <summary>
    /// A two-round, "adaptive" G-PTM-D workflow.
    ///
    /// Round 1 runs a default G-PTM-D task (which annotates the database with candidate
    /// modifications) followed by a search task. The search task's
    /// <see cref="ModificationAnalysisEngine"/> reports how many times each modification was
    /// confidently localized (<see cref="ModificationAnalysisResults.CountOfModsSeenAndLocalized"/>).
    ///
    /// A count threshold is then derived from the distribution of all of those localized-mod
    /// counts and only the modifications whose localized count is above the threshold are kept.
    /// A new G-PTM-D toml is written containing just the surviving modifications.
    ///
    /// Round 2 re-runs G-PTM-D with the reduced modification list and searches the resulting
    /// database. This keeps the (expensive) second-pass search focused on the modifications that
    /// were actually observed in the data, instead of the full default G-PTM-D mod list.
    /// </summary>
    public class GptmdCountThresholdTask
    {
        /// <summary>
        /// How the count threshold is computed from the distribution of all localized-mod counts.
        /// </summary>
        public enum ThresholdMethod
        {
            /// <summary>Arithmetic mean of all localized-mod counts.</summary>
            Mean,

            /// <summary>Median (50th percentile) of all localized-mod counts.</summary>
            Median,

            /// <summary>A caller-specified percentile of all localized-mod counts (see <see cref="ThresholdPercentile"/>).</summary>
            Percentile
        }

        public GptmdCountThresholdTask()
        {
            CommonParameters = new CommonParameters();
            GptmdParameters = new GptmdParameters();
            SearchParameters = new SearchParameters();
            ThresholdSelectionMethod = ThresholdMethod.Median;
            ThresholdPercentile = 50;
        }

        /// <summary>Common parameters shared by both the G-PTM-D and search sub-tasks.</summary>
        public CommonParameters CommonParameters { get; set; }

        /// <summary>The default G-PTM-D parameters used for the first (unfiltered) round.</summary>
        public GptmdParameters GptmdParameters { get; set; }

        /// <summary>The search parameters used for both search rounds.</summary>
        public SearchParameters SearchParameters { get; set; }

        /// <summary>The strategy used to turn the distribution of localized-mod counts into a threshold.</summary>
        public ThresholdMethod ThresholdSelectionMethod { get; set; }

        /// <summary>
        /// Percentile in [0, 100] used when <see cref="ThresholdSelectionMethod"/> is
        /// <see cref="ThresholdMethod.Percentile"/>. Lower values keep more modifications.
        /// </summary>
        public double ThresholdPercentile { get; set; }

        /// <summary>
        /// Runs the full two-round workflow.
        /// </summary>
        /// <param name="outputFolder">Root output folder. Each sub-task writes into its own child folder.</param>
        /// <param name="dbFilenameList">The starting protein/nucleic-acid databases.</param>
        /// <param name="currentRawFileList">The spectra files to search.</param>
        public GptmdCountThresholdResults Run(string outputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList)
        {
            if (!Directory.Exists(outputFolder))
                Directory.CreateDirectory(outputFolder);

            // ---------- Round 1: default G-PTM-D ----------
            var gptmd1Folder = Path.Combine(outputFolder, "Task1-GPTMD");
            Directory.CreateDirectory(gptmd1Folder);
            var gptmd1 = new GptmdTask
            {
                CommonParameters = CommonParameters,
                GptmdParameters = GptmdParameters
            };
            var gptmd1Results = gptmd1.RunTask(gptmd1Folder, dbFilenameList, currentRawFileList, "Task1-GPTMD");
            var round1Databases = gptmd1Results.NewDatabases;
            string round1GptmdDatabase = round1Databases
                .Where(d => !d.IsSpectralLibrary)
                .Select(d => d.FilePath)
                .FirstOrDefault();

            // ---------- Round 1: search the G-PTM-D database, capturing localized-mod counts ----------
            var search1Folder = Path.Combine(outputFolder, "Task2-Search");
            Directory.CreateDirectory(search1Folder);
            var search1 = new SearchTask
            {
                CommonParameters = CommonParameters,
                SearchParameters = SearchParameters
            };

            var localizedModCounts = new Dictionary<string, int>();
            EventHandler<SingleEngineFinishedEventArgs> captureHandler = (sender, args) =>
            {
                // The modification-analysis engine reports the localized-mod counts we threshold on.
                if (args.MyResults is ModificationAnalysisResults modResults && modResults.CountOfModsSeenAndLocalized != null)
                {
                    foreach (var kvp in modResults.CountOfModsSeenAndLocalized)
                        localizedModCounts[kvp.Key] = kvp.Value;
                }
            };

            MyTaskResults round1SearchResults;
            try
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler += captureHandler;
                round1SearchResults = search1.RunTask(search1Folder, round1Databases, currentRawFileList, "Task2-Search");
            }
            finally
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= captureHandler;
            }

            // ---------- Threshold + mod selection ----------
            double countThreshold = DetermineCountThreshold(localizedModCounts.Values.ToList(), ThresholdSelectionMethod, ThresholdPercentile);
            List<string> selectedModIds = SelectModsAboveThreshold(localizedModCounts, countThreshold);
            List<(string, string)> filteredModList = FilterGptmdModList(GptmdParameters.ListOfModsGptmd, selectedModIds);

            // ---------- Build + persist the filtered G-PTM-D task as a new toml ----------
            var gptmd2 = new GptmdTask
            {
                CommonParameters = CommonParameters,
                GptmdParameters = new GptmdParameters
                {
                    ListOfModsGptmd = filteredModList,
                    GptmdFilters = GptmdParameters.GptmdFilters,
                    WriteDecoys = GptmdParameters.WriteDecoys
                }
            };
            string filteredTomlPath = Path.Combine(outputFolder, "Task3-GPTMD-CountFilteredconfig.toml");
            Toml.WriteFile(gptmd2, filteredTomlPath, MetaMorpheusTask.tomlConfig);

            // ---------- Round 2: filtered G-PTM-D + search ----------
            var gptmd2Folder = Path.Combine(outputFolder, "Task3-GPTMD-CountFiltered");
            Directory.CreateDirectory(gptmd2Folder);
            var gptmd2Results = gptmd2.RunTask(gptmd2Folder, dbFilenameList, currentRawFileList, "Task3-GPTMD-CountFiltered");
            var round2Databases = gptmd2Results.NewDatabases;
            string round2GptmdDatabase = round2Databases
                .Where(d => !d.IsSpectralLibrary)
                .Select(d => d.FilePath)
                .FirstOrDefault();

            var search2Folder = Path.Combine(outputFolder, "Task4-Search");
            Directory.CreateDirectory(search2Folder);
            var search2 = new SearchTask
            {
                CommonParameters = CommonParameters,
                SearchParameters = SearchParameters
            };
            var round2SearchResults = search2.RunTask(search2Folder, round2Databases, currentRawFileList, "Task4-Search");

            return new GptmdCountThresholdResults
            {
                Round1GptmdDatabase = round1GptmdDatabase,
                Round1SearchResults = round1SearchResults,
                LocalizedModCounts = localizedModCounts,
                CountThreshold = countThreshold,
                SelectedModIdsWithMotif = selectedModIds,
                FilteredGptmdModList = filteredModList,
                FilteredGptmdTomlPath = filteredTomlPath,
                Round2GptmdDatabase = round2GptmdDatabase,
                Round2SearchResults = round2SearchResults
            };
        }

        /// <summary>
        /// Computes a count threshold from the distribution of all localized-mod counts.
        /// Returns 0 when there are no counts (so nothing survives selection).
        /// </summary>
        public static double DetermineCountThreshold(IReadOnlyCollection<int> counts, ThresholdMethod method, double percentile)
        {
            if (counts == null || counts.Count == 0)
                return 0;

            switch (method)
            {
                case ThresholdMethod.Mean:
                    return counts.Average();
                case ThresholdMethod.Median:
                    return GetPercentile(counts, 50);
                case ThresholdMethod.Percentile:
                    return GetPercentile(counts, percentile);
                default:
                    throw new ArgumentOutOfRangeException(nameof(method), method, "Unknown threshold method");
            }
        }

        /// <summary>
        /// Returns the ids (with motif, e.g. "Phosphorylation on S") of every modification whose
        /// localized count is strictly greater than <paramref name="threshold"/>.
        /// </summary>
        public static List<string> SelectModsAboveThreshold(IReadOnlyDictionary<string, int> localizedModCounts, double threshold)
        {
            if (localizedModCounts == null)
                return new List<string>();

            return localizedModCounts
                .Where(kv => kv.Value > threshold)
                .Select(kv => kv.Key)
                .ToList();
        }

        /// <summary>
        /// Restricts a G-PTM-D modification list to those whose IdWithMotif is in
        /// <paramref name="keptModIdsWithMotif"/>. The (ModificationType, IdWithMotif) tuple shape is
        /// preserved so the result can be assigned straight back to <see cref="GptmdParameters.ListOfModsGptmd"/>.
        /// </summary>
        public static List<(string, string)> FilterGptmdModList(IEnumerable<(string, string)> originalGptmdMods, ICollection<string> keptModIdsWithMotif)
        {
            var kept = new HashSet<string>(keptModIdsWithMotif ?? Array.Empty<string>());
            return (originalGptmdMods ?? Enumerable.Empty<(string, string)>())
                .Where(m => kept.Contains(m.Item2))
                .ToList();
        }

        /// <summary>
        /// Linear-interpolation percentile (matching Excel's PERCENTILE.INC / NIST method) of a set of counts.
        /// </summary>
        private static double GetPercentile(IReadOnlyCollection<int> counts, double percentile)
        {
            var sorted = counts.OrderBy(c => c).ToList();
            if (sorted.Count == 1)
                return sorted[0];

            double clamped = Math.Clamp(percentile, 0, 100);
            double rank = clamped / 100.0 * (sorted.Count - 1);
            int lowerIndex = (int)Math.Floor(rank);
            int upperIndex = (int)Math.Ceiling(rank);
            double fraction = rank - lowerIndex;
            return sorted[lowerIndex] + fraction * (sorted[upperIndex] - sorted[lowerIndex]);
        }
    }

    /// <summary>
    /// The artifacts produced by a <see cref="GptmdCountThresholdTask"/> run.
    /// </summary>
    public class GptmdCountThresholdResults
    {
        /// <summary>Path to the G-PTM-D database written in round 1 (the unfiltered pass).</summary>
        public string Round1GptmdDatabase { get; init; }

        /// <summary>Task results from the round-1 search.</summary>
        public MyTaskResults Round1SearchResults { get; init; }

        /// <summary>Localized-mod counts (IdWithMotif -> count) observed in the round-1 search.</summary>
        public IReadOnlyDictionary<string, int> LocalizedModCounts { get; init; }

        /// <summary>The count threshold derived from <see cref="LocalizedModCounts"/>.</summary>
        public double CountThreshold { get; init; }

        /// <summary>Ids (with motif) of the modifications kept because their count exceeded the threshold.</summary>
        public IReadOnlyList<string> SelectedModIdsWithMotif { get; init; }

        /// <summary>The reduced G-PTM-D modification list (subset of the original) used in round 2.</summary>
        public IReadOnlyList<(string, string)> FilteredGptmdModList { get; init; }

        /// <summary>Path to the newly generated, count-filtered G-PTM-D toml.</summary>
        public string FilteredGptmdTomlPath { get; init; }

        /// <summary>Path to the G-PTM-D database written in round 2 (the filtered pass).</summary>
        public string Round2GptmdDatabase { get; init; }

        /// <summary>Task results from the round-2 search.</summary>
        public MyTaskResults Round2SearchResults { get; init; }
    }
}
