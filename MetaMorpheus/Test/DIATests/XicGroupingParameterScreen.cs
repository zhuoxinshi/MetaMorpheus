using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.DIATests
{
    /// <summary>
    /// Parameter screen over <see cref="XicGroupingEngine"/> settings for the ISD pseudo-MS2 pipeline.
    ///
    /// Everything except the grouping engine is identical to
    /// <see cref="IsdMsAlignExportTest.SearchAllPseudoScansMgf"/>: the same HCD-YBC ISD60-80-100 raw files,
    /// the same NeutralMassXicConstructor settings for the MS1 and MS2 sides, the same search / GPTMD tomls
    /// with the same overrides, and the same yeast + human XML databases. Only ApexRTTolerance,
    /// CorrelationThreshold and FragmentRankThreshold vary. Each set writes to its own folder
    /// neutralMass1, neutralMass2, ... under MM\ISD_search.
    ///
    /// Which sets run is controlled by the XIC_SCREEN_SET environment variable:
    ///   unset / "all"  -> every set, in order
    ///   "3"            -> only set 3
    ///   "1,2,4"        -> those sets, in that order
    ///
    /// A set whose folder already contains allResults.txt is treated as finished and skipped, so the
    /// driver script can be re-run without redoing completed work.
    /// </summary>
    [TestFixture]
    public static class XicGroupingParameterScreen
    {
        // Sets 1-4 are the requested screen. Sets 5-8 are the single-factor contrasts the first four
        // confound: 5 vs 2, 7 vs 1 and 8 vs 3 isolate the fragment rank cap at fixed RT/correlation,
        // and 2 / 6 / 4 form a correlation ladder (0.75 -> 0.5 -> 0) at apexRT 0.1 with cap 1000.
        private static readonly (int Index, float ApexRt, double Overlap, double Correlation, int FragmentRank)[] AllSets =
        {
            (1, 0.2f, 0.0, 0.75, 200),
            (2, 0.1f, 0.0, 0.75, 1000),
            (3, 0.3f, 0.0, 0.50, 200),
            (4, 0.1f, 0.0, 0.00, 1000),
            (5, 0.1f, 0.0, 0.75, 200),
            (6, 0.1f, 0.0, 0.50, 1000),
            (7, 0.2f, 0.0, 0.75, 1000),
            (8, 0.3f, 0.0, 0.50, 1000),
            // Added 2026-08-26 after sets 1-4 landed. Sets 4, 6, 8 and this one form a 2x2 factorial,
            // correlation {0, 0.5} x apexRT {0.1, 0.3}, all at cap 1000 - the only clean way to separate
            // the two knobs, which every earlier set moves together. This is the missing corner:
            // the loosest condition tested, combining both directions that improved PSM counts.
            (9, 0.3f, 0.0, 0.00, 1000),
        };

        private const string DataDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
        private const string OutDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\MM\ISD_search";

        [Test]
        public static void RunNeutralMassXicGroupingScreen()
        {
            if (!Directory.Exists(OutDir)) Directory.CreateDirectory(OutDir);

            var fileList = Directory.GetFiles(DataDir, "*ISD60-80-100*")
                .Where(f => f.Contains(".raw") & f.Contains("HCD-YBC")).ToList();

            Log($"{fileList.Count} spectra files matched in {DataDir}");
            foreach (var f in fileList) Log($"   {Path.GetFileName(f)}");
            Assert.That(fileList.Count, Is.GreaterThan(0), "no raw files matched the HCD-YBC ISD60-80-100 filter");

            foreach (var set in ResolveRequestedSets())
            {
                RunOneSet(set, fileList);
            }
        }

        private static void RunOneSet(
            (int Index, float ApexRt, double Overlap, double Correlation, int FragmentRank) set,
            List<string> fileList)
        {
            var outFolder = Path.Combine(OutDir, "neutralMass" + set.Index);
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

            string label = $"set {set.Index} (apexRT {set.ApexRt}, overlap {set.Overlap}, " +
                           $"corr {set.Correlation}, fragRank {set.FragmentRank})";

            if (File.Exists(Path.Combine(outFolder, "allResults.txt")))
            {
                Log($"SKIP {label} - allResults.txt already present in {outFolder}");
                return;
            }

            var sw = Stopwatch.StartNew();
            Log($"START {label} -> {outFolder}");

            // ---- pseudo-MS2 construction. Only the XicGroupingEngine differs between sets. ----
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                    new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3,
                    xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                    new ClassicDeconvolutionParameters(1, 20, 4, 3),
                    xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                new XicGroupingEngine(set.ApexRt, set.Overlap, set.Correlation,
                    maxThreadsForGrouping: 10, fragmentRankThreshold: set.FragmentRank),
                PseudoMs2ConstructionType.Mass, combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            diaParams.WriteSettingsToFile(Path.Combine(outFolder, "DIAparameters.txt"));

            foreach (var file in fileList)
            {
                var fileSw = Stopwatch.StartNew();
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();
                var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();

                var stem = Path.GetFileNameWithoutExtension(dataFile.FilePath);
                ISDEngine.WriteMsAlignFile(Path.Combine(outFolder, stem + "_pseudo_ms2.msalign"), pseudoScans);
                IsdMsAlignExporter.WriteMgf(pseudoScans, Path.Combine(outFolder, stem + "_pseudo_ms2.mgf"));

                var fragCounts = pseudoScans.Select(s => s.ExperimentalFragments.Length).OrderBy(n => n).ToArray();
                string fragSummary = fragCounts.Length > 0
                    ? $"frags/precursor mean={fragCounts.Average():F1} median={fragCounts[fragCounts.Length / 2]} " +
                      $"max={fragCounts[^1]} total={fragCounts.Sum():N0}"
                    : "no pseudo scans";
                Log($"   [set {set.Index}] {stem}: {pseudoScans.Length} pseudo-MS2; {fragSummary}; " +
                    $"{fileSw.Elapsed:hh\\:mm\\:ss}");
            }

            // ---- search + GPTMD + post-GPTMD search, unchanged from SearchAllPseudoScansMgf ----
            string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
            var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
            searchTask.CommonParameters.DoPrecursorDeconvolution = false;   // already deconvoluted via FromFile
            searchTask.CommonParameters.UseProvidedPrecursorInfo = true;    // FromFile assembly provides the precursor info
            searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
            searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;

            var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = searchTask.CommonParameters;

            string yeast_xml = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.xml";
            string human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
            var databases = new List<DbForTask> { new DbForTask(yeast_xml, false), new DbForTask(human_xml, false) };

            var searchFileList = Directory.GetFiles(outFolder, "*_pseudo_ms2.mgf").ToList();
            Log($"   [set {set.Index}] searching {searchFileList.Count} mgf files");

            var taskList = new List<(string, MetaMorpheusTask)>
            {
                ("search", searchTask),
                ("GPTMD", gptmdTask),
                ("postGPTMD-search", searchTask)
            };
            var engine = new EverythingRunnerEngine(taskList, searchFileList, databases, outFolder);
            engine.Run();

            Log($"DONE {label} in {sw.Elapsed:hh\\:mm\\:ss}");
            AppendSummary(set, outFolder, sw.Elapsed);
        }

        /// <summary>
        /// Pull the final (post-GPTMD) target counts out of allResults.txt and append one row per set to
        /// screen_summary.tsv, so the whole screen is readable in one place without opening each folder.
        /// </summary>
        private static void AppendSummary(
            (int Index, float ApexRt, double Overlap, double Correlation, int FragmentRank) set,
            string outFolder, TimeSpan elapsed)
        {
            string allResults = Path.Combine(outFolder, "allResults.txt");
            string psms = "NA", proteoforms = "NA", proteinGroups = "NA";

            if (File.Exists(allResults))
            {
                // later lines win: the last block in the file is the post-GPTMD search
                foreach (var line in File.ReadAllLines(allResults))
                {
                    if (line.StartsWith("All target PSMs with q-value <= 0.01:"))
                        psms = line.Split(':').Last().Trim();
                    else if (line.StartsWith("All target proteoforms with q-value <= 0.01:"))
                        proteoforms = line.Split(':').Last().Trim();
                    else if (line.StartsWith("All target protein groups with q-value <= 0.01"))
                        proteinGroups = line.Split(':').Last().Trim();
                }
            }

            string summaryPath = Path.Combine(OutDir, "screen_summary.tsv");
            try
            {
                if (!File.Exists(summaryPath))
                {
                    File.AppendAllText(summaryPath,
                        "Folder\tApexRTTolerance\tOverlapThreshold\tCorrelationThreshold\tFragmentRankThreshold\t" +
                        "TargetPSMs\tTargetProteoforms\tTargetProteinGroups\tElapsed\n");
                }
                File.AppendAllText(summaryPath,
                    $"neutralMass{set.Index}\t{set.ApexRt}\t{set.Overlap}\t{set.Correlation}\t{set.FragmentRank}\t" +
                    $"{psms}\t{proteoforms}\t{proteinGroups}\t{elapsed:hh\\:mm\\:ss}\n");
            }
            catch (Exception e)
            {
                Log($"   could not write screen_summary.tsv: {e.Message}");
            }

            Log($"   [set {set.Index}] PSMs={psms} proteoforms={proteoforms} proteinGroups={proteinGroups}");
        }

        private static IEnumerable<(int Index, float ApexRt, double Overlap, double Correlation, int FragmentRank)> ResolveRequestedSets()
        {
            var raw = Environment.GetEnvironmentVariable("XIC_SCREEN_SET");
            if (string.IsNullOrWhiteSpace(raw) || raw.Trim().Equals("all", StringComparison.OrdinalIgnoreCase))
            {
                Log("XIC_SCREEN_SET not set - running every parameter set in order");
                return AllSets;
            }

            var selected = raw.Split(',', StringSplitOptions.RemoveEmptyEntries)
                .Select(s => int.Parse(s.Trim()))
                .Select(i => AllSets.First(s => s.Index == i))
                .ToArray();
            Log($"XIC_SCREEN_SET={raw} - running {selected.Length} set(s)");
            return selected;
        }

        private static void Log(string message)
        {
            string line = $"[{DateTime.Now:yyyy-MM-dd HH:mm:ss}] {message}";
            TestContext.Progress.WriteLine(line);
            Console.WriteLine(line);
            try
            {
                if (!Directory.Exists(OutDir)) Directory.CreateDirectory(OutDir);
                File.AppendAllText(Path.Combine(OutDir, "screen_progress.log"), line + Environment.NewLine);
            }
            catch
            {
                // logging must never fail the run
            }
        }
    }
}
