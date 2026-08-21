using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using MassSpectrometry;
using Microsoft.Win32;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;
using UsefulProteomicsDatabases;
using static iText.StyledXmlParser.Jsoup.Select.Evaluator;

namespace Test.DIATests
{
    [TestFixture]
    public static class IsdMsAlignExportTest
    {
        /// <summary>
        /// End-to-end: take a real ISD mzML, build pseudo-MS2 scans using CONSENSUS FEATURE TRACING
        /// (<see cref="ConsensusMassXicConstructor"/> = mzLib MassTraceBuilder→TraceCorrector→MassFeatureBuilder)
        /// as the XIC source, combined with apex-RT/Pearson precursor-fragment grouping and pseudo-MS2 generation,
        /// then export them to ms1/ms2.msalign and confirm the ms2.msalign round-trips through the MsAlign reader
        /// (i.e. it is a valid file TopPIC could consume).
        /// </summary>
        [Test]
        public static void ExportIsdPseudoScansAndSearchMgf()
        {
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HC_YC";
            var fileList = Directory.GetFiles(dataDir, "*ISD60-80-100*").Where(f => f.Contains(".raw") & f.Contains("HC-YC")).ToList();

            string yeast_xml = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.xml";
            var human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
            var yeast_fasta = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.fasta";
            string pierceDb = @"E:\ISD Project\071626_ISD\Quant\Pierce\Pierce_6pro.fasta";
            string human_fasta = @"E:\Databases\Human_9606.fasta";
            var databases = new List<DbForTask> { new DbForTask(yeast_xml, false), new DbForTask(human_xml, false) };

            var outDir = @"E:\ISD Project\071626_ISD\Quant\HC_YC\ISD_search";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);

            bool group = true;
            bool write = true;
            bool search = true;
            foreach (var file in fileList)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();

                // Consensus feature tracing as the XIC front-end for both the intact (precursor) and fragment channels.
                var diaParams = new DIAparameters(
                    DIAanalysisType.ISD,
                    //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), traceToleranceDa: 1.5, minMass: 4000, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    //new XicGroupingEngine(0.2f, 0.3, 0.5, maxThreadsForGrouping: 10),
                    new XicGroupingEngine(0.2f, 0, 0.5, maxThreadsForGrouping: 10, fragmentRankThreshold: 200), PseudoMs2ConstructionType.Mass, combineFragments: true);
                var commonParams = new CommonParameters { DIAparameters = diaParams };

                Ms2ScanWithSpecificMass[] pseudoScans = null;
                if (group)
                {
                    pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();
                }
                
                var outFolder = Path.Combine(outDir, Path.GetFileNameWithoutExtension(dataFile.FilePath), "neutralMass1");
                if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);
                var msAlignOutPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.msalign");
                var mgfOutPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.mgf");

                if (write)
                {
                    diaParams.WriteSettingsToFile(Path.Combine(outFolder, "DIAparameters.txt"));
                    ISDEngine.WriteMsAlignFile(msAlignOutPath, pseudoScans);
                    IsdMsAlignExporter.WriteMgf(pseudoScans, mgfOutPath);
                }
                
                if (search)
                {
                    string calibrationToml = @"E:\ISD Project\071626_ISD\Quant\HC_YC\HC_YC_DDA_xml-cali-xml\Task Settings\Task2-CalibrateTaskconfig.toml";
                    var calibrationTask = Toml.ReadFile<CalibrationTask>(calibrationToml, MetaMorpheusTask.tomlConfig);
                    string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
                    var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
                    searchTask.CommonParameters.DoPrecursorDeconvolution = false; // already deconvoluted via FromFile
                    searchTask.CommonParameters.UseProvidedPrecursorInfo = true; // FromFile assembly provides the precursor info
                    searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                    searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;
                    //searchTask.SearchParameters.MinAllowedInternalFragmentLength = 0;

                    var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                    var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                    var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);
                    gptmdTask.CommonParameters = searchTask.CommonParameters;

                    var taskList = new List<(string, MetaMorpheusTask)> { ("calibration", searchTask), ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) }; //("GPTMD", gptmdTask)
                    var engine = new EverythingRunnerEngine(taskList, new List<string> { mgfOutPath }, databases, outFolder);
                    engine.Run();
                }
            }
        }

        [Test]
        public static void SearchAllPseudoScansMgf()
        {
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HC_YC";
            //var fileList = Directory.GetFiles(dataDir, "*ISD60-80-100*").Where(f => f.Contains(".raw") & f.Contains("HC-YC")).ToList();
            var fileList = new List<string>
            {
                @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1.raw"
            };

            var outDir = @"E:\ISD Project\071626_ISD\Quant\HC_YC\ISD_search";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);
            var outFolder = Path.Combine(outDir, "neutralMass3");
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

            bool group = true;
            bool write = true;
            bool search = true;

            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), traceToleranceDa: 1.5, minMass: 4000, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new XicGroupingEngine(0.2f, 0.3, 0.5, maxThreadsForGrouping: 10),
                new XicGroupingEngine(0.2f, 0, 0.6, maxThreadsForGrouping: 10, fragmentRankThreshold: 200), PseudoMs2ConstructionType.Mass, combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            if (group)
            {
                foreach(var file in fileList)
                {
                    var dataFile = MsDataFileReader.GetDataFile(file);
                    dataFile.LoadAllStaticData();
                    var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();
                    var msAlignOutPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.msalign");
                    var mgfOutPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.mgf");
                    if (write)
                    {
                        diaParams.WriteSettingsToFile(Path.Combine(outFolder, "DIAparameters.txt"));
                        ISDEngine.WriteMsAlignFile(msAlignOutPath, pseudoScans);
                        IsdMsAlignExporter.WriteMgf(pseudoScans, mgfOutPath);
                    }
                }
            }

            if (search)
            {
                string calibrationToml = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\HC_YC_DDA_xml-cali-xml\Task Settings\Task2-CalibrateTaskconfig.toml";
                var calibrationTask = Toml.ReadFile<CalibrationTask>(calibrationToml, MetaMorpheusTask.tomlConfig);
                string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
                var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DoPrecursorDeconvolution = false; // already deconvoluted via FromFile
                searchTask.CommonParameters.UseProvidedPrecursorInfo = true; // FromFile assembly provides the precursor info
                searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;
                //searchTask.SearchParameters.MinAllowedInternalFragmentLength = 0;

                var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters;

                string yeast_xml = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.xml";
                var human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
                var yeast_fasta = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.fasta";
                string pierceDb = @"E:\ISD Project\071626_ISD\Quant\Pierce\Pierce_6pro.fasta";
                string human_fasta = @"E:\Databases\Human_9606.fasta";
                var databases = new List<DbForTask> { new DbForTask(yeast_xml, false), new DbForTask(human_xml, false) };

                var searchFileList = Directory.GetFiles(outFolder, "*_pseudo_ms2.mgf").ToList();
                var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) }; //("GPTMD", gptmdTask)
                var engine = new EverythingRunnerEngine(taskList, searchFileList, databases, outFolder);
                engine.Run();
            }
        }

        /// <summary>
        /// Same pipeline as <see cref="SearchAllPseudoScansMgf"/> - ISD pseudo-MS2 construction, msAlign/mgf
        /// export, then MetaMorpheus search + GPTMD + post-GPTMD search - but the two XIC constructors read
        /// EXTERNAL deconvolution results (TopFD msAlign) instead of running MetaMorpheus's classic
        /// deconvolution. Everything downstream of the neutral masses is unchanged.
        ///
        /// TopFD is run per ISD channel, so each raw file has one msAlign per channel. The MS1-channel msAlign
        /// (*_ms1_ms1.msalign) drives the precursor XICs; the source-voltage channels (*_isd60/80/100_ms1.msalign)
        /// drive the fragment XICs. The join back to the raw file is by scan number, which the channel splits
        /// preserve (channel k of cycle n is raw scan 4n+k).
        /// </summary>
        [Test]
        public static void SearchAllPseudoScansMgf_ExternalDecon()
        {
            // Each entry is a spectra file paired with the TopFD output base name for that run, i.e. the prefix
            // its "<base>_feature.xml" / "<base>_ms2.msalign" companions share.
            //
            // REQUIRES TopFD 1.8 OR LATER. The MS1 side reads the per-charge <envelope_list> blocks out of
            // feature.xml, and those were only added in 1.8 - a 1.7.x feature.xml stops at <envelope_num> and
            // carries no per-scan points at all (197 KB vs 105 MB for the same kind of run), so this path yields
            // zero precursor XICs on it. The YB runs deconvoluted with 1.7.4 therefore cannot be used here; they
            // still work with the msAlign-tracing front end (see FeatureGroupedXicConstructor).
            //
            // The spectra file is the labelCorrected mzML TopFD actually saw, so scan numbers are guaranteed to
            // line up with what it wrote.
            var topFdRoot = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            var fileList = new List<(string SpectraFile, string TopFdBase)>
            {
                (Path.Combine(topFdRoot, "07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected.mzML"),
                 Path.Combine(topFdRoot, "07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected")),
                (Path.Combine(topFdRoot, "07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep2_labelCorrected.mzML"),
                 Path.Combine(topFdRoot, "07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep2_labelCorrected")),
            };

            var outDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\ISD_search_externalDecon";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);
            var outFolder = Path.Combine(outDir, "topfdDecon1");
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

            bool group = true;
            bool write = false; // counting fragments per precursor only; flip on to write the msalign/mgf
            // the goal here is the pseudo ms2.msalign itself; flip to true to run the MetaMorpheus search on it
            bool search = false;

            if (group)
            {
                foreach (var (spectraFile, topFdBase) in fileList)
                {
                    string stem = Path.GetFileNameWithoutExtension(spectraFile);
                    // TopFD writes both companions beside the spectra file it was run on. The feature xml covers
                    // MS1 scans ONLY - every envelope in it sits on an MS1 scan - which is why the fragment
                    // channels below stay on the msAlign path.
                    string precursorFeatureXml = topFdBase + "_feature.xml";
                    string combinedMs2MsAlign = topFdBase + "_ms2.msalign";

                    // the data file has to be loaded up front here: splitting the combined ms2.msalign into ISD
                    // channels needs the scan filters, so the constructors depend on the scans
                    var dataFile = MsDataFileReader.GetDataFile(spectraFile);
                    dataFile.LoadAllStaticData();
                    var allScans = dataFile.GetAllScansList().ToArray();

                    TestContext.Progress.WriteLine($"{stem}: {allScans.Length} scans; precursor feature xml " +
                        $"{Path.GetFileName(precursorFeatureXml)}; combined fragment msAlign " +
                        $"{Path.GetFileName(combinedMs2MsAlign)}");

                    // ---- the change from SearchAllPseudoScansMgf: masses come from TopFD, not ClassicDecon ----
                    // MS1: feature XICs read STRAIGHT out of TopFD's feature xml - its per-charge envelope lists
                    //      already give scan, RT and intensity for every point of every feature, so there is no
                    //      tracing at all. That sidesteps the two things that crippled tracing the ms1.msalign on
                    //      this data: masses recurring in only ~2 scans per feature, and TopFD flipping a
                    //      species' monoisotope between scans (M vs M+1.003 Da). On a matched slice this yields
                    //      ~3x the precursor XICs, with a median length of 6 points rather than the 3-point floor.
                    // was: new FeatureGroupedXicConstructor(...precursorMsAlign, precursorFeatureFile...)
                    // was: new NeutralMassXicConstructor(..., new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3, ...)
                    //
                    // The tolerance is unused here (nothing is matched by mass) but is still passed so the
                    // constructor stays interchangeable with the others and the value is recorded in the settings
                    // dump. minEcScore is now the quality lever: TopFD writes this file at ECScore >= 0.1, which
                    // is permissive, so raise it to drop low-confidence features.
                    var ms1Constructor = new FeatureXmlXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                        precursorFeatureXml, minMass: 4000, maxMass: 70000, minEcScore: 0,
                        xicSpline: new XicCubicSpline(0.05, 1, 1, true));

                    // MS2: still the msAlign path, because the feature xml holds MS1 scans only. This run was
                    //      deconvoluted in one go rather than per channel split, so all three source voltages sit
                    //      in ONE ms2.msalign and are told apart purely by scan number - FromCombinedMsAlign reads
                    //      it once and splits it using the sid= filters of the scans.
                    // was: new NeutralMassXicConstructor(..., new ClassicDeconvolutionParameters(1, 20, 4, 3), ...)
                    var ms2Constructor = ExternalDeconXicConstructor.FromCombinedMsAlign(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                        combinedMs2MsAlign, allScans, maxMass: 70000, xicSpline: new XicCubicSpline(0.05, 1, 1, true));

                    // apex RT tolerance 0.1 min, no overlap requirement, NO correlation filter. The fragment rank
                    // cap is left off deliberately: with correlationThreshold 0 nothing else limits the fragment
                    // list, so a cap would just report itself back rather than the true assignment count.
                    var diaParams = new DIAparameters(
                        DIAanalysisType.ISD, ms1Constructor, ms2Constructor,
                        new XicGroupingEngine(0.1f, 0, 0, maxThreadsForGrouping: 10),
                        PseudoMs2ConstructionType.Mass, combineFragments: true);
                    var commonParams = new CommonParameters { DIAparameters = diaParams };

                    var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();
                    TestContext.Progress.WriteLine($"{stem}: {pseudoScans.Length} pseudo-MS2 scans; MS1 " +
                        $"{ms1Constructor.FeatureCountBeforeFilters} features -> {ms1Constructor.FeatureCount} after filters " +
                        $"-> {ms1Constructor.LastXicCount} feature XICs " +
                        $"({ms1Constructor.LastFeaturesBelowMinPeaksCount} below minNumberOfPeaks)");

                    // fragments assigned per precursor - the number the grouping thresholds actually control
                    var fragmentCounts = pseudoScans.Select(s => s.ExperimentalFragments.Length).OrderBy(n => n).ToArray();
                    if (fragmentCounts.Length > 0)
                    {
                        TestContext.Progress.WriteLine($"   fragments/precursor: mean={fragmentCounts.Average():F1} " +
                            $"median={fragmentCounts[fragmentCounts.Length / 2]} " +
                            $"min={fragmentCounts[0]} max={fragmentCounts[^1]} " +
                            $"p25={fragmentCounts[fragmentCounts.Length / 4]} " +
                            $"p75={fragmentCounts[3 * fragmentCounts.Length / 4]} " +
                            $"p90={fragmentCounts[9 * fragmentCounts.Length / 10]}; " +
                            $"total assignments={fragmentCounts.Sum():N0}");
                    }

                    var msAlignOutPath = Path.Combine(outFolder, stem + "_pseudo_ms2.msalign");
                    var mgfOutPath = Path.Combine(outFolder, stem + "_pseudo_ms2.mgf");
                    if (write)
                    {
                        diaParams.WriteSettingsToFile(Path.Combine(outFolder, "DIAparameters.txt"));
                        ISDEngine.WriteMsAlignFile(msAlignOutPath, pseudoScans);
                        IsdMsAlignExporter.WriteMgf(pseudoScans, mgfOutPath);
                    }
                }
            }

            if (search)
            {
                string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
                var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DoPrecursorDeconvolution = false; // already deconvoluted via FromFile
                searchTask.CommonParameters.UseProvidedPrecursorInfo = true; // FromFile assembly provides the precursor info
                searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;

                var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters;

                // YB is yeast; add the human xml back if the data files are switched to the HC_YC mixes
                string yeast_xml = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.xml";
                var human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
                var databases = new List<DbForTask> { new DbForTask(yeast_xml, false) };

                var searchFileList = Directory.GetFiles(outFolder, "*_pseudo_ms2.mgf").ToList();
                var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) };
                var engine = new EverythingRunnerEngine(taskList, searchFileList, databases, outFolder);
                engine.Run();
            }
        }

        [Test]
        public static void SearchMgf()
        {
            string filePath = @"E:\ISD Project\Paper\Tentitative\Lysate_id\fullWindow\YB\Consensus_MM\rep1_xicGrouping_0.5_outsideFilter+2\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_pseudo_ms2.mgf";
            string dbPath = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            string combined_fasta = @"E:\Databases\human9606_yeast559292_reviewed_merged_2026_07_24.fasta";

            string outDir = @"E:\ISD Project\Paper\Tentitative\Lysate_id\fullWindow\YB\Consensus_MM\rep1_xicGrouping_0.5_outsideFilter+2";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);

            string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
            var searchTomlOpen = @"E:\ISD Project\ISD_250428\new_topDown_openSearch_toml\Task Settings\Task1-SearchTaskconfig.toml";
            var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
            searchTask.CommonParameters.DoPrecursorDeconvolution = false; // already deconvoluted via FromFile
            searchTask.CommonParameters.UseProvidedPrecursorInfo = true; // FromFile assembly provides the precursor info
            searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
            searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;
            //searchTask.SearchParameters.MinAllowedInternalFragmentLength = 0;

            var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_noFilter_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = searchTask.CommonParameters;

            var taskList = new List<(string, MetaMorpheusTask)> { ("search_YB", searchTask) }; //("GPTMD", gptmdTask)
            var engine = new EverythingRunnerEngine(taskList, new List<string> { filePath }, new List<DbForTask> { new DbForTask(dbPath, false) }, outDir);
            engine.Run();
        }

        [Test]
        public static void TestISD()
        {
            string file1 = @"E:\ISD Project\071626_ISD\07-16-26_0919YD_81min_ISD60-80-100_full400-1600.raw";
            var file2 = @"E:\ISD Project\071626_ISD\07-16-26_0919YD_81min_ISD60-80-100_preFilter800-1000-1200.raw";
            var file3 = @"E:\ISD Project\071626_ISD\07-16-26_0919YD_81min_ISD100_preFilter1200.raw";

            var fileList1 = new List<string> { file3 };
            var outputFolder = @"E:\Proteomics_software\TopPIC\toppic-windows-1.7.4\ISD\ISD_vs_DDA\YB_ISD\Consensus\rep1_umpire_MMsearch";
            var outDir = @"E:\ISD Project\071626_ISD\Test\0919YD\MM1";
            if (!Directory.Exists(outDir))
            {
                Directory.CreateDirectory(outDir);
            }

            string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
            var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
            var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_noFilter_toml, MetaMorpheusTask.tomlConfig);

            //DIA parameters
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), traceToleranceDa: 1.5, minMass: 4000, xicSpline: new XicCubicSpline(0.2, 1, 1, true)),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new XicGroupingEngine(0.2f, 0.3, 0.5, maxThreadsForGrouping: 10),
                new XicGroupingEngine(0.3f, 0, 0.5, maxThreadsForGrouping: 10), PseudoMs2ConstructionType.Mass, combineFragments: true);

            searchTask.CommonParameters = new CommonParameters { DIAparameters = diaParams };
            searchTask.CommonParameters.PrecursorMassTolerance = new PpmTolerance(10);
            gptmdTask.CommonParameters = searchTask.CommonParameters;

            var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) }; //("GPTMD", gptmdTask)
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";

            var engine = new EverythingRunnerEngine(taskList, fileList1, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outDir);
            engine.Run();
        }


        [Test]
        public static void SearchISDConsensus()
        {
            var path1 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1.raw";
            string dbPath = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            string outDir = @"E:\ISD Project\Paper\Tentitative\Lysate_id\YB\Consensus_MM\rep1_gptmd";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);

            string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
            var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
            var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_noFilter_toml, MetaMorpheusTask.tomlConfig);

            //Figure out the correct traceToleranceDa, default is 0.02 Da, should use something similar from consensus paper
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, xicSpline: new Bspline(2, 150), traceToleranceDa: 1.5),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), aggregateCharges: false, xicSpline: new Bspline(2, 150), traceToleranceDa: 1.5),
                new UmpirePfGroupingEngine(150, 0.2f, 0.3, 0.7, maxThreadsForGrouping: 10, fragmentRankThreshold: 200),
                PseudoMs2ConstructionType.Mass, combineFragments: true);
            searchTask.CommonParameters.DIAparameters = diaParams;

            var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", searchTask)}; //("GPTMD", gptmdTask)
            var engine = new EverythingRunnerEngine(taskList, new List<string> { path1 }, new List<DbForTask> { new DbForTask(dbPath, false) }, outDir);
            engine.Run();
        }

        [Test]
        public static void RNA()
        {
            var filePath = @"R:\Zhuoxin\2026_08_14-mix_134min_ISD40-60-80_filter.raw";
            string AllStandards_fasta = @"E:\ISD Project\RNA\AllStandards.fasta";
            var databases = new List<DbForTask> { new DbForTask(AllStandards_fasta, false) };
            var fileList = new List<string> { filePath };

            var searchToml = @"E:\ISD Project\RNA\2026-08-15-10-48-07\Task Settings\Task1-RnaSearchTaskconfig.toml";
            var searchTask = Toml.ReadFile<SearchTask>(searchToml, MetaMorpheusTask.tomlConfig);
            var outDir = @"E:\ISD Project\RNA\ISD_search";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);

            bool group = true;
            bool write = true;
            bool search = true;
            foreach (var file in fileList)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();

                // Consensus feature tracing as the XIC front-end for both the intact (precursor) and fragment channels.
                var diaParams = new DIAparameters(
                    DIAanalysisType.ISD,
                    //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, searchTask.CommonParameters.PrecursorDeconvolutionParameters, traceToleranceDa: 1.5, minMass: 4000,xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 2, 5, searchTask.CommonParameters.PrecursorDeconvolutionParameters, minMass: 500, minCharge: -40, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, searchTask.CommonParameters.ProductDeconvolutionParameters, minCharge: -10, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    //new XicGroupingEngine(0.2f, 0.3, 0.5, maxThreadsForGrouping: 10),
                    new XicGroupingEngine(0.2f, 0, 0, maxThreadsForGrouping: 10), PseudoMs2ConstructionType.Mass, combineFragments: true);
                searchTask.CommonParameters.DIAparameters = diaParams;

                var outFolder = Path.Combine(outDir, Path.GetFileNameWithoutExtension(dataFile.FilePath), "neutralMass3");
                if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

                var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) }; //("GPTMD", gptmdTask)
                var engine = new EverythingRunnerEngine(taskList, fileList, databases, outFolder);
                engine.Run();
            }
        }

        /// <summary>
        /// Environment-driven harness: generate consensus-path pseudo-MS2 MGFs at several precursor-fragment
        /// correlation thresholds, reusing one consensus-XIC build. Set ISD_MZML (input), ISD_OUTDIR (output),
        /// and optionally ISD_CORRS (comma list, default "0,0.4,0.6,0.7"). Writes consensus_corr&lt;t&gt;.mgf per
        /// threshold; each MGF is then searched externally with the top-down td_pseudoMS2 config to count IDs.
        /// Ignored unless ISD_MZML is set, so it never runs in normal CI.
        /// </summary>
        //[Test]
        //public static void ConsensusIdSweep_FromEnv()
        //{
        //    string mzml = Environment.GetEnvironmentVariable("ISD_MZML");
        //    if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set ISD_MZML to run this harness"); return; }
        //    string outDir = Environment.GetEnvironmentVariable("ISD_OUTDIR") ?? TestContext.CurrentContext.TestDirectory;
        //    Directory.CreateDirectory(outDir);
        //    var corrs = (Environment.GetEnvironmentVariable("ISD_CORRS") ?? "0,0.4,0.6,0.7")
        //        .Split(',').Select(s => double.Parse(s, System.Globalization.CultureInfo.InvariantCulture)).ToArray();
        //    float apex = 0.2f;
        //    double overlap = 0.5;
        //    var decon = new ClassicDeconvolutionParameters(1, 30, 4, 3);

        //    var dataFile = MsDataFileReader.GetDataFile(mzml);
        //    dataFile.LoadAllStaticData();
        //    var allScans = dataFile.GetAllScansList().ToArray();
        //    var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
        //    ISDEngine.ReLabelIsdScans(isdMap, allScans);

        //    // consensus feature tracing -> XICs, built ONCE and reused across thresholds.
        //    // PRECURSOR channel: keep only intact proteoforms (>=3 kDa, >=3 charge states) so small noise
        //    // features do not become precursors. FRAGMENT channel: unfiltered (fragments are low mass/charge).
        //    double precMinMass = double.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINMASS") ?? "3000",
        //        System.Globalization.CultureInfo.InvariantCulture);
        //    int precMinCharge = int.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINCHARGE") ?? "3");
        //    var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
        //        minMass: precMinMass, minChargeCount: precMinCharge);
        //    var ms2Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon);
        //    // ISD_FRAG_AGG=false -> fragments use mass tracing only (one XIC per charge-locked trace, no
        //    // cross-charge aggregation); precursors always use full feature tracing.
        //    bool fragAgg = (Environment.GetEnvironmentVariable("ISD_FRAG_AGG") ?? "true").ToLowerInvariant() != "false";
        //    ms2Ctor.AggregateCharges = fragAgg;
        //    var ms1Xics = ms1Ctor.GetAllXics(ms1Scans);
        //    var ms2Xics = isdMap.Values.SelectMany(v => ms2Ctor.GetAllXics(v.ToArray())).ToList();
        //    TestContext.WriteLine($"consensus XICs: {ms1Xics.Count} precursor (>= {precMinMass} Da, >= {precMinCharge} charges), {ms2Xics.Count} fragment (chargeAggregation={fragAgg})");

        //    foreach (var corr in corrs)
        //    {
        //        var grouper = new XicGroupingEngine(apex, overlap, corr, Environment.ProcessorCount);
        //        var groups = grouper.PrecursorFragmentGrouping(ms1Xics, ms2Xics).ToList();
        //        int idx = 0;
        //        var scans = new List<Ms2ScanWithSpecificMass>();
        //        foreach (var g in groups)
        //        {
        //            idx++;
        //            g.PFgroupIndex = idx;
        //            scans.Add(g.GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.Mass, new CommonParameters(), mzml));
        //        }
        //        string outMgf = Path.Combine(outDir, $"consensus_corr{corr.ToString(System.Globalization.CultureInfo.InvariantCulture)}.mgf");
        //        IsdMsAlignExporter.WriteMgf(scans, outMgf);
        //        TestContext.WriteLine($"corr={corr}: {scans.Count} pseudo scans -> {outMgf}");
        //    }
        //}

    }
}
