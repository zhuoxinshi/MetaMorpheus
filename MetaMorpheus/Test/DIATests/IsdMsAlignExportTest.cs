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
using System.IO;
using System.Linq;
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
            string file1 = @"E:\ISD Project\071626_ISD\Quant\PierceYB\07-30-26_PierceYB_0-005_ISD_81min_rep1.raw";
            var file2 = @"E:\ISD Project\071626_ISD\Quant\PierceYB\07-30-26_PierceYB_0-005_ISD_81min_rep2.raw";
            var file3 = @"E:\ISD Project\071626_ISD\Quant\PierceYB\07-30-26_PierceYB_0-005_ISD_81min_rep3.raw";
            var file4 = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\07-25-26_HD_81min_ISD_1-4.raw";

            var fastYB = @"E:\ISD Project\ISD_250428\05-02-25_PEPPI-YB_60min_ISD60-80-100_preFilter700-900-1100_RF.raw";
            var fileList = new List<string> { file1, file2, file3 };
            var outDir = @"E:\ISD Project\071626_ISD\Quant\PierceYB";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);

            foreach (var file in fileList)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();

                // Consensus feature tracing as the XIC front-end for both the intact (precursor) and fragment channels.
                var diaParams = new DIAparameters(
                    DIAanalysisType.ISD,
                    new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), traceToleranceDa: 1.5, minMass: 4000, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    //new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                    //new XicGroupingEngine(0.2f, 0.3, 0.5, maxThreadsForGrouping: 10),
                    new XicGroupingEngine(0.3f, 0, 0.5, maxThreadsForGrouping: 10), PseudoMs2ConstructionType.Mass, combineFragments: true);
                var commonParams = new CommonParameters { DIAparameters = diaParams };

                var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();
                var outFolder = Path.Combine(outDir, Path.GetFileNameWithoutExtension(dataFile.FilePath), "consensus_neutralMass2");
                if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

                // record the exact DIA/ISD engine settings used to build these pseudo-scans (the MGF search below
                // runs without DIAparameters, so the auto-dump in SearchTask does not fire for this flow)
                diaParams.WriteSettingsToFile(Path.Combine(outFolder, "DIAparameters.txt"));

                var pseudoScansFileName = Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.msalign";

                var msAlignOutPath = Path.Combine(outFolder, pseudoScansFileName);
                ISDEngine.WriteMsAlignFile(msAlignOutPath, pseudoScans);
                var mgfOutPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.mgf");
                IsdMsAlignExporter.WriteMgf(pseudoScans, mgfOutPath);

                string dbPath = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
                string pierceDb = @"E:\ISD Project\071626_ISD\Quant\Pierce\Pierce_6pro.fasta";
                string combined_fasta = @"E:\Databases\human9606_yeast559292_reviewed_merged_2026_07_24.fasta";
                string human_fasta = @"E:\Databases\Human_9606.fasta";
                var human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";

                string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
                var searchTomlOpen = @"E:\ISD Project\ISD_250428\new_topDown_openSearch_toml\Task Settings\Task1-SearchTaskconfig.toml";
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

                var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) }; //("GPTMD", gptmdTask)
                var engine = new EverythingRunnerEngine(taskList, new List<string> { mgfOutPath }, new List<DbForTask> { new DbForTask(dbPath, false) }, outFolder);
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
