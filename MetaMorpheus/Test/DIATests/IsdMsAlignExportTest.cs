using Chemistry;
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
using System.IO;
using System.Linq;
using TaskLayer;

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
        public static void ExportIsdPseudoScansToMsAlign()
        {
            string mzmlPath = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1.raw";
            Assert.That(File.Exists(mzmlPath), $"missing test fixture: {mzmlPath}");

            var dataFile = MsDataFileReader.GetDataFile(mzmlPath);
            dataFile.LoadAllStaticData();

            // Consensus feature tracing as the XIC front-end for both the intact (precursor) and fragment channels.
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), minMass: 4000, xicSpline: new Bspline(2, 150)),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new Bspline(2, 150)),
                new XicGroupingEngine(0.2f, 0.3, 0.7, maxThreadsForGrouping: 10),
                PseudoMs2ConstructionType.Mass,
                combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, mzmlPath, commonParams).ToArray();
            var outDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.7.4\ISD\ISD_vs_DDA\YB_ISD\Consensus\rep1_xicGrouping";
            if (!Directory.Exists(outDir))
            {
                Directory.CreateDirectory(outDir);
            }
            var pseudoScansFileName = Path.GetFileNameWithoutExtension(mzmlPath) + "_pseudo_ms2.msalign";
            var outPath = Path.Combine(outDir, pseudoScansFileName);
            ISDEngine.WriteMsAlignFile(outPath, pseudoScans);
        }

        [Test]
        public static void TestISD()
        {
            var path1 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1.raw";

            var fileList1 = new List<string> { path1 };
            var outputFolder = @"E:\Proteomics_software\TopPIC\toppic-windows-1.7.4\ISD\ISD_vs_DDA\YB_ISD\Consensus\rep1_umpire_MMsearch";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string tomlDir = @"E:\GitClones\MetaMorpheus\isd-search-reproduction\tomls";
            var searchTask = Toml.ReadFile<SearchTask>(Path.Combine(tomlDir, "td_pseudoMS2_search.toml"), MetaMorpheusTask.tomlConfig);
            var gptmdTask = Toml.ReadFile<GptmdTask>(Path.Combine(tomlDir, "gptmd_isd.toml"), MetaMorpheusTask.tomlConfig);

            //DIA parameters
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), minMass: 4000, xicSpline: new Bspline(2, 150)),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new Bspline(2, 150)),
                new UmpirePfGroupingEngine(150, 0.2f, 0.3, 0.7, maxThreadsForGrouping: 10),
                PseudoMs2ConstructionType.Mass,
                combineFragments: true);
            searchTask.CommonParameters = new CommonParameters { DIAparameters = diaParams };
            searchTask.CommonParameters.PrecursorMassTolerance = new PpmTolerance(10);
            gptmdTask.CommonParameters = searchTask.CommonParameters;

            var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) }; //("GPTMD", gptmdTask)
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";

            var engine = new EverythingRunnerEngine(taskList, fileList1, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder);
            engine.Run();
        }

        /// <summary>
        /// Environment-driven harness: generate consensus-path pseudo-MS2 MGFs at several precursor-fragment
        /// correlation thresholds, reusing one consensus-XIC build. Set ISD_MZML (input), ISD_OUTDIR (output),
        /// and optionally ISD_CORRS (comma list, default "0,0.4,0.6,0.7"). Writes consensus_corr&lt;t&gt;.mgf per
        /// threshold; each MGF is then searched externally with the top-down td_pseudoMS2 config to count IDs.
        /// Ignored unless ISD_MZML is set, so it never runs in normal CI.
        /// </summary>
        [Test]
        public static void ConsensusIdSweep_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("ISD_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set ISD_MZML to run this harness"); return; }
            string outDir = Environment.GetEnvironmentVariable("ISD_OUTDIR") ?? TestContext.CurrentContext.TestDirectory;
            Directory.CreateDirectory(outDir);
            var corrs = (Environment.GetEnvironmentVariable("ISD_CORRS") ?? "0,0.4,0.6,0.7")
                .Split(',').Select(s => double.Parse(s, System.Globalization.CultureInfo.InvariantCulture)).ToArray();
            float apex = 0.2f;
            double overlap = 0.5;
            var decon = new ClassicDeconvolutionParameters(1, 30, 4, 3);

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // consensus feature tracing -> XICs, built ONCE and reused across thresholds.
            // PRECURSOR channel: keep only intact proteoforms (>=3 kDa, >=3 charge states) so small noise
            // features do not become precursors. FRAGMENT channel: unfiltered (fragments are low mass/charge).
            double precMinMass = double.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINMASS") ?? "3000",
                System.Globalization.CultureInfo.InvariantCulture);
            int precMinCharge = int.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINCHARGE") ?? "3");
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
                minMass: precMinMass, minChargeCount: precMinCharge);
            var ms2Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon);
            // ISD_FRAG_AGG=false -> fragments use mass tracing only (one XIC per charge-locked trace, no
            // cross-charge aggregation); precursors always use full feature tracing.
            bool fragAgg = (Environment.GetEnvironmentVariable("ISD_FRAG_AGG") ?? "true").ToLowerInvariant() != "false";
            ms2Ctor.AggregateCharges = fragAgg;
            var ms1Xics = ms1Ctor.GetAllXics(ms1Scans);
            var ms2Xics = isdMap.Values.SelectMany(v => ms2Ctor.GetAllXics(v.ToArray())).ToList();
            TestContext.WriteLine($"consensus XICs: {ms1Xics.Count} precursor (>= {precMinMass} Da, >= {precMinCharge} charges), {ms2Xics.Count} fragment (chargeAggregation={fragAgg})");

            foreach (var corr in corrs)
            {
                var grouper = new XicGroupingEngine(apex, overlap, corr, Environment.ProcessorCount);
                var groups = grouper.PrecursorFragmentGrouping(ms1Xics, ms2Xics).ToList();
                int idx = 0;
                var scans = new List<Ms2ScanWithSpecificMass>();
                foreach (var g in groups)
                {
                    idx++;
                    g.PFgroupIndex = idx;
                    scans.Add(g.GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.Mass, new CommonParameters(), mzml));
                }
                string outMgf = Path.Combine(outDir, $"consensus_corr{corr.ToString(System.Globalization.CultureInfo.InvariantCulture)}.mgf");
                IsdMsAlignExporter.WriteMgf(scans, outMgf);
                TestContext.WriteLine($"corr={corr}: {scans.Count} pseudo scans -> {outMgf}");
            }
        }

        /// <summary>
        /// Search DDA the CONSENSUS-PAPER way: precursor feature tracing + the built-in "FromFile" precursor
        /// assembly. Consensus-trace the DDA MS1 into features, write them as a `.ms1.feature` file, then use
        /// mzLib's real <c>FromFileDeconvolutionParameters</c> join — `ms2.GetIsolatedMassesAndCharges(ms1, ff)`
        /// — to assemble each MS2's precursor exactly as MetaMorpheus's precursor pipeline does (matching
        /// features to the isolation window + precursor RT), then attach the scan's real fragments. This is the
        /// mechanism from the consensus-deconvolution paper's ExternalMs1Features / FromFile search, as opposed
        /// to the hand-rolled nearest-feature matching in SearchDdaWithConsensusPrecursors_FromEnv.
        /// Env-driven: DDAFF_MZML (input), DDAFF_OUT (output MGF). Ignored unless DDAFF_MZML is set.
        /// </summary>
        [Test]
        public static void SearchDdaWithFromFileConsensusFeatures_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("DDAFF_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set DDAFF_MZML to run this harness"); return; }
            string outMgf = Environment.GetEnvironmentVariable("DDAFF_OUT")
                ?? Path.Combine(TestContext.CurrentContext.TestDirectory, "dda_fromfile.mgf");
            string featPath = Path.Combine(Path.GetDirectoryName(outMgf) ?? ".",
                Path.GetFileNameWithoutExtension(outMgf) + "_ms1.feature"); // reader keys on the "_ms1.feature" suffix
            var decon = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            var fragParams = new CommonParameters(productDeconParams: new ClassicDeconvolutionParameters(1, 20, 4, 3));

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2Scans = allScans.Where(s => s.MsnOrder == 2).ToArray();

            // (A) precursor feature tracing (consensus) -> external _ms1.feature file (intact precursors only)
            var features = ConsensusMassXicConstructor.TraceFeatures(ms1Scans, decon);
            IsdMsAlignExporter.WriteMs1FeatureFile(features, featPath, minMass: 3000, minChargeCount: 3);

            // (B) assemble Ms2ScanWithSpecificMass via mzLib's FromFile join (the consensus-paper mechanism)
            var ff = new Readers.FromFileDeconvolutionParameters(featPath, 1, 60);
            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            int assigned = 0;
            foreach (var ms2 in ms2Scans)
            {
                if (!ms2.OneBasedPrecursorScanNumber.HasValue) continue;
                var ms1 = dataFile.GetOneBasedScan(ms2.OneBasedPrecursorScanNumber.Value);
                IsotopicEnvelope[] frags = null;
                foreach (var env in ms2.GetIsolatedMassesAndCharges(ms1, ff))
                {
                    frags ??= Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2, fragParams);
                    assigned++;
                    double precMz = env.MonoisotopicMass.ToMz(env.Charge);
                    pseudoScans.Add(new Ms2ScanWithSpecificMass(ms2, precMz, env.Charge, mzml, fragParams, frags,
                        env.TotalIntensity, env.Peaks.Count));
                }
            }
            IsdMsAlignExporter.WriteMgf(pseudoScans, outMgf);
            TestContext.WriteLine($"FromFile consensus: {features.Count} features, {assigned} precursor assignments over {ms2Scans.Length} MS2 -> {outMgf}");
            Assert.That(assigned, Is.GreaterThan(0), "no MS2 scans got a FromFile precursor");
        }
    }
}
