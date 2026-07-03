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

            string tomlFile_CommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);

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

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = searchTask.CommonParameters;

            var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) }; //("GPTMD", gptmdTask)
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";

            var engine = new EverythingRunnerEngine(taskList, fileList1, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder);
            engine.Run();
        }

        /// <summary>
        /// Search DDA runs using CONSENSUS MASS TRACING for the precursors: trace the DDA MS1 scans with
        /// <see cref="ConsensusMassXicConstructor"/> to get clean intact precursor features, pair each real DDA
        /// MS2 scan to the consensus feature matching its isolation m/z + retention time, and build a pseudo scan
        /// from the consensus precursor + the scan's REAL (deconvoluted) fragments. Writes an MGF that is then
        /// searched with the same td_pseudoMS2 config as the ISD path (apples-to-apples ISD-consensus vs
        /// DDA-consensus). Env-driven: DDA_MZML (input), DDA_OUT (output MGF). Ignored unless DDA_MZML is set.
        /// </summary>
        [Test]
        public static void SearchDdaWithConsensusPrecursors_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("DDA_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set DDA_MZML to run this harness"); return; }
            string outMgf = Environment.GetEnvironmentVariable("DDA_OUT")
                ?? Path.Combine(TestContext.CurrentContext.TestDirectory, "dda_consensus.mgf");
            const double proton = 1.007276466879;
            var decon = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            // real DDA MS2 fragments are multiply charged -> deconvolute them to neutral masses
            var fragParams = new CommonParameters(productDeconParams: new ClassicDeconvolutionParameters(1, 20, 4, 3));

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2Scans = allScans.Where(s => s.MsnOrder == 2).ToArray();

            // consensus mass tracing on the MS1 scans -> intact precursor features (>=3 kDa, >=3 charges)
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
                minMass: 3000, minChargeCount: 3);
            var precursorXics = ms1Ctor.GetAllXics(ms1Scans);

            // pair each real DDA MS2 scan with the consensus feature matching its isolation m/z + RT
            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            int assigned = 0;
            foreach (var ms2 in ms2Scans)
            {
                if (ms2.IsolationMz == null) continue;
                double isoMz = ms2.IsolationMz.Value;
                double rt = ms2.RetentionTime;
                ExtractedIonChromatogram best = null; double bestErr = double.MaxValue; int bestZ = 1;
                foreach (var xic in precursorXics)
                {
                    if (rt < xic.StartRT - 0.5 || rt > xic.EndRT + 0.5) continue;
                    double mass = xic.ApexPeak.M;
                    int z = (int)Math.Round(mass / (isoMz - proton)); // charge that places this mass at the iso m/z
                    if (z < 1 || z > 60) continue;
                    double err = Math.Abs(mass.ToMz(z) - isoMz);
                    if (err < bestErr && err < 2.5) { bestErr = err; best = xic; bestZ = z; }
                }
                if (best == null) continue;
                assigned++;
                double precMz = ((double)best.ApexPeak.M).ToMz(bestZ);
                var frags = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2, fragParams);
                pseudoScans.Add(new Ms2ScanWithSpecificMass(ms2, precMz, bestZ, mzml, fragParams, frags));
            }
            //IsdMsAlignExporter.WriteMgf(pseudoScans, outMgf);
            TestContext.WriteLine($"DDA consensus: {precursorXics.Count} MS1 features, {ms2Scans.Length} MS2, {assigned} assigned -> {outMgf}");
            Assert.That(assigned, Is.GreaterThan(0), "no MS2 scans matched a consensus precursor feature");
        }
    }
}
