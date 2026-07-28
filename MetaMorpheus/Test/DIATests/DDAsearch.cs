using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using EngineLayer;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;
using Chemistry;

namespace Test.DIATests
{
    internal class DDAsearch
    {
        /// <summary>
        /// COMPLETE, self-contained DDA search via the consensus-paper method — no external toml or CMD needed.
        /// Does the whole pipeline in-process: consensus-trace the DDA MS1 → write a `_ms1.feature` → assemble
        /// each MS2's precursor with mzLib's FromFile join → write an MGF → run a real MetaMorpheus top-down
        /// <see cref="SearchTask"/> (config built in code below) against a protein database → report proteoform
        /// and PSM counts at 1% FDR, exactly like a normal MetaMorpheus search.
        /// Env-driven: DDAFF_MZML (DDA input), DDAFF_DB (protein .xml/.fasta), DDAFF_OUT (output dir, optional).
        /// Ignored unless both DDAFF_MZML and DDAFF_DB are set.
        /// </summary>
        [Test]
        public static void SearchDdaWithFromFileConsensusFeatures_Complete_FromEnv()
        {
            string filePath = @"E:\ISD Project\071626_ISD\07-18-26_YBCDE_81min_DDA_HCD30.raw";
            string dbPath = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";

            var fastYB = @"E:\ISD Project\ISD_250428\05-02-25_PEPPI-YB_60min_DDA.raw";

            var fileList = new List<string> { filePath };
            string outDir = @"E:\ISD Project\071626_ISD\YBCDE";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);

            foreach(var file in fileList)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();

                var outFolder = Path.Combine(outDir, Path.GetFileNameWithoutExtension(dataFile.FilePath), "consensus");
                if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

                string outMgf = Path.Combine(outFolder, "dda_fromfile.mgf");
                string featPath = Path.Combine(outFolder, "dda_fromfile_ms1.feature");

                // ---- (1) precursor feature tracing + FromFile assembly -> MGF ----
                var decon = new ClassicDeconvolutionParameters(1, 60, 4, 3);

                var fragParams = new CommonParameters(productDeconParams: new ClassicDeconvolutionParameters(1, 20, 4, 3), assumeOrphanPeaksAreZ1Fragments: false);
                var allScans = dataFile.GetAllScansList().ToArray();
                var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
                var ms2Scans = allScans.Where(s => s.MsnOrder == 2).ToArray();
                var features = ConsensusMassXicConstructor.TraceFeatures(ms1Scans, decon, traceToleranceDa: 1.5);
                IsdMsAlignExporter.WriteMs1FeatureFile(features, featPath, minMass: 4000, minChargeCount: 3);
                var ff = new FromFileDeconvolutionParameters(featPath, 1, 60);
                var pseudoScans = new List<Ms2ScanWithSpecificMass>();
                foreach (var ms2 in ms2Scans)
                {
                    if (!ms2.OneBasedPrecursorScanNumber.HasValue) continue;
                    var ms1 = dataFile.GetOneBasedScan(ms2.OneBasedPrecursorScanNumber.Value);
                    IsotopicEnvelope[] frags = null;
                    foreach (var env in ms2.GetIsolatedMassesAndCharges(ms1, ff))
                    {
                        frags ??= Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2, fragParams);
                        pseudoScans.Add(new Ms2ScanWithSpecificMass(ms2, env.MonoisotopicMass.ToMz(env.Charge),
                            env.Charge, filePath, fragParams, frags, env.TotalIntensity, env.Peaks.Count));
                    }
                }
                IsdMsAlignExporter.WriteMgf(pseudoScans, outMgf);
                var outMsAlign = Path.Combine(outFolder, "dda_fromfile_ms2.msalign");
                ISDEngine.WriteMsAlignFile(outMsAlign, pseudoScans);

                // ---- (2) a complete top-down MetaMorpheus search (config in code — the "toml" as C#) ----
                string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
                var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DoPrecursorDeconvolution = false; // already deconvoluted via FromFile
                searchTask.CommonParameters.UseProvidedPrecursorInfo = true; // FromFile assembly provides the precursor info
                searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;

                var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters;

                var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) }; //("GPTMD", gptmdTask)
                var engine = new EverythingRunnerEngine(taskList, new List<string> { outMgf }, new List<DbForTask> { new DbForTask(dbPath, false) }, outFolder);
                engine.Run();
            }
        }
    }
}
