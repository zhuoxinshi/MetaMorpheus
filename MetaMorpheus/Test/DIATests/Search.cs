using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;
using TaskLayer;
using System.IO;
using Chemistry;
using Readers;
using MathNet.Numerics.Interpolation;
using System.Drawing.Imaging;
using Nett;

namespace Test.DIATests
{
    public class Search
    {
        [Test]
        public static void SearchCE()
        {
            var filePath1 = @"E:\CE\250730_CE\08-19-25_CE_ammon-acet_PEPPI-YB_500nL-pHjunction_ISD60-80_preFilter700-900-1100.raw";
            var fileList = new List<string> { filePath1};
            var outputFolder = @"E:\ISD Project\TestSearch\DIAupdate\test1_ISD";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";

            SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);
            searchTask.CommonParameters.PrecursorMassTolerance = new PpmTolerance(10);
            var ms1XicConstructor = new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 2, 2), 2, 0.5, 3, searchTask.CommonParameters.PrecursorDeconvolutionParameters, 4000, 4, new XicLinearSpline(0.05, numberOfPeaksToAdd: 1));
            var ms2XicConstructor = new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 2, 2), 2, 0.5, 3, searchTask.CommonParameters.ProductDeconvolutionParameters, 0, 1, new XicLinearSpline(0.05, numberOfPeaksToAdd: 1));
            var pfGroupingEngine = new XicGroupingEngine(0.15f, 0.2, 0.5, 10, 10);
            searchTask.CommonParameters.DIAparameters = new DIAparameters(AnalysisType.ISD, ms1XicConstructor, ms2XicConstructor, pfGroupingEngine, PseudoMs2ConstructionType.Mass);

            //match all charge fragment ions
            searchTask.SearchParameters.WriteSpectralLibrary = true;

            //Use IsoDec
            //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
            //task.CommonParameters.ProductDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
            //task.CommonParameters.ProductDeconvolutionParameters.MaxAssumedChargeState = 20;

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var GPTMD_max2 = @"E:\ISD Project\ISD_250428\GPTMD_max2\Task Settings\Task1-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = searchTask.CommonParameters;

            var taskList = new List<(string, MetaMorpheusTask)> {("search", searchTask) }; //("GPTMD", gptmdTask), ("Calibration", cali_task),

            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var yeast_fasta = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";

            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder);
            engine.Run();
        }

        [Test]
        public static void TestBottomUpDIASearch()
        {
            string tomlFile = @"E:\Aneuploidy\searchToml_commonFixedVariable_noTrim_writeLib\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            searchTask.CommonParameters.PrecursorMassTolerance = new PpmTolerance(10);
            string outputFolder = @"E:\DIA\TestSearch\bottomUp_update\oldData\umpire_try_25-300_afterMerge";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            var ms1XicConstructor = new NeutralMassXicConstructor(new PpmTolerance(10), 1, 0.5, 3, searchTask.CommonParameters.PrecursorDeconvolutionParameters, 0, 1, new Bspline(2, 150));
            var ms2XicConstructor = new MzPeakXicConstructor(new PpmTolerance(20), 1, 0.5, 3, new Bspline(2, 150));
            //var xicGroupingEngine = new XicGroupingEngine(0.2f, 0.2, 0.5, 10, 0, precursorRankThreshold: 10, fragmentRankThreshold: 200);
            var umpireGroupingEngine = new UmpirePfGroupingEngine(150, 0.3f, 0.2, 0.5, 10, 0, 25, 300);
            searchTask.CommonParameters.DIAparameters = new DIAparameters(AnalysisType.DIA, ms1XicConstructor, ms2XicConstructor, umpireGroupingEngine, PseudoMs2ConstructionType.MzPeak);

            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string umpireFile = @"E:\DIA\DIA-Umpire data\18300_REP2_500ng_HumanLysate_SWATH_1.mzML";
            string benchmarkDb = @"E:\REF_EColi_K12_UPS1_combined.fasta";
            string humanDb = @"E:\ISD Project\Claire's human data\Human_9606.fasta";
            searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(humanDb, false) }, new List<string> { DIAfile }, "test");
        }
    }
}
