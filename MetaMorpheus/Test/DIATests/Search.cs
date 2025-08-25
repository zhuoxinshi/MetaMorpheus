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

            SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);
            searchTask.CommonParameters.PrecursorMassTolerance = new PpmTolerance(10);
            var ms1XicConstructor = new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 2, 2), 2, 0.5, 3, searchTask.CommonParameters.PrecursorDeconvolutionParameters, 4000, 4, new XicLinearSpline(0.05, numberOfPeaksToAdd: 1));
            var ms2XicConstructor = new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 2, 2), 2, 0.5, 3, searchTask.CommonParameters.ProductDeconvolutionParameters, 0, 1, new XicLinearSpline(0.05, numberOfPeaksToAdd: 1));
            var pfGroupingEngine = new XicGroupingEngine(0.15f, 0.2, 0.5, 10, 10);
            searchTask.CommonParameters.DIAparameters = new DIAparameters(ms1XicConstructor, ms2XicConstructor, pfGroupingEngine, PseudoMs2ConstructionType.Mass);

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
    }
}
