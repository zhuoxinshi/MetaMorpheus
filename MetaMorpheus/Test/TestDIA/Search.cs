using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using EngineLayer;
using NUnit.Framework;
using System.IO;
using TaskLayer;
using EngineLayer.DIA;
using static System.Net.WebRequestMethods;
using MzLibUtil;
using Plotly.NET;
using Nett;
using System.IO.Compression;
using static Plotly.NET.StyleParam.Range;
using System.Windows;
using FlashLFQ;
using static iText.StyledXmlParser.Jsoup.Select.Evaluator;

namespace Test.TestDIA
{
    public class Search
    {
        [Test]
        public static void TestBottomUpDIASearch()
        {
            string tomlFile = @"E:\DIA\Data\DIA_bu_250114\bu-DIA_toml\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.TrimMsMsPeaks = false;
            task.CommonParameters.TrimMs1Peaks = false;
            string outputFolder = @"E:\DIA\TestSearch\DIAEngine_static_deconHighestPeak-Peak_Umpire_corr0.5_apex0.2_overlap0.3_frag500_pre25_maxRT0.6_num2";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 
                maxNumMissedScan:1, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance:0.2, 
                fragmentRankCutOff: 500, precursorRankCutOff: 25, maxRTrangeMS1: 0.6, maxRTrangeMS2: 0.6, highCorrThreshold: 0.5, numHighCorrFragments: 0, 
                precursorIntensityCutOff: 0, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.05f, minMS1Mass: 0, maxMass: 10000, type: "DIA", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, cutMs1Peaks: false, cutMs2Peaks: false, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.DIAEngine_static, ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, 
                splineRtInterval: 0.005, ms1NumPeaksThreshold: 2);
            // Use reflection to set max threads
            //task.CommonParameters.GetType().GetProperty("MaxThreadsToUsePerFile").SetMethod.Invoke(task.CommonParameters, new object[] { 1 });
            //var type = task.CommonParameters.GetType();
            //var property = type.GetProperty("MaxThreadsToUsePerFile");
            //property.SetMethod.Invoke(task.CommonParameters, new object[] { 1 });

            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string benchmarkDb = @"E:\REF_EColi_K12_UPS1_combined.fasta";
            string humanDb = @"E:\ISD Project\Claire's human data\Human_9606.fasta";
            string yeast_xml = @"E:\ISD Project\ISD_240812\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(humanDb, false) }, new List<string> { DIAfile }, "test");
        }

        [Test]
        public static void TestTopDownDIASearch()
        {
            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_noFixedMods = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_withFixedMods = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_withFixedMods, MetaMorpheusTask.tomlConfig);
            string outputFolder = @"E:\DIA\TestSearch\050425YB\Ms1SpaceCubicSpline_GPTMD_20ppm_corr0.75_apexRT0.5_overlap0.2_num4";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string DIAfile1 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep1.raw";
            string DIAfile2 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep2.raw";
            string DIAfile3 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep3.raw";
            string DIAfile4 = @"E:\ISD Project\ISD_250428\0504YB_rep_DIA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep1-calib-averaged.mzML";
            string DIAfile5 = @"E:\ISD Project\ISD_250428\0504YB_rep_DIA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep2-calib-averaged.mzML";
            string DIAfile6 = @"E:\ISD Project\ISD_250428\0504YB_rep_DIA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep3-calib-averaged.mzML";
            string DIAfile7 = @"E:\ISD Project\ISD_250428\05-01-25_PEPPI-YC_105min_td-DIA_25mz-overlap2_HCD21-23-25_AGC1e6.raw";
            string DIAfile8 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep1_Ms2Averaged.mzML";
            string topDIAfile = @"E:\DIA\TopDIA\20231117_DIA_720_800_rep2.mzML";
            var fileList = new List<string> {DIAfile8}; //DIAfile2, DIAfile3, DIAfile4, DIAfile5, DIAfile6 

            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.75, apexRtTolerance: 0.5,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 4000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.DIAEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.CubicSpline, ms2SplineType: SplineType.Ms1SpaceCubicSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 4, ms2NumPeaksThreshold: 4,
                combineFragments: false,rankFilter: false, minPFpairCount: 10);

           //set threads
           task.CommonParameters.MaxThreadsToUsePerFile = 15;
            //match all charge fragment ions
            task.SearchParameters.WriteSpectralLibrary = true;

            //Use IsoDec
            //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
            //task.CommonParameters.ProductDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
            //task.CommonParameters.ProductDeconvolutionParameters.MaxAssumedChargeState = 20;

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var GPTMD_noFixedMod = @"E:\DIA\TopDIA\2025-05-10-14-22-53\Task Settings\Task1-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = task.CommonParameters.Clone();
            var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", task) };//("GPTMD", gptmdTask),

            string YC_gptmd = @"E:\CE\250318_CE\YC_cali-avged_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string gptmdDb = @"E:\ISD Project\ISD_250428\0428YB_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            string ecoli_fasta = @"E:\DIA\TopDIA\Ecoli.fasta";
            string gptmd_ecoli = @"E:\DIA\TopDIA\2025-05-10-14-22-53\Task1-GPTMDTask\EcoliGPTMD.xml";
            string YB_rep_DDA_prunedDb = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task4-SearchTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMDpruned.xml";
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder);
            engine.Run();
        }

        [Test]
        public static void TestTopDIA()
        {
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            string outputFolder = @"E:\DIA\TopDIA\SearchResultsMM\test2.0-TD_corr0.5_highestPeakXIC_ms1Tol5ppm_cubic_apexRT0.2_maxMissed2_overlap0.3_FragRank2000_maxRT1.5_300000";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.2,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 1.5, maxRTrangeMS2: 1.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA", 
                apexCycleTolerance:2, scanCycleSplineInterval: 0.05);
            string DIAfile = @"E:\DIA\TopDIA\20231117_DIA_720_800_rep2.mzML";
            string dataBase = @"E:\DIA\TopDIA\Ecoli.fasta";
            var myFileManager = new MyFileManager(true);
            var myMsDataFile = myFileManager.LoadFile(DIAfile, task.CommonParameters);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(dataBase, false) }, new List<string> { DIAfile }, "test");
        }

        [Test]
        public static void TestBottomUpISD()
        {
            var task = new SearchTask();
            task.CommonParameters.TrimMsMsPeaks = false;
            task.CommonParameters.TrimMs1Peaks = false;
            string outputFolder = @"E:\ISD Project\TestSearchBottomUp\ISDEngine_static_12-18-24_bu-ISD100_yeast_RT36-37_5ppm_DeconHighestMS1&PeakMS2_corr0.5_MaxRT0.5-NoSpline_apexRT0.2_overlap0.3";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 1, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.2,
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.05f, minMS1Mass: 0, maxMass: 10000, 
                ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: true, cutMs2Peaks: true, 
                ms1SplineType:SplineType.SavgolSmoothed, ms2SplineType: SplineType.SavgolSmoothed, sgFilterWindowSize: 9);
            string DIAfile = @"E:\ISD Project\ISD_bu\01-14-25_bu_Yeast_SP3_1ug_rep1_ISD100_labeled.mzML";
            string dataBase = @"E:\ISD Project\ISD_250128\2025-02-01-13-52-03\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            var myFileManager = new MyFileManager(true);
            var myMsDataFile = myFileManager.LoadFile(DIAfile, task.CommonParameters);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(dataBase, false) }, new List<string> { DIAfile }, "test");
        }

        [Test]
        public static void TestTopDownISDSearch()
        {
            var filePath1 = @"E:\ISD Project\ISD_250128\01-31-25_td-ISD_PEPPI-YD_105min_ISD60-80-100_micro1_labelCorrected.mzML";
            var filePath2 = @"E:\ISD Project\ISD_250128\01-28-25_td-ISD_PEPPI-YB_105min_ISD60-80-100_60k_micro1_labelCorrected.mzML";
            var filePath3 = @"E:\ISD Project\ISD_250128\01-31-25_td-ISD_PEPPI-YD_105min_gradient5_ISD60-80-100_micro3_labelCorrected.mzML";
            var filePath4 = @"E:\ISD Project\ISD_250128\01-30-25_td-ISD_PEPPI-YB_105min_ISD60-80-100_120k_micro1_labelCorrected.mzML";
            var filePath5 = @"E:\ISD Project\ISD_250128\01-30-25_td-ISD_PEPPI-YB_105min_ISD60-80-100_120k_micro4_2ug_labelCorrected.mzML";
            var filePath6 = @"E:\ISD Project\ISD_250128\01-28-25_td-ISD_PEPPI-YD_105min_ISD60-80-100_120k_micro1_labelCorrected.mzML";
            var filePath7 = @"E:\DIA\FW-DIA data\20230727 RPLC ribosomal protein CV20-70V scan time 5s S1_filtered.mzML";
            var filePath8 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YB_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected.mzML";
            var filePath9 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YB_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected-calib-averaged.mzML";
            var filePath10 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YD_105min_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected.mzML";
            var filePath11 = @"E:\ISD Project\ISD_250428\04-30-25_PEPPI-YD_105min_gradient5_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected.mzML";
            var filePath12 = @"E:\ISD Project\ISD_250428\04-30-25_PEPPI-YD_105min_gradient5_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";
            var filePath13 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";
            var filePath14 = @"E:\ISD Project\ISD_250428\05-02-25_PEPPI-YB_60min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected.mzML";
            var filePath15 = @"E:\ISD Project\ISD_250428\05-03-25_PEPPI-YD_60min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";
            var filePath16 = @"E:\ISD Project\ISD_250428\0428YB_cali-avg-gptmd-xml\Task2-AveragingTask\04-29-25_PEPPI-YB_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected-calib-averaged.mzML";
            var filePath17 = @"E:\ISD Project\ISD_250428\0429YD_DDA&ISD_cali-avg-gptmd-xml\Task2-AveragingTask\04-29-25_PEPPI-YD_105min_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected-calib-averaged.mzML";
            var filePath18 = @"E:\ISD Project\ISD_250428\05-03-25_PEPPI-YD_81min_ISD60-80-100_400-1200_300mz-overlap50_rep1_labelCorrected.mzML";
            var filePath19 = @"E:\ISD Project\ISD_250428\05-03-25_PEPPI-YD_81min_ISD60-80-100_400-1200_300mz-overlap50_rep2_labelCorrected.mzML";
            var filePath20 = @"E:\ISD Project\ISD_250428\05-03-25_PEPPI-YD_81min_ISD60-80-100_400-1200_300mz-overlap50_rep3_labelCorrected.mzML";
            var filePath21 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected.mzML";
            var filePath22 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep2_labelCorrected.mzML";
            var filePath23 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep3_labelCorrected.mzML";
            var filePath24 = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected-calib-averaged.mzML";
            var filePath25 = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep2_labelCorrected-calib-averaged.mzML";
            var filePath26 = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep3_labelCorrected-calib-averaged.mzML";
            var filePath27 = @"E:\ISD Project\ISD_250428\0503YD_rep_ISD&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-03-25_PEPPI-YD_81min_ISD60-80-100_400-1200_300mz-overlap50_rep1_labelCorrected-calib-averaged.mzML";
            var filePath28 = @"E:\ISD Project\ISD_250428\0503YD_rep_ISD&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-03-25_PEPPI-YD_81min_ISD60-80-100_400-1200_300mz-overlap50_rep2_labelCorrected-calib-averaged.mzML";
            var filePath29 = @"E:\ISD Project\ISD_250428\0503YD_rep_ISD&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-03-25_PEPPI-YD_81min_ISD60-80-100_400-1200_300mz-overlap50_rep3_labelCorrected-calib-averaged.mzML";
            var filePath30 = @"E:\ISD Project\ISD_250428\05-03-25_MixS4_60min_ISD100_preFilter1200_labelCorrected.mzML";
            var filePath31 = @"E:\ISD Project\ISD_250428\05-03-25_MixS5_60min_ISD100_preFilter1200_labelCorrected.mzML";
            var filePath32 = @"E:\ISD Project\ISD_250428\05-03-25_MixS6_60min_ISD100_preFilter1200_labelCorrected.mzML";
            var filePath33 = @"E:\ISD Project\ISD_250428\05-03-25_MixS7_60min_ISD100_preFilter1200_labelCorrected.mzML";
            var filePath34 = @"E:\ISD Project\ISD_250428\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected.mzML";
            var filePath35 = @"E:\ISD Project\ISD_250428\05-01-25_PEPPI-YC_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected.mzML";
            var filePath36 = @"E:\ISD Project\ISD_250428\05-01-25_PEPPI-YC_105min_gradient5_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected.mzML";
            var filePath37 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_gradient5_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected-calib-averaged.mzML";
            var filePath38 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected-calib-averaged.mzML";
            var filePath39 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected-calib-averaged.mzML";
            var filePath40 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected.mzML";
            var filePath41 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";
            var filePath42 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_averaged_labelCorrected.mzML";
            var filePath43 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep2_averaged_labelCorrected.mzML";
            var filePath44 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep3_averaged_labelCorrected.mzML";
            var filePath45 = @"E:\ISD Project\ISD_250428\0429YD_DDA&ISD_cali-avg-gptmd-xml\Task2-AveragingTask\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";
            var filePath46 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_gradient5_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected-calib-averaged.mzML";
            var filePath47 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected-calib-averaged.mzML";
            var filePath48 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected-calib-averaged.mzML";
            var filePath49 = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100.mzML";
            var filePath50 = @"E:\ISD Project\ISD_250428\YD_ISD_cali-avg_1.0.8\Task2-AveragingTask\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";
            

            var fileList = new List<string> { filePath24, filePath25, filePath26, };//filePath21, filePath22, filePath23,
            var outputFolder = @"E:\ISD Project\TestSearch\random\b_rep_cali_GPTMD";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_neutralLossSearch = @"E:\CE\250318_CE\YB_seq400-1100_300mz_100overlap_neutralLossSearch\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            string tomlFile_noMods = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_variableMods = @"E:\ISD Project\ISD_250428\variableMods\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly_noInternal = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_gptmd-xml_prunedDb_noInternal\Task Settings\Task2-SearchTaskconfig.toml";

            SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);
            searchTask.CommonParameters.PrecursorMassTolerance = new PpmTolerance(10);
            searchTask.CommonParameters.DIAparameters = new DIAparameters(new PpmToleranceWithNotch(20, 2), new PpmToleranceWithNotch(20, 2),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.75, apexRtTolerance: 0.15,
                fragmentRankCutOff: 150, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 4000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 21, ms1NumPeaksThreshold: 2, ms2NumPeaksThreshold:2, combineFragments: true,
                rankFilter: false, minPFpairCount: 10, sharedXICCutOff: 0.5) ;

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

            var cali_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task1-CalibrateTaskconfig.toml";
            var cali_task = Toml.ReadFile<CalibrationTask>(cali_toml, MetaMorpheusTask.tomlConfig);
            cali_task.CommonParameters = searchTask.CommonParameters.Clone();

            var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", searchTask), }; //("GPTMD", gptmdTask), ("Calibration", cali_task),

            string YC_gptmd = @"E:\CE\250318_CE\YC_cali-avged_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string gptmdDb = @"E:\ISD Project\ISD_250428\0428YB_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var yeast_fasta = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            string YB_gptmd = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            string ecoli_fasta = @"E:\DIA\FW-DIA data\uniprotkb_taxonomy_id_469008_2025_04_24.fasta";
            string YB_rep_DDA_prunedDb = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task4-SearchTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMDpruned.xml";
            string rep_gptmdDb = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DDA_gptmd-xml_prunedDb\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string dda_gptmdDb = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_gptmd-xml_prunedDb\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string all_prunedDb = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DDA_gptmd-xml_prunedDb\Task2-SearchTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMDpruned.xml";
            string D_gptmdDb = @"E:\ISD Project\ISD_250428\0429YD_DDA&ISD_gptmd-xml_noInternal\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder);
            engine.Run();

            var fileList2 = new List<string> { filePath34, filePath35, filePath36, filePath46, filePath47, filePath48 };
            var outputFolder2 = @"E:\ISD Project\TestSearch\0501YC\all+cali-avg_Umpire_GPTMD_ClassicDecon_20ppm_apex0.3_corr0.5_overlap0_maxRT0.5_preSN0.01_num2";
            if (!Directory.Exists(outputFolder2))
            {
                Directory.CreateDirectory(outputFolder2);
            }
            searchTask.CommonParameters.DIAparameters = new DIAparameters(new PpmToleranceWithNotch(50, 2), new PpmToleranceWithNotch(50,2),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, combineFragments: false,
                rankFilter: false, minPFpairCount: 10);

            var engine2 = new EverythingRunnerEngine(taskList, fileList2, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder2);
            engine2.Run();

            var outputFolder3 = @"E:\ISD Project\TestSearch\0501YC\all+cali-avg_Umpire_GPTMD_ClassicDecon_20ppm_apex0.15_corr0.75_overlap0_maxRT0.5_preSN0.01_num2";
            if (!Directory.Exists(outputFolder3))
            {
                Directory.CreateDirectory(outputFolder3);
            }
            searchTask.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.75, apexRtTolerance: 0.15,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, combineFragments: false,
                rankFilter: false, minPFpairCount: 10);

            var engine3 = new EverythingRunnerEngine(taskList, fileList2, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder3);
            engine3.Run();
        }
        //remove minimum number of points required for correlation?

        [Test]
        public static void DDAQuant()
        {
            var psmFile2 = @"E:\ISD Project\ISD_250428\0501YC_DDA_cali-avg-gptmd-xml_1.0.8\Task4-SearchTask\AllProteoforms.psmtsv";
            var psmFile1 = @"E:\ISD Project\ISD_250428\0429YB_DDA_cali-avg-gptmd-xml_1.0.8\Task4-SearchTask\AllProteoforms.psmtsv";
            var psmFile3 = @"E:\ISD Project\ISD_250428\0429YD_DDA_cali-avg-gptmd-xml_1.0.8\Task4-SearchTask\AllProteoforms.psmtsv";
            var psmFile4 = @"E:\ISD Project\ISD_250428\0429YE_DDA_cali-avg-gptmd-xml_1.0.8\Task4-SearchTask\AllProteoforms.psmtsv";
            var psmFile5 = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_gptmd-xml\Task2-SearchTask\Individual File Results\05-04-25_PEPPI-YB_81min_DDA_rep1_Proteoforms.psmtsv";
            var psmFile6 = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_gptmd-xml\Task2-SearchTask\Individual File Results\05-04-25_PEPPI-YB_81min_DDA_rep2_Proteoforms.psmtsv";
            var psmFile7 = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_gptmd-xml\Task2-SearchTask\Individual File Results\05-04-25_PEPPI-YB_81min_DDA_rep3_Proteoforms.psmtsv";
            var psmFile8 = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task4-SearchTask\Individual File Results\05-04-25_PEPPI-YB_81min_DDA_rep1-calib-averaged_Proteoforms.psmtsv";
            var psmFile9 = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task4-SearchTask\Individual File Results\05-04-25_PEPPI-YB_81min_DDA_rep2-calib-averaged_Proteoforms.psmtsv";
            var psmFile10 = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task4-SearchTask\Individual File Results\05-04-25_PEPPI-YB_81min_DDA_rep3-calib-averaged_Proteoforms.psmtsv";

            var filePath2 = @"E:\ISD Project\ISD_250428\0501YC_DDA_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_DDA-calib-averaged.mzML";
            var filePath1 = @"E:\ISD Project\ISD_250428\0429YB_DDA_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\04-29-25_PEPPI-YB_105min_DDA-calib-averaged.mzML";
            var filePath3 = @"E:\ISD Project\ISD_250428\0429YD_DDA_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\04-29-25_PEPPI-YD_105min_DDA-calib-averaged.mzML";
            var filePath4 = @"E:\ISD Project\ISD_250428\0429YE_DDA_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\04-29-25_PEPPI-YE_105min_gradient5_DDA-calib-averaged.mzML";
            var filePath5 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_DDA_rep1.raw";
            var filePath6 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_DDA_rep2.raw";
            var filePath7 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_DDA_rep3.raw";
            var filePath8 = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_DDA_rep1-calib-averaged.mzML";
            var filePath9 = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_DDA_rep2-calib-averaged.mzML";
            var filePath10 = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_DDA_rep3-calib-averaged.mzML";

            var psmFileList = new List<string> {  psmFile8, psmFile9, psmFile10 };
            var fileList = new List<string> { filePath8, filePath9, filePath10 };
            var outputFolder = @"E:\ISD Project\TestSearch\random\DDAQuant_cali-avg"; 
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_neutralLossSearch = @"E:\CE\250318_CE\YB_seq400-1100_300mz_100overlap_neutralLossSearch\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_noMods = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";

            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_noMods, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmToleranceWithNotch(20,2), new PpmToleranceWithNotch(20,2),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.25, apexRtTolerance: 0.5,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 4000, minMS1Charge: 3, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.DDAQuant, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.NoSpline, ms2SplineType: SplineType.NoSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 4, combineFragments: false,
                rankFilter: false, minPFpairCount: 10, numScanPerCycle: 6);

            //Use IsoDec
            //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
            //task.CommonParameters.ProductDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
            //task.CommonParameters.ProductDeconvolutionParameters.MaxAssumedChargeState = 20;

            var combinedResults = new List<DIAProteoformQuant>();
            for (int i = 0; i < fileList.Count; i++)
            {
                var myFileManager = new MyFileManager(true);
                var dataFile = myFileManager.LoadFile(fileList[i], task.CommonParameters);
                var quantFile = DIAProteoformQuantFile.DDAQuantFromExistingSearch(dataFile, psmFileList[i], task.CommonParameters);
                var fileName = Path.GetFileNameWithoutExtension(fileList[i]);
                var resultPath = Path.Combine(outputFolder, fileName + "_DDAQuant.tsv");
                combinedResults.AddRange(quantFile.Results);
                quantFile.WriteResults(resultPath);
            }
            var combinedFilePath = Path.Combine(outputFolder, "Combined_DDAQuant.tsv");
            var combinedQuantFile = new DIAProteoformQuantFile { FilePath = combinedFilePath, Results = combinedResults };
            combinedQuantFile.WriteResults(combinedFilePath);
        }

        [Test]
        public static void SearchCEISD()
        {
            var filePath1 = @"E:\CE\250318_CE\03-22-25_CE_PEPPI-YC_90min_ISD60-80-100_seq400-1350_350mz-50overlap_labelCorrected.mzML";
            var filePath2 = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100.mzML";
            var filePath3 = @"E:\CE\250318_CE\03-21-25_CE_PEPPI-YB_90min_ISD60-80-100_seq400-1100_300mz-100overlap_redo_labelCorrected.mzML";
            var filePath4 = @"E:\CE\250318_CE\YC_cali-xml-avg-xml\Task3-AveragingTask\03-22-25_CE_PEPPI-YC_90min_ISD60-80-100_seq400-1350_350mz-50overlap_labelCorrected-calib-averaged.mzML";
            var filePath5 = @"E:\CE\250318_CE\03-22-25_CE_PEPPI-YC_90min_ISD60-80-100_seq400-1100_300mz-100overlap_labelCorrected.mzML";
            var filePath6 = @"E:\DIA\Data\250313_DIA\03-18-25_CE_5pro_90min_ISD60-80-100_normal_AGC3e6_labelCorrected.mzML";
            var filePath7 = @"E:\DIA\Data\250313_DIA\03-18-25_CE_5pro_90min_ISD60-80-100_equal_labelCorrected.mzML";
            var filePath8 = @"E:\DIA\Data\250313_DIA\03-18-25_CE_5pro_90min_ISD60-80-100_overlap-equal_AGC3e6_labelCorrected.mzML";
            var filePath9 = @"E:\DIA\Data\250313_DIA\03-18-25_CE_5pro_90min_ISD80_4subRange_400-1200_labelCorrected.mzML";
            var filePath10 = @"E:\ISD Project\ISD_250128\01-28-25_td-DIA_PEPPI-YB_105min_50mz_21-23-25HCD_AGC1e6_200ms.raw";
            var filePath11 = @"E:\DIA\Data\250313_DIA\03-17-25_CE_5pro_90min_ISD60-80-100_overlap-equal_labelCorrected.mzML";
            var filePath12 = @"E:\CE\250318_CE\2025-06-18-11-12-21\Task2-AveragingTask\03-21-25_CE_PEPPI-YB_90min_ISD60-80-100_seq400-1100_300mz-100overlap_redo_labelCorrected-averaged.mzML";

            var fileList = new List<string> { filePath3 };
            var outputFolder = @"E:\DIA\TestSearch\YB_seq400-1100_300mz-100overlap_redo\test";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_neutralLossSearch = @"E:\CE\250318_CE\YB_seq400-1100_300mz_100overlap_neutralLossSearch\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmToleranceWithNotch(20, 2), new PpmToleranceWithNotch(20, 2),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.15,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 4000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.ExtendedCycleSpline, ms2SplineType: SplineType.ExtendedCycleSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2,
                ms2NumPeaksThreshold: 2, combineFragments: true);

            //match all charge fragment ions
            task.SearchParameters.WriteSpectralLibrary = true;

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = task.CommonParameters.Clone();
            var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", task), }; //("GPTMD", gptmdTask),

            string YC_gptmd = @"E:\CE\250318_CE\YC_cali-avged_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string yeast_xml = @"E:\ISD Project\ISD_250128\2025-02-01-13-52-03\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder);
            engine.Run();
        }

        [Test]
        public static void SearchLoop()
        {
            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_neutralLossSearch = @"E:\CE\250318_CE\YB_seq400-1100_300mz_100overlap_neutralLossSearch\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            string tomlFile_noMods = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var GPTMD_max2 = @"E:\ISD Project\ISD_250428\GPTMD_max2\Task Settings\Task1-GPTMDTaskconfig.toml";

            string gptmdDb = @"E:\ISD Project\ISD_250428\0428YB_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var yeast_fasta = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            string YB_gptmd = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            string ecoli_fasta = @"E:\DIA\FW-DIA data\uniprotkb_taxonomy_id_469008_2025_04_24.fasta";
            string YB_rep_DDA_prunedDb = @"E:\ISD Project\ISD_250428\0504YB_rep_DDA_cali-avg-gptmd-xml_1.0.8_library-prunedDb\Task4-SearchTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMDpruned.xml";

            var filePath1 = @"E:\ISD Project\ISD_250428\0504YB_ISD_rep_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected-calib-averaged.mzML";
            var filePath2 = @"E:\ISD Project\ISD_250428\0504YB_ISD_rep_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep2_labelCorrected-calib-averaged.mzML";
            var filePath3 = @"E:\ISD Project\ISD_250428\0504YB_ISD_rep_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep3_labelCorrected-calib-averaged.mzML";
            var filePath4 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_ISD60-80-100_400-1100_300mz-overlap100_RF_labelCorrected-calib-averaged.mzML";
            var filePath5 = @"E:\ISD Project\ISD_250428\0501YC_ISD&DIA&DDA_cali-avg-gptmd-xml\Task2-AveragingTask\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected-calib-averaged.mzML";
            
            var outPath1 = @"E:\ISD Project\TestSearch\0501YC\cali-avged-seq_GPTMD_Umpire_20ppm_apex0.3_corr0.75_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10";
            var outPath2 = @"E:\ISD Project\TestSearch\0501YC\cali-avged-preFilter_GPTMD_Umpire_20ppm_apex0.3_corr0.75_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10";
            var outPath3 = @"E:\ISD Project\TestSearch\0501YC\cali-avged-rep3_GPTMD_Umpire_20ppm_apex0.3_corr0.75_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10";
            var outPath4 = @"E:\ISD Project\TestSearch\0501YC\cali-avged-seq_GPTMD_Umpire_20ppm_apex0.3_corr0.75_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_combineFrags";
            var outPath5 = @"E:\ISD Project\TestSearch\0501YC\cali-avged-preFilter_GPTMD_Umpire_20ppm_apex0.3_corr0.75_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_combineFrags";
            
            var fileList = new List<string> {  };
            var outFolderList = new List<string> { outPath4, outPath5};

            for(int i = 0; i < fileList.Count; i++)
            {
                var searchFileList = new List<string> { fileList[i] };
                if (!Directory.Exists(outFolderList[i]))
                {
                    Directory.CreateDirectory(outFolderList[i]);
                }
                SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                    maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.5,
                    fragmentRankCutOff: 150, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                    precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                    apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 4000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
            ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                    pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                    ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, ms2NumPeaksThreshold: 2, combineFragments: true,
                    rankFilter: false, minPFpairCount: 10, sharedXICCutOff: 0.5);

                //match all charge fragment ions
                searchTask.SearchParameters.WriteSpectralLibrary = true;

                //Use IsoDec
                //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
                //task.CommonParameters.ProductDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
                //task.CommonParameters.ProductDeconvolutionParameters.MaxAssumedChargeState = 20;

                var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters.Clone();

                var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", searchTask), }; //("GPTMD", gptmdTask), ("Calibration", cali_task),

                //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

                var engine = new EverythingRunnerEngine(taskList, searchFileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outFolderList[i]);
                engine.Run();

            }

            var filePath6 = @"E:\ISD Project\ISD_250428\0429YD_ISD_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\04-29-25_PEPPI-YD_105min_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected-calib-averaged.mzML";
            var filePath7 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YD_105min_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected.mzML";
            var filePath8 = @"E:\ISD Project\ISD_250428\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected.mzML";
            var filePath9 = @"E:\ISD Project\ISD_250428\allYE_cali-avg-gptmd-xml\Task2-AveragingTask\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";
            var filePath10 = @"E:\ISD Project\ISD_250428\allYE_cali-avg-gptmd-xml\Task2-AveragingTask\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_400-1200_300mz-overlap50_RF_labelCorrected-calib-averaged.mzML";
            var filePath11 = @"E:\ISD Project\ISD_250428\0429YD_ISD_cali-avg-gptmd-xml_1.0.8\Task2-AveragingTask\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";

            var outPath6 = @"E:\ISD Project\TestSearch\0429_YD_seq\cali-avged_GPTMD_Umpire_20ppm_apex0.3_corr0.5_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_new";
            var outPath7 = @"E:\ISD Project\TestSearch\0429_YD_seq\GPTMD_Umpire_20ppm_apex0.3_corr0.5_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_new_combineFrags";
            var outPath8 = @"E:\ISD Project\TestSearch\0429YE\seq_GPTMD_Umpire_20ppm_apex0.3_corr0.75_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_new";
            var outPath9 = @"E:\ISD Project\TestSearch\0429YE\cali-avged-preFilter_GPTMD_Umpire_20ppm_apex0.3_corr0.75_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_new";
            var outPath10 = @"E:\ISD Project\TestSearch\0429YE\cali-avged-seq_GPTMD_sharedXIC0.25_apex0.3";
            var outPath11 = @"E:\ISD Project\TestSearch\0429_YD_preFilter\cali-avged-GPTMD_sharedXIC0.25_apex0.3";

            var fileList2 = new List<string> { filePath11, filePath10};
            var outFolderList2 = new List<string> { outPath11, outPath10};
            for (int i = 0; i < fileList2.Count; i++)
            {
                var searchFileList = new List<string> { fileList2[i] };
                if (!Directory.Exists(outFolderList2[i]))
                {
                    Directory.CreateDirectory(outFolderList2[i]);
                }
                SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                    maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                    fragmentRankCutOff: 150, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                    precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                    apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 8000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
            ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                    pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                    ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, ms2NumPeaksThreshold: 2, combineFragments: false,
                    rankFilter: false, minPFpairCount: 10, sharedXICCutOff: 0.25);

                //match all charge fragment ions
                searchTask.SearchParameters.WriteSpectralLibrary = true;

                //Use IsoDec
                //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
                //task.CommonParameters.ProductDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
                //task.CommonParameters.ProductDeconvolutionParameters.MaxAssumedChargeState = 20;

                var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters.Clone();

                var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", searchTask), }; //("GPTMD", gptmdTask), ("Calibration", cali_task),

                //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

                var engine = new EverythingRunnerEngine(taskList, searchFileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outFolderList2[i]);
                engine.Run();
            }

            //var outPath10 = @"E:\ISD Project\TestSearch\0429_YD_preFilter\cali-avged_GPTMD_Umpire_20ppm_apex0.3_corr0.5_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_new";
            //var outPath11 = @"E:\ISD Project\TestSearch\0429_YD_seq\cali-avged_GPTMD_Umpire_20ppm_apex0.3_corr0.5_overlap0_maxRT0.5_preIntensity0.01_num2_minCount10_new_combineFrags";
            var outPath12 = @"E:\ISD Project\TestSearch\0429YE\cali-avged-seq_GPTMD_Umpire_20ppm_apex0.3_corr0.5_overlap0_maxRT0.5_preIntensity0.01_num2_minCount10_new_combineFrags";
            var outPath13 = @"E:\ISD Project\TestSearch\0429YE\cali-avged-preFilter_GPTMD_Umpire_20ppm_apex0.3_corr0.5_overlap0.2_maxRT0.5_preIntensity0.01_num2_minCount10_new";

            var outFolderList3 = new List<string> {  outPath11, outPath12 };
            for (int i = 0; i < fileList2.Count; i++)
            {
                var searchFileList = new List<string> { fileList2[i] };
                if (!Directory.Exists(outFolderList3[i]))
                {
                    Directory.CreateDirectory(outFolderList3[i]);
                }
                SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                    maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                    fragmentRankCutOff: 150, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                    precursorIntensityCutOff: 0.01, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                    apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 8000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
            ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                    pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                    ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, ms2NumPeaksThreshold: 2, combineFragments: false,
                    rankFilter: false, minPFpairCount: 10, sharedXICCutOff: 0.5);

                //match all charge fragment ions
                searchTask.SearchParameters.WriteSpectralLibrary = true;

                //Use IsoDec
                //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
                //task.CommonParameters.ProductDeconvolutionParameters = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);
                //task.CommonParameters.ProductDeconvolutionParameters.MaxAssumedChargeState = 20;

                var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters.Clone();

                var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", searchTask), }; //("GPTMD", gptmdTask), ("Calibration", cali_task),

                //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

                var engine = new EverythingRunnerEngine(taskList, searchFileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outFolderList3[i]);
                engine.Run();
            }
        }

        [Test]
        public static void ASMS_individual()
        {
            var filePathB = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YB_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected.mzML";
            var filePathC = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected.mzML";
            var filePathD = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";
            var filePathE = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";

            var fileB_cali_avg = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YB_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected-calib-averaged.mzML";
            var fileC_cali_avg = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected-calib-averaged.mzML";
            var fileD_cali_avg = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";
            var fileE_cali_avg = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";

            var rep1 = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected.mzML";
            var rep2 = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep2_labelCorrected.mzML";
            var rep3 = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep3_labelCorrected.mzML";

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = searchTask.CommonParameters.Clone();

            var fileList = new List<string> { fileB_cali_avg, fileC_cali_avg, fileD_cali_avg, fileE_cali_avg };
            var outFolder = @"E:\ISD Project\TestSearch\ASMS\Individual_cali";
            if (!Directory.Exists(outFolder))
            {
                Directory.CreateDirectory(outFolder);
            }

            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var yeast_fasta = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";

            var fileList2 = new List<string> { filePathB, filePathC, filePathD, filePathE };
            var outFolder2 = @"E:\ISD Project\TestSearch\ASMS\ind_nocali";
            if (!Directory.Exists(outFolder))
            {
                Directory.CreateDirectory(outFolder);
            }
            var corrList = new List<double> { 0.25, 0.5, 0.75 };
            var apexList = new List<double> { 0.15, 0.3, 0.5 };

            foreach(double corr in corrList)
            {
                foreach(double apex in apexList)
                {
                    foreach (var file in fileList2)
                    {
                        var searchList = new List<string> { file };
                        SearchTask searchTask1 = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);
                        var resultFolder = Path.Combine(outFolder2, Path.GetFileNameWithoutExtension(file), $"{corr}, {apex}");
                        searchTask1.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                       maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: corr, apexRtTolerance: apex,
                       fragmentRankCutOff: 150, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                       precursorIntensityCutOff: 0.01, splineTimeInterval: 0.005f, type: "DIA", scanCycleSplineInterval: 0.05, minMS1Mass: 4000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
                       ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                       pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                       ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, ms2NumPeaksThreshold: 2, combineFragments: false,
                       rankFilter: false, minPFpairCount: 10, sharedXICCutOff: 0.5);
                        if (!Directory.Exists(resultFolder))
                        {
                            Directory.CreateDirectory(resultFolder);
                        }
                        gptmdTask.CommonParameters = searchTask1.CommonParameters;
                        var taskList = new List<(string, MetaMorpheusTask)> {  ("search", searchTask1) }; //("GPTMD", gptmdTask), ("Calibration", cali_task),
                        var engine = new EverythingRunnerEngine(taskList, searchList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, resultFolder);
                        engine.Run();
                    }
                }
            }
            
        }

        [Test]
        public static void ASMS_all()
        {
            var filePathB = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YB_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected.mzML";
            var filePathC = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected.mzML";
            var filePathD = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";
            var filePathE = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected.mzML";

            var fileC_cali_avg = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-01-25_PEPPI-YC_105min_ISD60-80-100_preFilter700-900-1100_RF_labelCorrected-calib-averaged.mzML";
            var fileD_cali_avg = @"E:\ISD Project\ISD_250428\YD_ISD_cali-avg_1.0.8\Task2-AveragingTask\04-29-25_PEPPI-YD_105min_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";
            var fileE_cali_avg = @"E:\ISD Project\Presentations\ASMS 2025 Data\04-29-25_PEPPI-YE_105min_gradient5_ISD60-80-100_preFilter800-1000-1200_RF_labelCorrected-calib-averaged.mzML";

            var rep1 = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1_labelCorrected.mzML";
            var rep2 = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep2_labelCorrected.mzML";
            var rep3 = @"E:\ISD Project\Presentations\ASMS 2025 Data\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep3_labelCorrected.mzML";

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask searchTask = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = searchTask.CommonParameters.Clone();

            var fileList = new List<string> { filePathB, filePathC, filePathD, filePathE };
            var outFolder = @"E:\ISD Project\TestSearch\ASMS\all_nocali\0.75,0.15";
            if (!Directory.Exists(outFolder))
            {
                Directory.CreateDirectory(outFolder);
            }

            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var yeast_fasta = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";

            searchTask.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
               maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.2, correlationCutOff: 0.75, apexRtTolerance: 0.15,
               fragmentRankCutOff: 150, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
               precursorIntensityCutOff: 0.01, splineTimeInterval: 0.005f, type: "DIA", scanCycleSplineInterval: 0.05, minMS1Mass: 4000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
               ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
               pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
               ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, ms2NumPeaksThreshold: 2, combineFragments: false,
               rankFilter: false, minPFpairCount: 10, sharedXICCutOff: 0.5);

            var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", searchTask), }; //("GPTMD", gptmdTask), ("Calibration", cali_task),
            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outFolder);
            engine.Run();
        }
    }
}
