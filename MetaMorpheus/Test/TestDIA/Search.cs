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
            string outputFolder = @"E:\DIA\TestSearch\TopDIA\Umpire_GPTMD_20ppm_corr0.5_apexRT0.3_overlap0_num2";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string DIAfile1 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep1.raw";
            string DIAfile2 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep2.raw";
            string DIAfile3 = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_td-DIA_700-1000_25mz-2overlap_21-23-25HCD_rep3.raw";
            string topDIAfile = @"E:\DIA\TopDIA\20231117_DIA_720_800_rep2.mzML";
            var fileList = new List<string> { topDIAfile };

            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 30000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.DIAEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.SimpleGaussian, ms2SplineType: SplineType.SimpleGaussian, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2, combineFragments: false,
                rankFilter: false, minPFpairCount: 1);

            //match all charge fragment ions
            task.SearchParameters.WriteSpectralLibrary = true;

            //Use IsoDec
            var isoDecDeconParamMS1 = new IsoDecDeconvolutionParameters();
            var isoDecDeconParamMS2 = new IsoDecDeconvolutionParameters();
            isoDecDeconParamMS2.MaxAssumedChargeState = 20;
            task.CommonParameters.PrecursorDeconvolutionParameters = isoDecDeconParamMS1;
            task.CommonParameters.ProductDeconvolutionParameters = isoDecDeconParamMS2;

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = task.CommonParameters.Clone();
            var taskList = new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask), ("search", task) };//("GPTMD", gptmdTask),

            string YC_gptmd = @"E:\CE\250318_CE\YC_cali-avged_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string gptmdDb = @"E:\ISD Project\ISD_250428\0428YB_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            string ecoli_fasta = @"E:\DIA\TopDIA\Ecoli.fasta";
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(ecoli_fasta, false) }, outputFolder);
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

            var fileList = new List<string> { filePath33 };
            var outputFolder = @"E:\ISD Project\TestSearch\Mix\debug_MixS7_NoSpline_20ppm_apex0.15_corr0.75_overlap0_maxRT0.5_preIntensity1e5_num4_minCount10";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_neutralLossSearch = @"E:\CE\250318_CE\YB_seq400-1100_300mz_100overlap_neutralLossSearch\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            string tomlFile_noMods = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";

            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_noMods, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.75, apexRtTolerance: 0.15,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 100000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.NoSpline, ms2SplineType: SplineType.NoSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 4, combineFragments: false, 
                rankFilter:false, minPFpairCount: 10);

            //match all charge fragment ions
            task.SearchParameters.WriteSpectralLibrary = true;

            //Use IsoDec
            //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters();
            //task.CommonParameters.ProductDeconvolutionParameters = new IsoDecDeconvolutionParameters();
            //task.CommonParameters.ProductDeconvolutionParameters.MaxAssumedChargeState = 20;

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var GPTMD_max2 = @"E:\ISD Project\ISD_250428\GPTMD_max2\Task Settings\Task1-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = task.CommonParameters.Clone();
            var taskList = new List<(string, MetaMorpheusTask)> {  ("search", task), }; //("GPTMD", gptmdTask), 

            string YC_gptmd = @"E:\CE\250318_CE\YC_cali-avged_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string gptmdDb = @"E:\ISD Project\ISD_250428\0428YB_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string yeast_xml = @"E:\ISD Project\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            string YB_gptmd = @"E:\ISD Project\ISD_250428\0504YB_rep_ISD&DIA&DDA_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            string ecoli_fasta = @"E:\DIA\FW-DIA data\uniprotkb_taxonomy_id_469008_2025_04_24.fasta";
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(standard_xml, false) }, outputFolder);
            engine.Run();

            var outputFolder2 = @"E:\ISD Project\TestSearch\0504YB_rep_preFilter\rep123_CubicSpline_GPTMD_20ppm_apex0.15_corr0.5_overlap0.5_maxRT0.5_preIntensity1e5_num4_minCount10";
            if (!Directory.Exists(outputFolder2))
            {
                Directory.CreateDirectory(outputFolder2);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.15,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 100000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.CubicSpline, ms2SplineType: SplineType.CubicSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 4, combineFragments: false,
                rankFilter: false, minPFpairCount: 10);

            var engine2 = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder2);
            engine2.Run();

            var outputFolder3 = @"E:\ISD Project\TestSearch\0504YB_rep_preFilter\rep123_CubicSpline_GPTMD_20ppm_apex0.15_corr0.75_overlap0.5_maxRT0.5_preIntensity1e5_num4_minCount10";
            if (!Directory.Exists(outputFolder3))
            {
                Directory.CreateDirectory(outputFolder3);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0, correlationCutOff: 0.75, apexRtTolerance: 0.15,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 100000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.CubicSpline, ms2SplineType: SplineType.CubicSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 4, combineFragments: false,
                rankFilter: false, minPFpairCount: 10);

            var engine3 = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder3);
            engine3.Run();
        }
        //remove minimum number of points required for correlation?

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

            var fileList = new List<string> { filePath1 };
            var outputFolder = @"E:\DIA\TestSearch\test";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_neutralLossSearch = @"E:\CE\250318_CE\YB_seq400-1100_300mz_100overlap_neutralLossSearch\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.2, correlationCutOff: 0.25, apexRtTolerance: 0.2,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 100000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 3000, minMS1Charge: 4, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.NoSpline, ms2SplineType: SplineType.NoSpline, sgFilterWindowSize: 7, ms1NumPeaksThreshold: 2);

            //match all charge fragment ions
            task.SearchParameters.WriteSpectralLibrary = true;

            var lessGPTMD_toml = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task3-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(lessGPTMD_toml, MetaMorpheusTask.tomlConfig);
            gptmdTask.CommonParameters = task.CommonParameters.Clone();
            var taskList = new List<(string, MetaMorpheusTask)> {  ("search", task), };

            string YC_gptmd = @"E:\CE\250318_CE\YC_cali-avged_gptmd-xml\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string yeast_xml = @"E:\ISD Project\ISD_250128\2025-02-01-13-52-03\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(YC_gptmd, false) }, fileList, "test");

            var engine = new EverythingRunnerEngine(taskList, fileList, new List<DbForTask> { new DbForTask(yeast_xml, false) }, outputFolder);
            engine.Run();
        }
    }
}
