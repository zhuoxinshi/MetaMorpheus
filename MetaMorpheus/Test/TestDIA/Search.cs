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
                splineRtInterval: 0.005, numPeaksThreshold: 2);
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
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_CommonFixedVariable, MetaMorpheusTask.tomlConfig);
            string outputFolder = @"E:\DIA\TestSearch\0128-YB-50mz_DIAEngine_static_20ppm_Umpire_massCurve_corr0_apexRT0.3_maxMissed2_overlap0.3_maxRT1_num2_combinePfGroup_f200_p20";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.3, correlationCutOff: 0, apexRtTolerance: 0.3,
                fragmentRankCutOff:200, precursorRankCutOff: 20, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA", minMS1Mass: 5000, minMS1Charge: 4,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, cutMs1Peaks: false, cutMs2Peaks: false, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.massCurve, analysisType: AnalysisType.DIAEngine_static, ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline,
                splineRtInterval: 0.005, numPeaksThreshold: 2);
            string DIAfile1 = @"E:\ISD Project\ISD_250128\01-28-25_td-DIA_PEPPI-YB_105min_50mz_21-23-25HCD_AGC1e6_200ms.raw"; 
            string DIAfile2 = @"E:\DIA\Data\250313_DIA\03-14-25_td-DIA_PEPPI-YC_575-1300_25mz_nce22-25-27_AGC1e6_500ms_micro1_charge10.raw";
            string myDatabase = @"E:\ISD Project\ISD_240812\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { DIAfile1 }, "test");
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
            var filePath1 = @"E:\DIA\Data\250313_DIA\03-16-25_PEPPI-YD_ISD40-60-80-100_normal_120k_250ms_micro1_labelCorrected.mzML";
            var filePath2 = @"E:\ISD Project\TestIsdDataAnalysis\02-01-25_td-ISD_PEPPI-YC_105min_ISD60-80-100_micro1_labelCorrected-calib-averaged.mzML";
            string tomlFile_noFixedMods = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_searchOnly = @"E:\ISD Project\CE_241118\1122_DDA&ISD-rep123\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_CommonFixedVariable = @"E:\CE\250318_CE\0322_YC_SearchOnly\Task Settings\Task1-SearchTaskconfig.toml";
            string tomlFile_FixedOnly = @"E:\ISD Project\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile_FixedOnly, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.MaxThreadsToUsePerFile = 10;
            //task.CommonParameters.PrecursorDeconvolutionParameters = new IsoDecDeconvolutionParameters();
            var outputFolder = @"E:\ISD Project\TestSearch\massCurve-03-16-25_PEPPI-YD_ISD40-60-80-100_normal_20ppmMS1MS2_Umpire_corr0.5_apexRT0.3_overlap0.3_num2";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(20), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 1, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA",
                apexCycleTolerance: 2, scanCycleSplineInterval: 0.05, minMS1Mass: 5000, minMS1Charge: 5, minMS2Charge: 1, minMS2Mass: 0, splineRtInterval: 0.005,
        ms1XICType: XICType.MassCurve, ms2XICType: XICType.MassCurve, pfGroupingType: PFGroupingType.Umpire,
                pseudoMs2Type: PseudoMs2ConstructionType.massCurve, analysisType: AnalysisType.ISDEngine_static, cutMs1Peaks: false, cutMs2Peaks: false,
                ms1SplineType: SplineType.UmpireBSpline, ms2SplineType: SplineType.UmpireBSpline, sgFilterWindowSize: 7, numPeaksThreshold: 2);
            string yeast_xml = @"E:\ISD Project\ISD_250128\2025-02-01-13-52-03\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string standard_xml = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(yeast_xml, false) }, new List<string> { filePath1 }, "test");
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

            var fileList = new List<string> { filePath10 };
            var outputFolder = @"E:\DIA\TestSearch\DIAEngine_static_YB_50mz_20ppmMS1MS2_NoSpline_corr0.25_apexRT0.2_overlap0.2_maxRT1_num2_preIntensity1e5";
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
                ms1SplineType: SplineType.NoSpline, ms2SplineType: SplineType.NoSpline, sgFilterWindowSize: 7, numPeaksThreshold: 2);

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
