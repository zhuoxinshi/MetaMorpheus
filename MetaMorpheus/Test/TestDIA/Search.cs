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

namespace Test.TestDIA
{
    public class Search
    {
        [Test]
        public static void TestDIASearch()
        {
            var task = new SearchTask();
            task.CommonParameters.TrimMsMsPeaks = false;
            task.CommonParameters.TrimMs1Peaks = false;
            string outputFolder = @"E:\DIA\TestSearch\DIA_scanBased_corr0.5_highestPeakXIC_ms1Tol5ppmMs2Tol10ppm_RetentionTime_apexRT0.25_overlap0.2_corr0.5_maxRT0.5_Frank100_Prank10_miss1";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(10), 
                maxNumMissedScan:1, overlapRatioCutOff: 0.2, correlationCutOff: 0.5, apexRtTolerance:0.25, 
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0, 
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 0, maxMass: 10000, type: "DIA", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.DIA_scanBased, correlationType: CorrelationType.CubicSpline_RT);
            // Use reflection to set max threads
            //task.CommonParameters.GetType().GetProperty("MaxThreadsToUsePerFile").SetMethod.Invoke(task.CommonParameters, new object[] { 1 });
            //var type = task.CommonParameters.GetType();
            //var property = type.GetProperty("MaxThreadsToUsePerFile");
            //property.SetMethod.Invoke(task.CommonParameters, new object[] { 1 });

            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string myDatabase = @"E:\ISD Project\Claire's human data\Human_9606.fasta";
            string myDatabase_xml = @"E:\ISD Project\ISD_240812\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { DIAfile }, "test");
        }

        [Test]
        public static void TestTopDown()
        {
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.TrimMsMsPeaks = false;
            task.CommonParameters.TrimMs1Peaks = false;
            string outputFolder = @"E:\DIA\TestSearch\test2.0-TD-5pro_corr0.5_highestPeakXIC_ms1Tol5ppm_cubic_apexRT0.2_maxMissed2_overlap0.3_FragRank2000_maxRT1.5_300000";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.2,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 1.5, maxRTrangeMS2: 1.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, type: "DIA", minMass: 8000, minCharge: 4);
            string DIAfile = @"E:\DIA\Data\DIA_241108\11-08-24_td-DIA_5pro-sample10_105min_30mz_25stpHCD_AGC400_500ms_micro1_RT42.76-45.9.mzML";
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { DIAfile }, "test");
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

    }
}
