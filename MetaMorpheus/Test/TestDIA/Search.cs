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
            string outputFolder = @"E:\DIA\TestSearch\scanBased2_pfGroup_corr0.5_highestPeakXIC_ms1Tol5ppm_cubicSpline0.05_apexRT0.25_noPeakTrim_maxMissed1_overlap0.2_Frank150_Prank10_maxRT0.5";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20), 1, 100, 0.2, 0.5, 0.25, 150, 10, maxRTrange: 0.5);
            //string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var fm = new MyFileManager(false);
            var dataFile = fm.LoadFile(DIAfile, task.CommonParameters);
            var ms1scans = dataFile.GetMS1Scans().ToList();
            string myDatabase = @"E:\ISD Project\Claire's human data\Human_9606.fasta";
            string myDatabase_xml = @"E:\ISD Project\ISD_240812\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { DIAfile }, "test");
        }

        [Test]
        public static void TestISDSearchStandardMix()
        {
            string filePath1 = @"E:\ISD Project\TestIsdDataAnalysis\data\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string filePath2 = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD_RT32.52-35.65.mzML";

            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            //var file = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD.mzML";

            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\isdScanBased_corr0.5_highestPeakXIC_ms1Tol10ppm_apexRT0.3_maxMissed2_overlap0.5_maxRT2_300000_averageMs2_FragIntFilter";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.5, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrange:2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, minMass: 6000);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { filePath2 }, "test");
        }

        [Test]
        public static void TestISDyeast()
        {
            string file = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT31.39-36.17.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240812\FB_FD_NoGPTMD\Task Settings\Task3-SearchTaskconfig.toml";
            string dataBase = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task3-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            string outputFolder = @"E:\ISD Project\TestSearch\isdScanBased_FractionD-GPTMDdb_orbiMS1_ISD60-80-100_RT31.39-36.17_corr0.5_highestPeakXIC_ms1Tol10ppm_apexRT0.3_maxMissed2_overlap0.3_maxRT0.5_100000_averageMS1&2_FragIntFilter";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrange: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 100000, minMass: 12000);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(dataBase, false) }, new List<string> { file }, "test");
        }
    }
}
