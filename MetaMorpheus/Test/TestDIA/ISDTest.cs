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
using MzLibUtil;
using Plotly.NET;
using UsefulProteomicsDatabases;
using Nett;
using System.IO.Compression;

namespace Test.TestDIA
{
    public class ISDTest
    {
        [Test]
        public static void TestConstructMs2Groups()
        {
            string filePath1 = @"E:\ISD Project\TestIsdDataAnalysis\data\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string filePath2 = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task2-AveragingTask\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100-calib-averaged_RT28.83-32.69.mzML";

            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            //var file = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD.mzML";

            string myDatabase = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task3-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\isdEngine_cali-avg-FD-RT28.83-32.69_corr0.3-cubic_highestPeakXIC_ms1Tol10ppm_apexRT0.3_maxMissed2_overlap0.5_maxRT1_300000";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.5, correlationCutOff: 0.3, apexRtTolerance: 0.3,
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD");
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { filePath2 }, "test");
        }
    }
}
