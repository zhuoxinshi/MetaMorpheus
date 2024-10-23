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
            string filePath2 = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD_RT32.52-35.65.mzML";

            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            //var file = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD.mzML";

            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\isdEngine_corr0.5_highestPeakXIC_ms1Tol5ppm_apexRT0.5_maxMissed2_overlap0.3_maxRT3_300000";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(10),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.5, correlationCutOff: 0.5, apexRtTolerance: 0.5,
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrangeMS1: 3, maxRTrangeMS2: 3, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 8000, type: "ISD");
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { filePath2 }, "test");
        }
    }
}
