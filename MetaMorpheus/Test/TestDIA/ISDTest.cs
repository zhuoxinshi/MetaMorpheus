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
using Chemistry;

namespace Test.TestDIA
{
    public class ISDTest
    {
        [Test]
        public static void TestConstructMs2Groups()
        {
            string filePath1 = @"E:\ISD Project\TestIsdDataAnalysis\data\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string filePath2 = @"E:\ISD Project\Data for test\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100-calib-averaged.mzML";
            string snip = @"E:\ISD Project\ISD_240606\06-09-24_mix_sample10_5uL_ISD_shorter-gradient_RT27.6-32.44.mzML";

            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            //var file = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD.mzML";

            string myDatabase = @"E:\ISD Project\ISD_240927\60only-RT29-33\Task1-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\isdEngine_FD60only-RT29.03-33.52_corr0.5_highestPeakXIC_ms1Tol10ppm_apexRT0.3_maxMissed2_overlap0.3";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.05, cutPeaks: true);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { snip }, "test");
        }

        [Test]
        public static void TestEnvelopeXIC()
        {
            string snip = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT28.01-32.69.mzML";
            var myFileManagers = new MyFileManager(true);
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD");
            var myMsDataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            var isdEngine = new ISDEngine(myMsDataFile, task.CommonParameters, task.CommonParameters.DIAparameters);
            isdEngine.Ms1PeakIndexing();
            isdEngine.ConstructMs2Group();
            isdEngine.GetMs1PeakEnvelopeCurve();  
            //isdEngine.GetMs1PeakCurves();
            //isdEngine.GetPrecursorEnvelopeXIC2();
            //isdEngine.GetEnvelopeXICForMs2();
            //isdEngine.PrecursorFragmentPairing();
        }

        [Test]
        public static void TestDecon()
        {
            string file = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100.mzML";
            var myFileManagers = new MyFileManager(true);
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myMsDataFile = myFileManagers.LoadFile(file, task.CommonParameters);
            //var scan = myMsDataFile.GetOneBasedScan(1330);
            //var envelopes = Deconvoluter.Deconvolute(scan, task.CommonParameters.PrecursorDeconvolutionParameters).OrderBy(e => e.MostAbundantObservedIsotopicMass.ToMz(e.Charge)).ToList();

        }
    }
}
