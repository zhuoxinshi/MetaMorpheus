using NUnit.Framework;
using EngineLayer;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using MzLibUtil;
using Newtonsoft.Json.Linq;
using Proteomics.ProteolyticDigestion;
using Chemistry;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;
using System.IO;
using System.Text.RegularExpressions;
using System.Collections.Concurrent;
using EngineLayer.ISD;
using MathNet.Numerics.Distributions;
using Omics.Fragmentation;


namespace Test.TestISD
{
    public class TestSearchISD
    {
        [Test]
        public static void TestIsdSearchSnip()
        {
            string filePath1 = @"E:\ISD Project\TestIsdDataAnalysis\data\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string filePath2 = @"E:\ISD Project\TestIsdDataAnalysis\data\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            //task.SearchParameters.WriteSpectralLibrary = true;
            var myMsDataFile = myFileManager.LoadFile(filePath1, task.CommonParameters);
            //var testScans = MetaMorpheusTask._GetMs2Scans_DeconvoluteMs2(myMsDataFile, filePath, task.CommonParameters);
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\Search results\1pmol_5uL_ISD_RT32.16-35.59_DeconvoluteMs2_corr0.5_no-toleranceForXIC";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            //DbForTask lib = new DbForTask(library, false);
            task.RunTask(outputFolder, new List<DbForTask> { db }, new List<string> { filePath2 }, "normal");

        }

        [Test]
        public static void TestIsdSearchWhole()
        {
            string filePath1 = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";
            string filePath2 = @"E:\ISD Project\ISD_240606\06-11-24_mix_sample10_2uL_ISD.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            task.SearchParameters.WriteSpectralLibrary = true;
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\Search results\1pmol_5uL_ISD_FilterPrecursors_NormalMS2GroupingPeaks_corr0.9_allXIC_noTolerance_allprecursors_noDeconMs2";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            //DbForTask lib = new DbForTask(library, false);
            task.RunTask(outputFolder, new List<DbForTask> { db }, new List<string> { filePath1 }, "normal");

        }

        [Test]
        public static void TestIsdSearchSnip_LFQmethod()
        {
            string filePath = @"E:\ISD Project\TestIsdDataAnalysis\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT45.01-48.09_XIC100_LFQmethod_corr0.9_NoMassFilterDecon";
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            DbForTask lib = new DbForTask(library, false);
            task.RunTask(outputFolder, new List<DbForTask> { db, lib }, new List<string> { filePath }, "normal");

        }

        [Test]
        public static void TestAllMethods()
        {
            string filePath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\Search results\06-07-24_mix_1pmol_5uL_ISD_Ms1peakGroup";
            task.RunTask(outputFolder, new List<DbForTask> { db}, new List<string> { filePath }, "normal");
        }
    }
}
