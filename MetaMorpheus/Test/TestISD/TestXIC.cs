using NUnit.Framework;
using EngineLayer;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using MzLibUtil;
using Newtonsoft.Json.Linq;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;
using Nett;
using System.IO;
using static System.Net.WebRequestMethods;
using System.Text.RegularExpressions;
using System.Collections.Concurrent;
using FlashLFQ;
using EngineLayer.ISD;

namespace Test.TestISD
{
    public static class TestXIC
    {
        [Test]
        public static void TestPeakCurveOnISD()
        {
            string filePath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";
            MyFileManager myFileManager = new MyFileManager(true);
            var digestionParam = new DigestionParams(protease: "top-down");
            CommonParameters isdCommonParameters = new CommonParameters(digestionParams: digestionParam);
            var myMsDataFile = myFileManager.LoadFile(filePath, isdCommonParameters);
            var allMs1Scans = myMsDataFile.Scans.Where(s => s.MsnOrder == 1).ToList();
            var testMs1 = allMs1Scans.Where(p => p.OneBasedScanNumber == 1085).First();
            var filteredPeaksIndices = testMs1.MassSpectrum.YArray.Select((value, index) => new { Index = index, Value = value })
                .Where(y => y.Value / testMs1.MassSpectrum.SumOfAllY > 0.01)
                .Select(y => y.Index);
            var filteredMzs = testMs1.MassSpectrum.XArray.Select((value, index) => new { Value = value, Index = index })
                .Where(x => filteredPeaksIndices.Contains(x.Index))
                .Select(x => x.Value);
            var testPeakCurves = new List<PeakCurve>();
            var XICs = new List<(double[] RT, double[] intensity)>();
            foreach (double mz in filteredMzs)
            {
                PeakCurve p = PeakCurve.GetPeakCurve(allMs1Scans.ToArray(), testMs1, mz, isdCommonParameters);
                testPeakCurves.Add(p);
                double[] rt = p.Peaks.Select(p => p.RT).ToArray();
                double[] intensities = p.Peaks.Select(p => p.Intensity).ToArray();
                XICs.Add((rt, intensities));
            }
        }

        [Test]
        public static void TestIsdSearchSnip()
        {
            string filePath = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT45.67-51.56.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\1pmol_ISD_RT45.67-51.56_XIC100";
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            DbForTask db = new DbForTask(myDatabase, false);
            task.RunTask(outputFolder, new List<DbForTask> { db }, new List<string> { filePath }, "normal");

        }

        [Test]
        public static void TestIsdSearchWhole()
        {
            string filePath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\Whole";
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            DbForTask db = new DbForTask(myDatabase, false);
            task.RunTask(outputFolder, new List<DbForTask> { db }, new List<string> { filePath }, "normal");

        }

        [Test]
        public static void TestConvertingAllMs2Peaks()
        {
            string filePath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            MyFileManager myFileManager = new MyFileManager(true);
            var digestionParam = new DigestionParams(protease: "top-down");
            CommonParameters isdCommonParameters = new CommonParameters(digestionParams: digestionParam, trimMsMsPeaks: false);
            var myMsDataFile = myFileManager.LoadFile(filePath, isdCommonParameters);
            var allMs2Scans = myMsDataFile.Scans.Where(s => s.MsnOrder == 2).ToList();
            var allPeaks = new List<Peak>();
            int index = 0;
            foreach (var scan in allMs2Scans)
            {
                var spectrum = scan.MassSpectrum;
                for (int i = 0; i < spectrum.XArray.Length; i++)
                {
                    Peak newPeak = new Peak(spectrum.XArray[i], scan.RetentionTime, spectrum.YArray[i], scan.MsnOrder, scan.OneBasedScanNumber, index);
                    allPeaks.Add(newPeak);
                    index++;
                }
            }

            Parallel.ForEach(Partitioner.Create(0, allPeaks.Count), new ParallelOptions { MaxDegreeOfParallelism = 18 }, //max number of threads modified to use locally
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        allPeaks[i].XIC = PeakCurve.GetXIC(allPeaks[i], allPeaks, isdCommonParameters);
                    }
                });
            var sorted = allPeaks.OrderByDescending(p => p.XIC.Peaks.Count()).ToList();
            var count = sorted.Select(p => p.XIC.Peaks.Count).ToList();

            int k = 0;
        }

        [Test]
        public static void TestXICFromLFQ()
        {
            string filePath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";
            MyFileManager myFileManager = new MyFileManager(true);
            var digestionParam = new DigestionParams(protease: "top-down");
            CommonParameters isdCommonParameters = new CommonParameters(digestionParams: digestionParam, trimMsMsPeaks: false);
            var myMsDataFile = myFileManager.LoadFile(filePath, isdCommonParameters);
            var allMs2Scans = myMsDataFile.Scans.Where(s => s.MsnOrder == 2).ToList();
            var allPeaks = new List<Peak>();
            int index = 0;
            foreach (var scan in allMs2Scans)
            {
                var spectrum = scan.MassSpectrum;
                for (int i = 0; i < spectrum.XArray.Length; i++)
                {
                    Peak newPeak = new Peak(spectrum.XArray[i], scan.RetentionTime, spectrum.YArray[i], scan.MsnOrder, scan.OneBasedScanNumber, index);
                    allPeaks.Add(newPeak);
                    index++;
                }
            }
            var table = XICfromLFQ.GetXICTable(allPeaks, 1000);
            var allRTs = new List<double>[allPeaks.Count];
            for (int i = 0; i < allPeaks.Count; i++)
            {
                var RTs = allPeaks[i].otherPeaks.Select(p => p.RT).ToList();
                allRTs[i] = RTs;
            }
        }
    }
}
