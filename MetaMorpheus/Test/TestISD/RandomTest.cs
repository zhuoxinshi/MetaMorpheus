using Easy.Common.Extensions;
using EngineLayer.HistogramAnalysis;
using EngineLayer.ISD;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.RootFinding;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using MzLibUtil;
using Nett;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Data;
using System.IO;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using TaskLayer;
using EngineLayer;
using Omics.Fragmentation;
using Proteomics;
using System.Threading.Tasks;
using EngineLayer.ClassicSearch;
using SpectralAveraging;

namespace Test.TestISD
{
    public class RandomTest
    {
        [Test]
        public static void TestPeakNum()
        {
            string snipPath = @"E:\ISD Project\TestIsdDataAnalysis\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string rawPath = @"E:\ISD Project\ISD_240606\06-11-24_mix_sample1_2uL_ISD.raw";
            string mzmlPath = @"E:\ISD Project\ISD_240606\06-11-24_mix_sample1_2uL_ISD.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var snipFile = myFileManager.LoadFile(snipPath, task.CommonParameters);
            var rawFile = myFileManager.LoadFile(rawPath, task.CommonParameters);
            var mzmlFile = myFileManager.LoadFile(mzmlPath, task.CommonParameters);
            var peakNum = new List<(double, double, double)>();
            foreach(var scan in snipFile.Scans)
            {
                var rawNum = rawFile.Scans.First(s => s.RetentionTime == scan.RetentionTime).MassSpectrum.Size;
                var mzmlNum = mzmlFile.Scans.First(s => s.RetentionTime == scan.RetentionTime).MassSpectrum.Size;
                peakNum.Add((scan.MassSpectrum.Size, rawNum, mzmlNum));
            }
        }

        [Test]
        public static void TestXICUsingWholeFile()
        {
            string mzmlPath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var mzmlFile = myFileManager.LoadFile(mzmlPath, task.CommonParameters);
            var allScans = mzmlFile.Scans.Where(s => s.RetentionTime >32 && s.RetentionTime < 36).ToArray();
            int binSize = 100;
            var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2Scans = allScans.Where(p => p.MsnOrder == 2).ToArray();
            var ms2Peaks = Peak.GetAllPeaks(ms2Scans);
            var ms1Peaks = Peak.GetAllPeaks(ms1Scans);
            var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
            var ms2Table = XICfromLFQ.GetXICTable(ms2Peaks, binSize);

            var tol = new PpmTolerance(5);
            double mz = 725.92407;
            int roundedMz = 72592;
            var peaks = ms2Table[roundedMz];
            if (ms2Table[roundedMz - 1] != null)
            {
                peaks.AddRange(ms2Table[roundedMz - 1].Where(p => tol.Within(p.Mz, mz)));
            }
            if (ms2Table[roundedMz + 1] != null)
            {
                peaks.AddRange(ms2Table[roundedMz + 1].Where(p => tol.Within(p.Mz, mz)));
            }

            List<double> rtShift = new List<double>();
            for (int i = 0; i < allScans.Length - 1; i++)
            {
                double rtDiff = allScans[i].RetentionTime - allScans[i + 1].RetentionTime;
                rtShift.Add(rtDiff);
            }
            double meanRTdiff = rtShift.Average();
            var scansWithPrecursors = MetaMorpheusTask._GetMs2ScansForDIA_XIC(mzmlFile, mzmlPath, task.CommonParameters);
            int ms2ScanNum = 1084;
            int index = ms2ScanNum / 2 - 1;
            //var scan = scansWithPrecursors[index].Where(p => p.PrecursorCharge == 12 && Math.Abs(p.PrecursorMass - 8559.5) < 0.1).First();
            //XICfromLFQ.CompareMatchedWithFilteredFragmentIons(scan, ms1Table, ms2Table, task.CommonParameters, 100, meanRTdiff, ms2ScanNum);
        }

        //[Test]
        //public static void TestLFQmethod()
        //{
        //    string mzmlPath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";
        //    string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
        //    MyFileManager myFileManager = new MyFileManager(true);
        //    SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
        //    task.CommonParameters.DoDIA = true;
        //    var mzmlFile = myFileManager.LoadFile(mzmlPath, task.CommonParameters);
        //    var allScans = mzmlFile.Scans.Where(s => s.RetentionTime > 45 && s.RetentionTime < 48.1).ToArray();
        //    int binSize = 100;
        //    var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
        //    var ms2Scans = allScans.Where(p => p.MsnOrder == 2).ToArray();
        //    var ms2Peaks = Peak.GetAllPeaks(ms2Scans);
        //    var ms1Peaks = Peak.GetAllPeaks(ms1Scans);
        //    var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
        //    var ms2Table = XICfromLFQ.GetXICTable(ms2Peaks, binSize);
        //    double meanRTdiff = XICfromLFQ.GetRTshift(allScans);

        //    var scansWithPrecursors = MetaMorpheusTask._GetMs2ScansForDIA_XIC(mzmlFile, mzmlPath, task.CommonParameters);
        //    int ms2ScanNum = 1084;
        //    int index = ms2ScanNum / 2 - 1;
        //    var scan = scansWithPrecursors[index].Where(p => p.PrecursorCharge == 12 && Math.Abs(p.PrecursorMass - 8559.5) < 0.1).First();

        //    var allXICs = new List<(double rt, double intensity, double mz, double corr)>();
        //    var preXIC = XICfromLFQ.GetXIC(scan.MostAbundantPrePeak, ms1Scans, ms1Table, task.CommonParameters.PrecursorMassTolerance, binSize);
        //    foreach(var peak in preXIC)
        //    {
        //        allXICs.Add((peak.RT, peak.Intensity, peak.Mz, 1));
        //    }
        //    for (int i = 0; i < scan.TheScan.MassSpectrum.XArray.Length; i++)
        //    {
        //        var ms2XIC = XICfromLFQ.GetXIC(scan.TheScan.MassSpectrum.XArray[i], ms2Scans, ms2Table, task.CommonParameters.ProductMassTolerance, binSize);
        //        double corr = XICfromLFQ.CalculatePearsonCorr(preXIC, ms2XIC, meanRTdiff);
        //        foreach(var peak in ms2XIC)
        //        {
        //            allXICs.Add((peak.RT, peak.Intensity, peak.Mz, corr));
        //        }
        //    }
        //    string outputDirectory = @"E:\ISD Project\TestIsdDataAnalysis\XIC_visualization";
        //    string outputPath = Path.Combine(outputDirectory, $"XICs_{ms2ScanNum}_wholeFile_LFQmethod_apexFilter.csv");
        //    XICfromLFQ.VisualizeXICs(allXICs, outputPath);
        //    bool stop = true;
        //    bool stop2 = true;

        //    //alternative way to generate XICs for ms2 peaks:
        //    //Instead of generating ms2Table at the begining for all the ms2 peaks of the entire data file, when grouping peaks for each MS2WithSpecificMass,
        //    //generate small ms2Tables only in the precursor elution window

        //    //sample data: sample 1 45-55min, myo & beta-globulin
        //}

        //[Test]
        //public static void TestISDfragmentIonMatching()
        //{
        //    string filePath = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
        //    string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
        //    MyFileManager myFileManager = new MyFileManager(true);
        //    SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
        //    var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
        //    string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\1pmol_5uL_ISD_RT32.16-35.59_RT_debug";
        //    string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
        //    string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
        //    DbForTask db = new DbForTask(myDatabase, false);
        //    DbForTask lib = new DbForTask(library, false);
        //    task.RunTask(outputFolder, new List<DbForTask> { db, lib }, new List<string> { filePath }, "normal");
        //    var ms2Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 2).ToArray();
        //    var scansWithPrecursors = MetaMorpheusTask._GetMs2ScansForDIA_XIC(myMsDataFile, filePath, task.CommonParameters);

        //    int ms2ScanNum = 50;
        //    int index = ms2ScanNum / 2 - 1;
        //    var scan = scansWithPrecursors[index].Where(p => Math.Abs(p.PrecursorMass - 8559.5) < 0.1).First();
        //    var scanWithPre = new Ms2ScanWithSpecificMass[] { scan };
        //    var targetFragments = new Dictionary<DissociationType, List<Product>>();
        //    Protein ubiquitin = new Protein("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG", "P0CG48");
        //    var proteinList = new List<Protein>();
        //    proteinList.Add(ubiquitin);
        //    var fileSpecificPsms = new SpectralMatch[1];
        //    MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(task.CommonParameters.PrecursorMassTolerance, task.SearchParameters.MassDiffAcceptorType, 
        //        task.SearchParameters.CustomMdac);
        //    var newClassicSearchEngine = new ClassicSearchEngine(fileSpecificPsms, scanWithPre, null, null, task.SearchParameters.SilacLabels,
        //               task.SearchParameters.StartTurnoverLabel, task.SearchParameters.EndTurnoverLabel, proteinList, massDiffAcceptor, task.CommonParameters, 
        //               new List<(string FileName, CommonParameters Parameters)> { ("name", task.CommonParameters) }, 
        //               null, new List<string> { "id" }, task.SearchParameters.WriteSpectralLibrary);

        //    newClassicSearchEngine.Run();
        //}

        [Test]
        public static void TestGroupingMs1Peaks()
        {
            SearchTask task = new SearchTask();
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPsmCount2");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            Directory.CreateDirectory(outputFolder);
            MyFileManager myFileManager = new MyFileManager(true);
            var myMsDataFile = myFileManager.LoadFile(myFile, task.CommonParameters);
            var ms2WithPre = MetaMorpheusTask._GetMs2Scans(myMsDataFile, myFile, task.CommonParameters);
            var scans = myMsDataFile.GetAllScansList();
            var ms1scans = scans.Where(x => x.MsnOrder == 1).ToArray();
            var ms1Peaks = Peak.GetAllPeaks(ms1scans);
            var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, 100);
            var allXICs = EngineLayer.ISD.XIC.GetAllXICs(ms1Table);
            var allXICsForGrouping = new List<XIC>();
            foreach (var xic in allXICs)
            {
                allXICsForGrouping.Add(xic);
            }
            //Check peaks
            var xic1 = allXICs.Where(x => x.RoundedMz == 54980).First();
            var xic2 = allXICs.Where(x => x.RoundedMz == 55030).First();
            var xic3 = allXICs.Where(x => x.RoundedMz == 55080).First();
            double testCorr = XICgroup.GetCorr(xic1.XICpeaks, xic2.XICpeaks, 0);

            var allXICgroups = XICgroup.GroupWithApex(allXICsForGrouping, 0.5);
            var XICgroupsForDecon = allXICgroups.Where(g => g.XIClist.Count > 2).ToArray();

            var allSpectrum = new List<MzSpectrum>();
            foreach (var group in XICgroupsForDecon)
            {
                var mz = group.XIClist.OrderBy(x => x.AveragedMz).Select(x => x.AveragedMz).ToArray();
                var intensities = group.XIClist.OrderBy(x => x.AveragedMz).Select(x => x.AveragedIntensity).ToArray();
                MzSpectrum ms1spectrum = new MzSpectrum(mz, intensities, false);
                allSpectrum.Add(ms1spectrum);
            }
            var results = new List<IsotopicEnvelope>[allSpectrum.Count];
            var masses = new List<double>();
            for (int i = 0; i < allSpectrum.Count; i++)
            {
                results[i] = Deconvoluter.Deconvolute(allSpectrum[i], task.CommonParameters.PrecursorDeconvolutionParameters).ToList();
                masses.AddRange(results[i].Select(x => x.MonoisotopicMass));
            }
            
            //normal decon comparison
            var scansWithPre = MetaMorpheusTask.GetMs2Scans(myMsDataFile, myFile, task.CommonParameters);
            var allMasses = scansWithPre.Select(scan => scan.PrecursorMass).ToList();

            //For debug
            var peakNum = ms1scans.Sum(s => s.MassSpectrum.Size);
            var newPeakNum = allSpectrum.Sum(s => s.Size);
            var countXICinGroups = allXICgroups.Sum(g => g.XIClist.Count);
            var allPeaksInXICgroups = allXICgroups.Select(g => g.XIClist).SelectMany(xic => xic).Sum(xic => xic.XICpeaks.Count);
            var peakTableCount = ms1Table.Where(x => x != null).Sum(x => x.Count);
            var peakTableCount2 = ms1Table.Where(x => x != null).ToArray().Length;
            var peakCountAllXICs = allXICs.Sum(xic => xic.XICpeaks.Count);
        }

        [Test]
        public static void TestXIC()
        {
            SearchTask task = new SearchTask();
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPsmCount2");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            Directory.CreateDirectory(outputFolder);
            MyFileManager myFileManager = new MyFileManager(true);
            var myMsDataFile = myFileManager.LoadFile(myFile, task.CommonParameters);
            var scans = myMsDataFile.GetAllScansList();
            var ms1scans = scans.Where(x => x.MsnOrder == 1).ToArray();
            var ms1Peaks = Peak.GetAllPeaks(ms1scans);
            var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, 100);

            List<double> mzs = new List<double> { 549.797, 550.29947, 550.802 };
            foreach(double mz in mzs)
            {
                int roundedMz = (int)Math.Round(mz * 100, 0);
                var peaks = ms1Table[roundedMz];
            }
            
        }

        [Test]
        public static void search()
        {
            SearchTask task = new SearchTask();
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\Search results\smallCalibrateTest";
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            Directory.CreateDirectory(outputFolder);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
        }

        [Test]
        public static void TestAverageSpectra()
        {
            string filePath = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
            var allScans = myMsDataFile.GetAllScansList().ToArray();
            var msNScans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms1Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
            var ms2RawScans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
            var averageParam = new SpectralAveragingParameters()
            {
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectraFileAveragingType = SpectraFileAveragingType.AverageEverynScansWithOverlap,
                NumberOfScansToAverage = 5,
                ScanOverlap = 4,
                NormalizationType = NormalizationType.RelativeToTics,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly
            };

            var ms2Scans = SpectraFileAveraging.AverageSpectraFile(ms2RawScans.ToList(), averageParam);
            var ms1 = SpectraFileAveraging.AverageSpectraFile(ms1Scans.ToList(), averageParam);
            var all = SpectraFileAveraging.AverageSpectraFile(allScans.ToList(), averageParam);
            var ms2 = new MsDataScan[ms2RawScans.Length];
            for (int i = 0; i < ms2RawScans.Length; i++)
            {

            }
        }

        public static MsDataScan GetAveragedDataScanFromAveragedSpectrum(MzSpectrum averagedSpectrum,
        MsDataScan centralScan)
        {
            MsDataScan averagedScan = new(averagedSpectrum,
                centralScan.OneBasedScanNumber,
                centralScan.MsnOrder,
                centralScan.IsCentroid,
                centralScan.Polarity,
                centralScan.RetentionTime,
                averagedSpectrum.Range, null,
                centralScan.MzAnalyzer,
                averagedSpectrum.SumOfAllY,
                centralScan.InjectionTime,
                centralScan.NoiseData,
                centralScan.NativeId,
                centralScan.SelectedIonMZ,
                centralScan.SelectedIonChargeStateGuess,
                centralScan.SelectedIonIntensity,
                centralScan.IsolationMz,
                centralScan.IsolationWidth,
                centralScan.DissociationType,
                centralScan.OneBasedPrecursorScanNumber,
                centralScan.SelectedIonMonoisotopicGuessIntensity);
            var newNativeId =
                averagedScan.NativeId.Replace(averagedScan.NativeId.Split("=").Last(), centralScan.OneBasedScanNumber.ToString());
            averagedScan.SetNativeID(newNativeId);
            return averagedScan;
        }

        [Test]
        public static void TestGroupMs1()
        {
            string filePath = @"E:\ISD Project\TestIsdDataAnalysis\data\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            task.SearchParameters.WriteSpectralLibrary = true;
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
            string outputFolder = @"E:\ISD Project\TestIsdDataAnalysis\Search results\sample1_2uL_ISD_RT45.01-48.09_GroupXIC";
            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";

            int binSize = 100;
            var ms1Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
            var ms1Peaks = Peak.GetAllPeaks(ms1Scans);
            var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
            var allXICsTest = XIC.GetAllXICs_LFQ(ms1Peaks, ms1Scans, ms1Table, new PpmTolerance(20), binSize);
            var sortedXICs = allXICsTest.OrderBy(x => x.ApexRT).ToArray();
            var filteredXICs = sortedXICs.Where(x => x.XICpeaks.Count >= 10).OrderBy(x => x.AveragedMz).ToList();
            var allGroups = XICgroup.FindAllGroups2(filteredXICs.ToList(), 8, 0.7);

            var scan = ms1Scans.First();
            var filteredGroups = allGroups.Where(g => g.XIClist.Count >= 5);
            var filteredGroups2 = allGroups.Where(g => g.XIClist.Count >= 7);
            var filteredGroups3 = allGroups.Where(g => g.XIClist.Count >= 10);
            var filteredGroups4 = allGroups.Where(g => g.XIClist.Count >= 20);
            var deconResults = XICgroup.DeconvoluteNewMs1Scans(filteredGroups.ToArray(), task.CommonParameters);
            var sortedResults = deconResults.SelectMany(x => x).OrderByDescending(x => x.MonoisotopicMass).ToList();
            var groupsToLook = filteredXICs.Where(x => x.AveragedMz >= 1224.8 && x.AveragedMz <= 1225.8).OrderBy(x => x.AveragedMz).ToList();

            var peaksToLook = scan.MassSpectrum.XArray.Where(x => x >= 1218.9 && x <= 1220.16).ToArray();
            var XICtoLook = filteredXICs.Where(x => x.AveragedMz >= 854.4 && x.AveragedMz <= 854.9);
            var highestXIC = XICtoLook.OrderByDescending(x => x.AveragedIntensity).First();
            var corrlist = new List<double>();
            foreach(var xic in XICtoLook)
            {
                corrlist.Add(XICgroup.GetCorr(highestXIC.XICpeaks, xic.XICpeaks, 0));
            }

            var highestPeak = peaksToLook[8];
            
            var allXICs = new List<(double rt, double intensity, double mz, double corr)>();
            
            var higestPeakXIC = XICfromLFQ.GetXIC(highestPeak, ms1Scans, ms1Table, new PpmTolerance(10), binSize);
            
            //var higestPeakXIC = XIC.GetXICPeaksWithMz(ms1Table, highestPeak, binSize);

            foreach (var peak in higestPeakXIC)
            {
                allXICs.Add((peak.RT, peak.Intensity, peak.Mz, 1));
            }
            foreach (var xic in XICtoLook)
            {
                double corr = XIC.XICCorr(higestPeakXIC, xic.XICpeaks, 0);
                foreach(var peak in xic.XICpeaks)
                {
                    allXICs.Add((peak.RT, peak.Intensity, peak.Mz, corr));
                }
            }
            string output = @"E:\ISD Project\TestIsdDataAnalysis\XIC_visualization\TestMS1Group_New.csv";
            XICfromLFQ.VisualizeXICs(allXICs, output);

            foreach (double mz in peaksToLook)
            {
                var xic = XICfromLFQ.GetXIC(mz, ms1Scans, ms1Table, task.CommonParameters.DeconvolutionMassTolerance, binSize);
                //var xic = XIC.GetXICPeaksWithMz(ms1Table, mz, binSize);
                var corr = XICfromLFQ.CalculatePearsonCorr(higestPeakXIC, xic, 0);
                foreach (var peak in xic)
                {
                    allXICs.Add((peak.RT, peak.Intensity, peak.Mz, corr));
                }
            }

            string output2 = @"E:\ISD Project\TestIsdDataAnalysis\XIC_visualization\TestMS1Group_NoTol.csv";
            XICfromLFQ.VisualizeXICs(allXICs, output2);
        }
    }
}
