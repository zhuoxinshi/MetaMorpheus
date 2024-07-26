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
using Readers;
using static Nett.TomlObjectFactory;
using System.Diagnostics;

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
        public static void TestConvertingAllMs2Peaks()
        {
            string filePath = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            MyFileManager myFileManager = new MyFileManager(true);
            var digestionParam = new DigestionParams(protease: "top-down");
            CommonParameters isdCommonParameters = new CommonParameters(digestionParams: digestionParam, trimMsMsPeaks: false);
            var myMsDataFile = myFileManager.LoadFile(filePath, isdCommonParameters);
            var allMs2Scans = myMsDataFile.Scans.Where(s => s.MsnOrder == 2).ToList();
            var allPeaks = new List<EngineLayer.Peak>();
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
                        allPeaks[i].XICforDIA = PeakCurve.GetXIC(allPeaks[i], allPeaks, isdCommonParameters);
                    }
                });
            var sorted = allPeaks.OrderByDescending(p => p.XICforDIA.Peaks.Count()).ToList();
            var count = sorted.Select(p => p.XICforDIA.Peaks.Count).ToList();

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
                var RTs = allPeaks[i].XICpeaks.Select(p => p.RT).ToList();
                allRTs[i] = RTs;
            }
        }

        [Test]
        public static void TestIndividualScan()
        {
            string filePath = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);
            var allScans = myMsDataFile.Scans.OrderBy(s => s.OneBasedScanNumber).ToArray();
            List<double> rtShift = new List<double>();
            for (int i = 0; i < allScans.Length - 1; i++)
            {
                double rtDiff = allScans[i].RetentionTime - allScans[i + 1].RetentionTime;
                rtShift.Add(rtDiff);
            }
            double meanRTdiff = rtShift.Average();
            int binSize = 100;
            var msNScans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms1Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
            var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
            var ms2Peaks = Peak.GetAllPeaks(ms2Scans);
            var ms1Peaks = Peak.GetAllPeaks(ms1Scans);
            var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
            var ms2Table = XICfromLFQ.GetXICTable(ms2Peaks, binSize);

            //Test 2 ms1 peaks
            //883.578 ms2
            //var corr = XICfromLFQ.CalculateCorrelation(714.63904, 857.36583, ms1Table, ms1Table, 100, 0);
            //var corr2 = XICfromLFQ.CalculateCorrelation(714.63904, 883.578, ms1Table, ms2Table, 100, meanRTdiff);
            MsDataScan ms1scan = myMsDataFile.Scans.Where(s => s.OneBasedScanNumber == 113).First();
            MsDataScan ms2scan = myMsDataFile.Scans.Where(s => s.OneBasedScanNumber == 114).First();
            List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[1];
            List<(double, int, double, double)> precursors = new List<(double, int, double, double)>();

            var massList = new List<double>();
            foreach (MassSpectrometry.IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
                                    ms1scan.MassSpectrum, task.CommonParameters.PrecursorDeconvolutionParameters))
            {
                if (envelope.Charge == 1)
                {
                    continue;
                }
                if (massList.Contains(envelope.MonoisotopicMass))
                {
                    continue;
                }
                massList.Add(envelope.MonoisotopicMass);
                double monoPeakMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
                var highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).First().mz;
                precursors.Add((monoPeakMz, envelope.Charge, ms1scan.RetentionTime, highestPeakMz));
            }
            scansWithPrecursors[0] = new List<Ms2ScanWithSpecificMass>();
            IsotopicEnvelope[] neutralExperimentalFragments = null;
            if (task.CommonParameters.DissociationType != DissociationType.LowCID)
            {
                neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scan, task.CommonParameters);
            }

            foreach (var precursor in precursors)
            {
                // assign precursor for this MS2 scan
                var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.Item1,
                    precursor.Item2, filePath, task.CommonParameters, neutralExperimentalFragments, precursor.Item3, mostAbundantPrePeak: precursor.Item4);

                scansWithPrecursors[0].Add(scan);
            }
            //scansWithPrecursors = XICfromLFQ.GroupFragmentIonsXIC(scansWithPrecursors, ms1Table, ms2Table, task.CommonParameters, binSize, meanRTdiff);
            SpectralMatch[] fileSpecificPsms = new PeptideSpectralMatch[scansWithPrecursors[0].Count];
        }

        [Test]
        public static void CheckCorrelationOfMatchedIons()
        {
            string psmFileDDA = @"E:\ISD Project\ISD_240606\TestSnipDDA\Task1-SearchTask\AllPSMs.psmtsv";
            var psmsDDA = PsmTsvReader.ReadTsv(psmFileDDA, out List<string> warnings);
            string psmFileISD = @"E:\ISD Project\TestIsdDataAnalysis\1pmol_5uL_ISD_RT32.16-35.59_RT_XIC100_2\AllPSMs.psmtsv";
            var psmsISD = PsmTsvReader.ReadTsv(psmFileISD, out List<string> warnings2);
            string psmFileISDnoGroup = @"E:\ISD Project\TestIsdDataAnalysis\ISD_noGrouping_libraryAdded\Task1-SearchTask\AllPSMs.psmtsv";
            var psmsISD_ng = PsmTsvReader.ReadTsv(psmFileISDnoGroup, out List<string> warnings3);

            List<PsmFromTsv>[] matchedPsms = new List<PsmFromTsv>[psmsDDA.Count];
            for (int i = 0; i < psmsDDA.Count; i++)
            {
                matchedPsms[i] = new List<PsmFromTsv>();
                foreach (var psmISD in psmsISD)
                {
                    if (psmsDDA[i].PrecursorCharge == psmISD.PrecursorCharge && psmsDDA[i].FullSequence == psmISD.FullSequence && Math.Abs(psmsDDA[i].PrecursorMass - psmISD.PrecursorMass) < 0.05
                        && Math.Abs(psmsDDA[i].PrecursorMz - psmISD.PrecursorMz) < 0.01)
                    {
                        matchedPsms[i].Add(psmISD);
                    }
                }
            }

            string filePath_ISD = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            MyFileManager myFileManager = new MyFileManager(true);
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myMsDataFile = myFileManager.LoadFile(filePath_ISD, task.CommonParameters);
            int binSize = 100;
            var allScans = myMsDataFile.Scans.ToArray();
            var msNScans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms1Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
            var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
            var ms2Peaks = Peak.GetAllPeaks(ms2Scans);
            var ms1Peaks = Peak.GetAllPeaks(ms1Scans);
            var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
            var ms2Table = XICfromLFQ.GetXICTable(ms2Peaks, binSize);
            List<double> rtShift = new List<double>();
            for (int i = 0; i < allScans.Length - 1; i++)
            {
                double rtDiff = allScans[i].RetentionTime - allScans[i + 1].RetentionTime;
                rtShift.Add(rtDiff);
            }
            double meanRTdiff = rtShift.Average();

            int ms2ScanNum = 64;
            string outputDirectory = @"E:\ISD Project\TestIsdDataAnalysis";
            string outputPath = Path.Combine(outputDirectory, $"matchedIons_{ms2ScanNum}.csv");
            var matchedIonsISD = psmsISD.Where(p => p.Ms2ScanNumber == ms2ScanNum).First().MatchedIons;
            var matchedIonsDDA = psmsDDA[9].MatchedIons;
            var matchedIonsISD_ng = psmsISD_ng.Where(p => p.Ms2ScanNumber == ms2ScanNum).First().MatchedIons;

            using (var sw = new StreamWriter(File.Create(outputPath)))
            {
                sw.WriteLine("Mz,Intensity");
                foreach (var ion in matchedIonsDDA)
                {
                    sw.WriteLine($"{ion.Mz},{ion.Intensity}");
                }
                sw.WriteLine("Mz,Intensity");
                foreach (var ion in matchedIonsISD)
                {
                    sw.WriteLine($"{ion.Mz},{ion.Intensity}");
                }
                sw.WriteLine("Mz,Intensity");
                foreach (var ion in matchedIonsISD_ng)
                {
                    sw.WriteLine($"{ion.Mz},{ion.Intensity}");
                }
            }

            var matchedXICs = new List<(double rt, double intensity, double mz, double corr)>();
            string outputPath2 = Path.Combine(outputDirectory, $"XICs_matchedIons_{ms2ScanNum}.csv");
            using (var sw2 = new StreamWriter(File.Create(outputPath2)))
            {
                double preMz = 714.72;
                int roundedPreMz = 71472;
                var prePeaks = ms1Table[roundedPreMz].OrderBy(p => p.RT);
                foreach (var peak in prePeaks)
                {
                    matchedXICs.Add((peak.RT, peak.Intensity, peak.Mz, 1));
                }

                foreach (var ion in matchedIonsISD)
                {
                    int roundedMz = (int)Math.Round(ion.Mz * binSize, 0);
                    if (ms2Table[roundedMz] != null)
                    {
                        var peaks = ms2Table[roundedMz].OrderBy(p => p.RT);
                        double corr = XICfromLFQ.CalculateCorrelation(ms1Peaks.ToList(), ms2Peaks, meanRTdiff);
                        foreach (var peak in peaks)
                        {
                            matchedXICs.Add((peak.RT, peak.Intensity, peak.Mz, corr));
                        }
                    }
                    else
                    {
                        matchedXICs.Add((double.NaN, double.NaN, ion.Mz, double.NaN));
                    }
                }
                sw2.WriteLine("Retention Time,Intensity,rounded_mz,corr");
                foreach (var peak in matchedXICs)
                {
                    sw2.WriteLine($"{peak.rt}, {peak.intensity}, {peak.mz}, {peak.corr}");
                }
            }

            string outputPath3 = Path.Combine(outputDirectory, $"XICs_FilteredIons_{ms2ScanNum}.csv");
            var matchedMzs = matchedIonsISD.Select(p => p.Mz).ToList();
            var filteredIons = matchedIonsISD_ng.Where(i => !matchedMzs.Contains(i.Mz)).ToList();
            using (var sw3 = new StreamWriter(File.Create(outputPath3)))
            {
                sw3.WriteLine("rt, intensity, rounded_mz, corr");
                double preMz = 714.72;
                int roundedPreMz = 71472;
                var prePeaks = ms1Table[roundedPreMz].OrderBy(p => p.RT);
                foreach (var peak in prePeaks)
                {
                    sw3.WriteLine($"{peak.RT}, {peak.Intensity}, {preMz}, {1}");
                }
                foreach (var ion in filteredIons)
                {
                    int roundedMz = (int)Math.Round(ion.Mz * binSize, 0);
                    if (ms2Table[roundedMz] != null)
                    {
                        var peaks = ms2Table[roundedMz].OrderBy(p => p.RT);
                        double corr = XICfromLFQ.CalculateCorrelation(ms1Peaks.ToList(), ms2Peaks, meanRTdiff);
                        foreach (var peak in peaks)
                        {
                            sw3.WriteLine($"{peak.RT}, {peak.Intensity}, {roundedMz / (double)binSize}, {corr}");
                        }
                    }
                    else
                    {
                        sw3.WriteLine($"{double.NaN}, {double.NaN}, {ion.Mz}, {double.NaN}");
                    }
                }            
            }

            var sameIonsDDA = new List<MatchedFragmentIon>();
            foreach(var ion in matchedIonsDDA)
            {
                foreach(var ionISD in matchedIonsISD)
                {
                    if(Math.Abs(ionISD.Mz - ion.Mz) < 0.01)
                    {
                        sameIonsDDA.Add(ion);
                    }
                }
            }
            string outputPath4 = Path.Combine(outputDirectory, $"XICs_sameIonsDDA_{ms2ScanNum}.csv");
            using (var sw = new StreamWriter(File.Create(outputPath2)))
            {
                sw.WriteLine("rt, intensity, rounded_mz, corr");
                double preMz = 714.72;
                int roundedPreMz = 71472;
                var prePeaks = ms1Table[roundedPreMz].OrderBy(p => p.RT);
                foreach (var peak in prePeaks)
                {
                    sw.WriteLine($"{peak.RT}, {peak.Intensity}, {preMz}, {1}");
                }

                foreach (var ion in sameIonsDDA)
                {
                    int roundedMz = (int)Math.Round(ion.Mz * binSize, 0);
                    
                    if (ms2Table[roundedMz] != null)
                    {
                        var peaks = ms2Table[roundedMz].OrderBy(p => p.RT);
                        double corr = XICfromLFQ.CalculateCorrelation(ms1Peaks.ToList(), ms2Peaks, meanRTdiff);
                        foreach (var peak in peaks)
                        {
                            sw.WriteLine($"{peak.RT}, {peak.Intensity}, {(double)roundedMz / binSize}, {corr}");
                        }
                    }
                    else
                    {
                        sw.WriteLine($"{double.NaN}, {double.NaN}, {ion.Mz}, {double.NaN}");
                    }
                }
            }
        }

        //[Test]
        //public static void TestCompareMatchedWithFilteredFragmentIons()
        //{
        //    string filePath_ISD = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
        //    MyFileManager myFileManager = new MyFileManager(true);
        //    string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
        //    SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
        //    var myMsDataFile = myFileManager.LoadFile(filePath_ISD, task.CommonParameters);
        //    int binSize = 100;
        //    var allScans = myMsDataFile.Scans.ToArray();
        //    var msNScans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
        //    var ms1Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
        //    var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
        //    var ms2Peaks = Peak.GetAllPeaks(ms2Scans);
        //    var ms1Peaks = Peak.GetAllPeaks(ms1Scans);
        //    var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
        //    var ms2Table = XICfromLFQ.GetXICTable(ms2Peaks, binSize);

        //    var tol = new PpmTolerance(5);
        //    int roundedMz = 107224;
        //    if (ms2Table[roundedMz - 1] != null)
        //    {
        //        var peaks = ms2Table[roundedMz - 1];
        //    }
        //    if (ms2Table[roundedMz + 1] != null)
        //    {
        //        var peaks2  = ms2Table[roundedMz + 1];
        //    }

        //    List<double> rtShift = new List<double>();
        //    for (int i = 0; i < allScans.Length - 1; i++)
        //    {
        //        double rtDiff = allScans[i].RetentionTime - allScans[i + 1].RetentionTime;
        //        rtShift.Add(rtDiff);
        //    }
        //    double meanRTdiff = rtShift.Average();
        //    var scansWithPrecursors = MetaMorpheusTask._GetMs2ScansForDIA_XIC(myMsDataFile, filePath_ISD, task.CommonParameters);

        //    int ms2ScanNum = 50;
        //    int index = ms2ScanNum / 2 - 1;
        //    var scan1 = scansWithPrecursors[index].Where(p => p.PrecursorCharge == 12 && Math.Abs(p.PrecursorMass - 8559.5) < 0.1).First();
        //    var scan2 = scansWithPrecursors[index].Where(p => Math.Abs(p.MostAbundantPrePeak - 855.5) < 0.1).First();
        //    var scanWithPre1 = XICfromLFQ.GroupFragmentIonsForOneScan(scan1, ms1Scans, ms2Scans, ms1Table, ms2Table, task.CommonParameters, binSize, meanRTdiff);
        //    var scanWithPre2 = XICfromLFQ.GroupFragmentIonsForOneScan(scan2, ms1Scans, ms2Scans, ms1Table, ms2Table, task.CommonParameters, binSize, meanRTdiff);
        //    var scansWithPre = new MsDataScan[] { scanWithPre1.TheScan, scanWithPre2.TheScan };

        //    //string outputDirectory = @"E:\ISD Project\TestIsdDataAnalysis\XIC_visualization";
        //    //string XICpath1 = Path.Combine(outputDirectory, $"XICs_{ms2ScanNum}_{scan1.MostAbundantPrePeak}_LFQmethod_apexFilter.csv");
        //    //XICfromLFQ.VisualizeCorrelations(scan1, ms1Scans, ms2Scans, ms1Table, ms2Table, task.CommonParameters.PrecursorMassTolerance, binSize, meanRTdiff, XICpath1);
        //    //string XICpath2 = Path.Combine(outputDirectory, $"XICs_{ms2ScanNum}_{scan2.MostAbundantPrePeak}_LFQmethod_apexFilter.csv");
        //    //XICfromLFQ.VisualizeCorrelations(scan2, ms1Scans, ms2Scans, ms1Table, ms2Table, task.CommonParameters.PrecursorMassTolerance, binSize, meanRTdiff, XICpath2);

        //    string outPath = @"E:\ISD Project\TestIsdDataAnalysis\scansAfterGrouping_0.9.mzML";
        //    SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
        //        null, null, null);
        //    GenericMsDataFile msFile = new GenericMsDataFile(scansWithPre, genericSourceFile);
        //    msFile.ExportAsMzML(outPath, false);
        //    //XICfromLFQ.CompareMatchedWithFilteredFragmentIons(scan1, ms1Table, ms2Table, task.CommonParameters, 100, meanRTdiff, ms2ScanNum);
        //}

        //[Test]
        //public static void TestXICwithMatchedFragmentIon()
        //{
        //    string filePath_ISD = @"E:\ISD Project\ISD_240606\06-11-24_mix_sample1_2uL_ISD.mzML";
        //    MyFileManager myFileManager = new MyFileManager(true);
        //    string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
        //    SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
        //    var myMsDataFile = myFileManager.LoadFile(filePath_ISD, task.CommonParameters);
        //    int binSize = 100;
        //    var allScans = myMsDataFile.Scans.Where(s => s.RetentionTime > 45 && s.RetentionTime < 48.1).ToArray();
        //    var ms1Scans = allScans.Where(x => x.MsnOrder == 1).ToArray();
        //    var ms2Scans = allScans.Where(p => p.MsnOrder == 2).ToArray();
        //    var ms2Peaks = Peak.GetAllPeaks(ms2Scans);
        //    var ms1Peaks = Peak.GetAllPeaks(ms1Scans);
        //    var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
        //    var ms2Table = XICfromLFQ.GetXICTable(ms2Peaks, binSize);
        //    double rtShift = XICfromLFQ.GetRTshift(allScans);
        //    var scansWithPrecursors = MetaMorpheusTask._GetMs2ScansForDIA_XIC(myMsDataFile, filePath_ISD, task.CommonParameters);

        //    int ms2ScanNum = 1520;
        //    int index = ms2ScanNum / 2 - 1;
        //    var scans = scansWithPrecursors[index];
        //    var scan = scans.Where(s => Math.Abs(s.PrecursorMass - 18351) < 1).First();

        //    string outputDirectory = @"E:\ISD Project\TestIsdDataAnalysis\XIC_visualization";
        //    string XICpath = Path.Combine(outputDirectory, $"XICs_{ms2ScanNum}_{Math.Round(scan.MostAbundantPrePeak, 2)}_LFQmethod_apexFilter_NoPeakFilter.csv");
        //    XICfromLFQ.VisualizeCorrelations(scan, ms1Scans, ms2Scans, ms1Table, ms2Table, task.CommonParameters.PrecursorMassTolerance, binSize, rtShift, XICpath);

        //    var ms1XIC = XICfromLFQ.GetXIC(scan.MostAbundantPrePeak, ms1Scans, ms1Table, task.CommonParameters.PrecursorMassTolerance, binSize);
        //    string psmsPath = @"E:\ISD Project\TestIsdDataAnalysis\sample1_2uL_RT45.01-48.09_normalSearch\Task1-SearchTask\AllPSMs.psmtsv";
        //    var psms = PsmTsvReader.ReadTsv(psmsPath, out List<string> warnings);
        //    var ionsCA = psms.Take(26).Select(p => p.MatchedIons).SelectMany(p => p).ToList();
        //    var ionMzs = ionsCA.GroupBy(i => Math.Round(i.Mz, 2)).Select(g => g.Key).ToList();
        //    var XICs = new List<(double rt, double intensity, double mz, double corr)>();
        //    foreach(var mz in ionMzs)
        //    {
        //        var XIC = XICfromLFQ.GetXIC(mz, ms2Scans, ms2Table, task.CommonParameters.ProductMassTolerance, binSize);
        //        double corr = XICfromLFQ.CalculatePearsonCorr(ms1XIC, XIC, rtShift);
        //        foreach(var xic in XIC)
        //        {
        //            XICs.Add((xic.RT, xic.Intensity, xic.Mz, corr));
        //        }
        //    }
        //}

        [Test]
        public static void TestDDAfragments()
        {
            string filePath_DDA = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_DDA_RT31.94-35.64.mzML";
            MyFileManager myFileManager = new MyFileManager(true);
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myMsDataFile = myFileManager.LoadFile(filePath_DDA, task.CommonParameters);
            var scansWithPrecursors = MetaMorpheusTask._GetMs2Scans(myMsDataFile, filePath_DDA, task.CommonParameters);
            for(int i = 0; i < scansWithPrecursors.Length; i++)
            {
                
            }

        }
    }
}
    
