using Chemistry;
using Easy.Common.Extensions;
using EngineLayer.HistogramAnalysis;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.RootFinding;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using MzLibUtil;
using Newtonsoft.Json.Linq;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Data;
using System.IO;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using ThermoFisher.CommonCore.Data;
using ThermoFisher.CommonCore.Data.Business;
using static System.Formats.Asn1.AsnWriter;

namespace EngineLayer.ISD
{
    public class XICfromLFQ
    {
        public static List<Peak>[] GetXICTable(List<Peak> allPeaks, int binsPerDalton)
        {
            var table = new List<Peak>[(int)Math.Ceiling(allPeaks.Max(p => p.Mz) * binsPerDalton) + 1];
            //var peaks = allPeaks.ToArray();
            for (int i = 0; i< allPeaks.Count; i++)
            {
                int roundedMz = (int)Math.Round(allPeaks[i].Mz * binsPerDalton, 0);

                if (table[roundedMz] == null)
                {
                    table[roundedMz] = new List<Peak>();
                }
                table[roundedMz].Add(allPeaks[i]);
                allPeaks[i].XICpeaks = table[roundedMz];
            }
            return table;
        }

        public static Peak GetPeakFromScan(double mz, MsDataScan scan, List<Peak>[] peakTable, Tolerance tolerance, int binSize)
        {
            Peak bestPeak = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(mz) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(mz) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < peakTable.Length && peakTable[j] != null)
                {
                    var bin = peakTable[j].OrderBy(p => p.ScanNumber).ToList();
                    int index = Array.BinarySearch(bin.Select(p => p.ScanNumber).ToArray(), scan.OneBasedScanNumber);
                    if (index < 0)
                    {
                        continue;
                    }

                    for (int i = index; i < bin.Count; i++)
                    {
                        Peak peak = bin[i];

                        if (peak.ScanNumber > scan.OneBasedScanNumber)
                        {
                            break;
                        }

                        if (tolerance.Within(peak.Mz, mz) && peak.ScanNumber == scan.OneBasedScanNumber
                            && (bestPeak == null || Math.Abs(peak.Mz - mz) < Math.Abs(bestPeak.Mz - mz)))
                        {
                            bestPeak = peak;
                        }
                    }
                }
            }
            return bestPeak;
        }
        public static List<Peak> GetXIC(double mz, MsDataScan[] scans, List<Peak>[] peakTable, Tolerance tolerance, int binSize)
        {
            List<Peak> XICpeaks = new List<Peak>();
            for(int i = 0; i < scans.Length; i++)
            {
                var peak = GetPeakFromScan(mz, scans[i], peakTable, tolerance, binSize);
                if(peak != null)
                {
                    XICpeaks.Add(peak);
                }
            }
            return XICpeaks;
        }

        public static double CalculatePearsonCorr(List<Peak> ms1XIC, List<Peak> ms2XIC, double rtShift)
        {
            if (ApexFilter(ms1XIC, ms2XIC, 0.3) == true)
            {
                return 0;
            }
            var RT_1 = ms1XIC.OrderBy(p => p.RT).Select(p => Math.Round(p.RT, 2)).ToArray();
            var RT_2 = ms2XIC.OrderBy(p => p.RT).Select(p => Math.Round((p.RT + rtShift), 2)).ToArray();

            if (RT_1 == null || RT_2 == null)
            {
                return double.NaN;
            }

            var ms1Intensity = new List<double>();
            var ms2Intensity = new List<double>();
            for (int i = 0; i < RT_1.Length; i++)
            {
                int index = Array.BinarySearch(RT_2, RT_1[i]);
                if (index >= 0)
                {
                    ms1Intensity.Add(ms1XIC[i].Intensity);
                    ms2Intensity.Add(ms2XIC[index].Intensity);
                }
            }
            if (ms1Intensity.Count >= 5 && ms2Intensity.Count >=5)
            {
                // Calculate Pearson correlation
                double correlation = Correlation.Pearson(ms1Intensity, ms2Intensity);
                return correlation;
            }
            else
            {
                return double.NaN;
            }
        } 

        public static List<Ms2ScanWithSpecificMass>[] GroupFragmentIonsXIC(List<Ms2ScanWithSpecificMass>[] scansWithPrecursors, MsDataScan[] ms1scans, MsDataScan[] ms2scans, List<Peak>[] ms1Table, 
            List<Peak>[] ms2Table, CommonParameters commonParameters, int binSize, double rtShift)
        {
            for (int i = 0; i < scansWithPrecursors.Length; i++)
            {
                if (scansWithPrecursors[i].Count > 0)
                {
                    for (int j = 0; j < scansWithPrecursors[i].Count; j++)
                    {
                        var targetScan = scansWithPrecursors[i][j];
                        scansWithPrecursors[i][j] = GroupFragmentIonsForOneScan(targetScan, ms1scans, ms2scans, ms1Table, ms2Table, commonParameters, binSize, rtShift);
                    }
                }
            }
            return scansWithPrecursors;
        }

        public static Ms2ScanWithSpecificMass GroupFragmentIonsForOneScan(Ms2ScanWithSpecificMass targetScan, MsDataScan[] ms1scans, MsDataScan[] ms2scans, List<Peak>[] ms1Table,
            List<Peak>[] ms2Table, CommonParameters commonParameters, int binSize, double rtShift)
        {
            List<double> diaMzs = new List<double>();
            List<double> diaIntensities = new List<double>();
            var allMzs = targetScan.TheScan.MassSpectrum.XArray;
            var ms1Peaks = GetXIC(targetScan.MostAbundantPrePeak, ms1scans, ms1Table, commonParameters.PrecursorMassTolerance, binSize);
            for (int k = 0; k < allMzs.Length; k++)
            {
                var ms2Peaks = GetXIC(allMzs[k], ms2scans, ms2Table, commonParameters.PrecursorMassTolerance, binSize);
                double correlation = CalculatePearsonCorr(ms1Peaks, ms2Peaks, rtShift);
                if (correlation > 0.8)
                {
                    diaMzs.Add(allMzs[k]);
                    diaIntensities.Add(targetScan.TheScan.MassSpectrum.YArray[k]);
                }
            }
            //change scan number!
            MzSpectrum diaSpectrum = new MzSpectrum(diaMzs.ToArray(), diaIntensities.ToArray(), false);
            MsDataScan newScan = new MsDataScan(diaSpectrum, targetScan.TheScan.OneBasedScanNumber, targetScan.TheScan.MsnOrder, targetScan.TheScan.IsCentroid,
                targetScan.TheScan.Polarity, targetScan.TheScan.RetentionTime, targetScan.TheScan.ScanWindowRange, targetScan.TheScan.ScanFilter, targetScan.TheScan.MzAnalyzer,
                targetScan.TheScan.TotalIonCurrent, targetScan.TheScan.InjectionTime, targetScan.TheScan.NoiseData, targetScan.TheScan.NativeId,
                targetScan.TheScan.SelectedIonMZ, targetScan.TheScan.SelectedIonChargeStateGuess, targetScan.TheScan.SelectedIonIntensity,
                targetScan.TheScan.IsolationMz, targetScan.TheScan.IsolationWidth, targetScan.TheScan.DissociationType, null,
                targetScan.TheScan.SelectedIonMonoisotopicGuessMz, targetScan.TheScan.HcdEnergy, targetScan.TheScan.ScanDescription);
            var scanWithPrecursor = new Ms2ScanWithSpecificMass(newScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath,
    commonParameters, null, targetScan.Pre_RT, mostAbundantPrePeak: targetScan.MostAbundantPrePeak);
            return scanWithPrecursor;
        }

        public static double CalculateCorrelation(List<Peak> ms1Peaks, List<Peak> ms2Peaks, double rtShift)
        {
            var matchedPeaks = MatchRTs(ms1Peaks, ms2Peaks, rtShift);
            var dataSet1 = matchedPeaks.Select(p => p.Item1).ToArray();
            var dataSet2 = matchedPeaks.Select(p => p.Item2).ToArray();
            if (matchedPeaks.Count >= 5)
            {
                // Calculate Pearson correlation
                double correlation = Correlation.Pearson(dataSet1, dataSet2);
                return correlation;
            }
            else
            {
                return double.NaN;
            }
        }

        public static bool ApexFilter(List<Peak> ms1Peaks, List<Peak> ms2Peaks, double threshold)
        {
            double ms1Apex = ms1Peaks.OrderByDescending(p => p.Intensity).First().RT;
            double ms2Apex = ms2Peaks.OrderByDescending(p => p.Intensity).First().RT;
            if (Math.Abs(ms1Apex - ms2Apex) < threshold)
            {
                return false;
            }
            return true;
        }
        public static List<(double, double)> MatchRTs(List<Peak> ms1Peaks, List<Peak> ms2Peaks, double rtShift)
        {
            var RT_1 = ms1Peaks.Select(p => Math.Round(p.RT, 2)).ToArray();
            var RT_2 = ms2Peaks.OrderBy(p => Math.Round(p.RT, 2)).Select(p => Math.Round((p.RT + rtShift), 2)).ToArray();
            if (ApexFilter(ms1Peaks.ToList(), ms2Peaks, 0.3) == true)
            {
                return new List<(double, double)> {(0,0)};
            }
            if (RT_1 == null || RT_2 == null)
            {
                return new List<(double, double)> { (-1.0, -1.0) };
            }

            List<(double, double)> list_intensity = new List<(double, double)>();
            List<(double, double)> list_rt = new List<(double, double)>();
            for (int i = 0; i < RT_1.Length; i++)
            {
                int index = Array.BinarySearch(RT_2, RT_1[i]);
                if(index >= 0)
                {
                    list_intensity.Add((ms1Peaks[i].Intensity, ms2Peaks[index].Intensity));
                    list_rt.Add((ms1Peaks[i].RT, ms2Peaks[index].RT));
                }
            }
            return list_intensity;
        }

        public static void CompareMatchedWithFilteredFragmentIons(Ms2ScanWithSpecificMass targetScan, List<Peak>[] ms1Table,
            List<Peak>[] ms2Table, CommonParameters commonParameters, int binSize, double rtShift, int ms2ScanNum)
        {
            List<double> diaMzs = new List<double>();
            List<double> diaIntensities = new List<double>();
            var matchedXICs = new List<(double rt, double intensity, double mz, double corr)>();
            var unmatchedXICs = new List<(double rt, double intensity, double mz, double corr)>();
            int roundedMs1mz = (int)Math.Round(targetScan.MostAbundantPrePeak * binSize, 0);
            var ms1Peaks = ms1Table[roundedMs1mz].OrderBy(p => p.RT).ToArray();
            foreach (var peak in ms1Peaks)
            {
                matchedXICs.Add((peak.RT, peak.Intensity, peak.Mz, 1));
                unmatchedXICs.Add((peak.RT, peak.Intensity, peak.Mz, 1));
            }

            var allMzs = targetScan.TheScan.MassSpectrum.XArray;
            for (int k = 0; k < allMzs.Length; k++)
            {
                int roundedMs2mz = (int)Math.Round(allMzs[k] * binSize, 0);
                var ms2Peaks = ms2Table[roundedMs2mz].ToList();
                if (ms2Table[roundedMs2mz - 1] != null)
                {
                    ms2Peaks.AddRange(ms2Table[roundedMs2mz - 1]);
                }
                if (ms2Table[roundedMs2mz + 1] != null)
                {
                    ms2Peaks.AddRange(ms2Table[roundedMs2mz + 1]);
                }
                double correlation = CalculateCorrelation(ms1Peaks.ToList(), ms2Peaks, rtShift);
                if (correlation > 0.9)
                {
                    diaMzs.Add(allMzs[k]);
                    diaIntensities.Add(targetScan.TheScan.MassSpectrum.YArray[k]);
                    foreach (var peak in ms2Peaks)
                    {
                        matchedXICs.Add((peak.RT, peak.Intensity, peak.Mz, correlation));
                    }
                }
                else
                {
                    foreach (var peak in ms2Peaks)
                    {
                        unmatchedXICs.Add((peak.RT, peak.Intensity, peak.Mz, correlation));
                    }
                }
            }

            //export to excel
            string outputDirectory = @"E:\ISD Project\TestIsdDataAnalysis\XIC_visualization";
            string outputPath1 = Path.Combine(outputDirectory, $"matchedIons_{ms2ScanNum}_addAdjacentXICs_wholeFile_apexFilter.csv");
            using (var sw1 = new StreamWriter(File.Create(outputPath1)))
            {
                sw1.WriteLine("Retention Time,Intensity,rounded_mz,corr");
                foreach (var xic in matchedXICs)
                {
                sw1.WriteLine($"{xic.rt},{xic.intensity},{xic.mz},{xic.corr}");
                }
            }
            string outputPath2 = Path.Combine(outputDirectory, $"unmatchedIons_{ms2ScanNum}_addAdjacentXICs_wholeFile_apexFilter.csv");
            using (var sw2 = new StreamWriter(File.Create(outputPath2)))
            {
                sw2.WriteLine("Retention Time,Intensity,rounded_mz,corr");
                foreach (var xic in unmatchedXICs)
                {
                    sw2.WriteLine($"{xic.rt},{xic.intensity},{xic.mz},{xic.corr}");
                }
            }
        }

        public static void VisualizeXICs(List<(double rt, double intensity, double mz, double corr)> XICs, string outputPath)
        {
            using (var sw = new StreamWriter(File.Create(outputPath)))
            {
                sw.WriteLine("Retention Time,Intensity,rounded_mz,corr");
                foreach (var xic in XICs)
                {
                    sw.WriteLine($"{xic.rt},{xic.intensity},{xic.mz},{xic.corr}");
                }
            }
        }

        public static void VisualizeCorrelations(Ms2ScanWithSpecificMass scan, MsDataScan[] ms1scans, MsDataScan[] ms2scans, List<Peak>[] ms1Table,
            List<Peak>[] ms2Table, Tolerance tolerance, int binSize, double rtShift, string outputPath)
        {
            var XICs = new List<(double rt, double intensity, double mz, double corr)>();
            var preXIC = GetXIC(scan.MostAbundantPrePeak, ms1scans, ms1Table, tolerance, binSize);
            foreach(var peak in preXIC)
            {
                XICs.Add((peak.RT, peak.Intensity, peak.Mz, 1));
            }
            foreach(double mz in scan.TheScan.MassSpectrum.XArray)
            {
                var XIC = GetXIC(mz, ms2scans, ms2Table, tolerance, binSize);
                double corr = CalculatePearsonCorr(preXIC, XIC, rtShift);
                foreach(var peak in XIC)
                {
                    XICs.Add((peak.RT, peak.Intensity, peak.Mz, corr));
                }
            }
            VisualizeXICs(XICs, outputPath);
        }

        public static double GetRTshift(MsDataScan[] allScans)
        {
            List<double> rtShift = new List<double>();
            for (int i = 0; i < allScans.Length - 1; i++)
            {
                double rtDiff = allScans[i].RetentionTime - allScans[i + 1].RetentionTime;
                rtShift.Add(rtDiff);
            }
            double meanRTdiff = rtShift.Average();
            return meanRTdiff;
        }

        }
              
    }

