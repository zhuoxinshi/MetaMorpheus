using Easy.Common.Extensions;
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
using System.Data;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;


namespace EngineLayer
{
    public class PeakCurve
    {
        public PeakCurve(List<Peak> peaks, int msLevel, MzRange isolationWindow = null)
        {
            Peaks = peaks;
            MsLevel = msLevel;
            IsolationWindow = isolationWindow;

            if (peaks != null && peaks.Count > 0)
            {
                Average_Mz = peaks.Average(p => p.Mz);
            }
            else
            {
                Average_Mz = 0;
            }
        }

        public double Average_Mz { get; set; }
        public List<Peak> Peaks { get; set; }
        public int MsLevel { get; set; }
        public MzRange IsolationWindow { get; set; }


        public static double CalPeakCorr(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            var matchedPeaks = RTPairs(peakCurve1, peakCurve2);
            var dataSet1 = matchedPeaks.Select(p => p.Item1).ToArray();
            var dataSet2 = matchedPeaks.Select(p => p.Item2).ToArray();
            if (matchedPeaks.Count >= 5)
            {
                // Calculate Pearson correlation
                double correlation = CalculatePearsonCorrelation(dataSet1, dataSet2);
                return correlation;
            }
            else
            {
                return double.NaN;
            }
        }
        public static double CalculatePearsonCorrelation(double[] x, double[] y)
        {
            int n = x.Length;
            double sumX = 0, sumY = 0, sumXy = 0, sumX2 = 0, sumY2 = 0;
            for (int i = 0; i < n; i++)
            {
                sumX += x[i];
                sumY += y[i];
                sumXy += x[i] * y[i];
                sumX2 += x[i] * x[i];
                sumY2 += y[i] * y[i];
            }
            double numerator = n * sumXy - sumX * sumY;
            double denominator = Math.Sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));

            return numerator / denominator;
        }

        private static List<(double, double)> RTPairs(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            double[] RT_1 = peakCurve1.Peaks.Select(p => p.RT).ToArray();
            double[] Intensity_1 = peakCurve1.Peaks.Select(p => p.Intensity).ToArray();
            double[] RT_2 = peakCurve2.Peaks.Select(p => p.RT).ToArray();
            double[] Intensity_2 = peakCurve2.Peaks.Select(p => p.Intensity).ToArray();

            if (RT_1 == null || RT_2 == null)
            {
                return new List<(double, double)> { (-1.0, -1.0) };
            }

            List<(double, double)> list = new List<(double, double)>();
            List<(double, double)> list2 = new List<(double, double)>();
            List<(double, double)> list3 = new List<(double, double)>();

            for (int j = 0; j < RT_1.Length; j++)
            {
                list2.Add((RT_1[j], Intensity_1[j]));
            }

            for (int k = 0; k < RT_2.Length; k++)
            {
                list3.Add((RT_2[k], Intensity_2[k]));
            }

            list2 = list2.OrderByDescending(((double, double) i) => i.Item2).ToList();
            list3 = list3.OrderByDescending(((double, double) i) => i.Item2).ToList();
            foreach (var item in list3)
            {
                int num = 0;
                while (list2.Count > 0 && num < list2.Count)
                {
                    if (Math.Abs(list2[num].Item1 - item.Item1) / (list2[num].Item1 + item.Item1) < 0.01)
                    {
                        list.Add((list2[num].Item2, item.Item2));
                        list2.RemoveAt(num);
                        num = -1;
                        break;
                    }

                    num++;
                }

                if (list2.Count == 0)
                {
                    num++;
                }

                if (num > 0)
                {
                    list.Add((0.0, item.Item2));
                }
            }

            return list;
        }

        public static PeakCurve GetPeakCurve(MsDataScan[] scans, MsDataScan targetScan, double mz, CommonParameters commonParameters)
        {
            //EngineLayer.Peak newPeak = new EngineLayer.Peak(mz, targetScan.RetentionTime, targetScan.MassSpectrum.YArray[index]);
            PeakCurve newPeakCurve = new PeakCurve(new List<Peak> { }, targetScan.MsnOrder);

            foreach (var scan in scans)
            {
                try
                {
                    int indexByMz = scan.MassSpectrum.XArray.ToList().BinarySearch(mz);
                    if (indexByMz < 0)
                    {
                        indexByMz = ~indexByMz;
                    }
                    if (indexByMz > 0 && Math.Abs(scan.MassSpectrum.XArray[indexByMz] - mz) > Math.Abs(mz - scan.MassSpectrum.XArray[indexByMz - 1]))
                    {
                        indexByMz = indexByMz - 1;
                    }
                    if (commonParameters.PrecursorMassTolerance.Within(scan.MassSpectrum.XArray[indexByMz], mz))
                    {
                        EngineLayer.Peak newPeakForAdded = new EngineLayer.Peak(scan.MassSpectrum.XArray[indexByMz], scan.RetentionTime, scan.MassSpectrum.YArray[indexByMz], scan.MsnOrder);
                        newPeakCurve.Peaks.Add(newPeakForAdded);
                    }
                }
                catch
                {

                }

            }
            return newPeakCurve;
        }

        
        //For a given target Ms2WithSpecificMass Scan, extract all the peaks that come from the precursor ion by scoring the correlation
        //between the XIC of all the fragment ions in the scan and the XIC of the precursor ion
        public static Ms2ScanWithSpecificMass GroupPrecursorPeakAndFragmentIons(Ms2ScanWithSpecificMass targetScan, List<Ms2ScanWithSpecificMass> scans, CommonParameters commonParameters)
        {
            List<EngineLayer.Peak> DIA_peaks = new List<EngineLayer.Peak> { };
            PeakCurve prePeakCurve = targetScan.PrecursurPeak;
            //PeakCurve prePeakCurve = new PeakCurve(scans.Select(p => p.PrecursorMonoisotopicPeakMz).ToList().Average(), scans.Select(p => p.PrecursurPeak).ToList(), 1);
            List<PeakCurve> ms2PeakCurves = new List<PeakCurve>();
            for (int i = 0; i < targetScan.TheScan.MassSpectrum.XArray.Length; i++)
            {
                double RT = targetScan.Pre_RT;
                EngineLayer.Peak targetPeak = new EngineLayer.Peak(targetScan.TheScan.MassSpectrum.XArray[i], RT, targetScan.TheScan.MassSpectrum.YArray[i], targetScan.TheScan.MsnOrder);
                PeakCurve newPeakCurve = new PeakCurve(new List<EngineLayer.Peak> { }, 2, targetScan.TheScan.IsolationRange);

                foreach (var ms2scan in scans)
                {
                    int index = ms2scan.TheScan.MassSpectrum.XArray.ToList().BinarySearch(targetPeak.Mz);
                    if (index < 0)
                    {
                        index = ~index;
                    }
                    try
                    {
                        if (index > 0)
                        {
                            if (Math.Abs(ms2scan.TheScan.MassSpectrum.XArray[index] - targetPeak.Mz) > Math.Abs(targetPeak.Mz - ms2scan.TheScan.MassSpectrum.XArray[index - 1]))
                            {
                                index = index - 1;
                            }
                        }
                        if (commonParameters.ProductMassTolerance.Within(ms2scan.TheScan.MassSpectrum.XArray[index], targetPeak.Mz))
                        {
                            EngineLayer.Peak newPeak = new EngineLayer.Peak(ms2scan.TheScan.MassSpectrum.XArray[index], ms2scan.Pre_RT, ms2scan.TheScan.MassSpectrum.YArray[i], targetScan.TheScan.MsnOrder);
                            newPeakCurve.Peaks.Add(newPeak);
                        }
                    }
                    catch
                    {

                    }


                }

                double score = PeakCurve.CalPeakCorr(prePeakCurve, newPeakCurve);
                if (score > 0.5)
                {
                    DIA_peaks.Add(targetPeak);
                }
            }
            targetScan = new Ms2ScanWithSpecificMass(targetScan.TheScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath, 
                commonParameters, null, targetScan.Pre_RT, targetScan.PrecursurPeak, DIA_peaks);
            return targetScan;
        }

        //For each Ms2ScanWithSpecificMass, find the subset of Ms2ScanWithSpecificMass where the precursor associated with that scan is in the 
        //isolation window, create XICs based on the range and group the precursor ion with fragment ions
        //This method searches all copies of MS2 Scans in all MS2WithSpecificMass, which needs to be fixed
        public static List<Ms2ScanWithSpecificMass>[] GroupPrecursorPeaksAndFragmentIonsForDIA(List<Ms2ScanWithSpecificMass>[] scansWithPrecursors, CommonParameters commonParameters)
        {
            for (int i = 0; i < scansWithPrecursors.Length; i++)
            {
                if (scansWithPrecursors[i].Count > 0)
                {
                    for (int j = 0; j < scansWithPrecursors[i].Count; j++)
                    {
                        var target_scan = scansWithPrecursors[i][j];
                        List<Ms2ScanWithSpecificMass> Ms2ScanWithSpecificMz = new List<Ms2ScanWithSpecificMass>();
                        foreach (var scan in scansWithPrecursors.Where(p => p.Count > 0))
                        {
                            var sortedList = scan.OrderBy(p => p.PrecursorMonoisotopicPeakMz).ToList();
                            int indexByMz = sortedList.Select(p => p.PrecursorMonoisotopicPeakMz).ToList().BinarySearch(target_scan.PrecursorMonoisotopicPeakMz);
                            if (indexByMz < 0)
                            {
                                indexByMz = ~indexByMz;
                                if (indexByMz == sortedList.Count)
                                {
                                    indexByMz = indexByMz - 1;
                                }

                            }
                            try
                            {
                                if (Math.Abs(sortedList[indexByMz].PrecursorMonoisotopicPeakMz - target_scan.PrecursorMonoisotopicPeakMz) > Math.Abs(target_scan.PrecursorMonoisotopicPeakMz - sortedList[indexByMz - 1].PrecursorMonoisotopicPeakMz))
                                {
                                    indexByMz = indexByMz - 1;
                                }
                            }
                            catch
                            {

                            }

                            //}
                            try
                            {
                                if (commonParameters.PrecursorMassTolerance.Within(sortedList[indexByMz].PrecursorMonoisotopicPeakMz, target_scan.PrecursorMonoisotopicPeakMz))
                                {
                                    Ms2ScanWithSpecificMz.Add(sortedList[indexByMz]);
                                }
                            }
                            catch
                            {

                            }

                        }
                        
                        if (Ms2ScanWithSpecificMz.Count >= 5)
                        {
                            Ms2ScanWithSpecificMass newList = GroupPrecursorPeakAndFragmentIons(target_scan, Ms2ScanWithSpecificMz, commonParameters);
                            scansWithPrecursors[i][j] = newList;
                        }

                    }
                }
            }
            return scansWithPrecursors;
        }

        
        public static List<Ms2ScanWithSpecificMass>[] GroupPrecursorPeaksAndFragmentIonsISD(List<Ms2ScanWithSpecificMass>[] scansWithPrecursors, List<Peak> ms2Peaks, CommonParameters commonParameters)
        {
            for (int i = 0; i < scansWithPrecursors.Length; i++)
            {
                if (scansWithPrecursors[i].Count > 0)
                {
                    for (int j = 0; j < scansWithPrecursors[i].Count; j++)
                    {
                        var targetScan = scansWithPrecursors[i][j];
                        List<EngineLayer.Peak> DIA_peaks = new List<EngineLayer.Peak> { };
                        var ms2PeaksToMatch = ms2Peaks.Where(p => p.ScanNumber == targetScan.OneBasedScanNumber).ToList();
                        List<Ms2ScanWithSpecificMass> Ms2ScanWithSpecificMz = new List<Ms2ScanWithSpecificMass>();

                        foreach(Peak peak in ms2PeaksToMatch)
                        {
                            if(peak.XIC != null)
                            {
                                double score = PeakCurve.CalPeakCorr(targetScan.PrecursurPeak, peak.XIC);
                                if (score > 0.5)
                                {
                                    DIA_peaks.Add(peak);
                                }
                            }
                        }
                        DIA_peaks = DIA_peaks.OrderBy(p => p.Mz).ToList();
                        MzSpectrum diaSpectrum = new MzSpectrum(DIA_peaks.Select(p => p.Mz).ToArray(), DIA_peaks.Select(p => p.Intensity).ToArray(), false);
                        MsDataScan newScan = new MsDataScan(diaSpectrum, targetScan.OneBasedScanNumber, targetScan.TheScan.MsnOrder, targetScan.TheScan.IsCentroid,
                            targetScan.TheScan.Polarity, targetScan.TheScan.RetentionTime, targetScan.TheScan.ScanWindowRange, targetScan.TheScan.ScanFilter, targetScan.TheScan.MzAnalyzer,
                            targetScan.TheScan.TotalIonCurrent, targetScan.TheScan.InjectionTime, targetScan.TheScan.NoiseData, targetScan.TheScan.NativeId);
                        scansWithPrecursors[i][j] = new Ms2ScanWithSpecificMass(newScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath,
                commonParameters, null, targetScan.Pre_RT, targetScan.PrecursurPeak, DIA_peaks);
                    }
                }
            }
            return scansWithPrecursors;
        }

        //public static PeakCurve GetXIC(Peak targetPeak, List<Peak> allPeaks, CommonParameters commonParameters)
        //{
        //    PeakCurve peakCurve = new PeakCurve(new List<Peak> { }, targetPeak.MsLevel);
        //    foreach (MsDataScan scan in scans)
        //    {
        //        int index = scan.MassSpectrum.XArray.ToList().BinarySearch(targetPeak.Mz);
        //        if (index < 0)
        //        {
        //            index = ~index;
        //        }
        //        try
        //        {
        //            if (index > 0)
        //            {
        //                if (Math.Abs(scan.MassSpectrum.XArray[index] - targetPeak.Mz) > Math.Abs(targetPeak.Mz - scan.MassSpectrum.XArray[index - 1]))
        //                {
        //                    index = index - 1;
        //                }
        //            }
        //            if (commonParameters.ProductMassTolerance.Within(scan.MassSpectrum.XArray[index], targetPeak.Mz))
        //            {
        //                Peak newPeak = new Peak(scan.MassSpectrum.XArray[index], scan.RetentionTime, scan.MassSpectrum.YArray[index]);
        //                peakCurve.Peaks.Add(newPeak);
        //            }
        //        }
        //        catch
        //        {

        //        }

        //    }
        //    return peakCurve;
        //}
        public static int Match(double mz, List<double> list, Tolerance tolerance)
        {
            int index = list.BinarySearch(mz);
            if (index < 0)
            {
                index = ~index;
                if (index == list.Count)
                {
                    index = index - 1;
                }
                if (index != 0)
                {
                    if (Math.Abs(list[index] - mz) > Math.Abs(mz - list[index - 1]))
                    {
                        index = index - 1;
                    }
                }
            }
            if (tolerance.Within(mz, list[index]))
            {
                return index;
            }
            else
            {
                return -1;
            }
        }

        public static PeakCurve GetXIC(Peak targetPeak, List<Peak> allPeaks, CommonParameters commonParameters)
        {
            if (targetPeak.XIC != null) 
            { 
                return targetPeak.XIC; 
            }
            PeakCurve peakCurve = new PeakCurve(new List<Peak> { }, targetPeak.MsLevel);
            var groups = allPeaks.GroupBy(peak => peak.ScanNumber);
            foreach (var group in groups)
            {
                List<Peak> sortedList = group.ToList().OrderBy(p => p.Mz).ToList();
                List<double> sortedMzs = sortedList.Select(p => p.Mz).ToList();
                int index = PeakCurve.Match(targetPeak.Mz, sortedMzs, commonParameters.PrecursorMassTolerance);
                if (index > 0)
                {
                    peakCurve.Peaks.Add(sortedList[index]);
                    sortedList[index].XIC = peakCurve;
                }
            }
            return peakCurve;
        }

        public static PeakCurve GetXIC_noRemove(Peak targetPeak, List<Peak> allPeaks, CommonParameters commonParameters)
        {
            PeakCurve peakCurve = new PeakCurve(new List<Peak> { }, targetPeak.MsLevel);
            var groups = allPeaks.GroupBy(peak => peak.ScanNumber);
            foreach (var group in groups)
            {
                List<Peak> sortedList = group.ToList().OrderBy(p => p.Mz).ToList();
                List<double> sortedMzs = sortedList.Select(p => p.Mz).ToList();
                int index = PeakCurve.Match(targetPeak.Mz, sortedMzs, commonParameters.PrecursorMassTolerance);
                if (index > 0)
                {
                    peakCurve.Peaks.Add(sortedList[index]);
                }
            }
            return peakCurve;
        }

        //public static List<PeakCurve> GetXICForAll(CommonParameters commonParameters, List<Peak> peaks, int msnOrder)
        //{
        //    //EngineLayer.Peak newPeak = new EngineLayer.Peak(mz, targetScan.RetentionTime, targetScan.MassSpectrum.YArray[index]);
        //    List<PeakCurve> allXICs = new List<PeakCurve>();
        //    peaks.OrderBy(p => p.Index);
        //    List<int> indices = peaks.Select(p => p.Index).ToList();

        //    foreach (int i in indices)
        //    {
        //        if (peaks[i].XIC != null)
        //        {
        //            continue;
        //        }

        //        PeakCurve newPeakCurve = new PeakCurve(new List<Peak> { }, msnOrder);
        //        var listOfMzs = peaks.Where(p => indices.Contains(p.Index)).Select(p => p.Mz).ToList();
        //        try
        //        {
        //            int indexByMz = listOfMzs.BinarySearch(peaks[i].Mz);
        //            if (indexByMz < 0)
        //            {
        //                indexByMz = ~indexByMz;
        //            }
        //            if (indexByMz > 0 && Math.Abs(listOfMzs[indexByMz] - peaks[i].Mz) > Math.Abs(peaks[i].Mz - listOfMzs[indexByMz - 1]))
        //            {
        //                indexByMz = indexByMz - 1;
        //            }
        //            if (commonParameters.PrecursorMassTolerance.Within(listOfMzs[indexByMz], peaks[i].Mz))
        //            {
        //                newPeakCurve.Peaks.Add(peaks[indexByMz]);
        //                indices.Remove(peaks[indexByMz].Index);
        //            }
        //        }
        //        catch
        //        {

        //        }


        //    }
        //    return allXICs;
        //}
    }
}
