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
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;


namespace EngineLayer
{
    public class PeakCurve
    {
        public PeakCurve(double average_mz, List<Peak> peaks, int msLevel, MzRange isolationWindow = null)
        {
            Average_Mz = average_mz;
            Peaks = peaks;
            MsLevel = msLevel;
            IsolationWindow = isolationWindow;
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

        public static PeakCurve GetPeakCurve(MsDataScan[] scans, MsDataScan targetScan, double mz, double monomz, int index, CommonParameters commonParameters)
        {
            //EngineLayer.Peak newPeak = new EngineLayer.Peak(mz, targetScan.RetentionTime, targetScan.MassSpectrum.YArray[index]);
            PeakCurve newPeakCurve = new PeakCurve(mz, new List<EngineLayer.Peak> { }, targetScan.MsnOrder);

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
                        EngineLayer.Peak newPeakForAdded = new EngineLayer.Peak(scan.MassSpectrum.XArray[indexByMz], scan.RetentionTime, scan.MassSpectrum.YArray[indexByMz]);
                        newPeakCurve.Peaks.Add(newPeakForAdded);
                    }
                }
                catch
                {

                }

            }
            return newPeakCurve;
        }
    }
}

        //public static Ms2ScanWithSpecificMass GroupPrecursorPeakAndFragmentIons(Ms2ScanWithSpecificMass targetScan, List<Ms2ScanWithSpecificMass> scans, CommonParameters commonParameters)
        //{
        //    List<EngineLayer.Peak> DIA_peaks = new List<EngineLayer.Peak> { };
        //    PeakCurve prePeakCurve = targetScan.PrecursurPeak;
        //    //PeakCurve prePeakCurve = new PeakCurve(scans.Select(p => p.PrecursorMonoisotopicPeakMz).ToList().Average(), scans.Select(p => p.PrecursurPeak).ToList(), 1);
        //    List<PeakCurve> ms2PeakCurves = new List<PeakCurve>();
        //    for (int i = 0; i < targetScan.TheScan.MassSpectrum.XArray.Length; i++)
        //    {
        //        double RT = targetScan.Pre_RT;
        //        EngineLayer.Peak targetPeak = new EngineLayer.Peak(targetScan.TheScan.MassSpectrum.XArray[i], RT, targetScan.TheScan.MassSpectrum.YArray[i]);
        //        PeakCurve newPeakCurve = new PeakCurve(targetScan.TheScan.MassSpectrum.XArray[i], new List<EngineLayer.Peak> { }, 2, targetScan.TheScan.IsolationRange);

        //        foreach (var ms2scan in scans)
        //        {
        //            int index = ms2scan.TheScan.MassSpectrum.XArray.ToList().BinarySearch(targetPeak.Mz);
        //            if (index < 0)
        //            {
        //                index = ~index;
        //            }
        //            try
        //            {
        //                if (index > 0)
        //                {
        //                    if (Math.Abs(ms2scan.TheScan.MassSpectrum.XArray[index] - targetPeak.Mz) > Math.Abs(targetPeak.Mz - ms2scan.TheScan.MassSpectrum.XArray[index - 1]))
        //                    {
        //                        index = index - 1;
        //                    }
        //                }
        //                if (commonParameters.ProductMassTolerance.Within(ms2scan.TheScan.MassSpectrum.XArray[index], targetPeak.Mz))
        //                {
        //                    EngineLayer.Peak newPeak = new EngineLayer.Peak(ms2scan.TheScan.MassSpectrum.XArray[index], ms2scan.Pre_RT, ms2scan.TheScan.MassSpectrum.YArray[i]);
        //                    newPeakCurve.Peaks.Add(newPeak);
        //                }
        //            }
        //            catch
        //            {

        //            }


        //        }

        //        double score = PeakCurve.CalPeakCorr(prePeakCurve, newPeakCurve);
        //        if (score > 0.5)
        //        {
        //            DIA_peaks.Add(targetPeak);
        //        }
        //    }
        //    targetScan = new Ms2ScanWithSpecificMass(targetScan.TheScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath, commonParameters, null, targetScan.Pre_RT, targetScan.PrecursurPeak, DIA_peaks);
        //    return targetScan;
        //}

        //public static List<Ms2ScanWithSpecificMass> GroupPrecursorPeaksAndFragmentIonsForDIA(Dictionary<MzRange, List<List<Ms2ScanWithSpecificMass>>> dic, List<Ms2ScanWithSpecificMass> scansWithPrecursors, CommonParameters commonParameters)
        //{
        //    var allScanIntheRange = dic[scansWithPrecursors.First().TheScan.IsolationRange];

        //    for (int i = 0; i < scansWithPrecursors.Count; i++)
        //    {
        //        var ms2ScanWithSpecificMass = scansWithPrecursors[i];
        //        List<Ms2ScanWithSpecificMass> Ms2ScanWithSpecificMz = new List<Ms2ScanWithSpecificMass>();
        //        foreach (var list in allScanIntheRange)
        //        {
        //            var sortedList = list.OrderByDescending(p => p.PrecursorMonoisotopicPeakMz).ToList();
        //            int indexByMz = sortedList.Select(p => p.PrecursorMonoisotopicPeakMz).ToList().BinarySearch(ms2ScanWithSpecificMass.PrecursorMonoisotopicPeakMz);
        //            if (indexByMz < 0)
        //            {
        //                indexByMz = ~indexByMz;
        //            }
        //            if (Math.Abs(sortedList[indexByMz].PrecursorMonoisotopicPeakMz - ms2ScanWithSpecificMass.PrecursorMonoisotopicPeakMz) > Math.Abs(ms2ScanWithSpecificMass.PrecursorMonoisotopicPeakMz - sortedList[indexByMz].PrecursorMonoisotopicPeakMz))
        //            {
        //                indexByMz = indexByMz - 1;
        //            }


        //        }
        //    }
