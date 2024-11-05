using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common;
using MathNet.Numerics.Statistics;
using ThermoFisher.CommonCore.Data.Interfaces;

namespace EngineLayer.DIA
{
    public class PrecursorFragmentPair
    {
        public PrecursorFragmentPair(PeakCurve pre, PeakCurve frag, double corr)
        {
            PrecursorPeakCurve = pre;
            FragmentPeakCurve = frag;
            Correlation = corr;
        }

        public PeakCurve PrecursorPeakCurve { get; set; }
        public PeakCurve FragmentPeakCurve { get; set; }    
        public double Correlation { get; set; }
        public int FragmentRank { get; set; }
        public int PrecursorRank { get; set; }

        public static double CalculateRTOverlapRatio(PeakCurve curve1, PeakCurve curve2)
        {
            double overlap = 0;
            var ms1rtrange = curve1.EndRT - curve1.StartRT;
            var ms2rtrange = curve2.EndRT - curve2.StartRT;
            if (curve1.StartRT >= curve2.StartRT && curve1.StartRT <= curve2.EndRT && curve1.EndRT >= curve2.EndRT)
            {
                overlap = (curve2.EndRT - curve1.StartRT) / ms1rtrange;
            }
            else if (curve1.EndRT >= curve2.StartRT && curve1.EndRT <= curve2.EndRT && curve1.StartRT <= curve2.StartRT)
            {
                overlap = (curve1.EndRT - curve2.StartRT) / ms1rtrange;
            }
            else if (curve1.StartRT <= curve2.StartRT && curve1.EndRT >= curve2.EndRT)
            {
                overlap = ms2rtrange / ms1rtrange;
            }
            else if (curve1.StartRT >= curve2.StartRT && curve1.EndRT <= curve2.EndRT)
            {
                overlap = 1;
            }
            return overlap;
        }

        public static double CalculateRTOverlapRatio_scanCycle(PeakCurve curve1, PeakCurve curve2)
        {
            double overlap = 0;
            var ms1rtrange = curve1.EndCycle - curve1.StartCycle;
            var ms2rtrange = curve2.EndCycle - curve2.StartCycle;
            if (curve1.StartCycle >= curve2.StartCycle && curve1.StartCycle <= curve2.EndCycle && curve1.EndCycle >= curve2.EndCycle)
            {
                overlap = (curve2.EndCycle - curve1.StartCycle) / ms1rtrange;
            }
            else if (curve1.EndCycle >= curve2.StartRT && curve1.EndCycle <= curve2.EndCycle && curve1.StartCycle <= curve2.StartRT)
            {
                overlap = (curve1.EndCycle - curve2.StartCycle) / ms1rtrange;
            }
            else if (curve1.StartCycle <= curve2.StartCycle && curve1.EndCycle >= curve2.EndCycle)
            {
                overlap = ms2rtrange / ms1rtrange;
            }
            else if (curve1.StartCycle >= curve2.StartCycle && curve1.EndCycle <= curve2.EndCycle)
            {
                overlap = 1;
            }
            return overlap;
        }

        public static double CalculateCorr_spline(PeakCurve peakCurve1, PeakCurve peakCurve2, string splineType, double timeInterval)
        {
            if (peakCurve2.Peaks.Count < 5 || peakCurve1.Peaks.Count < 5)
            {
                return 0;
            }
            //change
            //var startRT = Math.Max(peakCurve1.StartRT, peakCurve2.StartRT);
            //var endRT = Math.Min(peakCurve1.EndRT, peakCurve2.EndRT);
            var startRT = peakCurve1.StartRT;
            var endRT = peakCurve1.EndRT;

            if (splineType == "Bspline")
            {
                var rtSeqB = new List<float>();
                float interval = (float)timeInterval;
                for (float i = (float)startRT; i < (float)endRT; i += interval)
                {
                    rtSeqB.Add(i);
                }
                int numPt = rtSeqB.Count;
                var smoothedIntensities1 = peakCurve1.GetBsplineData(numPt, 2).Select(p => (double)p.Item2).ToArray();
                var smoothedIntensities2 = peakCurve2.GetBsplineData(numPt, 2).Select(p => (double)p.Item2).ToArray();

                double corr = MathNet.Numerics.Statistics.Correlation.Pearson(smoothedIntensities1, smoothedIntensities2);
                return corr;

                //float sumA = 0;
                //float sumB = 0;
                //float sumA2 = 0;
                //float sumB2 = 0;
                //float sumAB = 0;

                //for (int k = 0; k < numPt; k++)
                //{
                //    float a = smoothedIntensities1[k];
                //    float b = smoothedIntensities2[k];

                //    sumA += a;
                //    sumB += b;
                //    sumA2 += a * a;
                //    sumB2 += b * b;
                //    sumAB += a * b;
                //}

                //float numerator = numPt * sumAB - sumA * sumB;
                //var denominator = Math.Sqrt((numPt * sumA2 - sumA * sumA) * (numPt * sumB2 - sumB * sumB));

                //if (denominator == 0)
                //{
                //    return 0; // Avoid division by zero
                //}

                //return numerator / denominator;
            }
            var rtSeq = new List<double>();
            for (double i = startRT; i < endRT; i += timeInterval)
            {
                rtSeq.Add(i);
            }
            var intensities1 = new List<double>();
            var intensities2 = new List<double>();

            if (splineType == "linear")
            {
                if (peakCurve1.LinearSpline == null)
                {
                    peakCurve1.GetLinearSpline();
                }
                if (peakCurve2.LinearSpline == null)
                {
                    peakCurve2.GetLinearSpline();
                }
                foreach (var rt in rtSeq)
                {
                    intensities1.Add(peakCurve1.LinearSpline.Interpolate(rt));
                    intensities2.Add(peakCurve2.LinearSpline.Interpolate(rt));
                }
                double corr = MathNet.Numerics.Statistics.Correlation.Pearson(intensities1, intensities2);
                return corr;
            }

            if (splineType == "cubic")
            {
                if (peakCurve1.CubicSpline == null)
                {
                    peakCurve1.GetCubicSpline();
                }
                if (peakCurve2.CubicSpline == null)
                {
                    peakCurve2.GetCubicSpline();
                }
                foreach (var rt in rtSeq)
                {
                    intensities1.Add(peakCurve1.CubicSpline.Interpolate(rt));
                    intensities2.Add(peakCurve2.CubicSpline.Interpolate(rt));
                }
                double corr = MathNet.Numerics.Statistics.Correlation.Pearson(intensities1, intensities2);
                return corr;
            }
            return 0;
        }

        public static double CalculateCorr_diaUmpire(PeakCurve curve1, PeakCurve curve2, float timeInterval)
        {
            float start = Math.Max(curve1.SmoothedData[0].Item1, curve2.SmoothedData[0].Item1);
            float end = Math.Min(curve1.SmoothedData[curve1.SmoothedData.Count - 1].Item1, curve2.SmoothedData[curve2.SmoothedData.Count - 1].Item1);

            //int num = Math.Min(curve1.SmoothedData.Count, curve2.SmoothedData.Count) / 2;
            int num = (int)((end - start) / timeInterval + 1);
            if (num < 6)
            {
                return 0f;
            }
            float[] arrayA = new float[num];
            float[] arrayB = new float[num];

            int i = 0;
            float low = start;
            float up = start + timeInterval;

            for (int j = 0; j < curve1.SmoothedData.Count; j++)
            {
                while (curve1.SmoothedData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (curve1.SmoothedData[j].Item1 >= low && curve1.SmoothedData[j].Item1 < up)
                {
                    if (curve1.SmoothedData[j].Item2 > arrayA[i])
                    {
                        arrayA[i] = curve1.SmoothedData[j].Item2;
                    }
                }
            }

            i = 0;
            low = start;
            up = start + timeInterval;

            for (int j = 0; j < curve2.SmoothedData.Count; j++)
            {
                while (curve2.SmoothedData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (curve2.SmoothedData[j].Item1 >= low && curve2.SmoothedData[j].Item1 < up)
                {
                    if (curve2.SmoothedData[j].Item2 > arrayB[i])
                    {
                        arrayB[i] = curve2.SmoothedData[j].Item2;
                    }
                }
            }

            for (int idx = 1; idx < num - 1; idx++)
            {
                if (arrayA[idx] == 0f)
                {
                    arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
                }
                if (arrayB[idx] == 0f)
                {
                    arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
                }
            }

            float sumA = 0;
            float sumB = 0;
            float sumA2 = 0;
            float sumB2 = 0;
            float sumAB = 0;

            for (int k = 0; k < num; k++)
            {
                float a = arrayA[k];
                float b = arrayB[k];

                sumA += a;
                sumB += b;
                sumA2 += a * a;
                sumB2 += b * b;
                sumAB += a * b;
            }

            float numerator = num * sumAB - sumA * sumB;
            var denominator = Math.Sqrt((num * sumA2 - sumA * sumA) * (num * sumB2 - sumB * sumB));

            if (denominator == 0)
            {
                return 0; // Avoid division by zero
            }

            return numerator / denominator;
        }

        public static double CalculatePeakCurveCorr(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            var peakList1 = peakCurve1.Peaks.ToArray();
            var peakList2 = peakCurve2.Peaks.ToArray();
            var scanCycles1 = peakList1.Select(p => p.ZeroBasedScanIndex).ToArray();
            var scanCycles2 = peakList2.Select(p => p.ZeroBasedScanIndex).ToArray();

            var start = (int)Math.Max(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var end = (int)Math.Min(peakCurve1.EndCycle, peakCurve2.EndCycle);
            int numPoints = end - start  + 1;
            if (numPoints < 6)
            {
                return 0;
            }

            var list1 = new List<double>();
            var list2 = new List<double>();

            for (int i = start; i <= end; i++)
            {
                var index1 = Array.BinarySearch(scanCycles1, i);
                var index2 = Array.BinarySearch(scanCycles2, i);
                if (index1 < 0 || index2 < 0)
                {
                    continue;
                }
                else
                {
                    list1.Add(peakList1[index1].Intensity);
                    list2.Add(peakList2[index2].Intensity);
                }
            }
            if (list1.Count < 6)
            {
                return 0;
            }
            double corr = MathNet.Numerics.Statistics.Correlation.Pearson(list1, list2);
            return corr;
        }

        public static double CalculateCorr_spline_scanCycle(PeakCurve peakCurve1, PeakCurve peakCurve2, double interval)
        {
            var start = Math.Max(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var end = Math.Min(peakCurve1.EndCycle, peakCurve2.EndCycle);
            var numPoints = (int)((end - start) / interval) + 1;
            if (numPoints < 6)
            {
                return 0;
            }

            var effectivePoints1 = peakCurve1.ScanCycleSoomthedData.Where(p => p.Item1 > start - interval && p.Item1 < end + interval).ToArray();
            var effectivePoints2 = peakCurve2.ScanCycleSoomthedData.Where(p => p.Item1 > start - interval && p.Item1 < end + interval).ToArray();

            var list1 = new double[numPoints];
            var list2 = new double[numPoints];
            for(int i = 0; i< numPoints; i++)
            {
                list1[i] = effectivePoints1[i].Item2;
                list2[i] = effectivePoints2[i].Item2;
            }
            double corr = MathNet.Numerics.Statistics.Correlation.Pearson(list1, list2);

            return corr;
        }
    }
}
