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

        public static double CalculateCorr_spline(PeakCurve peakCurve1, PeakCurve peakCurve2, string splineType, double timeInterval)
        {
            if (peakCurve2.Peaks.Count < 5 || peakCurve1.Peaks.Count < 5)
            {
                return 0;
            }
            var startRT = Math.Max(peakCurve1.StartRT, peakCurve2.StartRT);
            var endRT = Math.Min(peakCurve1.EndRT, peakCurve2.EndRT);

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

            var intensityPair = new List<(double, double)>();
            //for plot
            var rtPair = new List<(double, double)>();

            if (peakList1.Length < peakList2.Length)
            {
                foreach (var peak in peakList1)
                {
                    var indexList = peakList2.Select(p => p.ZeroBasedScanIndex).ToArray();
                    var index = Array.BinarySearch(indexList, peak.ZeroBasedScanIndex);
                    if (index > 0)
                    {
                        intensityPair.Add((peak.Intensity, peakList2[index].Intensity));
                        rtPair.Add((peak.RetentionTime, peakList2[index].RetentionTime));
                    }
                }
            }
            else
            {
                foreach (var peak in peakList2)
                {
                    var indexList = peakList1.Select(p => p.ZeroBasedScanIndex).ToArray();
                    var index = Array.BinarySearch(indexList, peak.ZeroBasedScanIndex);
                    if (index > 0)
                    {
                        intensityPair.Add((peakList1[index].Intensity, peak.Intensity));
                        rtPair.Add((peakList1[index].RetentionTime, peak.RetentionTime));
                    }
                }
            }

            if (intensityPair.Count >= 5)
            {
                double corr = MathNet.Numerics.Statistics.Correlation.Pearson(intensityPair.Select(pair => pair.Item1), intensityPair.Select(pair => pair.Item2));

                //plot
                //if (intensityPair.Count >= 5 && corr > 0.7)
                //{
                //    var plot1 = Chart2D.Chart.Point<double, double, string>(
                //    x: rtPair.Select(pair => pair.Item1),
                //    y: intensityPair.Select(pair => pair.Item1/intensityPair.Sum(p => p.Item1))).WithTraceInfo("precursor").WithMarkerStyle(Color: Color.fromString("blue"));
                //    var plot2 = Chart2D.Chart.Point<double, double, string>(
                //        x: rtPair.Select(pair => pair.Item2),
                //        y: intensityPair.Select(pair => pair.Item2 / intensityPair.Sum(p => p.Item2))).WithTraceInfo("fragment").WithMarkerStyle(Color: Color.fromString("red"));
                //    var combinedPlot = Chart.Combine(new[] { plot1, plot2 }).WithTitle($"corr {corr}");
                //    combinedPlot.Show();
                //}

                return corr;
            }
            else
            {
                return 0;
            }
        }
    }
}
