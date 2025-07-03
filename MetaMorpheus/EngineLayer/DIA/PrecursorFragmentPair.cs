using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using Easy.Common;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Integration;
using MathNet.Numerics.Statistics;
using Plotly.NET;
using Readers;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.Interfaces;
using static Python.Runtime.TypeSpec;

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

        public PrecursorFragmentPair(PeakCurve pre, PeakCurve frag)
        {
            PrecursorPeakCurve = pre;
            FragmentPeakCurve = frag;
        }
        public PrecursorFragmentPair(PeakCurve pre, PeakCurve frag, double overlap, double corr, double sharedXIC = 0, string label = null)
        {
            PrecursorPeakCurve = pre;
            FragmentPeakCurve = frag;
            Overlap = overlap;
            Correlation = corr;
            SharedXIC = sharedXIC;
            Label = label;
        }
        public PeakCurve PrecursorPeakCurve { get; set; }
        public PeakCurve FragmentPeakCurve { get; set; }    
        public double Correlation { get; set; }
        public int FragmentRank { get; set; }
        public int PrecursorRank { get; set; }
        public double Overlap { get; set; }
        public double SharedXIC { get; set; } 
        public string Label { get; set; } 
        public int NormalizedIntensityRank { get; set; }

        public static double CalculateCorrelation((double, double)[] xy1, (double, double)[] xy2)
        {
            if (xy1 == null || xy2 == null)
            {
                return 0;
            }

            double start = Math.Max(xy1[0].Item1, xy2[0].Item1);
            double end = Math.Min(xy1[xy1.Length - 1].Item1, xy2[xy2.Length - 1].Item1);

            var validxy1 = xy1.Where(p => p.Item1 >= start && p.Item1 <= end).ToArray();
            var validxy2 = xy2.Where(p => p.Item1 >= start && p.Item1 <= end).ToArray();
            int numPoints = Math.Min(validxy1.Length, validxy2.Length);
            if (numPoints < 3)
            {
                return double.NaN;
            }
            var xy = validxy1.Take(numPoints).Zip(validxy2.Take(numPoints), (a, b) => (a.Item2, b.Item2)).ToArray();
            var y1 = xy.Select(p => p.Item1).ToArray();
            var y2 = xy.Select(p => p.Item2).ToArray();
            double corr = MathNet.Numerics.Statistics.Correlation.Pearson(y1, y2);

            return corr;
        }

        public static double CalculatePeakCurveCorrXYData(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            double corr = CalculateCorrelation(peakCurve1.XYData, peakCurve2.XYData);
            return corr;
        }

        public static double CalculatePeakCurveCorrXYData_Umpire(PeakCurve curve1, PeakCurve curve2, int NoPointPerInterval)
        {
            double start = Math.Max(curve1.UmpireBsplineData[0].Item1, curve2.UmpireBsplineData[0].Item1);

            //int num = Math.Min(curve1.SmoothedData.Count, curve2.SmoothedData.Count) / 2;
            int num = Math.Max(curve1.UmpireBsplineData.Count(), curve2.UmpireBsplineData.Count()) / 2;
            double timeInterval = (double) 2 / NoPointPerInterval;

            if (num < 6)
            {
                return 0f;
            }
            double[] arrayA = new double[num];
            double[] arrayB = new double[num];

            int i = 0;
            double low = start;
            double up = start + timeInterval;

            for (int j = 0; j < curve1.UmpireBsplineData.Count; j++)
            {
                while (curve1.UmpireBsplineData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (curve1.UmpireBsplineData[j].Item1 >= low && curve1.UmpireBsplineData[j].Item1 < up)
                {
                    if (curve1.UmpireBsplineData[j].Item2 > arrayA[i])
                    {
                        arrayA[i] = curve1.UmpireBsplineData[j].Item2;
                    }
                }
            }

            i = 0;
            low = start;
            up = start + timeInterval;

            for (int j = 0; j < curve2.UmpireBsplineData.Count; j++)
            {
                while (curve2.UmpireBsplineData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (curve2.UmpireBsplineData[j].Item1 >= low && curve2.UmpireBsplineData[j].Item1 < up)
                {
                    if (curve2.UmpireBsplineData[j].Item2 > arrayB[i])
                    {
                        arrayB[i] = curve2.UmpireBsplineData[j].Item2;
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

            double corr = MathNet.Numerics.Statistics.Correlation.Pearson(arrayA, arrayB);
            return corr;
        }

        //public static double CalculateRTOverlapRatio(PeakCurve curve1, PeakCurve curve2)
        //{
        //    double overlap = 0;
        //    var ms1rtrange = curve1.EndRT - curve1.StartRT;
        //    var ms2rtrange = curve2.EndRT - curve2.StartRT;
        //    if (curve1.StartRT >= curve2.StartRT && curve1.StartRT <= curve2.EndRT && curve1.EndRT >= curve2.EndRT)
        //    {
        //        overlap = (curve2.EndRT - curve1.StartRT) / ms1rtrange;
        //    }
        //    else if (curve1.EndRT >= curve2.StartRT && curve1.EndRT <= curve2.EndRT && curve1.StartRT <= curve2.StartRT)
        //    {
        //        overlap = (curve1.EndRT - curve2.StartRT) / ms1rtrange;
        //    }
        //    else if (curve1.StartRT <= curve2.StartRT && curve1.EndRT >= curve2.EndRT)
        //    {
        //        overlap = ms2rtrange / ms1rtrange;
        //    }
        //    else if (curve1.StartRT >= curve2.StartRT && curve1.EndRT <= curve2.EndRT)
        //    {
        //        overlap = 1;
        //    }
        //    return overlap;
        //}

        public static double CalculateRTOverlapRatio(PeakCurve curve1, PeakCurve curve2)
        {
            if (curve1.StartCycle <= curve2.StartCycle && curve1.EndCycle >= curve2.EndCycle)
            {
                return 1;
            }
            int maxStart = Math.Max(curve1.StartCycle, curve2.StartCycle);
            int minEnd = Math.Min(curve1.EndCycle, curve2.EndCycle);
            int minStart = Math.Min(curve1.StartCycle, curve2.StartCycle);
            int maxEnd = Math.Max(curve1.EndCycle, curve2.EndCycle);
            return (maxEnd - minStart) == 0 ? 0 : (minEnd - maxStart) / (double)(maxEnd - minStart);
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
            else if (curve1.EndCycle >= curve2.StartCycle && curve1.EndCycle <= curve2.EndCycle && curve1.StartCycle <= curve2.StartCycle)
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

        public static double CalculateRTOverlapRatio_scanCycle2(PeakCurve curve1, PeakCurve curve2)
        {
            var start = Math.Min(curve1.StartCycle, curve2.StartCycle);
            var end = Math.Max(curve1.EndCycle, curve2.EndCycle);
            var overlapStart = Math.Max(curve1.StartCycle, curve2.StartCycle);
            var overlapEnd = Math.Min(curve1.EndCycle, curve2.EndCycle);

            double overlap = (overlapEnd - overlapStart)/(end - start);
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
            //var startRT = peakCurve1.StartRT;
            //var endRT = peakCurve1.EndRT;

            if (splineType == "Bspline")
            {
                var rtSeqB = new List<float>();
                float interval = (float)timeInterval;
                for (float i = (float)startRT; i < (float)endRT; i += interval)
                {
                    rtSeqB.Add(i);
                }
                int numPt = rtSeqB.Count;
                var smoothedIntensities1 = peakCurve1.GetUmpireBSplineFloatData(numPt, 2).Select(p => (double)p.Item2).ToArray();
                var smoothedIntensities2 = peakCurve2.GetUmpireBSplineFloatData(numPt, 2).Select(p => (double)p.Item2).ToArray();

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
                foreach (var rt in rtSeq)
                {
                    intensities1.Add(peakCurve1.CubicSpline.Interpolate(rt));
                    intensities2.Add(peakCurve2.CubicSpline.Interpolate(rt));
                }
                double corr = MathNet.Numerics.Statistics.Correlation.Pearson(intensities1, intensities2);
                return corr;
            }

            if (splineType == "ms1cubic")
            {
                foreach (var rt in rtSeq)
                {
                    intensities1.Add(peakCurve1.CubicSpline.Interpolate(rt));
                    intensities2.Add(peakCurve2.Ms1SpaceSpline.Interpolate(rt));
                }
                double corr = MathNet.Numerics.Statistics.Correlation.Pearson(intensities1, intensities2);
                return corr;
            }
            return 0;
        }


        public static double CalculateCorr_diaUmpire(PeakCurve curve1, PeakCurve curve2, int NoPointPerInterval)
        {
            float start = Math.Max(curve1.BsplineSmoothedData[0].Item1, curve2.BsplineSmoothedData[0].Item1);
            float end = Math.Min(curve1.BsplineSmoothedData[curve1.BsplineSmoothedData.Count - 1].Item1, curve2.BsplineSmoothedData[curve2.BsplineSmoothedData.Count - 1].Item1);

            //int num = Math.Min(curve1.SmoothedData.Count, curve2.SmoothedData.Count) / 2;
            int num = Math.Max(curve1.BsplineSmoothedData.Count(), curve2.BsplineSmoothedData.Count()) / 2;
            float timeInterval = 2f / (float)NoPointPerInterval;

            if (num < 6)
            {
                return 0f;
            }
            float[] arrayA = new float[num];
            float[] arrayB = new float[num];

            int i = 0;
            float low = start;
            float up = start + timeInterval;

            for (int j = 0; j < curve1.BsplineSmoothedData.Count; j++)
            {
                while (curve1.BsplineSmoothedData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (curve1.BsplineSmoothedData[j].Item1 >= low && curve1.BsplineSmoothedData[j].Item1 < up)
                {
                    if (curve1.BsplineSmoothedData[j].Item2 > arrayA[i])
                    {
                        arrayA[i] = curve1.BsplineSmoothedData[j].Item2;
                    }
                }
            }

            i = 0;
            low = start;
            up = start + timeInterval;

            for (int j = 0; j < curve2.BsplineSmoothedData.Count; j++)
            {
                while (curve2.BsplineSmoothedData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (curve2.BsplineSmoothedData[j].Item1 >= low && curve2.BsplineSmoothedData[j].Item1 < up)
                {
                    if (curve2.BsplineSmoothedData[j].Item2 > arrayB[i])
                    {
                        arrayB[i] = curve2.BsplineSmoothedData[j].Item2;
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


            return 0;
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

        public static double CalculateCorr_scanCycleSpline_preCalculated(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            double interval = peakCurve1.ScanCycleSmoothedData[1].Item1 - peakCurve1.ScanCycleSmoothedData[0].Item1;
            var start = Math.Max(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var end = Math.Min(peakCurve1.EndCycle, peakCurve2.EndCycle);
            var numPoints = (int)((end - start) / interval) + 1;
            if (numPoints < 6)
            {
                return 0;
            }

            var effectivePoints1 = peakCurve1.ScanCycleSmoothedData.Where(p => p.Item1 > start - interval && p.Item1 < end + interval).ToArray();
            var effectivePoints2 = peakCurve2.ScanCycleSmoothedData.Where(p => p.Item1 > start - interval && p.Item1 < end + interval).ToArray();
            numPoints = Math.Min(effectivePoints1.Length, effectivePoints2.Length);

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


        public static double CalculateOverlapAreaRatio(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            var start = Math.Min(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var end = Math.Max(peakCurve1.EndCycle, peakCurve2.EndCycle);
            var overlapStart = Math.Max(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var overlapEnd = Math.Min(peakCurve1.EndCycle, peakCurve2.EndCycle);

            if (peakCurve1.NormalizedPeaks == null)
            {
                peakCurve1.GetNormalizedPeaks();
            }
            if (peakCurve2.NormalizedPeaks == null)
            {
                peakCurve2.GetNormalizedPeaks();
            }

            var overlapArea = new List<(double, double)>();
            var totalArea = new List<(double, double)>();
            var scanCycles1 = peakCurve1.NormalizedPeaks.Select(p => p.Item1).ToArray();
            var scanCycles2 = peakCurve2.NormalizedPeaks.Select(p => p.Item1).ToArray();

            for (int i = start; i <= end; i++)
            {
                var index1 = Array.BinarySearch(scanCycles1, i);
                var index2 = Array.BinarySearch(scanCycles2, i);

                if (index1 < 0 || index2 < 0)
                {
                    overlapArea.Add((i, 0));
                    if (index1 >= 0)
                    {
                        totalArea.Add((i, peakCurve1.NormalizedPeaks[index1].Item2));
                    } 
                    if (index2 >= 0)
                    {
                        totalArea.Add((i, peakCurve2.NormalizedPeaks[index2].Item2));
                    }
                }
                else
                {
                    totalArea.Add((i, Math.Max(peakCurve1.NormalizedPeaks[index1].Item2, peakCurve2.NormalizedPeaks[index2].Item2)));
                    overlapArea.Add((i, Math.Min(peakCurve1.NormalizedPeaks[index1].Item2, peakCurve2.NormalizedPeaks[index2].Item2)));
                }
            }
            double overlapAUC = CalculateArea(overlapArea);
            double totalAUC = CalculateArea(totalArea);
            double ratio = overlapAUC / totalAUC;
            return ratio;
        }

        public static double CalculateSharedXIC(PeakCurve peakCurve1, PeakCurve peakCurve2, bool visualize = false)
        {
            if (peakCurve1.EndCycle <= peakCurve2.StartCycle || peakCurve2.EndCycle <= peakCurve1.StartCycle)
            {
                return 0;
            }

            var overlapStart = Math.Max(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var overlapEnd = Math.Min(peakCurve1.EndCycle, peakCurve2.EndCycle);
            var maxLength = overlapEnd - overlapStart + 1;

            var overlapArea = new List<(double, double)>();
            var scanCycles1 = peakCurve1.NormalizedLinearSplinePeaks.Select(p => p.Item1).ToArray();
            var scanCycles2 = peakCurve2.NormalizedLinearSplinePeaks.Select(p => p.Item1).ToArray();
            var index1 = Array.BinarySearch(scanCycles1, overlapStart);
            var index2 = Array.BinarySearch(scanCycles2, overlapStart);

            for (int i = 0; i < maxLength - 1; i++)
            {
                double diff0 = peakCurve1.NormalizedLinearSplinePeaks[index1 + i].Item2 - peakCurve2.NormalizedLinearSplinePeaks[index2 + i].Item2;
                double diff1 = peakCurve1.NormalizedLinearSplinePeaks[index1 + i + 1].Item2 - peakCurve2.NormalizedLinearSplinePeaks[index2 + i + 1].Item2;
                
                overlapArea.Add((overlapStart +i, Math.Min(peakCurve1.NormalizedLinearSplinePeaks[index1 + i].Item2, peakCurve2.NormalizedLinearSplinePeaks[index2 + i].Item2))); 
                if (diff0 * diff1 < 0)
                {
                    double slope = peakCurve1.NormalizedLinearSplinePeaks[index1 + i + 1].Item2 - peakCurve1.NormalizedLinearSplinePeaks[index1 + i].Item2;
                    double y = peakCurve1.NormalizedLinearSplinePeaks[index1 + i].Item2 + slope * Math.Abs(diff0 / (diff1 - diff0));
                    overlapArea.Add((overlapStart + i + Math.Abs(diff0 / (diff1 - diff0)), y));
                }
            }
            double overlapAUC = CalculateNormalizedArea(overlapArea);

            if (visualize)
            {
                var overlapPlot = Chart2D.Chart.Line<double, double, string>(
                    x: overlapArea.Select(p => p.Item1),
                    y: overlapArea.Select(p => p.Item2)).WithMarkerStyle(Color: Color.fromString("green"));
                var plot1 = Chart2D.Chart.Line<double, double, string>(
                    x: peakCurve1.NormalizedLinearSplinePeaks.Select(p => (double)p.Item1),
                    y: peakCurve1.NormalizedLinearSplinePeaks.Select(p => p.Item2)).WithMarkerStyle(Color: Color.fromString("blue"));
                var plot2 = Chart2D.Chart.Line<double, double, string>(
                    x: peakCurve2.NormalizedLinearSplinePeaks.Select(p => (double)p.Item1),
                    y: peakCurve2.NormalizedLinearSplinePeaks.Select(p => p.Item2)).WithMarkerStyle(Color: Color.fromString("red"));
                var combinedPlot = Chart.Combine(new[] { overlapPlot, plot1, plot2 });
                combinedPlot.Show();
            }
            return overlapAUC;
        }


        public static double CalculateOverlapAreaRatio_spline(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            var start = Math.Min(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var end = Math.Max(peakCurve1.EndCycle, peakCurve2.EndCycle);
            var overlapStart = Math.Max(peakCurve1.StartCycle, peakCurve2.StartCycle);
            var overlapEnd = Math.Min(peakCurve1.EndCycle, peakCurve2.EndCycle);

            if (peakCurve1.NormalizedPeaks == null)
            {
                peakCurve1.GetNormalizedPeaks();
            }
            if (peakCurve2.NormalizedPeaks == null)
            {
                peakCurve2.GetNormalizedPeaks();
            }

            var overlapArea = new List<(double, double)>();
            var totalArea = new List<(double, double)>();
            var scanCycles1 = peakCurve1.NormalizedPeaks.Select(p => p.Item1).ToArray();
            var scanCycles2 = peakCurve2.NormalizedPeaks.Select(p => p.Item1).ToArray();

            for (int i = start; i <= end; i++)
            {
                var index1 = Array.BinarySearch(scanCycles1, i);
                var index2 = Array.BinarySearch(scanCycles2, i);

                if (index1 < 0 || index2 < 0)
                {
                    overlapArea.Add((i, 0));
                    if (index1 >= 0)
                    {
                        totalArea.Add((i, peakCurve1.NormalizedPeaks[index1].Item2));
                    }
                    if (index2 >= 0)
                    {
                        totalArea.Add((i, peakCurve2.NormalizedPeaks[index2].Item2));
                    }
                }
                else
                {
                    totalArea.Add(((double)i, Math.Max(peakCurve1.NormalizedPeaks[index1].Item2, peakCurve2.NormalizedPeaks[index2].Item2)));
                    overlapArea.Add(((double)i, Math.Min(peakCurve1.NormalizedPeaks[index1].Item2, peakCurve2.NormalizedPeaks[index2].Item2)));
                }
            }
            double overlapAUC = CalculateArea(overlapArea);
            double totalAUC = CalculateArea(totalArea);
            double ratio = overlapAUC / totalAUC;
            return ratio;
        }

        public static double CalculateArea(List<(double, double)> data)
        {
            double area = 0;
            for (int i = 1; i < data.Count; i++)
            {
                double x1 = data[i - 1].Item1;
                double x2 = data[i].Item1;
                double y1 = data[i - 1].Item2;
                double y2 = data[i].Item2;
                area += (x2 - x1) * (y1 + y2) / 2;
            }
            return area;
        }

        public static double CalculateNormalizedArea(List<(double, double)> data)
        {
            double area = data[0].Item2/2;
            for (int i = 1; i < data.Count; i++)
            {
                double x1 = data[i - 1].Item1;
                double x2 = data[i].Item1;
                double y1 = data[i - 1].Item2;
                double y2 = data[i].Item2;
                area += (x2 - x1) * (y1 + y2) / 2;
            }
            area += data[data.Count - 1].Item2 / 2;
            return area;
        }

        public static double CalculateOverlapRatioGaussian(Ms1Feature feature1, Ms1Feature feature2)
        {
            double start = Math.Min(feature1.RetentionTimeBegin, feature2.RetentionTimeBegin);
            double end = Math.Max(feature1.RetentionTimeEnd, feature2.RetentionTimeEnd);
            double overlapStart = Math.Max(feature1.RetentionTimeBegin, feature2.RetentionTimeBegin);
            double overlapEnd = Math.Min(feature1.RetentionTimeEnd, feature2.RetentionTimeEnd);

            double sd1 = (feature1.RetentionTimeEnd - feature1.RetentionTimeBegin)/4;
            double sd2 = (feature2.RetentionTimeEnd - feature2.RetentionTimeBegin)/4;
            double mean1 = feature1.RetentionTimeApex;
            double mean2 = feature2.RetentionTimeApex;
            Func<double, double> gaussian1 = x =>
            (1 / (Math.Sqrt(2 * Math.PI) * sd1)) * Math.Exp(-Math.Pow(x - mean1, 2) / (2 * Math.Pow(sd1, 2)));
            Func<double, double> gaussian2 = x =>
            (1 / (Math.Sqrt(2 * Math.PI) * sd2)) * Math.Exp(-Math.Pow(x - mean2, 2) / (2 * Math.Pow(sd2, 2)));

            Func<double, double> minPDF = x => Math.Min(gaussian1(x), gaussian2(x));
            Func<double, double> maxPDF = x => Math.Max(gaussian1(x), gaussian2(x));

            var overlapArea = GaussLegendreRule.Integrate(minPDF, overlapStart, overlapEnd, 1000);
            var totalArea = GaussLegendreRule.Integrate(maxPDF, start, end, 1000);
            double ratio = overlapArea/totalArea;

            return ratio;
        }
    }
}
