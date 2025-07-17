using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using Plotly.NET;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.Statistics;
using MathNet.Numerics.Distributions;
using TopDownProteomics;
using NWaves.Signals;
using NWaves.Filters.Polyphase;
using NWaves.Filters;
using Numerics.NET.Statistics;
using System.Xml.Linq;
using EngineLayer.DIA.CWT;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Optimization;
using static System.Runtime.InteropServices.JavaScript.JSType;
using static Python.Runtime.TypeSpec;

namespace EngineLayer.DIA
{
    public class PeakCurve
    {
        public PeakCurve(List<Peak> peaks, int msLevel, MzRange isolationRange, double mass = double.NaN, int charge = 0, double startMz = 0, double endMz = 0, int index = 0)
        {
            Peaks = peaks;
            MsLevel = msLevel;
            IsolationRange = isolationRange;
            MonoisotopicMass = mass;
            Charge = charge;
            StartMz = startMz;
            EndMz = endMz;
            Index = index;
            PFpairs = new List<PrecursorFragmentPair>();
            CwtParameters = new CwtParameters(2f, 150, 0.1f, 0.3f);
        }

        public PeakCurve()
        {
            Peaks = new List<Peak>();
            PFpairs = new List<PrecursorFragmentPair>();
        }

        public PeakCurve(List<Peak> peaks)
        {
            Peaks = peaks;
            MsLevel = peaks.First().MsLevel;
        }

        public List<Peak> Peaks { get; set; }
        public int MsLevel {  get; set; }
        public MzRange IsolationRange { get; set; }
        public double MonoisotopicMass { get; set; }
        public int Charge { get; set; }
        public double StartRT => Peaks.Select(p => p.RetentionTime).OrderBy(t => t).First();
        public double EndRT => Peaks.Select(p => p.RetentionTime).OrderByDescending(t => t).First();
        public int StartCycle => Peaks.Select(p => p.ZeroBasedScanIndex).OrderBy(t => t).First();
        public int EndCycle => Peaks.Select(p => p.ZeroBasedScanIndex).OrderByDescending(t => t).First();
        public double StartMz {  get; set; }
        public double EndMz { get; set; }
        public MzRange MzRange => new MzRange(StartMz, EndMz);
        public double ApexRT => Peaks.OrderByDescending(p => p.Intensity).First().RetentionTime;
        public int ApexCycle => Peaks.OrderByDescending(p => p.Intensity).First().ZeroBasedScanIndex;
        public double ApexSN => Peaks.OrderByDescending(p => p.Intensity).First().SN;
        public double ApexIntensity => Peaks.Max(p => p.Intensity);
        public double TotalIntensity => Peaks.Sum(p => p.Intensity);
        public virtual double AveragedMz => AverageMz();
        public double AveragedIntensity => AverageIntensity();
        public LinearSpline LinearSpline { get; set; }
        public CubicSpline CubicSpline { get; set; }    
        public CubicSpline Ms1SpaceSpline { get; set; }
        public int Index {  get; set; }
        public List<PrecursorFragmentPair> PFpairs { get; set; }
        public IsotopicEnvelope Envelope { get; set; }

        public List<PeakRidge> PeakRidgeList { get; set; }
        public List<(float rt1, float rt2, float rt3)> PeakRegionList { get; set; }
        public List<List<float>> NoRidgeRegion { get; set; }

        public WaveletMassDetector WaveletMassDetector { get; set; }
        public CwtParameters CwtParameters { get; set; }
        public List<(float, float)> BsplineSmoothedData { get; set; }
        public List<(double, double)> UmpireBsplineData { get; set; }
        public List<(double, double)> ScanCycleSmoothedData { get; set; }
        public (double, double)[] NormalizedPeaks { get; set; }
        public double NL { get; set; }
        public (double, double)[] XYData { get; set; }
        public (int, double)[] NormalizedLinearSplinePeaks { get; set; }

        public double AveragedMass => AverageMass();

        public virtual double AverageMz()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            double averagedMz = 0;
            foreach (var peak in Peaks)
            {
                double weight = peak.Intensity / sumIntensity;
                averagedMz += weight * peak.Mz;
            }
            return averagedMz;
        }

        public double AverageIntensity()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            double averagedIntensity = 0;
            foreach (var peak in Peaks)
            {
                double weight = peak.Intensity / sumIntensity;
                averagedIntensity += weight * peak.Intensity;
            }
            return averagedIntensity;
        }

        public double AverageMass()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity * p.Intensity);
            double averagedMass = 0;
            foreach (var mass in Peaks)
            {
                double weight = mass.Intensity * mass.Intensity / sumIntensity;
                averagedMass += weight * mass.MonoisotopicMass;
            }
            return averagedMass;
        }

        public void CalculateNL()
        {
            NL = Peaks.Min(p => p.Intensity);
        }


        public void GetNormalizedPeaks(string type = "rt")
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            NormalizedPeaks = new (double, double)[Peaks.Count];
            if (type == "rt")
            {
                for (int i = 0; i < Peaks.Count; i++)
                {
                    NormalizedPeaks[i] = (Peaks[i].RetentionTime, Peaks[i].Intensity / sumIntensity * 100);
                }
            }
            else
            {
                for (int i = 0; i < Peaks.Count; i++)
                {
                    NormalizedPeaks[i] = ((double)Peaks[i].ZeroBasedScanIndex, Peaks[i].Intensity / sumIntensity * 100);
                }
            }
        }

        public void GetNormalizedLinearSplinePeaks()
        {
            var cycleArray = Peaks.Select(p => (double)p.ZeroBasedScanIndex).ToArray();
            var intensityArray = Peaks.Select(p => p.Intensity).ToArray();
            var linearSpline = LinearSpline.InterpolateSorted(cycleArray, intensityArray);

            int maxLength = EndCycle - StartCycle + 1;
            var linearSplinePeaks = new (int, double)[maxLength];
            NormalizedLinearSplinePeaks = new (int, double)[maxLength];
            for (int i = 0; i < maxLength; i++)
            {
                int j = 0;
                if (Peaks[j].ZeroBasedScanIndex == StartCycle + i)
                {
                    linearSplinePeaks[i] = (Peaks[i].ZeroBasedScanIndex, Peaks[i].Intensity);
                    j++;
                }
                else
                {
                    linearSplinePeaks[i] = (StartCycle + i, linearSpline.Interpolate((double) (StartCycle + i)));
                }
            }
            for (int i = 0; i < maxLength; i++)
            {
                NormalizedLinearSplinePeaks[i] = (linearSplinePeaks[i].Item1, linearSplinePeaks[i].Item2 / linearSplinePeaks.Sum(p => p.Item2));
            }
        }

        private void GetLinearSpline()
        {
            var sortedPeaks = Peaks.OrderBy(p => p.RetentionTime).ToList();
            var rtArray = sortedPeaks.Select(p => p.RetentionTime).ToArray();
            var intensityArray = sortedPeaks.Select(p => p.Intensity).ToArray();
            var linearSpline = LinearSpline.InterpolateSorted(rtArray, intensityArray);
            this.LinearSpline = linearSpline;
        }

        private void GetCubicSpline()
        {
            var sortedPeaks = Peaks.OrderBy(p => p.RetentionTime).ToList();
            var rtArray = sortedPeaks.Select(p => p.RetentionTime).ToArray();
            var intensityArray = sortedPeaks.Select(p => p.Intensity).ToArray();
            var cubicSpline = CubicSpline.InterpolateAkima(rtArray, intensityArray);
            this.CubicSpline = cubicSpline;
        }

        private (double, double)[] CalculateSpline(double startRT, double endRT, double splineRtInterval, IInterpolation spline)
        {
            int numPoints = (int)Math.Floor((endRT - startRT) / splineRtInterval) + 1;
            var xyData = new (double, double)[numPoints];
            for (int i = 0; i < numPoints; i++)
            {
                var rt = startRT + i * splineRtInterval;
                var intensity = spline.Interpolate(rt);
                xyData[i] = (rt, intensity);
            }
            return xyData;
        }

        public void GetLinearSplineXYData(double splineRtInterval)
        {
            if (LinearSpline == null)
            {
                GetLinearSpline();
            }
            XYData = CalculateSpline(StartRT, EndRT, splineRtInterval, LinearSpline);
        }

        public void GetCubicSplineXYData(double splineRtInterval, Dictionary<int, double> rtIndexMap = null)
        {
            if (Peaks.Count < 5)
            {
                GetExtendedCubicSplineXYData(splineRtInterval, rtIndexMap);
                return;
            }
            if (CubicSpline == null)
            {
                GetCubicSpline();
            }
            XYData = CalculateSpline(StartRT, EndRT, splineRtInterval, CubicSpline);
        }

        public void GetExtendedCycleCubicSplineXYData( double splineRtInterval, int numberOfPeaksToAdd = 1)
        {
            AddPeaks(numberOfPeaksToAdd, out double[] newRtArray, out double[] newIntensityArray);
            var extendedCubicSpline = CubicSpline.InterpolateAkima(newRtArray.ToArray(), newIntensityArray.ToArray());
            XYData = CalculateSpline(newRtArray[0], newRtArray[newRtArray.Length - 1], splineRtInterval, extendedCubicSpline);
        }

        public void GetExtendedCycleCubicSplineSavgolSmoothedXYData(double splineRtInterval, int windowSize)
        {
            GetExtendedCycleCubicSplineXYData(splineRtInterval);
            var y = CalculateSavgolSmoothedData(XYData.Select(xy => xy.Item2).ToArray(), windowSize);
            for (int i = 0; i < XYData.Length; i++)
            {
                XYData[i] = (XYData[i].Item1, y[i]);
            }
        }

        public void GetExtendedCubicSplineXYData(double splineRtInterval, Dictionary<int, double> rtIndexMap)
        {
            var rtArray = Peaks.Select(p => p.RetentionTime).ToList();
            var intensityArray = Peaks.Select(p => p.Intensity).ToList();
            if (Peaks.First().ZeroBasedScanIndex != 0)
            {
                rtArray.Insert(0, rtIndexMap[Peaks.First().ZeroBasedScanIndex - 1]);
                intensityArray.Insert(0, 0);
            }
            if (Peaks.Last().ZeroBasedScanIndex != rtIndexMap.Count - 1)
            {
                rtArray.Add(rtIndexMap[Peaks.Last().ZeroBasedScanIndex + 1]);
                intensityArray.Add(0);
            }
            if (rtArray.Count < 5)
            {
                return;
            }
            var extendedCubicSpline = CubicSpline.InterpolateAkima(rtArray.ToArray(), intensityArray.ToArray());
            XYData = CalculateSpline(rtArray.First(), rtArray.Last(), splineRtInterval, extendedCubicSpline);
        }

        public void GetBSplineXYData(double splineRtInterval, int smoothDegree)
        {
            int PtNum = (int)Math.Floor((EndRT - StartRT) / splineRtInterval) + 1;
            var rawData = Peaks.Select(p => (p.RetentionTime, p.Intensity)).ToList();
            var smoothedData = new Bspline2().Run(rawData, PtNum, smoothDegree);
            XYData = smoothedData.ToArray();
        }

        public void GetMs1SpaceBSplineXYData(double splineRtInterval, int smoothDegree, Dictionary<double, double> rtMap)
        {
            int PtNum = (int)Math.Floor((EndRT - StartRT) / splineRtInterval) + 1;
            var peaks = Peaks.Select(p => (rtMap[p.RetentionTime], p.Intensity)).ToList();
            var splineData = new Bspline2().Run(peaks, PtNum, smoothDegree);
            XYData = splineData.ToArray();
        }

        public List<(float, float)> GetUmpireBSplineFloatData(int NoPeakPerMin, int smoothDegree)
        {
            var rawData = Peaks.Select(p => ((float)p.RetentionTime, (float)p.Intensity)).ToList();
            var smoothedData = new Bspline().Run(rawData, (int)Math.Max((EndRT - StartRT) * NoPeakPerMin, Peaks.Count), smoothDegree);
            BsplineSmoothedData = smoothedData;
            return smoothedData;
        }

        public void GetUmpireBSplineData(int NoPeakPerMin, int smoothDegree)
        {
            var rawData = Peaks.Select(p => (p.RetentionTime, p.Intensity)).ToList();
            var smoothedData = new Bspline2().Run(rawData, (int)Math.Max((EndRT - StartRT) * NoPeakPerMin, Peaks.Count), smoothDegree);
            UmpireBsplineData = smoothedData;
            XYData = UmpireBsplineData.ToArray();
        }

        public void GetRawXYData()
        {
            XYData = new (double, double)[Peaks.Count];
            for (int i = 0; i < Peaks.Count; i++)
            {
                XYData[i] = (Peaks[i].ZeroBasedScanIndex, Peaks[i].Intensity);
            }
        }

        public void NormalizeXYData()
        {
            if (XYData == null || XYData.Length == 0)
            {
                GetRawXYData();
            }
            double sumIntensity = XYData.Sum(p => p.Item2);
            for (int i = 0; i < XYData.Length; i++)
            {
                XYData[i] = (XYData[i].Item1, XYData[i].Item2 / sumIntensity * 100);
            }
        }

        public void GetSavgolSmoothedXYData(int windowSize)
        {
            if (CubicSpline == null)
            {
                GetCubicSpline();
            }
            double[] imputedY = ImputeMissingValues(CubicSpline, out double[] rtArray);
            double[] y = CalculateSavgolSmoothedData(imputedY, windowSize);
            XYData = new (double, double)[y.Length];
            double cycleInterval = (double)(EndCycle - StartCycle) / (y.Length - 1);
            for (int i = 0; i < y.Length; i++)
            {
                XYData[i] = (rtArray[i], y[i]);
            }
        }

        public void GetSavgolSmoothedCubicSplineXYData (int windowSize, double splineRtInterval)
        {
            GetSavgolSmoothedXYData(windowSize);
            double[] x = XYData.Select(p => p.Item1).ToArray();
            double[] y = XYData.Select(p => p.Item2).ToArray();
            var cubicSpline = CubicSpline.InterpolateAkima(x, y);
            XYData = CalculateSpline(StartRT, EndRT, splineRtInterval, cubicSpline);
        }

        public void GetExtendedSavgolSmoothedCycleCubicSplineXYData(int numberOfPeaksToAdd, int windowSize, double splineRtInterval)
        {
            AddPeaks(numberOfPeaksToAdd, out double[] newRtArray, out double[] newIntensityArray);
            double[] y = CalculateSavgolSmoothedData(newIntensityArray, windowSize);
            var cubicSpline = CubicSpline.InterpolateAkima(newRtArray, y);
            XYData = CalculateSpline(newRtArray[0], newRtArray[newRtArray.Length - 1], splineRtInterval, cubicSpline);
        }

        public void AddPeaks(int numberOfPeaksToAdd, out double[] newRtArray, out double[] newIntensityArray)
        {
            var rtArray = Peaks.Select(p => (double)p.ZeroBasedScanIndex).ToArray();
            var intensityArray = Peaks.Select(p => p.Intensity).ToArray();
            newRtArray = new double[Peaks.Count + numberOfPeaksToAdd * 2];
            newIntensityArray = new double[Peaks.Count + numberOfPeaksToAdd * 2];
            for (int i = 0; i < newRtArray.Length; i++)
            {
                if (i < numberOfPeaksToAdd)
                {
                    newRtArray[i] = rtArray[0] - (numberOfPeaksToAdd - i) * 1;
                    newIntensityArray[i] = 0;
                }
                else if (i >= rtArray.Length + numberOfPeaksToAdd)
                {
                    newRtArray[i] = newRtArray[i - 1] + 1;
                    newIntensityArray[i] = 0;
                }
                else
                {
                    newRtArray[i] = rtArray[i - numberOfPeaksToAdd];
                    newIntensityArray[i] = intensityArray[i - numberOfPeaksToAdd];
                }
            }
        }

        public void GetCubicSplineSavgolSmoothedXYData(int windowSize, double splineRtInterval)
        {
            GetCubicSplineXYData(splineRtInterval);
            double[] y = XYData.Select(p => p.Item2).ToArray();
            var smoothedY = CalculateSavgolSmoothedData(y, windowSize);
            for (int i = 0; i < XYData.Length; i++)
            {
                XYData[i] = (XYData[i].Item1, smoothedY[i]);
            }
        }

        public void GetScanCycleCubicSplineXYData(double splineCycleInterval)
        {
            int numPoints = (int)Math.Floor((EndCycle - StartCycle) / splineCycleInterval) + 1;
            double splineRtInterval = (EndRT - StartRT) / (numPoints - 1);
            GetCubicSplineXYData(splineRtInterval);
            for (int i = 0; i < XYData.Length; i++)
            {
                var cyclePoint = StartCycle + i * splineCycleInterval;
                XYData[i] = (cyclePoint, XYData[i].Item2);
            }
        }

        public void GetScanCycleCubicSplineSavgolSmoothedXYData(int windowSize, double splineRtInterval)
        {
            GetScanCycleCubicSplineXYData(splineRtInterval);
            double[] y = XYData.Select(p => p.Item2).ToArray();
            var smoothedY = CalculateSavgolSmoothedData(y, windowSize);
            for (int i = 0; i < XYData.Length; i++)
            {
                XYData[i] = (XYData[i].Item1, smoothedY[i]);
            }
        }

        public void GetMs1SpaceCubicSplineXYData(Dictionary<double, double> rtMap, double splineRTinterval)
        {
            var sortedPeaks = Peaks.OrderBy(p => p.RetentionTime).ToList();
            var rtArray = sortedPeaks.Select(p => rtMap[p.RetentionTime]).ToArray();
            var intensityArray = sortedPeaks.Select(p => p.Intensity).ToArray();
            var ms1SpaceSpline = CubicSpline.InterpolateAkima(rtArray, intensityArray);

            var startRT = rtArray[0];
            var endRT = rtArray[rtArray.Length - 1];
            XYData = CalculateSpline(startRT, endRT, splineRTinterval, ms1SpaceSpline);
        }

        public void GetMs1SpaceCubicSplineSavgolSmoothedXYData(Dictionary<double, double> rtMap, int windowSize, double splineRTinterval)
        {
            GetMs1SpaceCubicSplineXYData(rtMap, splineRTinterval);
            var y = XYData.Select(p => p.Item2).ToArray();
            var smoothedY = CalculateSavgolSmoothedData(y, windowSize);
            for (int i = 0 ; i < smoothedY.Length; i++)
            {
                XYData[i] = (XYData[i].Item1, smoothedY[i]);
            }
        }

        public void GetMs1SpaceSavgolSmoothedCubicSplineXYData(Dictionary<double, double> rtMap, int windowSize, double splineRTinterval)
        {
            var ms1SpaceRts = Peaks.Select(p => rtMap[p.RetentionTime]).ToArray();
            var ms1SpaceCubicSpline = CubicSpline.InterpolateAkima(ms1SpaceRts, Peaks.Select(p => p.Intensity).ToArray());

            double meanInterval = CalculateAverageInterval();
            var numPoints = EndCycle - StartCycle + 1;
            int[] cycles = Enumerable.Range(StartCycle, numPoints).ToArray();
            var x = new double[numPoints];
            var y = new double[numPoints];
            var indexList = Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            for (int i = 0; i < numPoints; i++)
            {
                var id = Array.BinarySearch(indexList.ToArray(), cycles[i]);
                if (id >= 0)
                {
                    x[i] = ms1SpaceRts[id];
                    y[i] = Peaks[id].Intensity;
                }
                else
                {
                    x[i] = x[i - 1] + meanInterval;
                    y[i] = ms1SpaceCubicSpline.Interpolate(x[i]);
                }
            }
            var smoothedY = CalculateSavgolSmoothedData(y, windowSize);
            var cubicSpline = CubicSpline.InterpolateAkima(x, smoothedY);
            XYData = CalculateSpline(x[0], x[x.Length - 1], splineRTinterval, cubicSpline);
        }

        public void GetSimpleGaussianXYData()
        {
            var x = Peaks.Select(p => p.RetentionTime).ToArray();
            var y = Peaks.Select(p => p.Intensity).ToArray();
            var normalizedY = y.Select(p => p / y.Max() * 100).ToArray();

            double A = y.Max();
            double mu = x[Array.IndexOf(y, A)];
            double sigma = (x.Max() - x.Min()) / 4;

            XYData = new (double, double)[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                var yFit = A * Math.Exp(-(x[i] - mu) * (x[i] - mu) / (2 * sigma * sigma));
                XYData[i] = (x[i], yFit);
            }
        }

        public void GetSimpleGaussianSplineXYData(double splineRTinterval)
        {
            var x = Peaks.Select(p => p.RetentionTime).ToArray();
            var y = Peaks.Select(p => p.Intensity).ToArray();
            var normalizedY = y.Select(p => p / y.Max() * 100).ToArray();

            double A = y.Max();
            double mu = x[Array.IndexOf(y, A)];
            double sigma = (x.Max() - x.Min()) / 4;

            int numPoints = (int)Math.Floor((EndRT - StartRT) / splineRTinterval) + 1;
            XYData = new (double, double)[numPoints];
            for (int i = 0; i < numPoints; i++)
            {
                var rt = StartRT + i * splineRTinterval;
                var intensity = A * Math.Exp(-(rt - mu) * (rt - mu) / (2 * sigma * sigma));
                XYData[i] = (rt, intensity);
            }
        }

        public void GetGaussianFitXYData()
        {
            var x = Peaks.Select(p => p.RetentionTime).ToArray();
            var y = Peaks.Select(p => p.Intensity).ToArray();
            var normalizedY = y.Select(p => p / y.Max() * 100).ToArray();

            var (A, mu, sigma) = FitGaussian(x, normalizedY);
            XYData = new (double, double)[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                var yFit = A * Math.Exp(-(x[i] - mu) * (x[i] - mu) / (2 * sigma * sigma));
                XYData[i] = (x[i], yFit);
            }
        }

        public static (double A, double mu, double sigma) FitGaussian(double[] x, double[] y)
        {
            double A_guess = y.Max();
            double mu_guess = x[Array.IndexOf(y, A_guess)];
            double sigma_guess = (x.Max() - x.Min()) / 4;

            Func<Vector<double>, double, double> model = (parameters, x) =>
            {
                double A = parameters[0];
                double mu = parameters[1];
                double sigma = parameters[2];
                return A * Math.Exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
            };

            var xdata = Vector<double>.Build.DenseOfArray(x);
            var ydata = Vector<double>.Build.DenseOfArray(y);

            // Define the objective function (sum of squared errors)
            var objective = ObjectiveFunction.NonlinearModel(
                model,
                observedX: xdata,
                observedY: ydata,
                weight: null,
                accuracyOrder: 2
            );

            var initialGuess = Vector<double>.Build.DenseOfArray(new double[] { A_guess, mu_guess, sigma_guess });
            var solver = new LevenbergMarquardtMinimizer();
            var result = solver.FindMinimum(objective, initialGuess);

            return (result.MinimizingPoint[0], result.MinimizingPoint[1], result.MinimizingPoint[2]);
        }

        //Test
        public void ScanCycleSpline(double interval)
        {
            var sortedIndex = Peaks.Select(p => (double)p.ZeroBasedScanIndex).OrderBy(t => t).ToArray();
            var sortedIntensity = Peaks.Select(p => p.Intensity).ToArray();
            var spline = CubicSpline.InterpolateAkima(sortedIndex, sortedIntensity);
            XYData = CalculateSpline(StartCycle, EndCycle, interval, spline);
            if (XYData.Select(xy => xy.Item2).Contains(double.NaN))
            {
                int stop = 0;
            }
        }


        private double CalculateAverageInterval()
        {
            var intervals = new List<double>();
            var rts = Peaks.Select(p => p.RetentionTime).ToArray();
            var indices = Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            for (int i = 0; i < rts.Length - 1; i++)
            {
                if (indices[i + 1] == indices[i] + 1)
                {
                    var rtInterval = rts[i + 1] - rts[i];
                    intervals.Add(rtInterval);
                }
            }
            return intervals.Mean();
        }

        private double[] ImputeMissingValues(IInterpolation spline, out double[] x)
        {
            double meanInterval = CalculateAverageInterval();
            var numPoints = EndCycle - StartCycle + 1;
            int[] cycles = Enumerable.Range(StartCycle, numPoints).ToArray();
            var rts = Peaks.Select(p => p.RetentionTime).ToArray();
            x = new double[numPoints];
            var y = new double[numPoints];
            var indexList = Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            for (int i = 0; i < numPoints; i++)
            {
                var id = Array.BinarySearch(indexList.ToArray(), cycles[i]);
                if (id >= 0)
                {
                    x[i] = rts[id];
                    y[i] = Peaks[id].Intensity;
                }
                else
                {
                    x[i] = x[i - 1] + meanInterval;
                    y[i] = spline.Interpolate(x[i]);
                }
            }
            return y;
        }
        private static double[] CalculateSavgolSmoothedData(double[] y, int windowSize)
        {
            float[] newY = y.Select(p => (float)p).ToArray();
            var smoother = new SavitzkyGolayFilter(windowSize);
            var signal = new DiscreteSignal(1, newY);
            var extra = smoother.Size / 2;
            var smoothedY = smoother.ApplyTo(signal).Samples.Skip(extra).Take(y.Length).ToArray();
            var yData = smoothedY.Select(p => (double)p).ToArray();

            return yData;
        }

        public void GetPrecursorRanks()
        {
            if (PFpairs == null || PFpairs.Count == 0)
            {
                return;
            }
            // Sort PFpairs by correlation in descending order and assign ranks
            var rankedPairs = PFpairs.OrderByDescending(p => p.Correlation)
                .Select((p, index) => new { PFpair = p, Rank = index + 1 }).ToList();

            // Update the PFpairs with their ranks
            foreach (var rankedPair in rankedPairs)
            {
                rankedPair.PFpair.PrecursorRank = rankedPair.Rank;
            }
        }

        public static Peak GetPeakFromScan(double targetMz, List<Peak>[] peakTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            Peak bestPeak = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(targetMz) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(targetMz) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < peakTable.Length && peakTable[j] != null)
                {
                    List<Peak> bin = peakTable[j];
                    int index = BinarySearchForIndexedPeak(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        Peak peak = bin[i];

                        if (peak.ZeroBasedScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        if (tolerance.Within(peak.Mz, targetMz) && peak.ZeroBasedScanIndex == zeroBasedScanIndex
                            && (bestPeak == null || Math.Abs(peak.Mz - targetMz) < Math.Abs(bestPeak.Mz - targetMz)))
                        {
                            bestPeak = peak;
                        }
                    }
                }
            }

            return bestPeak;
        }

        private static int BinarySearchForIndexedPeak(List<Peak> peakList, int zeroBasedScanIndex)
        {
            int m = 0;
            int l = 0;
            int r = peakList.Count - 1;

            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                {
                    break;
                }
                if (peakList[m].ZeroBasedScanIndex < zeroBasedScanIndex)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }

            for (int i = m; i >= 0; i--)
            {
                if (peakList[i].ZeroBasedScanIndex < zeroBasedScanIndex)
                {
                    break;
                }
                m--;
            }

            if (m < 0)
            {
                m = 0;
            }

            return m;
        }

        public static PeakCurve FindPeakCurve(Peak targetPeak, List<Peak>[] massTable, MsDataScan[] scans, MzRange isolationWindow, int maxMissedScans
            , Tolerance tolerance, int binSize, double maxRTrange)
        {
            var xic = new List<Peak>();
            xic.Add(targetPeak);
            PeakCurve newPeakCurve = new PeakCurve(xic, targetPeak.MsLevel, isolationWindow);
            targetPeak.PeakCurve = newPeakCurve;

            // go right
            int missedScans = 0;
            for (int t = targetPeak.ZeroBasedScanIndex + 1; t < massTable.Length; t++)
            {
                //Changed for test!! Remember to change back!!!
                var peak = GetPeakFromScan(targetPeak.Mz, massTable, t, tolerance, binSize);

                if (peak == null)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    //if(peak.PeakCurve == null || Math.Abs(peak.Mz - newPeakCurve.AveragedMz) < Math.Abs(peak.Mz - peak.PeakCurve.AveragedMz))
                    if(peak.PeakCurve == null)
                    {
                        missedScans = 0;
                        xic.Add(peak);
                        peak.PeakCurve = newPeakCurve;
                    }
                    else
                    {
                        missedScans++;
                    }
                }

                if (missedScans > maxMissedScans)
                {
                    break;
                }
                if (newPeakCurve.EndRT - newPeakCurve.ApexRT > maxRTrange)
                {
                    break;
                }
            }

            // go left
            missedScans = 0;
            for (int t = targetPeak.ZeroBasedScanIndex - 1; t >= 0; t--)
            {
                //Changed for test!! Remember to change back!!!
                var peak = GetPeakFromScan(targetPeak.Mz, massTable, t, tolerance, binSize);

                if (peak == null)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    //if (peak.PeakCurve == null || Math.Abs(peak.Mz - newPeakCurve.AveragedMz) < Math.Abs(peak.Mz - peak.PeakCurve.AveragedMz))
                    if (peak.PeakCurve == null)
                    {
                        missedScans = 0;
                        xic.Add(peak);
                        peak.PeakCurve = newPeakCurve;
                    }
                    else
                    {
                        missedScans++;
                    }
                }

                if (missedScans > maxMissedScans)
                {
                    break;
                }
                if (newPeakCurve.ApexRT - newPeakCurve.StartRT > maxRTrange)
                {
                    break;
                }
            }
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return newPeakCurve;
        }

        //public static PeakCurve FindPeakCurve_cutPeak(Peak targetPeak, List<Peak>[] peakTable, MsDataScan[] scans, MzRange isolationWindow, int maxMissedScans
        //    , Tolerance mzTolerance, int binSize, double maxRTrange = 2)
        //{
        //    var xic = new List<Peak>();
        //    var peakList = new List<Peak>();
        //    peakList.Add(targetPeak);
        //    PeakCurve newPeakCurve = new PeakCurve(xic, targetPeak.MsLevel, isolationWindow);
        //    targetPeak.PeakCurve = newPeakCurve;

        //    // go right
        //    int missedScans = 0;
        //    for (int t = targetPeak.ZeroBasedScanIndex + 1; t < scans.Length; t++)
        //    {
        //        var peak = GetPeakFromScan(newPeakCurve.AveragedMz, peakTable, t, mzTolerance, binSize);

        //        if (peak == null)
        //        {
        //            missedScans++;
        //        }
        //        else if (peak != null)
        //        {
        //            if (peak.RetentionTime - targetPeak.RetentionTime > maxRTrange)
        //            {
        //                break;
        //            }
        //            //if(peak.PeakCurve == null || Math.Abs(peak.Mz - newPeakCurve.AveragedMz) < Math.Abs(peak.Mz - peak.PeakCurve.AveragedMz))
        //            if (peak.PeakCurve == null)
        //            {
        //                missedScans = 0;
        //                peakList.Add(peak);
        //            }
        //            else
        //            {
        //                missedScans++;
        //            }
        //        }

        //        if (missedScans > maxMissedScans)
        //        {
        //            break;
        //        }
        //    }

        //    // go left
        //    missedScans = 0;
        //    for (int t = targetPeak.ZeroBasedScanIndex - 1; t >= 0; t--)
        //    {
        //        var peak = GetPeakFromScan(newPeakCurve.AveragedMz, peakTable, t, mzTolerance, binSize);

        //        if (peak == null)
        //        {
        //            missedScans++;
        //        }
        //        else if (peak != null)
        //        {
        //            if (peak.RetentionTime - targetPeak.RetentionTime > maxRTrange)
        //            {
        //                break;
        //            }
        //            //if (peak.PeakCurve == null || Math.Abs(peak.Mz - newPeakCurve.AveragedMz) < Math.Abs(peak.Mz - peak.PeakCurve.AveragedMz))
        //            if (peak.PeakCurve == null)
        //            {
        //                missedScans = 0;
        //                peakList.Add(peak);
        //            }
        //            else
        //            {
        //                missedScans++;
        //            }
        //        }

        //        if (missedScans > maxMissedScans)
        //        {
        //            break;
        //        }
        //    }

        //    if (peakList.Count > 5)
        //    {
        //        //DiscriminationFactorToCutPeak is a parameter, default 0.6 in FlashLFQ
        //        double DiscriminationFactorToCutPeak = 0.6;

        //        peakList = peakList.OrderBy(p => p.RetentionTime).ToList();
        //        var apexPeak = peakList.OrderByDescending(p => p.Intensity).First();
        //        int apexIndex = peakList.IndexOf(apexPeak);
        //        Peak valleyPeak = null;
        //        int indexOfValley = 0;

        //        //go left
        //        for (int i = apexIndex; i >= 0; i--)
        //        {
        //            Peak timepoint = peakList[i];

        //            if (valleyPeak == null || timepoint.Intensity < valleyPeak.Intensity)
        //            {
        //                valleyPeak = timepoint;
        //                indexOfValley = peakList.IndexOf(valleyPeak);
        //            }

        //            double discriminationFactor =
        //                (timepoint.Intensity - valleyPeak.Intensity) / timepoint.Intensity;

        //            if (discriminationFactor > DiscriminationFactorToCutPeak)
        //            {
        //                peakList.RemoveAll(p => p.RetentionTime < valleyPeak.RetentionTime);
        //                break;
        //            }
        //        }

        //        //go right
        //        valleyPeak = null;
        //        for (int i = apexIndex; i < peakList.Count; i++)
        //        {
        //            Peak timepoint = peakList[i];

        //            if (valleyPeak == null || timepoint.Intensity < valleyPeak.Intensity)
        //            {
        //                valleyPeak = timepoint;
        //                indexOfValley = peakList.IndexOf(valleyPeak);
        //            }

        //            double discriminationFactor =
        //                (timepoint.Intensity - valleyPeak.Intensity) / timepoint.Intensity;

        //            if (discriminationFactor > DiscriminationFactorToCutPeak)
        //            {
        //                peakList.RemoveAll(p => p.RetentionTime > valleyPeak.RetentionTime);
        //                break;
        //            }
        //        }
        //    }

        //    foreach(var peak in peakList)
        //    {
        //        xic.Add(peak);
        //        peak.PeakCurve = newPeakCurve;
        //    }

        //    xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

        //    return newPeakCurve;
        //}

        public GenericChart VisualizeBspline(out List<float> rtSeq)
        {
            var rawData = Peaks.Select(p => ((float)p.RetentionTime, (float)p.Intensity)).ToList();
            rtSeq = new List<float>();
            for (float i = (float)StartRT; i < (float)EndRT; i += 0.005f)
            {
                rtSeq.Add(i);
            }
            int numPt = rtSeq.Count;
            var smoothedData = new Bspline().Run(rawData, numPt, 2);
            var plot = Chart2D.Chart.Line<float, float, string>(
                x: smoothedData.Select(p => p.Item1),
                y: smoothedData.Select(p => p.Item2)).WithTraceInfo("spline").WithMarkerStyle(Color: Color.fromString("blue"));
            return plot;
        }

        public GenericChart VisualizeGeneral(string type, bool XYDataOnly = false)
        {
            GenericChart plot = null;
            if (XYDataOnly)
            {
                if (type == "rt")
                {
                    plot = Chart2D.Chart.Line<double, double, string>(
                            x: XYData.Select(p => p.Item1),
                            y: XYData.Select(p => p.Item2)).WithTraceInfo($"{Math.Round(AveragedMz, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
                    return plot;
                }
                if (type == "cycle")
                {
                    plot = Chart2D.Chart.Line<double, double, string>(
                            x: XYData.Select(p => p.Item1),
                            y: XYData.Select(p => p.Item2)).WithTraceInfo($"{Math.Round(AveragedMz, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
                    return plot;
                }
            }
            if (type == "rt")
            {
                plot = Chart2D.Chart.Line<double, double, string>(
                        x: Peaks.Select(p => p.RetentionTime),
                        y: Peaks.Select(p => p.Intensity)).WithTraceInfo($"{Math.Round(AveragedMz, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
            }
            if (type == "cycle")
            {
                plot = Chart2D.Chart.Line<int, double, string>(
                        x: Peaks.Select(p => p.ZeroBasedScanIndex),
                        y: Peaks.Select(p => p.Intensity)).WithTraceInfo($"{Math.Round(AveragedMz, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
            }
            if (XYData == null)
            {
                return plot;
            }

            var plot2 = Chart2D.Chart.Line<double, double, string>(
                        x: XYData.Select(xy => xy.Item1),
                        y: XYData.Select(xy => xy.Item2)).WithTraceInfo("XYData").WithMarkerStyle(Color: Color.fromString("blue"));
            var combinedPlot = Chart.Combine(new[] { plot, plot2 });
            return combinedPlot;
        }

        public static GenericChart VisualizeCombined(List<PeakCurve> pcs, string type, bool xyDataOnly = false)
        {
            var allPlots = new List<GenericChart>();
            foreach (var pc in pcs)
            {
                var plot = pc.VisualizeGeneral(type, xyDataOnly);
                if (plot != null)
                {
                    allPlots.Add(plot);
                }
            }
            var combinedPlot = Chart.Combine(allPlots);
            return combinedPlot;
        }

        public static GenericChart VisualizeCombinedNormalizedRaw(List<PeakCurve> pcs, string type = "rt")
        {
            var allPlots = new List<GenericChart>();
            foreach (var pc in pcs)
            {
                var plot = pc.VisualizeNormalizedRaw(type);
                if (plot != null)
                {
                    allPlots.Add(plot);
                }
            }
            var combinedPlot = Chart.Combine(allPlots);
            return combinedPlot;
        }

        public GenericChart VisualizeNormalizedRaw(string type = "rt")
        {
            if (NormalizedPeaks == null)
            {
                GetNormalizedPeaks(type);
            }
            var plot = Chart2D.Chart.Line<double, double, string>(
                        x: NormalizedPeaks.Select(p => p.Item1),
                        y: NormalizedPeaks.Select(p => p.Item2)).WithTraceInfo($"{Math.Round(AveragedMz, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("blue"));
            return plot;
        }

        public GenericChart VisualizeRaw()
        {
            var plot = Chart2D.Chart.Line<double, double, string>(
                        x: Peaks.Select(p => p.RetentionTime),
                        y: Peaks.Select(p => p.Intensity)).WithTraceInfo($"{Math.Round(AveragedMz, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("blue"));
            return plot;
        }
        public static GenericChart VisualizeCombinedRaw(List<PeakCurve> pcs)
        {
            var allPlots = new List<GenericChart>();
            foreach (var pc in pcs)
            {
                var plot = pc.VisualizeRaw();
                allPlots.Add(plot);
            }
            var combinedPlot = Chart.Combine(allPlots);
            return combinedPlot;
        }

        public GenericChart VisualizeUmpireBSplineData()
        {
            var plot1 = Chart2D.Chart.Line<double, double, string>(
                        x: Peaks.Select(p => p.RetentionTime),
                        y: Peaks.Select(p => p.Intensity)).WithTraceInfo($"{Math.Round(AveragedMz, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
            var plot2 = Chart2D.Chart.Line<double, double, string>(
                       x: UmpireBsplineData.Select(p => p.Item1),
                       y: UmpireBsplineData.Select(p => p.Item2)).WithTraceInfo("Umpire BSpline").WithMarkerStyle(Color: Color.fromString("blue"));
            var combinedPlot = Chart.Combine(new[] { plot1, plot2 });
            return combinedPlot;
        }

        public double IntegrateAreaUnderCurve(string type = "TotalIntensity")
        {
            var data = new List<(double, double)>();
            if (type == "TotalIntensity")
            {
                foreach (var peak in Peaks)
                {
                    data.Add((peak.RetentionTime, peak.TotalIntensity));
                }
            }
            if (type == "HighestPeakIntensity")
            {
                foreach (var peak in Peaks)
                {
                    data.Add((peak.RetentionTime, peak.HighestPeakIntensity));
                }
            }
            var area = PrecursorFragmentPair.CalculateArea(data);
            return area;
        }

        //public GenericChart VisualizeCubicSpline(float timeInterval)
        //{
        //    var smoothedData = Interpolte_cubic(timeInterval);
        //    var raw = Chart2D.Chart.Point<float, float, string>(
        //        x: Peaks.Select(p => (float)p.RetentionTime),
        //        y: Peaks.Select(p => (float)p.Intensity)).WithTraceInfo("raw").WithMarkerStyle(Color: Color.fromString("red"));
        //    var plot = Chart2D.Chart.Point<float, float, string>(
        //        x: smoothedData.Select(p => p.Item1),
        //        y: smoothedData.Select(p => p.Item2)).WithTraceInfo("spline").WithMarkerStyle(Color: Color.fromString("blue"));
        //    var combined = Chart.Combine(new[] { plot, raw });
        //    return combined;
        //}
        //public List<(float, float)> Interpolte_cubic(float timeInterval)
        //{
        //    if (CubicSpline == null)
        //    {
        //        GetCubicSpline();
        //    }
        //    var rtSeq = new List<float>();
        //    for (float i = (float)StartRT; i < (float)EndRT; i += timeInterval)
        //    {
        //        rtSeq.Add(i);
        //    }
        //    var smoothedData = new List<(float, float)>();
        //    for (int i = 0; i < rtSeq.Count; i++)
        //    {
        //        smoothedData.Add((rtSeq[i], (float)CubicSpline.Interpolate(rtSeq[i])));
        //    }
        //    return smoothedData;
        //}

        public void DetectPeakRegions()
        {
            PeakRidgeList = new List<PeakRidge>();
            PeakRegionList = new List<(float rt1, float rt2, float rt3)>();
            NoRidgeRegion = new List<List<float>>();

            //need a parameter
            if (EndRT - StartRT < CwtParameters.MinPeakWidth)
            {
                PeakRegionList.Add(((float)StartRT, (float)ApexRT, (float)EndRT));
                return;
            }

            //need to consider other spline types and time interval settings
            //if (CubicSpline == null)
            //{
            //    GetCubicSpline();
            //}

            var rtSeq = new List<float>();
            var RTwindow = EndRT - StartRT;
            for (float i = (float)StartRT; i < (float)EndRT; i += CwtParameters.SplineTimeInterval)
            {
                rtSeq.Add(i);
            }
            //for (int i = 0; i < rtSeq.Count; i++)
            //{
            //    peakArrayList[2 * i] = rtSeq[i];
            //    peakArrayList[2 * i + 1] = (float)CubicSpline.Interpolate(rtSeq[i]);
            //}
            
            var smoothedData = GetUmpireBSplineFloatData(rtSeq.Count, 2);
            rtSeq = smoothedData.Select(p => p.Item1).ToList();
            var peakArrayList = new float[rtSeq.Count * 2];
            for (int i = 0; i < rtSeq.Count; i++)
            {
                peakArrayList[2 * i] = smoothedData[i].Item1;
                peakArrayList[2 * i + 1] = smoothedData[i].Item2;
            }
            WaveletMassDetector = new WaveletMassDetector(peakArrayList, rtSeq.Count);
            WaveletMassDetector.Run();

            int maxScale = WaveletMassDetector.PeakRidge.Length - 1;

            float[] DisMatrixF = null;
                                      
            for (int i = maxScale; i >= 0; i--)//trace peak ridge from maximum wavelet scale to minimum scale
            {
                //Get peak ridge list (maximum RT points given a CWT scale
                var PeakRidgeArray = WaveletMassDetector.PeakRidge[i];

                if (PeakRidgeArray == null)
                {
                    maxScale = i;
                    continue;
                }
                if (PeakRidgeArray.Count == 0)
                {
                    continue;
                }

                //RT distance matrix between the existing peak riges and peak ridges extracted from current CWT scale
                int r = PeakRidgeList.Count(), c = PeakRidgeArray.Count();
                if (DisMatrixF == null || r * c > DisMatrixF.Length)
                {
                    DisMatrixF = new float[r * c * 2];
                }

                for (int k = 0; k < PeakRidgeList.Count(); k++)
                {   
                    ///For each existing peak ridge line
                    for (int l = 0; l < PeakRidgeArray.Count(); l++)
                    {
                        DisMatrixF[k * c + l] = (float)Math.Abs(PeakRidgeList[k].RT - PeakRidgeArray[l].rt);
                    }
                }

                var conti = true;
                var removedRidgeList = new List<(float rt, float intensity, int index)>();
                while (conti)
                {
                    //find the smallest value from the matrix first, if find the ridge it belongs to, add it then move to the second smallest; if not, stop looking
                    float closest = float.MaxValue;
                    int ExistingRideIdx = -1;
                    int PeakRidgeInx = -1;
                    for (int k = 0; k < PeakRidgeList.Count(); k++)
                    {
                        for (int l = 0; l < PeakRidgeArray.Count(); l++)
                        {
                            {
                                if (DisMatrixF[k * c + l] < closest)
                                {
                                    closest = DisMatrixF[k * c + l];
                                    ExistingRideIdx = k;
                                    PeakRidgeInx = l;
                                }
                            }
                        }
                    }
                    //if even closet is larger than MinRTRange (too far from the PeakRidge), stop looking => setting conti to false
                    if (closest < float.MaxValue && closest <= CwtParameters.MaxRTDiff)
                    {
                        PeakRidge ridge = PeakRidgeList[ExistingRideIdx]; //update the matched PeakRidge line
                        //PeakRidgeList.Remove(ridge); //remove the existing PeakRidge line from the list
                        ridge.lowScale = i; //update the lowest scale of the PeakRidge line
                        ridge.ContinuousLevel++; //update the continous level of the PeakRidge line
                        var nearestRidge = PeakRidgeArray[PeakRidgeInx]; //why updating the RT?
                        ridge.RT = nearestRidge.rt;
                        ridge.Intensity = peakArrayList[nearestRidge.index * 2 + 1];
                        //PeakRidgeList.Add(ridge); //re-add the updated PeakRidge to the list
                        removedRidgeList.Add(nearestRidge); //remove the potential PeakRidge from the array so we don't start a new PeakRidge line on it
                        for (int k = 0; k < PeakRidgeList.Count(); k++)
                        {
                            DisMatrixF[k * c + PeakRidgeInx] = float.MaxValue; //update the distance matrix so the potential PeakRidge point won't be matched again
                        }
                        for (int l = 0; l < PeakRidgeArray.Count(); l++)
                        {
                            DisMatrixF[ExistingRideIdx * c + l] = float.MaxValue; //update the distance matrix so the existing PeakRidge line won't be matched again
                        }
                    }
                    else
                    {
                        conti = false;
                    }
                }

                PeakRidgeArray.RemoveAll(x => removedRidgeList.Contains(x));
                removedRidgeList.Clear();
                removedRidgeList = null;

                //remove the existing PeakRidge line if it is too far from the current CWT scale and do not have enough continuous levels
                var removelist = new List<PeakRidge>();
                for (int k = 0; k < PeakRidgeList.Count(); k++)
                {
                    PeakRidge existridge = PeakRidgeList[k];
                    if (existridge.lowScale - i > 2 && existridge.ContinuousLevel < maxScale / 2)
                    {
                        removelist.Add(existridge);
                    }
                }
                PeakRidgeList.RemoveAll(x => removelist.Contains(x));
                removelist.Clear();
                removelist = null;

                //for those potential PeakRidge points that are not matched to any existing PeakRidge line, start a new PeakRidge line
                if (i > maxScale / 2)
                {
                    foreach (var ridge in PeakRidgeArray)
                    {
                        PeakRidge newRidge = new PeakRidge(ridge.rt, peakArrayList[ridge.index * 2 + 1], i);
                        newRidge.ContinuousLevel++;
                        //newRidge.intensity = SmoothData.GetPoinByXCloset(newRidge.RT).getY();
                        //don't understand why we need to find the point in the spline again; the intensity is already in the ridge
                        PeakRidgeList.Add(newRidge);
                    }
                }
                PeakRidgeArray.Clear();
                PeakRidgeArray = null;
            }

            //Added code
            //Remove the PeakRidge lines with intensity lower than NL * SNRThreshold
            CalculateNL();
            PeakRidgeList = PeakRidgeList.Where(r => r.Intensity > NL * CwtParameters.SNRThreshold).OrderBy(r => r.RT).ToList();

            //if PeakRidgeList is empty or only have one PeakRidge line, add the whole region as a peak region
            if (PeakRidgeList.Count() <= 1)
            {
                PeakRegionList.Add((rtSeq.First(), (float)ApexRT, rtSeq.Last()));
                var RidgeRTs = new List<float>();
                RidgeRTs.Add((float)ApexRT);
                NoRidgeRegion.Add(RidgeRTs);
            }

            //if we have more than one PeakRidge line, we can have valley points now (local minimum)
            if (PeakRidgeList.Count() > 1)
            {
                var ValleyPoints = new (float rt, float intensity)[PeakRidgeList.Count() + 1];
                ValleyPoints[0] = (peakArrayList[0], peakArrayList[1]);
                PeakRidge currentridge = PeakRidgeList.First();
                var localmin = (-1f, float.MaxValue);
                int startidx = rtSeq.IndexOf(currentridge.RT);//BUG!!

                //loop over all PeakRidge lines, find the local minimum between two PeakRidge lines next to each other
                for (int j = 1; j < PeakRidgeList.Count(); j++)
                {
                    PeakRidge nextridge = PeakRidgeList[j];
                    for (int i = startidx; i < peakArrayList.Count()/2; i++) //loop over all points in the spline
                    {
                        var point = (peakArrayList[2 * i], peakArrayList[2 * i + 1]);
                        if (point.Item1 > currentridge.RT && point.Item1 < nextridge.RT) //check if the point is between the two PeakRidge lines
                        {
                            if (localmin.Item2 > point.Item2) //if the point is lower than the current local minimum, update the local minimum
                            {
                                localmin = (point.Item1, point.Item2);
                            }
                        }
                        if (point.Item1 >= nextridge.RT) //when looping to the next PeakRidge line, stop the loop
                        {
                            startidx = i;
                            break;
                        }
                    }
                    ValleyPoints[j] = localmin; 
                    localmin = (-1f, float.MaxValue);
                    currentridge = nextridge;
                }
                ValleyPoints[PeakRidgeList.Count()] = (rtSeq.Last(), peakArrayList.Last());

                //Correct ridge rt and intensity
                startidx = 0;
                for (int i = 0; i < PeakRidgeList.Count(); i++)
                {
                    PeakRidge ridge = PeakRidgeList[i];
                    for (int j = startidx; j < rtSeq.Count(); j++)
                    {
                        var point = (rtSeq[j], peakArrayList[2*j +1]);
                        if (point.Item1 < ValleyPoints[i + 1].rt)
                        {
                            if (ridge.Intensity < point.Item2)
                            {
                                ridge.Intensity = point.Item2;
                                ridge.RT = point.Item1;
                            }
                        }
                        else
                        {
                            startidx = j;
                            break;
                        }
                    }
                }//go through each peak (between two valleypoints), check if local maximum is another point

                //Find split points to generate peak regions
                bool[] Splitpoints = new bool[PeakRidgeList.Count() - 1];
                int left = 0;
                int right = PeakRidgeList.Count() - 1;
                FindSplitPoint(left, right, ValleyPoints, Splitpoints);

                //check the intensity and width of the split region, if it is too small, merge it with the previous or next region
                //for (int i = 1; i < PeakRidgeList.Count() - 1; i++)
                //{
                //    if (PeakRidgeList[i].Intensity < NL * CwtParameters.SNRThreshold || ValleyPoints[i + 1].rt - ValleyPoints[i].rt < CwtParameters.MinPeakWidth)
                //    {
                //        if (i == 0)
                //        {
                //            Splitpoints[1] = false;
                //            continue;
                //        }
                //        if (i == PeakRidgeList.Count - 2)
                //        {
                //            Splitpoints[i] = false;
                //            continue;
                //        }
                //        if (PeakRidgeList[i + 1].Intensity < PeakRidgeList[i - 1].Intensity)
                //        {
                //            Splitpoints[i] = false;
                //        }
                //        else
                //        {
                //            Splitpoints[i + 1] = false;
                //        }
                //    }
                //}
                bool split = false;

                //VisualizeCubicSpline().Show();

                var RidgeRTs = new List<float>();
                startidx = 0;
                PeakRidge maxridge = PeakRidgeList.First();

                for (int i = 0; i < PeakRidgeList.Count() - 1; i++)
                {
                    RidgeRTs.Add(PeakRidgeList[i].RT);
                    if (PeakRidgeList[i].Intensity > maxridge.Intensity)
                    {
                        maxridge = PeakRidgeList[i]; //find the apex of this split region
                    }
                    if (Splitpoints[i])
                    {
                        //if (maxridge.Intensity < NL * CwtParameters.SNRThreshold)
                        //{
                        //    continue;
                        //}
                        //if (ValleyPoints[i + 1].Item1 - ValleyPoints[startidx].Item1 < CwtParameters.MinPeakWidth)
                        //{
                        //    continue;
                        //}
                        PeakRegionList.Add((ValleyPoints[startidx].Item1, maxridge.RT, ValleyPoints[i + 1].Item1));
                        NoRidgeRegion.Add(RidgeRTs);

                        maxridge = PeakRidgeList[i + 1];
                        RidgeRTs = new List<float>(); //ridgeRTs updates everytime there is a split, it contains all ridge RTs for a split region
                        startidx = i + 1;
                    }
                }
                RidgeRTs.Add(PeakRidgeList.Last().RT);
                if (PeakRidgeList.Last().Intensity > maxridge.Intensity)
                {
                    maxridge = PeakRidgeList[PeakRidgeList.Count() - 1];
                }
                PeakRegionList.Add((ValleyPoints[startidx].Item1, maxridge.RT, ValleyPoints[PeakRidgeList.Count()].Item1));
                //RidgeRTs.trimToSize();
                NoRidgeRegion.Add(RidgeRTs);
            }
            WaveletMassDetector = null;
            PeakRidgeList.Clear();
            PeakRidgeList = null;
        }

        private void FindSplitPoint(int left, int right, (float, float)[] ValleyPoints, bool[] splitpoints)
        {
            for (int i = left; i < right; i++)
            {
                if (ValidSplitPoint(left, right, i, ValleyPoints))
                {
                    splitpoints[i] = true;
                    FindSplitPoint(left, i, ValleyPoints, splitpoints);
                    FindSplitPoint(i + 1, right, ValleyPoints, splitpoints);
                    break;
                }
            }
        }

        private bool ValidSplitPoint(int left, int right, int cut, (float, float)[] ValleyPoints)
        {

            PeakRidge leftridge = PeakRidgeList[left];
            PeakRidge rightridge = PeakRidgeList[cut + 1];

            for (int i = left; i <= cut; i++)
            {
                if (PeakRidgeList[i].Intensity > leftridge.Intensity)
                {
                    leftridge = PeakRidgeList[i];
                }
            }
            for (int i = cut + 1; i <= right; i++)
            {
                if (PeakRidgeList[i].Intensity > rightridge.Intensity)
                {
                    rightridge = PeakRidgeList[i];
                }
            }
            return (Math.Abs(ValleyPoints[left].Item2 - ValleyPoints[cut + 1].Item2) / leftridge.Intensity < CwtParameters.SymThreshold
                && Math.Abs(ValleyPoints[cut + 1].Item2 - ValleyPoints[right + 1].Item2) / rightridge.Intensity < CwtParameters.SymThreshold);
        }

        public PeakCurve[] SeparatePeakByRegion()
        {
            var newPeakCurves = new PeakCurve[PeakRegionList.Count()];
            int i = 0;
            foreach (var region in PeakRegionList)
            {
                var newPeakCurve = new PeakCurve();
                newPeakCurve.Peaks = new List<Peak>();
                newPeakCurve.MsLevel = MsLevel;
                newPeakCurve.BsplineSmoothedData = new List<(float, float)>();
                newPeakCurves[i] = newPeakCurve;
                i++;
            }

            //check if any region is too wide
            //foreach (var region in PeakRegionList)
            //{
            //    var newPeakCurve = new PeakCurve();

            //    if (region.rt3 - region.rt1 > CwtParameters.MaxCurveRTRange)
            //    {
            //        var RTs = SmoothedData.Select(p => p.Item1).ToArray();
            //        int leftidx = BinarySearchLower(RTs, region.rt1);
            //        int rightidx = BinarySearchHigher(RTs, region.rt1);
            //        var left = SmoothedData[leftidx];
            //        var right = SmoothedData[rightidx];
            //        while ((right.Item1 - left.Item1) > CwtParameters.MaxCurveRTRange)
            //        {
            //            if (right.Item1 - region.rt2 <= CwtParameters.MaxCurveRTRange / 4f)
            //            {
            //                leftidx++;
            //            }
            //            else if (region.rt2 - left.Item1 <= CwtParameters.MaxCurveRTRange / 4f)
            //            {
            //                rightidx--;
            //            }
            //            else if (left.Item2 < right.Item2)
            //            {
            //                leftidx++;
            //            }
            //            else
            //            {
            //                rightidx--;
            //            }
            //            left = SmoothedData[leftidx];
            //            right = SmoothedData[rightidx];
            //        }
            //        var newRegion = (left.Item1, region.rt2, right.Item1);
            //        PeakRegionList.Add(newRegion);
            //        PeakRegionList.Remove(region);
            //    }
            //}

            // Add corresponding raw peaks
            foreach(var peak in Peaks)
            {
                for (int j = 0; j < PeakRegionList.Count; j++)
                {
                    if (peak.RetentionTime >= PeakRegionList[j].rt1 && peak.RetentionTime <= PeakRegionList[j].rt3)
                    {
                        newPeakCurves[j].Peaks.Add(peak);
                        break;
                    }
                }
            }

            //Add corresponding smoothed peaks
            //foreach(var point in SmoothedData)
            //{
            //    for (int j = 0; j < PeakRegionList.Count; j++)
            //    {
            //        if (point.Item1 >= PeakRegionList[j].rt1 && point.Item1 <= PeakRegionList[j].rt3)
            //        {
            //            newPeakCurves[j].SmoothedData.Add(point);
            //            break;
            //        }
            //    }
            //}

            return newPeakCurves;
        }

        public static int BinarySearchLower(float[] array, float value)
        {
            if (array.IsNotNullOrEmpty())
            {
                return 0;
            }
            int lower = 0;
            int upper = array.Length - 1;

            if (value - array[upper] >= 0)
            {
                return upper;
            }
            if (value - array[lower] <= 0)
            {
                return 0;
            }

            while (lower <= upper)
            {
                int middle = (lower + upper) / 2;
                float comparisonResult = value - array[middle];
                if (comparisonResult == 0)
                {
                    while (middle - 1 >= 0 && array[middle - 1] == value)
                    {
                        middle--;
                    }
                    return middle;
                }
                else if (comparisonResult < 0)
                {
                    upper = middle - 1;
                }
                else
                {
                    lower = middle + 1;
                }
            }
            if (upper < 0)
            {
                return 0;
            }
            while (upper > 0 && array[upper] >= value)
            {
                upper--;
            }
            return upper;
        }

        public static int BinarySearchHigher(float[] array, float value)
        {
            if (array.IsNotNullOrEmpty())
            {
                return 0;
            }
            int lower = 0;
            int upper = array.Length - 1;

            if (value - array[upper] >= 0)
            {
                return upper;
            }
            if (value - array[lower] <= 0)
            {
                return 0;
            }

            while (lower <= upper)
            {
                int middle = (lower + upper) / 2;
                float comparisonResult = value - array[middle];
                if (comparisonResult == 0)
                {
                    while (middle - 1 >= 0 && array[middle - 1] == value)
                    {
                        middle--;
                    }
                    return middle;
                }
                else if (comparisonResult < 0)
                {
                    upper = middle - 1;
                }
                else
                {
                    lower = middle + 1;
                }
            }
            if (lower > array.Length - 1)
            {
                return array.Length - 1;
            }
            while (lower < array.Length - 1 && array[lower] <= value)
            {
                lower++;
            }
            return lower;
        }

        public void VisualizePeakRegions()
        {
            if (PeakRegionList == null)
            {
                DetectPeakRegions();
            }
            var plot_spline = VisualizeBspline(out List<float> rtSeq);
            var markedRTs = new List<float> { PeakRegionList[0].rt1 };
            var smoothedData = GetUmpireBSplineFloatData(rtSeq.Count, 2);
            float yMin = smoothedData.Min(p => p.Item2);
            float yMax = smoothedData.Max(p => p.Item2);
            foreach (var region in PeakRegionList)
            {
                markedRTs.Add(region.rt3);
            }
            var verticalLines = new List<GenericChart>();
            foreach (var x in markedRTs)
            {
                var line = Chart2D.Chart.Line<double, double, string>(
                    new List<double> { x, x },
                    new List<double> { yMin, yMax } 
                ).WithLineStyle(Width: 2, Color: Color.fromString("red"));
                verticalLines.Add(line);
            }
            var combinedPlot = Chart.Combine(new[] { plot_spline }.Concat(verticalLines).ToArray());
            combinedPlot.Show();
        }

        public void CutPeak()
        {
            // find out if we need to split this peak by using the discrimination factor
            // this method assumes that the isotope envelopes in a chromatographic peak are already sorted by MS1 scan number
            bool cutThisPeak = false;

            if (Peaks.Count < 5)
            {
                return;
            }
            var peaksToRemove = new List<Peak>();

            //DiscriminationFactorToCutPeak is a parameter, default 0.6 in FlashLFQ
            double DiscriminationFactorToCutPeak = 0.6;

            var apexPeak = Peaks.First(p => p.RetentionTime == ApexRT);
            int apexIndex = Peaks.IndexOf(apexPeak);
            Peak valleyPeak = null;
            int indexOfValley = 0;

            //go left
            for (int i = apexIndex; i >= 0; i --)
            {
                Peak timepoint = Peaks[i];

                if (valleyPeak == null || timepoint.Intensity < valleyPeak.Intensity)
                {
                    valleyPeak = timepoint;
                    indexOfValley = Peaks.IndexOf(valleyPeak);
                }

                double discriminationFactor =
                    (timepoint.Intensity - valleyPeak.Intensity) / timepoint.Intensity;

                if (discriminationFactor > DiscriminationFactorToCutPeak)
                {
                    var secondValleyPoint = Peaks[indexOfValley - 1];
                    var discriminationFactor2 = (timepoint.Intensity - secondValleyPoint.Intensity) / valleyPeak.Intensity;
                    if (discriminationFactor2 > DiscriminationFactorToCutPeak)
                    {
                        //decide which peak the valley point belongs to 
                        var discriminationLeft = (Peaks[indexOfValley - 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;
                        var discriminationRight = (Peaks[indexOfValley + 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;

                        if (discriminationLeft > discriminationRight)
                        {
                            peaksToRemove.AddRange(Peaks.Where(p => p.RetentionTime <= valleyPeak.RetentionTime));
                        }
                        else
                        {
                            peaksToRemove.AddRange(Peaks.Where(p => p.RetentionTime < valleyPeak.RetentionTime));
                        }
                    }
                    break;
                }
            }

            //go right
            valleyPeak = null;
            for (int i = apexIndex; i < Peaks.Count; i++)
            {
                Peak timepoint = Peaks[i];

                if (valleyPeak == null || timepoint.Intensity < valleyPeak.Intensity)
                {
                    valleyPeak = timepoint;
                    indexOfValley = Peaks.IndexOf(valleyPeak);
                }

                double discriminationFactor =
                    (timepoint.Intensity - valleyPeak.Intensity) / timepoint.Intensity;

                if (discriminationFactor > DiscriminationFactorToCutPeak)
                {
                    var secondValleyPoint = Peaks[indexOfValley + 1];
                    var discriminationFactor2 = (timepoint.Intensity - secondValleyPoint.Intensity) / valleyPeak.Intensity;
                    if (discriminationFactor2 > DiscriminationFactorToCutPeak)
                    {
                        //decide which peak the valley point belongs to 
                        var discriminationLeft = (Peaks[indexOfValley - 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;
                        var discriminationRight = (Peaks[indexOfValley + 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;

                        if (discriminationLeft > discriminationRight)
                        {
                            peaksToRemove.AddRange(Peaks.Where(p => p.RetentionTime > valleyPeak.RetentionTime));
                        }
                        else
                        {
                            peaksToRemove.AddRange(Peaks.Where(p => p.RetentionTime >= valleyPeak.RetentionTime));
                        }
                    }
                    break;
                }
            }

            foreach (var peak in peaksToRemove)
            {
                peak.PeakCurve = null;
                Peaks.Remove(peak);
            }
        }

        //for scipy signal split
        //TODO
        public void ScipySplit(bool doSpline, bool doSmooth, double[] widths, IInterpolation spline = null, float splineRtInterval = 0.005f, int windowSize = 5)
        {
            var peakPoints = GetScipyPeakPoints(doSpline, doSmooth, widths, spline, splineRtInterval, windowSize);
            var valleyPointIndices = new List<int> { 0 };
            for (int i = 0; i < peakPoints.Count - 1; i++)
            {
                var valleyPointIndex = FindValleyPoints((float)peakPoints[i], (float)peakPoints[i + 1]);
                valleyPointIndices.Add(valleyPointIndex);
            }
            valleyPointIndices.Add(Peaks.Count - 1);
        }

        public int FindValleyPoints(float maxima1, float maxima2)
        {
            var rtArray = Peaks.Select(p => (float)p.RetentionTime).ToArray();
            var leftIndex = BinarySearchLower(rtArray, maxima1);
            var rightIndex = BinarySearchHigher(rtArray, maxima2);
            var valleyPointIndex = 0;
            for (int i = leftIndex; i <= rightIndex; i++)
            {
                if (Peaks[i].Intensity < Peaks[leftIndex].Intensity)
                {
                    valleyPointIndex = i;
                }
            }
            return valleyPointIndex;
        }

        public List<double> GetScipyPeakPoints(bool doSpline, bool doSmooth, double[] widths, IInterpolation spline = null, float splineRtInterval = 0.005f, int windowSize = 5)
        {
            var y = Peaks.Select(p => p.Intensity).ToArray();
            var x = Peaks.Select(p => p.RetentionTime).ToArray();
            if (doSpline)
            {
                var data = CalculateSpline(StartRT, EndRT, splineRtInterval, spline);
                y = data.Select(p => p.Item2).ToArray();
                x = data.Select(p => p.Item1).ToArray();
            }
            if (doSmooth)
            {
                double[] imputedY = ImputeMissingValues(spline, out double[] rtArray);
                y = CalculateSavgolSmoothedData(imputedY, windowSize);
            }
            var peaks = Scipy_signal.FindPeaks_cwt(y, widths);
            var peaksRT = peaks.Select(p => x[p]).ToList();
            return peaksRT;
        }

        public void CutPeak_smooth(int windowSize)
        {
            // find out if we need to split this peak by using the discrimination factor
            // this method assumes that the isotope envelopes in a chromatographic peak are already sorted by MS1 scan number
            bool cutThisPeak = false;

            if (Peaks.Count < 5)
            {
                return;
            }
            var peaksToRemove = new List<Peak>();

            //sg smooth the peak
            GetSavgolSmoothedXYData(windowSize);
            (double RetentionTime, double Intensity)[] smoothedData = XYData;

            //DiscriminationFactorToCutPeak is a parameter, default 0.6 in FlashLFQ
            double DiscriminationFactorToCutPeak = 0.6;

            var apexPeak = smoothedData.OrderByDescending(p => p.Intensity).First();
            int apexIndex = smoothedData.IndexOf(apexPeak);
            (double RetentionTime, double Intensity) valleyPeak = (0,0);
            int indexOfValley = 0;

            //go left
            for (int i = apexIndex; i >= 0; i--)
            {
                var timepoint = smoothedData[i];

                if (valleyPeak == (0,0) || timepoint.Intensity < valleyPeak.Intensity)
                {
                    valleyPeak = timepoint;
                    indexOfValley = smoothedData.IndexOf(valleyPeak);
                }

                double discriminationFactor =
                    (timepoint.Intensity - valleyPeak.Intensity) / timepoint.Intensity;

                if (discriminationFactor > DiscriminationFactorToCutPeak)
                {
                    //decide which peak the valley point belongs to 
                    var discriminationLeft = (smoothedData[indexOfValley - 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;
                    var discriminationRight = (smoothedData[indexOfValley + 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;

                    if (discriminationLeft > discriminationRight)
                    {
                        peaksToRemove = Peaks.Where(p => p.RetentionTime <= valleyPeak.RetentionTime).ToList();
                    }
                    else
                    {
                        peaksToRemove = Peaks.Where(p => p.RetentionTime < valleyPeak.RetentionTime).ToList();
                    }
                }
                break;
            }

            //go right
            valleyPeak = (0,0);
            for (int i = apexIndex; i < smoothedData.Length; i++)
            {
                var timepoint = smoothedData[i];

                if (valleyPeak == (0,0) || timepoint.Intensity < valleyPeak.Intensity)
                {
                    valleyPeak = timepoint;
                    indexOfValley = smoothedData.IndexOf(valleyPeak);
                }

                double discriminationFactor =
                    (timepoint.Intensity - valleyPeak.Intensity) / timepoint.Intensity;

                if (discriminationFactor > DiscriminationFactorToCutPeak)
                {
                    //decide which peak the valley point belongs to 
                    var discriminationLeft = (smoothedData[indexOfValley - 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;
                    var discriminationRight = (smoothedData[indexOfValley + 1].Intensity - valleyPeak.Intensity) / valleyPeak.Intensity;

                    if (discriminationLeft > discriminationRight)
                    {
                        peaksToRemove = Peaks.Where(p => p.RetentionTime > valleyPeak.RetentionTime).ToList();
                    }
                    else
                    {
                        peaksToRemove = Peaks.Where(p => p.RetentionTime >= valleyPeak.RetentionTime).ToList();
                    }
                }
                break;
            }
            foreach (var peak in peaksToRemove)
            {
                peak.PeakCurve = null;
                Peaks.Remove(peak);
            }
        }

        public static PeakCurve PeakTracing(double mz, int zeroBasedScanIndex, MsDataScan[] scans, Tolerance tolerance, int binSize, int maxMissedScans, double maxRTRange)
        {
            var allPeaks = Peak.GetAllPeaks(scans, binSize);
            var peakTable = Peak.GetPeakTable(allPeaks, binSize);
            var peak = GetPeakFromScan(mz, peakTable, zeroBasedScanIndex, tolerance, binSize);
            if (peak == null)
            {
                return null;
            }
            var peakCurve = FindPeakCurve(peak, peakTable, scans, null, maxMissedScans, tolerance, binSize, maxRTRange);
            return peakCurve;
        }


        public void Spline(SplineType splineType, DIAparameters diaParam, MsDataScan[] ms1Scans = null, MsDataScan[] ms2Scans = null)
        {
            //var rtIndexMap = GetRtIndexMap(ms1Scans);
            //var rtMap = GetRtMap(ms1Scans, ms2Scans);

            switch (splineType)
            {
                case SplineType.NoSpline:
                    GetRawXYData();
                    break;
                case SplineType.CubicSpline:
                    GetCubicSplineXYData(diaParam.SplineRtInterval);
                    break;
                case SplineType.BSpline:
                    GetBSplineXYData(diaParam.SplineRtInterval, 2);
                    break;
                case SplineType.UmpireBSpline:
                    GetUmpireBSplineData(diaParam.NoPointsPerMin, 2);
                    break;
                //case SplineType.Ms1SpaceBSpline:
                //    GetMs1SpaceBSplineXYData(diaParam.SplineRtInterval, 2, rtMap);
                //    break;
                case SplineType.ScanCycleCubicSpline:
                    GetScanCycleCubicSplineXYData(diaParam.ScanCycleSplineTimeInterval);
                    break;
                case SplineType.SavgolSmoothed:
                    GetSavgolSmoothedXYData(diaParam.SGfilterWindowSize);
                    break;
                case SplineType.CubicSplineSavgolSmoothed:
                    GetCubicSplineSavgolSmoothedXYData(diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                    break;
                case SplineType.ScanCycleCubicSplineSavgolSmoothed:
                    GetScanCycleCubicSplineSavgolSmoothedXYData(diaParam.SGfilterWindowSize, diaParam.ScanCycleSplineTimeInterval);
                    break;
                case SplineType.SavgolSmoothedCubicSpline:
                    GetSavgolSmoothedCubicSplineXYData(diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                    break;
                //case SplineType.Ms1SpaceCubicSpline:
                //    GetMs1SpaceCubicSplineXYData(rtMap, diaParam.SplineRtInterval);
                //    break;
                //case SplineType.Ms1SpaceCubicSplineSavgolSmoothed:
                //    GetMs1SpaceCubicSplineSavgolSmoothedXYData(rtMap, diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                //    break;
                //case SplineType.Ms1SpaceSavgolSmoothedCubicSpline:
                //    GetMs1SpaceSavgolSmoothedCubicSplineXYData(rtMap, diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                //    break;
            }
        }
    }
 
    }
