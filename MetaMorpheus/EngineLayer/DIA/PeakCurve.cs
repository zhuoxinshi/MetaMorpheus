using Chemistry;
using Easy.Common.Extensions;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using Plotly.NET;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
using System.Numerics;
using System.Threading;
using ThermoFisher.CommonCore.Data.Business;
using System.Collections.Concurrent;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PeakCurve
    {
        public PeakCurve(List<Peak> peaks, int msLevel, MzRange isolationRange, double mass = double.NaN, int charge = 0, double startMz = 0, double endMz = 0)
        {
            Peaks = peaks;
            MsLevel = msLevel;
            IsolationRange = isolationRange;
            MonoisotopicMass = mass;
            Charge = charge;
            StartMz = startMz;
            EndMz = endMz;
        }

        public List<Peak> Peaks { get; set; }
        public int MsLevel { get; set; }
        public MzRange IsolationRange { get; set; }
        public double MonoisotopicMass { get; set; }
        public int Charge { get; set; }
        public List<Peak> IsotopePeaks { get; set; }
        public double StartRT => Peaks.Select(p => p.RetentionTime).OrderBy(t => t).First();
        public double EndRT => Peaks.Select(p => p.RetentionTime).OrderByDescending(t => t).First();
        public double StartMz { get; set; }
        public double EndMz { get; set; }
        public MzRange MzRange => new MzRange(StartMz, EndMz);
        public double ApexRT => Peaks.OrderByDescending(p => p.Intensity).First().RetentionTime;
        public double TotalIntensity => Peaks.Sum(p => p.Intensity);
        public double AveragedMz => AverageMz();
        public double AveragedIntensity => AverageIntensity();
        public LinearSpline LinearSpline { get; set; }
        public CubicSpline CubicSpline { get; set; }

        public double AverageMz()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            double sumMz = Peaks.Sum(p => p.Mz);
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

        public void GetLinearSpline()
        {
            var sortedPeaks = Peaks.OrderBy(p => p.RetentionTime).ToList();
            var rtArray = sortedPeaks.Select(p => p.RetentionTime).ToArray();
            var intensityArray = sortedPeaks.Select(p => p.Intensity).ToArray();
            var linearSpline = LinearSpline.InterpolateSorted(rtArray, intensityArray);
            this.LinearSpline = linearSpline;

            //plot
            //var rtSeq = new List<double>();
            //var intensities = new List<double>();
            //for (double i = StartRT; i < EndRT; i += 0.0001)
            //{
            //    rtSeq.Add(i);
            //}
            //foreach(var rt in rtSeq)
            //{
            //    intensities.Add(linearSpline.Interpolate(rt));
            //}

            //var plot1 = Chart2D.Chart.Point<double, double, string>(
            //    x: Peaks.Select(p => p.RetentionTime),
            //    y: Peaks.Select(p => p.Intensity)).WithTraceInfo("original").WithMarkerStyle(Color: Color.fromString("red"));
            //var plot2 = Chart2D.Chart.Point<double, double, string>(
            //    x: rtSeq,
            //    y: intensities).WithTraceInfo("interpolate").WithMarkerStyle(Color: Color.fromString("blue"));
            //var combinedPlot = Chart.Combine(new[] { plot1, plot2 });
            //combinedPlot.Show();
        }

        public void GetCubicSpline()
        {
            var sortedPeaks = Peaks.OrderBy(p => p.RetentionTime).ToList();
            var rtArray = sortedPeaks.Select(p => p.RetentionTime).ToArray();
            var intensityArray = sortedPeaks.Select(p => p.Intensity).ToArray();
            var cubicSpline = CubicSpline.InterpolateAkima(rtArray, intensityArray);
            this.CubicSpline = cubicSpline;
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

        public static int BinarySearchForIndexedPeak(List<Peak> peakList, int zeroBasedScanIndex)
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

        public static PeakCurve FindPeakCurve(Peak targetPeak, List<Peak>[] peakTable, MsDataScan[] scans, MzRange isolationWindow, int maxMissedScans
            , Tolerance mzTolerance, int binSize)
        {
            var xic = new List<Peak>();
            xic.Add(targetPeak);
            PeakCurve newPeakCurve = new PeakCurve(xic, targetPeak.MsLevel, isolationWindow);
            targetPeak.PeakCurve = newPeakCurve;

            // go right
            int missedScans = 0;
            for (int t = targetPeak.ZeroBasedScanIndex + 1; t < scans.Length; t++)
            {
                var peak = GetPeakFromScan(newPeakCurve.AveragedMz, peakTable, t, mzTolerance, binSize);

                if (peak == null)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    //if(peak.PeakCurve == null || Math.Abs(peak.Mz - newPeakCurve.AveragedMz) < Math.Abs(peak.Mz - peak.PeakCurve.AveragedMz))
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
            }

            // go left
            missedScans = 0;
            for (int t = targetPeak.ZeroBasedScanIndex - 1; t >= 0; t--)
            {
                var peak = GetPeakFromScan(newPeakCurve.AveragedMz, peakTable, t, mzTolerance, binSize);

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
            }
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return newPeakCurve;
        }

        public static PeakCurve FindPeakCurve_mz(double targetMz, int zeroBasedScanIndex, List<Peak>[] peakTable, MsDataScan[] scans, MzRange isolationWindow, int maxMissedScans
            , Tolerance mzTolerance, int binSize)
        {
            var targetPeak = GetPeakFromScan(targetMz, peakTable, zeroBasedScanIndex, mzTolerance, binSize);
            var pc = FindPeakCurve(targetPeak, peakTable, scans, isolationWindow, maxMissedScans, mzTolerance, binSize);
            return pc;
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
                double corr = Correlation.Pearson(intensityPair.Select(pair => pair.Item1), intensityPair.Select(pair => pair.Item2));

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

        public static double CalculateCorr_spline(PeakCurve peakCurve1, PeakCurve peakCurve2, string splineType, double timeInterval)
        {
            if (peakCurve2.Peaks.Count < 5 || peakCurve1.Peaks.Count < 5)
            {
                return 0;
            }
            var startRT = Math.Max(peakCurve1.StartRT, peakCurve2.StartRT);
            var endRT = Math.Min(peakCurve1.EndRT, peakCurve2.EndRT);
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
                double corr = Correlation.Pearson(intensities1, intensities2);
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
                double corr = Correlation.Pearson(intensities1, intensities2);
                return corr;
            }
            return 0;
        }

        public static Dictionary<(double min, double max), List<MsDataScan>> ConstructMs2Group(MsDataScan[] ms2Scans)
        {
            var DIAScanWindowMap = new Dictionary<(double min, double max), List<MsDataScan>>();
            foreach (var ms2 in ms2Scans)
            {
                (double min, double max) range = new(ms2.IsolationRange.Minimum, ms2.IsolationRange.Maximum);
                if (!DIAScanWindowMap.ContainsKey(range))
                {
                    DIAScanWindowMap[range] = new List<MsDataScan>();
                    DIAScanWindowMap[range].Add(ms2);
                }
                else
                {
                    DIAScanWindowMap[range].Add(ms2);
                }
            }
            return DIAScanWindowMap;
        }

        public static Dictionary<MzRange, List<PeakCurve>> GetMs2PeakCurves(MsDataScan[] allMs2scans, List<Peak>[] allPeaks, DIAparameters DIAparam)
        {
            var ms2PeakCurves = new Dictionary<MzRange, List<PeakCurve>>();
            var DIAScanWindowMap = ConstructMs2Group(allMs2scans);
            foreach (var ms2Group in DIAScanWindowMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                MzRange range = ms2scans[0].IsolationRange;
                var allMs2Peaks = allPeaks.Where(v => v != null && v.FirstOrDefault().IsolationRange != null &&
                 v.FirstOrDefault().IsolationRange.Maximum == range.Maximum && v.FirstOrDefault().IsolationRange.Minimum == range.Minimum).ToArray();
                var rankedMs2Peaks = allMs2Peaks.SelectMany(p => p).OrderByDescending(p => p.Intensity).ToList();
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, DIAparam.PeakSearchBinSize, out Dictionary<int, int> scanIndexMap);
                ms2PeakCurves[range] = new List<PeakCurve>();
                foreach (var peak in rankedMs2Peaks)
                {
                    if (peak.PeakCurve == null)
                    {
                        var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms2PeakTable, ms2scans, ms2scans[0].IsolationRange,
                            DIAparam.MaxNumMissedScan, DIAparam.Ms2PeakFindingTolerance, DIAparam.PeakSearchBinSize);
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            ms2PeakCurves[range].Add(newPeakCurve);
                        }
                    }
                }
            }
            return ms2PeakCurves;
        }

        public static List<PeakCurve> GetAllMs1PeakCurves(MsDataScan[] ms1scans, List<Peak>[] allPeaks, List<Peak>[] ms1PeakTable, DIAparameters DIAparam)
        {
            var allMs1PeakCurves = new List<PeakCurve>();
            var allMs1Peaks = allPeaks.Where(v => v != null && v.FirstOrDefault().MsLevel == 1).ToArray();
            var rankedMs1Peaks = allMs1Peaks.SelectMany(p => p).OrderByDescending(p => p.Intensity).ToList();
            foreach (var peak in rankedMs1Peaks)
            {
                if (peak.PeakCurve == null)
                {
                    var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms1PeakTable, ms1scans, ms1scans[0].IsolationRange,
                        DIAparam.MaxNumMissedScan, DIAparam.Ms2PeakFindingTolerance, DIAparam.PeakSearchBinSize);
                    if (newPeakCurve.Peaks.Count > 4)
                    {
                        allMs1PeakCurves.Add(newPeakCurve);
                    }
                }
            }
            return allMs1PeakCurves;
        }

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

        public static List<Ms2ScanWithSpecificMass>[] GetPseudoMs2Scans(List<Ms2ScanWithSpecificMass>[] scansWithPrecursor, CommonParameters commonParam, DIAparameters DIAparam, List<Peak>[] allPeaks)
        {
            var newScansWithPre = new List<Ms2ScanWithSpecificMass>[scansWithPrecursor.Length];
            Parallel.ForEach(Partitioner.Create(0, scansWithPrecursor.Length), new ParallelOptions { MaxDegreeOfParallelism = 15 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        newScansWithPre[i] = new List<Ms2ScanWithSpecificMass>();
                        foreach (var scan in scansWithPrecursor[i])
                        {
                            var allms2PeaksForThisScan = allPeaks[scan.OneBasedScanNumber - 1];
                            var ms1curve = scan.PrecursorPeakCurve;
                            var DIApeaks = new List<Peak>();
                            var corrList = new List<double>();

                            if (ms1curve.Peaks.Count > 4)
                            {
                                foreach (var peak in allms2PeaksForThisScan)
                                {
                                    var ms2curve = peak.PeakCurve;
                                    if (ms2curve.ApexRT >= ms1curve.StartRT && ms2curve.ApexRT <= ms1curve.EndRT)
                                    {
                                        if (Math.Abs(ms2curve.ApexRT - ms1curve.ApexRT) < DIAparam.ApexRtTolerance)
                                        {
                                            var overlap = PeakCurve.CalculateRTOverlapRatio(ms1curve, ms2curve);
                                            if (overlap >= DIAparam.OverlapRatioCutOff)
                                            {
                                                double corr = PeakCurve.CalculateCorr_spline(ms1curve, ms2curve, "cubic", 0.01);
                                                if (corr > DIAparam.CorrelationCutOff)
                                                {
                                                    DIApeaks.Add(peak);
                                                    corrList.Add(corr);
                                                }
                                            }
                                        }
                                    }
                                }
                                if (corrList.Count > DIAparam.FragmentRankCutOff)
                                {
                                    var indices = corrList.Select((value, index) => new { Value = value, Index = index })
                                                  .OrderByDescending(x => x.Value).Take(DIAparam.FragmentRankCutOff).Select(x => x.Index).ToList();
                                    DIApeaks = DIApeaks.Select((value, index) => new {Value = value, Index = index}).Where(p => indices.Contains(p.Index)).
                                    OrderBy(p => p.Value.Mz).Select(p => p.Value).ToList();
                                }
                                if (DIApeaks.Count > 0)
                                {
                                    var newSpectrum = new MzSpectrum(DIApeaks.Select(p => p.Mz).ToArray(), DIApeaks.Select(p => p.Intensity).ToArray(), false);
                                    MsDataScan newScan = new MsDataScan(newSpectrum, scan.TheScan.OneBasedScanNumber, scan.TheScan.MsnOrder, scan.TheScan.IsCentroid,
                                                            scan.TheScan.Polarity, scan.TheScan.RetentionTime, scan.TheScan.ScanWindowRange, scan.TheScan.ScanFilter, scan.TheScan.MzAnalyzer,
                                                            scan.TheScan.TotalIonCurrent, scan.TheScan.InjectionTime, scan.TheScan.NoiseData, scan.TheScan.NativeId,
                                                            scan.TheScan.SelectedIonMZ, scan.TheScan.SelectedIonChargeStateGuess, scan.TheScan.SelectedIonIntensity,
                                                            scan.TheScan.IsolationMz, scan.TheScan.IsolationWidth, scan.TheScan.DissociationType, null,
                                                            scan.TheScan.SelectedIonMonoisotopicGuessMz, scan.TheScan.HcdEnergy, scan.TheScan.ScanDescription);
                                    var newScanWithPrecursor = new Ms2ScanWithSpecificMass(newScan, scan.PrecursorMonoisotopicPeakMz, scan.PrecursorCharge, scan.FullFilePath, commonParam);
                                    newScansWithPre[i].Add(newScanWithPrecursor);
                                }
                            }
                        }
                    }
                });
            return newScansWithPre;
        }

        public static void GetPrecursorPeakCurve(List<Ms2ScanWithSpecificMass>[] scansWithPrecursor, MsDataScan[] ms1scans, List<Peak>[] allPeaks,
            List<Peak>[] ms1PeakTable, DIAparameters DIAparam, Dictionary<int, int> scanIndexMap)
        {
            for (int i = 0; i < scansWithPrecursor.Length; i++)
            {
                foreach (var scan in scansWithPrecursor[i])
                {
                    var preScan = ms1scans.Where(s => s.OneBasedScanNumber == scan.OneBasedPrecursorScanNumber).First();
                    var precursorPeak = GetPeakFromScan(scan.HighestPeakMz, ms1PeakTable, scanIndexMap[preScan.OneBasedScanNumber], new PpmTolerance(0), DIAparam.PeakSearchBinSize);
                    if (precursorPeak.PeakCurve == null)
                    {
                        scan.PrecursorPeakCurve = FindPeakCurve(precursorPeak, ms1PeakTable, ms1scans, null, DIAparam.MaxNumMissedScan, DIAparam.Ms1PeakFindingTolerance, DIAparam.PeakSearchBinSize);
                    }
                    else
                    {
                        scan.PrecursorPeakCurve = precursorPeak.PeakCurve;
                    }
                }
            }
        }
    }
}
