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
using System.Numerics;
using System.Threading;
using System.Collections.Concurrent;
using System.Threading.Tasks;
using System.Text.RegularExpressions;
using SpectralAveraging;

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
        public int Index { get; set; }
        public List<PrecursorFragmentPair> PFpairs { get; set; }
        public IsotopicEnvelope Envelope { get; set; }

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

        public void GetPrecursorRanks()
        {
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
            , Tolerance mzTolerance, int binSize, double maxRTrange = 2)
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
                if (newPeakCurve.EndRT - newPeakCurve.ApexRT > maxRTrange)
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
                if (newPeakCurve.ApexRT - newPeakCurve.StartRT > maxRTrange)
                {
                    break;
                }
            }
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return newPeakCurve;
        }

        public static PeakCurve FindPeakCurve_inOrder(Peak targetPeak, List<Peak>[] peakTable, MsDataScan[] scans, MzRange isolationWindow, int maxMissedScans
            , Tolerance mzTolerance, int binSize, double maxRTrange = 2)
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
                if (newPeakCurve.EndRT - newPeakCurve.StartRT > maxRTrange)
                {
                    break;
                }
            }

            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return newPeakCurve;
        }
        //the old method searches all the scans: no limit, not just from left to right

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

        public static double CalculatePeakCurveCorr_spline(PeakCurve peakCurve1, PeakCurve peakCurve2, string splineType, double timeInterval)
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

        public static Dictionary<double, List<MsDataScan>> ConstructMs2Group(MsDataScan[] ms2Scans)
        {
            var DIAScanWindowMap = new Dictionary<double, List<MsDataScan>>();
            string pattern = $@"sid=(\d+)";
            foreach (var ms2 in ms2Scans)
            {
                var match = Regex.Match(ms2.ScanFilter, pattern);
                double voltage = double.Parse(match.Groups[1].Value);
                if (!DIAScanWindowMap.ContainsKey(voltage))
                {
                    DIAScanWindowMap[voltage] = new List<MsDataScan>();
                    DIAScanWindowMap[voltage].Add(ms2);
                }
                else
                {
                    DIAScanWindowMap[voltage].Add(ms2);
                }
            }
            return DIAScanWindowMap;
        }

        public static Dictionary<double, List<PeakCurve>> GetMs2PeakCurves(Dictionary<double, List<MsDataScan>> ISDScanVoltageMap, List<Peak>[] allPeaks, DIAparameters DIAparam)
        {
            var ms2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach (var ms2Group in ISDScanVoltageMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                Peak.GetAllPeaks(allPeaks, ms2scans);
                var allMs2Peaks = allPeaks.Where(v => v != null && v.First().Voltage == ms2Group.Key).ToArray();
                var rankedMs2Peaks = allMs2Peaks.SelectMany(p => p).OrderByDescending(p => p.Intensity).ToList();
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, DIAparam.PeakSearchBinSize, out Dictionary<int, int> scanIndexMap);
                ms2PeakCurves[ms2Group.Key] = new List<PeakCurve>();
                foreach (var peak in rankedMs2Peaks)
                {
                    if (peak.PeakCurve == null)
                    {
                        var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms2PeakTable, ms2scans, ms2scans[0].IsolationRange,
                            DIAparam.MaxNumMissedScan, DIAparam.Ms2PeakFindingTolerance, DIAparam.PeakSearchBinSize, DIAparam.MaxRTRange);
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            ms2PeakCurves[ms2Group.Key].Add(newPeakCurve);
                        }
                    }
                }
            }
            return ms2PeakCurves;
        }

        public static Dictionary<double, List<PeakCurve>> GetMs2PeakCurves_decon(Dictionary<double, List<MsDataScan>> ISDScanVoltageMap, List<Peak>[] allPeaks, CommonParameters commonParam)
        {
            var ms2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach (var ms2Group in ISDScanVoltageMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                var allMasses = new List<DeconvolutedMass>();
                for (int i = 0; i < ms2scans.Length; i++)
                {
                    var envelopes = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scans[i], commonParam);

                    foreach (var envelope in envelopes)
                    {
                        var charge = envelope.Charge;
                        double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                        double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                        var mass = new DeconvolutedMass(envelope, charge, ms2scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                            ms2scans[i].OneBasedScanNumber, i);
                        allMasses.Add(mass);
                    }
                }
                allMasses = allMasses.OrderByDescending(p => p.HighestPeakIntensity).ToList();

                Peak.GetAllPeaks(allPeaks, ms2scans);
                var allMs2Peaks = allPeaks.Where(v => v != null && v.First().Voltage == ms2Group.Key).ToArray();
                var rankedMs2Peaks = allMs2Peaks.SelectMany(p => p).OrderByDescending(p => p.Intensity).ToList();
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, commonParam.DIAparameters.PeakSearchBinSize, out Dictionary<int, int> scanIndexMap);
                ms2PeakCurves[ms2Group.Key] = new List<PeakCurve>();
                foreach (var mass in allMasses)
                {
                    var peak = PeakCurve.GetPeakFromScan(mass.HighestPeakMz, ms2PeakTable, mass.ZeroBasedScanIndex, new PpmTolerance(0),
                        commonParam.DIAparameters.PeakSearchBinSize);
                    if (peak.PeakCurve == null)
                    {
                        var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms2PeakTable, ms2scans, null, commonParam.DIAparameters.MaxNumMissedScan,
                            commonParam.DIAparameters.Ms1PeakFindingTolerance, commonParam.DIAparameters.PeakSearchBinSize, commonParam.DIAparameters.MaxRTRange);
                        newPeakCurve.MonoisotopicMass = mass.MonoisotopicMass;
                        newPeakCurve.Charge = mass.Charge;
                        newPeakCurve.Envelope = mass.Envelope;
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            ms2PeakCurves[ms2Group.Key].Add(newPeakCurve);
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
                        DIAparam.MaxNumMissedScan, DIAparam.Ms2PeakFindingTolerance, DIAparam.PeakSearchBinSize, DIAparam.MaxRTRange);
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
                            var allms2PeaksForThisScan = allPeaks[scan.OneBasedScanNumber];
                            var ms1curve = scan.PrecursorPeakCurve;
                            var DIApeaks = new List<Peak>();
                            //var corrList = new List<double>();

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
                                                double corr = PeakCurve.CalculatePeakCurveCorr(ms1curve, ms2curve);
                                                if (corr > DIAparam.CorrelationCutOff)
                                                {
                                                    DIApeaks.Add(peak);
                                                    //corrList.Add(corr);
                                                }
                                            }
                                        }
                                    }
                                }
                                //if (corrList.Count > DIAparam.FragmentRankCutOff)
                                //{
                                //    var indices = corrList.Select((value, index) => new { Value = value, Index = index })
                                //                  .OrderByDescending(x => x.Value).Take(DIAparam.FragmentRankCutOff).Select(x => x.Index).ToList();
                                //    DIApeaks = DIApeaks.Select((value, index) => new {Value = value, Index = index}).Where(p => indices.Contains(p.Index)).
                                //    OrderBy(p => p.Value.Mz).Select(p => p.Value).ToList();
                                //}
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

        public static List<Ms2ScanWithSpecificMass>[] GetPseudoMs2Scans_decon(List<Ms2ScanWithSpecificMass>[] scansWithPrecursor, CommonParameters commonParam, DIAparameters DIAparam, List<Peak>[] allPeaks)
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
                            var allms2PeaksForThisScan = allPeaks[scan.OneBasedScanNumber];
                            var ms1curve = scan.PrecursorPeakCurve;
                            var DIApeaks = new List<Peak>();

                            if (ms1curve.Peaks.Count > 4)
                            {
                                foreach (var peak in allms2PeaksForThisScan)
                                {
                                    var ms2curve = peak.PeakCurve;
                                    if (ms2curve == null)
                                    {
                                        continue;
                                    }
                                    if (ms2curve.ApexRT >= ms1curve.StartRT && ms2curve.ApexRT <= ms1curve.EndRT)
                                    {
                                        if (Math.Abs(ms2curve.ApexRT - ms1curve.ApexRT) < DIAparam.ApexRtTolerance)
                                        {
                                            var overlap = PeakCurve.CalculateRTOverlapRatio(ms1curve, ms2curve);
                                            if (overlap >= DIAparam.OverlapRatioCutOff)
                                            {
                                                double corr = PeakCurve.CalculatePeakCurveCorr_spline(ms1curve, ms2curve, "cubic", 0.005);
                                                if (corr > DIAparam.CorrelationCutOff)
                                                {
                                                    DIApeaks.Add(peak);
                                                }
                                            }
                                        }
                                    }
                                }
                                if (DIApeaks.Count > 0)
                                {
                                    var mzs = new double[] { 1 };
                                    var intensities = new double[] { double.MaxValue };
                                    var spectrum = new MzSpectrum(mzs, intensities, false);
                                    var neutralExperimentalFragments = DIApeaks.Select(p => p.PeakCurve.Envelope).ToArray();
                                    //MsDataScan newScan = new MsDataScan(spectrum, scan.TheScan.OneBasedScanNumber, scan.TheScan.MsnOrder, scan.TheScan.IsCentroid,
                                    //                        scan.TheScan.Polarity, scan.TheScan.RetentionTime, scan.TheScan.ScanWindowRange, scan.TheScan.ScanFilter, scan.TheScan.MzAnalyzer,
                                    //                        scan.TheScan.TotalIonCurrent, scan.TheScan.InjectionTime, scan.TheScan.NoiseData, scan.TheScan.NativeId,
                                    //                        scan.TheScan.SelectedIonMZ, scan.TheScan.SelectedIonChargeStateGuess, scan.TheScan.SelectedIonIntensity,
                                    //                        scan.TheScan.IsolationMz, scan.TheScan.IsolationWidth, scan.TheScan.DissociationType, null,
                                    //                        scan.TheScan.SelectedIonMonoisotopicGuessMz, scan.TheScan.HcdEnergy, scan.TheScan.ScanDescription);
                                    var newScan = new MsDataScan(spectrum, scan.TheScan.OneBasedScanNumber, 2, true, Polarity.Positive, scan.TheScan.RetentionTime, new MzRange(mzs.Min(), mzs.Max()),
                                        scan.TheScan.ScanFilter,MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
                                    var newScanWithPrecursor = new Ms2ScanWithSpecificMass(newScan, scan.PrecursorMonoisotopicPeakMz, scan.PrecursorCharge, scan.FullFilePath, commonParam,
                                        neutralExperimentalFragments);
                                    newScansWithPre[i].Add(newScanWithPrecursor);
                                }
                            }
                        }
                    }
                });
            return newScansWithPre;
        }

        public static List<Ms2ScanWithSpecificMass>[] GetPseudoMs2Scans_PFgroup(List<Ms2ScanWithSpecificMass>[] scansWithPrecursor, CommonParameters commonParam, DIAparameters DIAparam, 
            List<Peak>[] allPeaks, Dictionary<MzRange, List<PeakCurve>> allMs2PeakCurves, MsDataFile myMSDataFile)
        {
            var newScansWithPre = new List<Ms2ScanWithSpecificMass>[scansWithPrecursor.Length];
            var pfGroups = new List<PrecursorFragmentsGroup>[scansWithPrecursor.Length];
            Parallel.ForEach(Partitioner.Create(0, scansWithPrecursor.Length), new ParallelOptions { MaxDegreeOfParallelism = 15 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        pfGroups[i] = new List<PrecursorFragmentsGroup>();
                        foreach (var scan in scansWithPrecursor[i])
                        {
                            var allms2PeaksForThisScan = allPeaks[scan.OneBasedScanNumber - 1];
                            var ms1curve = scan.PrecursorPeakCurve;
                            var preFragGroup = new PrecursorFragmentsGroup(ms1curve);
                            preFragGroup.MonoPeakMz = scan.PrecursorMonoisotopicPeakMz;
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
                                                double corr = PeakCurve.CalculatePeakCurveCorr(ms1curve, ms2curve);
                                                if (corr > DIAparam.CorrelationCutOff)
                                                {
                                                    var PFpair = new PrecursorFragmentPair(ms1curve, ms2curve, corr);
                                                    lock (ms2curve.PFpairs)
                                                    {
                                                        ms2curve.PFpairs.Add(PFpair);
                                                    }
                                                    preFragGroup.PFpairs.Add(PFpair);
                                                }
                                            }
                                        }
                                    }
                                }
                                //if (preFragGroup.PFpairs.Count > DIAparam.FragmentRankCutOff)
                                //{
                                //    var filtered = preFragGroup.PFpairs.OrderByDescending(pair => pair.Correlation).Take(DIAparam.FragmentRankCutOff);
                                //    preFragGroup.PFpairs = filtered.ToList();
                                //}
                                if (preFragGroup.PFpairs.Count > 0)
                                {
                                    preFragGroup.PFpairs = preFragGroup.PFpairs.OrderBy(pair => pair.FragmentPeakCurve.AveragedMz).ToList();
                                    pfGroups[i].Add(preFragGroup);
                                }
                            }
                        }
                    }
                });

            var allMatchedMs2 = allMs2PeakCurves.SelectMany(p => p.Value).Where(pc => pc.PFpairs.Count > 0).ToList();
            foreach (var ms2curve in allMatchedMs2)
            {
                ms2curve.GetPrecursorRanks();
            }
            int oneBasedScanNum = 1;
            for (var i = 0; i < pfGroups.Length; i++)
            {
                foreach (var group in pfGroups[i])
                {
                    group.PFpairs = group.PFpairs.Where(pf => pf.PrecursorRank <= DIAparam.PrecursorRankCutOff).ToList();
                    group.GetNumberOfHighCorrFragments(DIAparam);
                }
                pfGroups[i] = pfGroups[i].Where(pf => pf.PFpairs.Count > 0 && pf.NumHighCorrFragments >= DIAparam.NumHighCorrFragments).ToList();
                newScansWithPre[i] = new List<Ms2ScanWithSpecificMass>();
                foreach (var pfGroup in pfGroups[i])
                {
                    var mzs = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedMz).ToArray();
                    var intensities = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedIntensity).ToArray();
                    var spectrum = new MzSpectrum(mzs, intensities, false);
                    var newMs2Scan = new MsDataScan(spectrum, oneBasedScanNum, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                                MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
                    var neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, commonParam);
                    var charge = pfGroup.PrecursorPeakCurve.Charge;
                    var highestPeakMz = pfGroup.PrecursorPeakCurve.AveragedMz;
                    Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, pfGroup.MonoPeakMz, charge
                        , myMSDataFile.FilePath, commonParam, neutralExperimentalFragments);
                    newScansWithPre[i].Add(scanWithprecursor);
                    oneBasedScanNum++;
                }
            }
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
                        scan.PrecursorPeakCurve = FindPeakCurve(precursorPeak, ms1PeakTable, ms1scans, null, DIAparam.MaxNumMissedScan, DIAparam.Ms1PeakFindingTolerance, 
                            DIAparam.PeakSearchBinSize, DIAparam.MaxRTRange);
                        scan.PrecursorPeakCurve.MonoisotopicMass = scan.PrecursorMass;
                        scan.PrecursorPeakCurve.Charge = scan.PrecursorCharge;
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
