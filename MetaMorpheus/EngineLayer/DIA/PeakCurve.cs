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
        public int MsLevel {  get; set; }
        public MzRange IsolationRange { get; set; }
        public double MonoisotopicMass { get; set; }
        public int Charge { get; set; }
        public List<Peak> IsotopePeaks { get; set; }
        public List<double> IsotopePeaksMz { get; set; }
        public double StartRT => Peaks.Select(p => p.RetentionTime).OrderBy(t => t).First();
        public double EndRT => Peaks.Select(p => p.RetentionTime).OrderByDescending(t => t).First();
        public double StartMz {  get; set; }
        public double EndMz { get; set; }
        public MzRange MzRange => new MzRange(StartMz, EndMz);
        public double ApexRT => Peaks.OrderByDescending(p => p.Intensity).First().RetentionTime;
        public double TotalIntensity => Peaks.Sum(p => p.Intensity);
        public double AveragedMz => AverageMz();
        public double AveragedIntensity => AverageIntensity();
        public LinearSpline LinearSpline { get; set; }
        public CubicSpline CubicSpline { get; set; }    
        public int Index {  get; set; }
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
                var peak = GetPeakFromScan(newPeakCurve.AveragedMz, peakTable, t, mzTolerance, binSize);

                if (peak == null)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    //if (peak.PeakCurve == null || Math.Abs(peak.Mz - newPeakCurve.AveragedMz) < Math.Abs(peak.Mz - peak.PeakCurve.AveragedMz))
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
                if (newPeakCurve.ApexRT - newPeakCurve.StartRT > maxRTrange)
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
            var targetPeak = GetPeakFromScan(targetMz,  peakTable, zeroBasedScanIndex, mzTolerance, binSize);
            var pc = FindPeakCurve(targetPeak, peakTable, scans, isolationWindow, maxMissedScans, mzTolerance, binSize);
            return pc;
        }

        public static double CalculatePeakCurveCorr(PeakCurve peakCurve1,  PeakCurve peakCurve2)
        {
            var peakList1 = peakCurve1.Peaks.ToArray();
            var peakList2 = peakCurve2.Peaks.ToArray();

            var intensityPair = new List<(double, double)>();
            //for plot
            var rtPair = new List<(double, double)>();

            if(peakList1.Length < peakList2.Length)
            {
                foreach(var peak in peakList1)
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
            if (peakCurve2.Peaks.Count <5 || peakCurve1.Peaks.Count < 5)
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

        public static List<PeakCurve> GetMs1PeakCurves(MsDataScan[] allMs1Scans, List<Peak>[] ms1PeakTable, DIAparameters DIAparameters, CommonParameters commonParameters)
        {
            //Get all precursors
            var allPrecursors = new List<DeconvolutedMass>();
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], commonParameters.PrecursorDeconvolutionParameters);
                foreach (var envelope in envelopes)
                {
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
                    allPrecursors.Add(precursor);
                }
            }
            //sort precursors by envelope totalintensity
            var preGroups = allPrecursors.GroupBy(p => new { p.MonoisotopicMass, p.Charge}).ToList();
            var referencePrecursors = preGroups.Select(g => g.OrderByDescending(p => p.Envelope.TotalIntensity).First()).ToList();
            var allMs1PeakCurves = new List<PeakCurve>();

            //debug
            var peak2 = GetPeakFromScan(626.997, ms1PeakTable, 1391, new PpmTolerance(5), 100);
            var pc2 = FindPeakCurve(peak2, ms1PeakTable, allMs1Scans, null, 2, new PpmTolerance(10), 100);
            var peak3 = GetPeakFromScan(627.331, ms1PeakTable, 1391, new PpmTolerance(5), 100);
            var pc3 = FindPeakCurve(peak3, ms1PeakTable, allMs1Scans, null, 2, new PpmTolerance(5), 100);
            var peak4 = GetPeakFromScan(627.666, ms1PeakTable, 1391, new PpmTolerance(5), 100);
            var pc4 = FindPeakCurve(peak4, ms1PeakTable, allMs1Scans, null, 2, new PpmTolerance(5), 100);

            //Find precursor XIC
            foreach (var precursor in allPrecursors)
            {
                var highestPeak = PeakCurve.GetPeakFromScan(precursor.HighestPeakMz, ms1PeakTable, precursor.ZeroBasedScanIndex, new PpmTolerance(0), DIAparameters.PeakSearchBinSize);
                if (Math.Abs(highestPeak.Mz - 626.997) < 0.005 && highestPeak.ScanNumber > 4170 && highestPeak.ScanNumber < 4190)
                {
                    int stop = 0;
                    var peak = GetPeakFromScan(626.663, ms1PeakTable, 1391, new PpmTolerance(5), 100);
                    var pc = FindPeakCurve(peak, ms1PeakTable, allMs1Scans, null, 2, new PpmTolerance(10), 100);
                    
                }
                if (highestPeak.PeakCurve == null)
                {
                    //TODO: label the peaks that have been included in a PeakCurve
                    var newPeakCurve = PeakCurve.FindPeakCurve(highestPeak, ms1PeakTable, allMs1Scans, null, DIAparameters.MaxNumMissedScan,
                        DIAparameters.Ms1PeakFindingTolerance, DIAparameters.PeakSearchBinSize);
                    newPeakCurve.MonoisotopicMass = precursor.MonoisotopicMass;
                    newPeakCurve.Charge = precursor.Charge;
                    newPeakCurve.StartMz = precursor.Envelope.Peaks.Min(p => p.mz);
                    newPeakCurve.EndMz = precursor.Envelope.Peaks.Max(p => p.mz);
                    allMs1PeakCurves.Add(newPeakCurve);
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

        public static List<GenericChart> VisualizePFgroups(MsDataFile dataFile, List<PsmFromTsv> psms, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var diaEngine2 = new DIAEngine2(dataFile, commonParameters, diaParam);
            diaEngine2.Ms1PeakIndexing();
            diaEngine2.ConstructMs2Group();
            diaEngine2.GetMs1PeakCurves();
            diaEngine2.GetMs2PeakCurves();

            var chartList = new List<GenericChart>();
            foreach (var psm in psms)
            {
                var precursorPeakCurve = diaEngine2.Ms1PeakCurves.Values.SelectMany(v => v).Where(pc => pc.Index == psm.Ms2ScanNumber).First();
                var ms2curves = diaEngine2.Ms2PeakCurves.Where(d => precursorPeakCurve.AveragedMz +1 > d.Key.min && precursorPeakCurve.AveragedMz -1 < d.Key.max)
                .First().Value;
                var pfgroup = DIAEngine2.GroupPrecursorFragments(precursorPeakCurve, ms2curves, diaParam);
                var plots = new List<GenericChart>();
                var normalizedIntensity = precursorPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                var precursorPlot = Chart2D.Chart.Line<double, double, string>(
                    x: precursorPeakCurve.Peaks.Select(p => p.RetentionTime),
                    y: normalizedIntensity).WithTraceInfo($"precursor_{Math.Round(precursorPeakCurve.AveragedMz, 3)}_{psm.DecoyContamTarget}").WithMarkerStyle(Color: Color.fromString("red"));
                plots.Add(precursorPlot);

                var matchedIonPeakCurves = new List<PeakCurve>();
                var matchedIonMzs = psm.MatchedIons.Select(i => i.Mz);
                foreach (double mz in matchedIonMzs)
                {
                    try
                    {
                        var fragmentPeakCurve = pfgroup.PFpairs.Where(p => Math.Abs(p.FragmentPeakCurve.AveragedMz - mz) < 0.005).First().FragmentPeakCurve;
                        matchedIonPeakCurves.Add(fragmentPeakCurve);
                    }
                    catch
                    {

                    }
                }
                foreach (var pf in pfgroup.PFpairs)
                {
                    if (!matchedIonPeakCurves.Contains(pf.FragmentPeakCurve))
                    {
                        var norm2 = pf.FragmentPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                        var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                                                       x: pf.FragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                                                       y: norm2).WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{pf.Correlation}_{pf.FragmentRank}")
                                                       .WithMarkerStyle(Color: Color.fromString("green"));
                        plots.Add(fragmentPlot2);
                    }
                    else
                    {
                        var norm2 = pf.FragmentPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                        var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                                                       x: pf.FragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                                                       y: norm2).WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{pf.Correlation}_{pf.FragmentRank}")
                                                       .WithMarkerStyle(Color: Color.fromString("blue"));
                        plots.Add(fragmentPlot2);
                    }
                }
                var combinedPlot = Chart.Combine(plots).WithLayout(Layout.init<string>(Title: Title.init($"{psm.Ms2ScanNumber}_{psm.Score}_{psm.DecoyContamTarget}")));
                chartList.Add(combinedPlot);
            }   
            return chartList;
        }
    }
}
