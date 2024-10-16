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
using Plotly.NET.TraceObjects;

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

        public List<PeakRidge> PeakRidgeList { get; set; }
        public List<(float rt1, float rt2, float rt3)> PeakRegionList { get; set; }
        public List<List<float>> NoRidgeRegion { get; set; }

        public WaveletMassDetector WaveletMassDetector { get; set; }
        public CwtParameters CwtParameters { get; set; }
        public List<(float, float)> SmoothedData { get; set; }
        public double NL { get; set; }

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

        public void CalculateNL()
        {
            NL = Peaks.Min(p => p.Intensity);
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

        public void Interpolte_cubic(int NoPeakPerMin)
        {
            if (CubicSpline == null)
            {
                GetCubicSpline();
            }
            var rtSeq = new List<float>();
            var RTwindow = EndRT - StartRT;
            var totalNumPoints = (int)(RTwindow * NoPeakPerMin);
            float timeInterval = (float)(RTwindow / (totalNumPoints - 1));
            for (float i = (float)StartRT; i < (float)EndRT; i += timeInterval)
            {
                rtSeq.Add(i);
            }
            SmoothedData = new List<(float, float)>();
            for (int i = 0; i < rtSeq.Count; i++)
            {
                SmoothedData.Add((rtSeq[i], (float)CubicSpline.Interpolate(rtSeq[i])));
            }
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

        public void DetectPeakRegions()
        {
            PeakRidgeList = new List<PeakRidge>();
            PeakRegionList = new List<(float rt1, float rt2, float rt3)>();
            NoRidgeRegion = new List<List<float>>();

            //need a parameter
            if (EndRT - StartRT < CwtParameters.MinPeakWidth)
            {
                return;
            }

            //need to consider other spline types and time interval settings
            if (CubicSpline == null)
            {
                GetCubicSpline();
            }
            if (SmoothedData == null)
            {
                Interpolte_cubic(150);
            }
 
            var rtSeq = SmoothedData.Select(p => p.Item1).ToList();
            var peakArrayList = new float[rtSeq.Count * 2];
            for (int i = 0; i < SmoothedData.Count; i++)
            {
                peakArrayList[2 * i] = SmoothedData[i].Item1;
                peakArrayList[2 * i + 1] = SmoothedData[i].Item2;
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
                        ridge.Intensity = nearestRidge.intensity;
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

            //if PeakRidgeList is empty or only have one PeakRidge line, add the whole region as a peak region
            if (PeakRidgeList.Count() <= 1)
            {
                PeakRegionList.Add((rtSeq.First(), (float)ApexRT, rtSeq.Last()));
                var RidgeRTs = new List<float>();
                RidgeRTs.Add((float)ApexRT);
                NoRidgeRegion.Add(RidgeRTs);
            }

            //Added code
            //Remove the PeakRidge lines with intensity lower than NL * SNRThreshold
            CalculateNL();
            PeakRidgeList = PeakRidgeList.Where(r => r.Intensity > NL * CwtParameters.SNRThreshold).ToList();

            //if we have more than one PeakRidge line, we can have valley points now (local minimum)
            if (PeakRidgeList.Count() > 1)
            {
                var ValleyPoints = new (float rt, float intensity)[PeakRidgeList.Count() + 1];
                ValleyPoints[0] = (peakArrayList[0], peakArrayList[1]);
                PeakRidge currentridge = PeakRidgeList.First();
                var localmin = (-1f, float.MaxValue);
                int startidx = rtSeq.IndexOf(currentridge.RT);

                //loop over all PeakRidge lines, find the local minimum between two PeakRidge lines next to each other
                for (int j = 1; j < PeakRidgeList.Count(); j++)
                {
                    PeakRidge nextridge = PeakRidgeList[j];
                    for (int i = startidx; i < rtSeq.Count(); i++) //loop over all points in the spline
                    {
                        var point = (rtSeq[i], peakArrayList[2 * i + 1]);
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

        public PeakCurve[] SeparatePeakByRegion(float SN)
        {
            var newPeakCurves = new PeakCurve[PeakRegionList.Count()];
            int i = 0;
            foreach (var region in PeakRegionList)
            {
                var newPeakCurve = new PeakCurve();
                newPeakCurve.Peaks = new List<Peak>();
                newPeakCurve.MsLevel = MsLevel;
                newPeakCurve.SmoothedData = new List<(float, float)>();
                newPeakCurves[i] = newPeakCurve;
                i++;
            }

            //check if any region is too wide
            foreach (var region in PeakRegionList)
            {
                var newPeakCurve = new PeakCurve();

                if (region.rt3 - region.rt1 > CwtParameters.MaxCurveRTRange)
                {
                    var RTs = SmoothedData.Select(p => p.Item1).ToArray();
                    int leftidx = BinarySearchLower(RTs, region.rt1);
                    int rightidx = BinarySearchHigher(RTs, region.rt1);
                    var left = SmoothedData[leftidx];
                    var right = SmoothedData[rightidx];
                    while ((right.Item1 - left.Item1) > CwtParameters.MaxCurveRTRange)
                    {
                        if (right.Item1 - region.rt2 <= CwtParameters.MaxCurveRTRange / 4f)
                        {
                            leftidx++;
                        }
                        else if (region.rt2 - left.Item1 <= CwtParameters.MaxCurveRTRange / 4f)
                        {
                            rightidx--;
                        }
                        else if (left.Item2 < right.Item2)
                        {
                            leftidx++;
                        }
                        else
                        {
                            rightidx--;
                        }
                        left = SmoothedData[leftidx];
                        right = SmoothedData[rightidx];
                    }
                    var newRegion = (left.Item1, region.rt2, right.Item1);
                    PeakRegionList.Add(newRegion);
                    PeakRegionList.Remove(region);
                }
            }

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
            foreach(var point in SmoothedData)
            {
                for (int j = 0; j < PeakRegionList.Count; j++)
                {
                    if (point.Item1 >= PeakRegionList[j].rt1 && point.Item1 <= PeakRegionList[j].rt3)
                    {
                        newPeakCurves[j].SmoothedData.Add(point);
                        break;
                    }
                }
            }

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

        public void VisualizeRaw(string chartType)
        {
            if (chartType == "line")
            {
                var plot = Chart2D.Chart.Line<double, double, string>(
                        x: Peaks.Select(p => p.RetentionTime),
                        y: Peaks.Select(p => p.Intensity)).WithTraceInfo("raw").WithMarkerStyle(Color: Color.fromString("red"));
                plot.Show();
            }
            if (chartType == "point")
            {
                var plot = Chart2D.Chart.Point<double, double, string>(
                        x: Peaks.Select(p => p.RetentionTime),
                        y: Peaks.Select(p => p.Intensity)).WithTraceInfo("raw").WithMarkerStyle(Color: Color.fromString("red"));
                plot.Show();
            }
        }

        public GenericChart VisualizeCubicSpline()
        {
            if (CubicSpline == null)
            {
                GetCubicSpline();
            }
            if (SmoothedData == null)
            {
                Interpolte_cubic(150);
            }
            var plot = Chart2D.Chart.Point<float, float, string>(
                x: SmoothedData.Select(p => p.Item1),
                y: SmoothedData.Select(p => p.Item2)).WithTraceInfo("spline").WithMarkerStyle(Color: Color.fromString("blue"));
            return plot;
        }

        public void VisualizePeakRegions()
        {
            if (PeakRegionList == null)
            {
                DetectPeakRegions();
            }
            var plot_spline = VisualizeCubicSpline();
            var markedRTs = new List<float> { PeakRegionList[0].rt1 };
            float yMin = SmoothedData.Min(p => p.Item2);
            float yMax = SmoothedData.Max(p => p.Item2);
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
    }
 
    }
