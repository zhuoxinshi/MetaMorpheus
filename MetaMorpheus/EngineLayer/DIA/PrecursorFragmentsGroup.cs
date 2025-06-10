using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Plotly.NET;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.DIA
{
    public class PrecursorFragmentsGroup
    {
        public PrecursorFragmentsGroup(PeakCurve precursorPeakCurve)
        {
            PrecursorPeakCurve = precursorPeakCurve;
            PFpairs = new List<PrecursorFragmentPair>();
        }

        public PrecursorFragmentsGroup(PeakCurve precursorPeakCurve, List<PeakCurve> fragments)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach(var fragment in fragments)
            {
                var pair = new PrecursorFragmentPair(precursorPeakCurve, fragment);
                pfPairs.Add(pair);
            }
            PrecursorPeakCurve = precursorPeakCurve;
            PFpairs = pfPairs;
        }

        public PrecursorFragmentsGroup(PeakCurve precursorPeakCurve, List<PeakCurve> fragments, int index)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragment in fragments)
            {
                var pair = new PrecursorFragmentPair(precursorPeakCurve, fragment);
                pfPairs.Add(pair);
            }
            PrecursorPeakCurve = precursorPeakCurve;
            PFpairs = pfPairs;
            PFgroupIndex = index;
        }

        public PrecursorFragmentsGroup(PeakCurve precursorPeakCurve, List<PrecursorFragmentPair> pfPairs, int index)
        {
            PrecursorPeakCurve = precursorPeakCurve;
            PFpairs = pfPairs;
            PFgroupIndex = index;
        }

        public PeakCurve PrecursorPeakCurve { get; set; }
        public List<PrecursorFragmentPair> PFpairs { get; set; }
        public int PFgroupIndex;
        public int NumHighCorrFragments { get; set; }

        //TODO: finish this
        public static PrecursorFragmentsGroup GroupPF(PeakCurve prePeakCurve, List<PeakCurve> fragPeakCurves)
        {
            return new PrecursorFragmentsGroup(prePeakCurve);
        }
        
        public void GetNumberOfHighCorrFragments(DIAparameters diaParam)
        {
            NumHighCorrFragments = PFpairs.Where(p => p.Correlation >= diaParam.HighCorrThreshold).Count();
        }

        public static PrecursorFragmentsGroup GroupPrecursorFragments(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters DIAparameters)
        {
            var preFragGroup = new PrecursorFragmentsGroup(precursor);
            foreach (var ms2curve in ms2curves)
            {
                //if (ms2curve.ApexIntensity > precursor.ApexIntensity)
                //{
                //    continue;
                //}
                if (Math.Abs(ms2curve.ApexRT - precursor.ApexRT) <= DIAparameters.ApexRtTolerance)
                {
                    var overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursor, ms2curve);
                    if (overlap > DIAparameters.OverlapRatioCutOff)
                    {
                        double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursor, ms2curve);
                        if (corr > DIAparameters.CorrelationCutOff)
                        {
                            var PFpair = new PrecursorFragmentPair(precursor, ms2curve, overlap, corr);
                            lock (ms2curve.PFpairs)
                            {
                                ms2curve.PFpairs.Add(PFpair);
                            }
                            preFragGroup.PFpairs.Add(PFpair);
                        }
                    }
                }

            }
            //if (preFragGroup.PFpairs.Count > DIAparameters.FragmentRankCutOff)
            //{
            //    var filtered = preFragGroup.PFpairs.OrderByDescending(pair => pair.Correlation).Take(DIAparameters.FragmentRankCutOff);
            //    preFragGroup.PFpairs = filtered.ToList();
            //}
            if (preFragGroup.PFpairs.Count > 0)
            {
                preFragGroup.PFpairs = preFragGroup.PFpairs.OrderBy(pair => pair.FragmentPeakCurve.AveragedMz).ToList();
                return preFragGroup;
            }
            else
            {
                return null;
            }
        }

        public static PrecursorFragmentsGroup UmpireGrouping(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters DIAparameters)
        {
            //if (Math.Abs(precursor.MonoisotopicMass - 9461) < 1 && precursor.Charge == 12)
            //{
            //    //debug
            //    int stop = 0;
            //}
            var preFragGroup = new PrecursorFragmentsGroup(precursor);
            foreach (var ms2curve in ms2curves)
            {
                if (Math.Abs(ms2curve.ApexRT - precursor.ApexRT) <= DIAparameters.ApexRtTolerance)
                {
                    var overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursor, ms2curve);
                    if (overlap > DIAparameters.OverlapRatioCutOff)
                    {
                        double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData_Umpire(precursor, ms2curve, DIAparameters.NoPointsPerMin);
                        if (corr > DIAparameters.CorrelationCutOff)
                        {
                            var PFpair = new PrecursorFragmentPair(precursor, ms2curve, overlap, corr);
                            lock (ms2curve.PFpairs)
                            {
                                ms2curve.PFpairs.Add(PFpair);
                            }
                            preFragGroup.PFpairs.Add(PFpair);
                        }
                    }
                }
            }
            //if (preFragGroup.PFpairs.Count > DIAparameters.FragmentRankCutOff)
            //{
            //    var filtered = preFragGroup.PFpairs.OrderByDescending(pair => pair.Correlation).Take(DIAparameters.FragmentRankCutOff);
            //    preFragGroup.PFpairs = filtered.ToList();
            //}

            //debug
            var allFragMzs = preFragGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.Peaks.First().HighestPeakMz).OrderBy(p => p).ToList();
            if (preFragGroup.PFpairs.Count > 0)
            {
                preFragGroup.PFpairs = preFragGroup.PFpairs.OrderBy(pair => pair.FragmentPeakCurve.AveragedMz).ToList();
                return preFragGroup;
            }
            else
            {
                return null;
            }
        }

        public static PrecursorFragmentsGroup GroupPrecursorFragments_overlapFirst(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters DIAparameters)
        {
            var preFragGroup = new PrecursorFragmentsGroup(precursor);

            //Get all ms2 XICs in range
            var ms2curvesInRange = ms2curves.Where(p => (p.StartCycle >= precursor.StartCycle && p.StartCycle <= precursor.EndCycle)
            || (p.EndCycle >= precursor.StartCycle && p.EndCycle <= precursor.EndCycle)).ToList();

            foreach (var ms2curve in ms2curvesInRange)
            {
                var ms2peaks = ms2curve.Peaks.Where(p => p.ZeroBasedScanIndex >= precursor.StartCycle && p.ZeroBasedScanIndex <= precursor.EndCycle).ToList();
                if (ms2peaks.Count < 5)
                {
                    continue;
                }
                var newMs2curve = new PeakCurve(ms2peaks);
                if (Math.Abs(newMs2curve.ApexRT - precursor.ApexRT) <= DIAparameters.ApexRtTolerance)
                {
                    var overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursor, newMs2curve);
                    if (overlap > DIAparameters.OverlapRatioCutOff)
                    {
                        double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursor, ms2curve);
                        if (corr > DIAparameters.CorrelationCutOff)
                        {
                            var PFpair = new PrecursorFragmentPair(precursor, ms2curve, corr);
                            lock (ms2curve.PFpairs)
                            {
                                ms2curve.PFpairs.Add(PFpair);
                            }
                            preFragGroup.PFpairs.Add(PFpair);
                        }
                    }
                }
            }
            if (preFragGroup.PFpairs.Count > 0)
            {
                preFragGroup.PFpairs = preFragGroup.PFpairs.OrderBy(pair => pair.FragmentPeakCurve.AveragedMz).ToList();
                return preFragGroup;
            }
            else
            {
                return null;
            }
        }


        public static PrecursorFragmentsGroup SharedXICGrouping(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters diaParam)
        {
            var preFragGroup = new PrecursorFragmentsGroup(precursor);
            foreach (var ms2curve in ms2curves)
            {
                //if (ms2curve.ApexIntensity > precursor.ApexIntensity)
                //{
                //    continue;
                //}
                if (Math.Abs(ms2curve.ApexRT - precursor.ApexRT) <= diaParam.ApexRtTolerance)
                {
                    var sharedXIC = PrecursorFragmentPair.CalculateSharedXIC(precursor, ms2curve);
                    if (sharedXIC >= diaParam.SharedXICCutOff)
                    {
                        var PFpair = new PrecursorFragmentPair(precursor, ms2curve, sharedXIC);
                        lock (ms2curve.PFpairs)
                        {
                            ms2curve.PFpairs.Add(PFpair);
                        }
                        preFragGroup.PFpairs.Add(PFpair);
                    }
                }

            }
            //if (preFragGroup.PFpairs.Count > DIAparameters.FragmentRankCutOff)
            //{
            //    var filtered = preFragGroup.PFpairs.OrderByDescending(pair => pair.Correlation).Take(DIAparameters.FragmentRankCutOff);
            //    preFragGroup.PFpairs = filtered.ToList();
            //}
            if (preFragGroup.PFpairs.Count > 0)
            {
                preFragGroup.PFpairs = preFragGroup.PFpairs.OrderBy(pair => pair.FragmentPeakCurve.AveragedMz).ToList();
                return preFragGroup;
            }
            else
            {
                return null;
            }
        }
        
        public GenericChart VisualizeLog10()
        {
            var normalizedIntensity = PrecursorPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
            var precursorPlot = Chart2D.Chart.Line<double, double, string>(
                    x: PrecursorPeakCurve.Peaks.Select(p => p.RetentionTime),
                    y: normalizedIntensity).WithTraceInfo($"precursor_{Math.Round(PrecursorPeakCurve.AveragedMz, 3)}").WithMarkerStyle(Color: Color.fromString("red"));
            var plots = new List<GenericChart> { precursorPlot };
            foreach (var pf in PFpairs)
            {
                var norm2 = pf.FragmentPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                        x: pf.FragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                        y: norm2).WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{Math.Round(pf.Correlation, 2)}_{pf.FragmentRank}")
                        .WithMarkerStyle(Color: Color.fromString("blue"));
                plots.Add(fragmentPlot2);
            }
            var combinedPlot = Chart.Combine(plots);
            return combinedPlot;
        }

        public GenericChart VisualizeNormalized()
        {
            var maxPreIntensity = PrecursorPeakCurve.Peaks.Max(p => p.Intensity);
            var normalizedIntensity = PrecursorPeakCurve.Peaks.Select(p => p.Intensity / maxPreIntensity);
            var precursorPlot = Chart2D.Chart.Line<double, double, string>(
                    x: PrecursorPeakCurve.Peaks.Select(p => p.RetentionTime),
                    y: normalizedIntensity, FillColor: Color.fromString("white")).WithTraceInfo($"precursor_{Math.Round(PrecursorPeakCurve.AveragedMz, 3)}").WithMarkerStyle(Color: Color.fromString("red"));
            var plots = new List<GenericChart> { };
            foreach (var pf in PFpairs)
            {
                var maxFragIntensity = pf.FragmentPeakCurve.Peaks.Max(p => p.Intensity);
                var norm2 = pf.FragmentPeakCurve.Peaks.Select(p => p.Intensity / maxFragIntensity);
                var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                        x: pf.FragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                        y: norm2).WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{Math.Round(pf.Correlation, 2)}_{pf.FragmentRank}")
                        .WithMarkerStyle(Color: Color.fromString("blue"));
                plots.Add(fragmentPlot2);
            }
            plots.Add(precursorPlot);
            var combinedPlot = Chart.Combine(plots);
            return combinedPlot;
        }

        public GenericChart VisualizeXYData(int numFrags = 0)
        {
            if (PrecursorPeakCurve.XYData == null)
            {
                return null;
            }
            var plots = new List<GenericChart> { };

            var maxPreIntensity = PrecursorPeakCurve.XYData.Max(p => p.Item2);
            var normalizedIntensity = PrecursorPeakCurve.XYData.Select(p => p.Item2 / maxPreIntensity);
            var precursorPlot = Chart2D.Chart.Line<double, double, string>(
                       x: PrecursorPeakCurve.XYData.Select(p => p.Item1),
                       y: normalizedIntensity).WithLayout(Layout.init<string>(PlotBGColor: Color.fromString("white")))
                       .WithTraceInfo($"precursor_{Math.Round(PrecursorPeakCurve.AveragedMz, 3)}")
                       .WithMarkerStyle(Color: Color.fromString("red")).WithLineStyle(Width: 8);

            var pairs = PFpairs;
            if (numFrags > 0 && PFpairs.Count > numFrags)
            {
                pairs = PFpairs.OrderByDescending(p => p.Correlation).Take(numFrags).ToList();
            }
            var maxFragIntensity = pairs.Max(pf => pf.FragmentPeakCurve.XYData.Max(p => p.Item2));
            foreach (var pf in pairs)
            {
                var norm2 = pf.FragmentPeakCurve.XYData.Select(p => p.Item2 / maxFragIntensity);
                var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                        x: pf.FragmentPeakCurve.XYData.Select(p => p.Item1),
                        y: norm2).WithLayout(Layout.init<string>(PlotBGColor: Color.fromString("white")))
                        .WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{Math.Round(pf.Correlation, 2)}_{pf.FragmentRank}")
                        .WithMarkerStyle(Color: Color.fromString("RoyalBlue")).WithLineStyle(Width: 4);
                plots.Add(fragmentPlot2);
            }
            plots.Add(precursorPlot);
            var combinedPlot = Chart.Combine(plots);
            return combinedPlot;
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
                var ms2curves = diaEngine2.Ms2PeakCurves.Where(d => precursorPeakCurve.AveragedMz + 1 > d.Key.min && precursorPeakCurve.AveragedMz - 1 < d.Key.max)
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

        public PeakCurve SearchFragment(double mz, Tolerance tolerance)
        {
            var mzs = PFpairs.Select(pair => pair.FragmentPeakCurve.AveragedMz).ToArray();
            int index = Array.BinarySearch(mzs, mz);
            if (index >= 0)
            {
                return PFpairs[index].FragmentPeakCurve;
            }
            index = ~index;
            if (index == 0)
            {
                if (tolerance.Within(mzs[index], mz))
                {
                    return PFpairs[index].FragmentPeakCurve;
                }
                else
                {
                    return null;
                }
            }
            else if (index == mzs.Length - 1)
            {
                if (tolerance.Within(mzs[index], mz))
                {
                    return PFpairs[index].FragmentPeakCurve;
                }
                else
                {
                    return null;
                }
            }
            else 
            {
                if (Math.Abs(mzs[index] - mz) < Math.Abs(mzs[index + 1] - mz))
                {
                    if (tolerance.Within(mzs[index], mz))
                    {
                        return PFpairs[index].FragmentPeakCurve;
                    }
                    else
                    {
                        return null;
                    }
                }
                else
                {
                    if (tolerance.Within(mzs[index + 1], mz))
                    {
                        return PFpairs[index + 1].FragmentPeakCurve;
                    }
                    else
                    {
                        return null;
                    }
                }
            }
        }

        public static PrecursorFragmentsGroup FindPfGroupBasedOnPsm(PsmFromTsv psmTsv, List<DeconvolutedMass>[] ms1PeakTable, Dictionary<double, List<DeconvolutedMass>[]> allMs2MassTables,
            DIAparameters diaParam)
        {
            var cycle = (psmTsv.Ms2ScanNumber - 1) % 4;
            List<DeconvolutedMass>[] ms2PeakTable = null;
            switch (cycle)
            {
                case 1:
                    ms2PeakTable = allMs2MassTables[60];
                    break;
                case 2:
                    ms2PeakTable = allMs2MassTables[80];
                    break;
                case 3:
                    ms2PeakTable = allMs2MassTables[100];
                    break;
            }
            var precursorMass = MassCurve.GetMassFromScan(psmTsv.PrecursorMass, psmTsv.PrecursorCharge, ms1PeakTable, (psmTsv.PrecursorScanNum - 1) / 4,
                diaParam.Ms1PeakFindingTolerance, diaParam.PeakSearchBinSize);
            if (precursorMass == null || precursorMass.PeakCurve == null)
            {
                return null;
            }
            var pfGroup = new PrecursorFragmentsGroup(precursorMass.PeakCurve);
            precursorMass.PeakCurve.GetRawXYData();

            foreach (var matchedIon in psmTsv.MatchedIons)
            {
                var fragmentMass = MassCurve.GetMassFromScan(matchedIon.Mz.ToMass(matchedIon.Charge), matchedIon.Charge, ms2PeakTable, (psmTsv.Ms2ScanNumber - 1) / 4,
                                       diaParam.Ms2PeakFindingTolerance, diaParam.PeakSearchBinSize);
                if ( fragmentMass.PeakCurve != null)//fragmentMass == null ||
                {
                    double overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    fragmentMass.PeakCurve.GetRawXYData();
                    double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    var pfPair = new PrecursorFragmentPair(precursorMass.PeakCurve, fragmentMass.PeakCurve, overlap, corr, psmTsv.DecoyContamTarget);
                    pfGroup.PFpairs.Add(pfPair);
                }
            }
            return pfGroup;
        }

    }
}
