using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.Modifications;
using Plotly.NET.CSharp;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using TopDownProteomics;

namespace EngineLayer.DIA
{
    public class ISDEngine_static
    {
        public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();

            //calculate number of scans per cycle
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            var ms2Scans = dataFile.Where(s => s.MsnOrder == 2).ToArray();
            int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;
            diaParam.NumScansPerCycle = scansPerCycle;

            //Get ms1 XICs
            var allMs1PeakCurves = GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance, diaParam.MaxRTRangeMS1,
                out List<Peak>[] peaksByScan, diaParam.CutMs1Peaks, null, diaParam.MinMS1Mass, diaParam.MinMS1Charge, diaParam.Ms1NumPeaksThreshold);
            PeakCurveSpline(allMs1PeakCurves.Where(p => p.Peaks.Count > 4).ToList(), diaParam.Ms1SplineType, diaParam, ms1Scans, ms2Scans);
            PeakCurveSpline(allMs1PeakCurves.Where(p => p.Peaks.Count <= 4).ToList(), diaParam.Ms1SplineType, diaParam, ms1Scans, ms2Scans);

            //Get ms2 XICs
            var isdScanVoltageMap = ConstructMs2Groups(ms2Scans);
            var allMs2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach(var ms2Group in isdScanVoltageMap)
            {
                allMs2PeakCurves[ms2Group.Key] = GetAllPeakCurves(ms2Group.Value.ToArray(), commonParameters, diaParam, diaParam.Ms2XICType,
                    diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out List<Peak>[] peaksByScan2, diaParam.CutMs2Peaks, null, diaParam.MinMS2Mass, diaParam.MinMS2Charge, diaParam.Ms2NumPeaksThreshold);
                PeakCurveSpline(allMs2PeakCurves[ms2Group.Key].Where(p => p.Peaks.Count > 4).ToList(), diaParam.Ms2SplineType, diaParam, ms1Scans, ms2Scans);
                PeakCurveSpline(allMs2PeakCurves[ms2Group.Key].Where(p => p.Peaks.Count <= 4).ToList(), diaParam.Ms2SplineType, diaParam, ms1Scans, ms2Scans);
            }

            //Group precursors
            //var highestIntensityPrecursors = allMs1PeakCurves
            //    .GroupBy(p => (Math.Round(p.MonoisotopicMass, 1), p.ApexRT))
            //    .Select(g => g.OrderByDescending(p => p.ApexIntensity).First())
            //    .ToList();

            //precursor fragment grouping
            var pfGroups = new List<PrecursorFragmentsGroup>();
            var allFragments = allMs2PeakCurves.Values.SelectMany(p => p).ToList();
            Parallel.ForEach(Partitioner.Create(0, allMs1PeakCurves.Count), new ParallelOptions { MaxDegreeOfParallelism = 15 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = allMs1PeakCurves[i];
                        if (precursor.ApexSN < diaParam.PrecursorSNCutOff)
                        {
                            continue;
                        }

                        if (diaParam.CombineFragments == false)
                        {
                            foreach (var fragments in allMs2PeakCurves.Values)
                            {
                                var pfGroup = PFgrouping(precursor, fragments, diaParam);
                                if (pfGroup != null && pfGroup.PFpairs.Count > diaParam.MinPFpairCount)
                                {
                                    lock (pfGroups)
                                        pfGroups.Add(pfGroup);
                                    //pfGroup.VisualizeXYData().Show();
                                }
                            }
                        } else
                        {
                            var pfGroup = PFgrouping(precursor, allFragments, diaParam);
                            if (pfGroup != null && pfGroup.PFpairs.Count > diaParam.MinPFpairCount)
                            {
                                lock (pfGroups)
                                    pfGroups.Add(pfGroup);
                            }
                        }
                    }
                });

            //precursor fragment rank filtering
            if (diaParam.RankFilter)
            {
                foreach (var ms2curve in allMs2PeakCurves.Values.SelectMany(p => p))
                {
                    ms2curve.GetPrecursorRanks();
                }
                foreach (var group in pfGroups)
                {
                    group.PFpairs = group.PFpairs.OrderByDescending(pf => pf.Correlation).Take(diaParam.FragmentRankCutOff).ToList();
                    group.PFpairs = group.PFpairs.Where(pf => pf.PrecursorRank <= diaParam.PrecursorRankCutOff).ToList();
                    group.GetNumberOfHighCorrFragments(diaParam);
                }
                pfGroups = pfGroups.Where(pf => pf.PFpairs.Count > 0 && pf.NumHighCorrFragments >= diaParam.NumHighCorrFragments).ToList();
            }

            //pfGroups grouping
            //var groupedPFgroups = pfGroups.GroupBy(pf => Math.Round(pf.PrecursorPeakCurve.MonoisotopicMass, 0)).ToList();
            //var filteredGroups = new List<PrecursorFragmentsGroup>();
            //foreach (var group in groupedPFgroups)
            //{
            //    var highestGroup = group.OrderByDescending(g => g.PrecursorPeakCurve.AveragedIntensity).FirstOrDefault();
            //    filteredGroups.Add(highestGroup);
            //}
            if (diaParam.PFgroupsDictionary == null)
            {
                diaParam.PFgroupsDictionary = new Dictionary<string, List<PrecursorFragmentsGroup>>();
            }
            diaParam.PFgroupsDictionary[dataFile.FilePath] = pfGroups;
            if (diaParam.PeakCurveDictionary == null)
            {
                diaParam.PeakCurveDictionary = new Dictionary<string, List<PeakCurve>>();
            }
            diaParam.PeakCurveDictionary[dataFile.FilePath] = allMs1PeakCurves.Concat(allMs2PeakCurves.SelectMany(p => p.Value)).ToList();

            //construct new ms2Scans
            int pfGroupIndex = 1;
            foreach (var pfGroup in pfGroups)
            {
                pfGroup.PFgroupIndex = pfGroupIndex;
                pfGroupIndex++;
                var newScans = ConstructNewMs2Scans(pfGroup, commonParameters, diaParam.PseudoMs2ConstructionType, dataFile.FilePath);
                pseudoMs2Scans.Add(newScans);
            }

            //debug
            var sortedGroups = pfGroups.OrderByDescending(p => p.PrecursorPeakCurve.MonoisotopicMass).ToList();

            return pseudoMs2Scans;
        }

        public static Dictionary<double, List<MsDataScan>> ConstructMs2Groups(MsDataScan[] ms2Scans)
        {
            var isdScanVoltageMap = new Dictionary<double, List<MsDataScan>>();
            string pattern = $@"sid=(\d+)";
            foreach (var ms2 in ms2Scans)
            {
                var match = Regex.Match(ms2.ScanFilter, pattern);
                double voltage = double.Parse(match.Groups[1].Value);
                if (!isdScanVoltageMap.ContainsKey(voltage))
                {
                    isdScanVoltageMap[voltage] = new List<MsDataScan>();
                    isdScanVoltageMap[voltage].Add(ms2);
                }
                else
                {
                    isdScanVoltageMap[voltage].Add(ms2);
                }
            }
            return isdScanVoltageMap;
        }

        public static List<PeakCurve> GetAllPeakCurves(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam, XICType xicType,
            Tolerance peakFindingTolerance, double maxRTRange, out List<Peak>[] allPeaksByScan, bool cutPeak = false, MzRange isolationWindow = null, double minMass = 0, int minCharge = 1
            , int numPeakThreshold = 2)
        {
            var peakCurves = new List<PeakCurve>();
            switch (xicType)
            {
                case XICType.Peak:
                    peakCurves = GetAllPeakCurves_Peak(scans, diaParam, peakFindingTolerance, maxRTRange, out allPeaksByScan, cutPeak);
                    break;

                case XICType.DeconHighestPeak:
                    peakCurves = GetAllPeakCurves_DeconHighestPeak(scans, commonParameters, diaParam, peakFindingTolerance, maxRTRange, out allPeaksByScan, cutPeak, isolationWindow,
                        minCharge, minMass);
                    break;

                //case XICType.isoEnvelopeTotal:
                //    peakCurves = GetAllPeakCurves_isoEnvelopeTotal(scans, commonParameters, diaParam, peakFindingTolerance, maxRTRange, out allPeaksByScan);
                //    break;

                //case XICType.Peak_cutPeak:
                //    peakCurves = GetAllPeakCurves_cutPeak(scans, diaParam, peakFindingTolerance, maxRTRange, out allPeaksByScan);
                //    break;
                case XICType.MassCurve:
                    peakCurves = MassCurve.GetAllMassCurves(scans, commonParameters, diaParam, peakFindingTolerance, maxRTRange, minMass, minCharge, out allPeaksByScan, cutPeak, isolationWindow);
                    break;
                default: throw new MzLibException("XICType");
            }
            peakCurves = peakCurves.Where(pc => pc.Peaks.Count > numPeakThreshold).ToList();
            int index = 1;
            foreach (var pc in peakCurves)
            {
                pc.Index = index;
                index++;
            }
            return peakCurves;
        }

        public static void PeakCurveSpline(List<PeakCurve> allPeakCurves, SplineType splineType, DIAparameters diaParam, MsDataScan[] ms1Scans, MsDataScan[] ms2Scans)
        {
            var rtIndexMap = GetRtIndexMap(ms1Scans);
            var rtMap = GetRtMap(ms1Scans, ms2Scans);

            switch (splineType)
            {
                case SplineType.NoSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetRawXYData();
                    break;
                case SplineType.CubicSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetCubicSplineXYData(diaParam.SplineRtInterval);
                    break;
                case SplineType.BSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetBSplineXYData(diaParam.SplineRtInterval, 2);
                    break;
                case SplineType.UmpireBSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetUmpireBSplineData(diaParam.NoPointsPerMin, 2);
                    break;
                case SplineType.Ms1SpaceBSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetMs1SpaceBSplineXYData(diaParam.SplineRtInterval, 2, rtMap);
                    break;
                case SplineType.ScanCycleCubicSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetScanCycleCubicSplineXYData(diaParam.ScanCycleSplineTimeInterval);
                    break;
                case SplineType.SavgolSmoothed:
                    foreach(var pc in allPeakCurves)
                        pc.GetSavgolSmoothedXYData(diaParam.SGfilterWindowSize);
                    break;
                case SplineType.CubicSplineSavgolSmoothed:
                    foreach (var pc in allPeakCurves)
                        pc.GetCubicSplineSavgolSmoothedXYData(diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                    break;
                case SplineType.ScanCycleCubicSplineSavgolSmoothed:
                    foreach (var pc in allPeakCurves)
                        pc.GetScanCycleCubicSplineSavgolSmoothedXYData(diaParam.SGfilterWindowSize, diaParam.ScanCycleSplineTimeInterval);
                    break;
                case SplineType.SavgolSmoothedCubicSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetSavgolSmoothedCubicSplineXYData(diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                    break;
                case SplineType.Ms1SpaceCubicSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetMs1SpaceCubicSplineXYData(rtMap, diaParam.SplineRtInterval);
                    break;
                case SplineType.Ms1SpaceCubicSplineSavgolSmoothed:
                    foreach(var pc in allPeakCurves)
                        pc.GetMs1SpaceCubicSplineSavgolSmoothedXYData(rtMap, diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                    break;
                case SplineType.Ms1SpaceSavgolSmoothedCubicSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetMs1SpaceSavgolSmoothedCubicSplineXYData(rtMap, diaParam.SGfilterWindowSize, diaParam.SplineRtInterval);
                    break;
                case SplineType.GaussianFit:
                    //    Parallel.ForEach(Partitioner.Create(0, allPeakCurves.Count), new ParallelOptions { MaxDegreeOfParallelism = 15 },
                    //(partitionRange, loopState) =>
                    //{
                    //    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    //    {
                    //        allPeakCurves[i].GetGaussianFitXYData();
                    //    }
                    //});
                    foreach (var pc in allPeakCurves)
                        pc.GetGaussianFitXYData();
                    break;
                case SplineType.SimpleGaussian:
                    foreach (var pc in allPeakCurves)
                        pc.GetSimpleGaussianXYData();
                    break;
                case SplineType.SimpleGaussianSpline:
                    foreach(var pc in allPeakCurves)
                        pc.GetSimpleGaussianSplineXYData(diaParam.SplineRtInterval);
                    break;
                case SplineType.NormalizedLinearSpline:
                    foreach (var pc in allPeakCurves)
                        pc.GetNormalizedLinearSplinePeaks();
                    break;
            }
        }

        public static List<PeakCurve> GetAllPeakCurves_Peak(MsDataScan[] scans, DIAparameters diaParam, Tolerance peakFindingTolerance, double maxRTRange
            , out List<Peak>[] allPeaksByScan, bool cutPeak = false)
        {
            var allPeakCurves = new List<PeakCurve>();
            if (diaParam.NumScansPerCycle != 0)
            {
                allPeaksByScan = Peak.GetAllPeaksByScan(scans, diaParam.NumScansPerCycle);
            }
            else
            {
                allPeaksByScan = Peak.GetAllPeaksByScan(scans);
            }

            var allPeaks = allPeaksByScan.Where(v => v != null).SelectMany(p => p).ToList();
            var rankedPeaks = allPeaks.OrderByDescending(p => p.Intensity).ToList();
            var peakTable = Peak.GetPeakTable(allPeaks, diaParam.PeakSearchBinSize);
            foreach (var peak in rankedPeaks)
            {
                if (peak.PeakCurve != null)
                    continue;

                var newPeakCurve = PeakCurve.FindPeakCurve(peak, peakTable, scans, scans[0].IsolationRange,
                    diaParam.MaxNumMissedScan, peakFindingTolerance, diaParam.PeakSearchBinSize, maxRTRange);
                if (cutPeak)
                {
                    newPeakCurve.CutPeak();
                }
                allPeakCurves.Add(newPeakCurve);
            }
            return allPeakCurves;
        }

        public static List<PeakCurve> GetAllPeakCurves_DeconHighestPeak(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam, 
            Tolerance peakFindingTolerance, double maxRTRange, out List<Peak>[] allPeaksByScan, bool cutPeak = false, MzRange isolationWindow = null, int minCharge = 1,
                        double minMass = 0)
        {
            var allPeakCurves = new List<PeakCurve>();
            if (diaParam.NumScansPerCycle != 0)
            {
                allPeaksByScan = Peak.GetAllPeaksByScan(scans, diaParam.NumScansPerCycle);
            }
            else
            {
                allPeaksByScan = Peak.GetAllPeaksByScan(scans);
            }
            var allPeaks = allPeaksByScan.Where(v => v != null).SelectMany(p => p).ToList();
            var peakTable = Peak.GetPeakTable(allPeaks, diaParam.PeakSearchBinSize);
            int index = 1;
            var allPrecursors = new List<DeconvolutedMass>();
            for (int i = 0; i < scans.Length; i++)
            {
                var envelopes = Deconvoluter.Deconvolute(scans[i], commonParameters.PrecursorDeconvolutionParameters, isolationWindow).OrderByDescending(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.Charge < minCharge || envelope.MonoisotopicMass < minMass || envelope.MonoisotopicMass > diaParam.MaxMass)
                    {
                        continue;
                    }
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    int zeroBasedScanIndex = i;
                    if (diaParam.NumScansPerCycle != 0)
                    {
                        zeroBasedScanIndex = (scans[i].OneBasedScanNumber - 1) / diaParam.NumScansPerCycle;
                    }
                    var precursor = new DeconvolutedMass(envelope, scans[i].RetentionTime, 1, scans[i].OneBasedScanNumber, zeroBasedScanIndex);
                    allPrecursors.Add(precursor);
                }
            }
            allPrecursors = allPrecursors.OrderByDescending(p => p.HighestPeakIntensity).ToList();
            //var precursorGroups = allPrecursors.GroupBy(p => new { mass = Math.Round(p.MonoisotopicMass, 0), p.Charge}).ToList();
            foreach (var precursor in allPrecursors)
            {
                //var precursor = group.OrderByDescending(p => p.HighestPeakIntensity).FirstOrDefault();
                //if (precursor.HighestPeakIntensity < DIAparameters.PrecursorIntensityCutOff)
                //{
                //    continue;
                //}
                var peak = PeakCurve.GetPeakFromScan(precursor.HighestPeakMz, peakTable, precursor.ZeroBasedScanIndex, new PpmTolerance(0),
                    diaParam.PeakSearchBinSize);
                if (peak.PeakCurve == null && peak.Intensity >= diaParam.PrecursorSNCutOff)
                {
                    var newPeakCurve = PeakCurve.FindPeakCurve(peak, peakTable, scans, null, diaParam.MaxNumMissedScan,
                    peakFindingTolerance, diaParam.PeakSearchBinSize, maxRTRange);
                    newPeakCurve.MonoisotopicMass = precursor.MonoisotopicMass;
                    newPeakCurve.Charge = precursor.Charge;
                    newPeakCurve.Envelope = precursor.Envelope;
                    if (cutPeak)
                    {
                        newPeakCurve.CutPeak();
                    }
                    allPeakCurves.Add(newPeakCurve);
                }
            }
            return allPeakCurves;
        }


        //public static List<PeakCurve> GetAllPeakCurves_isoEnvelopeTotal(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam, 
        //    Tolerance peakFindingTolerance, double maxRTRange, out List<Peak>[] allPeaksByScan)
        //{
        //    var allPeakCurves = new List<PeakCurve>();
        //    if (diaParam.NumScansPerCycle != 0)
        //    {
        //        allPeaksByScan = Peak.GetAllPeaksByScan(scans, diaParam.NumScansPerCycle);
        //    }
        //    else
        //    {
        //        allPeaksByScan = Peak.GetAllPeaksByScan(scans);
        //    }
        //    var allPeaks = allPeaksByScan.Where(v => v != null).SelectMany(p => p).ToList();
        //    var peakTable = Peak.GetPeakTable(allPeaks, diaParam.PeakSearchBinSize);
        //    var allPrecursors = new List<DeconvolutedMass>[scans.Length];
        //    for (int i = 0; i < scans.Length; i++)
        //    {
        //        allPrecursors[i] = new List<DeconvolutedMass>();
        //        var envelopes = Deconvoluter.Deconvolute(scans[i], commonParameters.PrecursorDeconvolutionParameters).OrderByDescending(E => E.MonoisotopicMass);
        //        foreach (var envelope in envelopes)
        //        {
        //            if (envelope.MonoisotopicMass < diaParam.MinMass || envelope.Charge < 5)
        //            {
        //                continue;
        //            }
        //            var charge = envelope.Charge;
        //            double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
        //            double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
        //            var precursor = new DeconvolutedMass(envelope, charge, scans[i].RetentionTime, 1, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
        //                scans[i].OneBasedScanNumber, i);
        //            allPrecursors[i].Add(precursor);
        //        }
        //        allPrecursors[i] = allPrecursors[i].OrderBy(p => p.MonoisotopicMass).ToList();
        //    }
        //    var sortedAllPrecursors = allPrecursors.SelectMany(p => p).OrderByDescending(p => p.Envelope.Peaks.Count).ThenByDescending(p => p.HighestPeakIntensity).ToList();

        //    int index = 1;
        //    var allPeakEnvelopeCurves = new List<PeakEnvelopeCurve>();
        //    foreach (var precursor in sortedAllPrecursors)
        //    {
        //        var targetPeaks = precursor.Envelope.Peaks.OrderBy(p => p.mz).ToArray();
        //        var theorIntensityRatio = targetPeaks.Select(p => p.intensity / targetPeaks.Sum(p => p.intensity)).ToArray();
        //        var highestPeak = targetPeaks.OrderByDescending(p => p.intensity).FirstOrDefault();
        //        var highestPeakIndex = Array.IndexOf(targetPeaks, highestPeak);
        //        var newPEC = PeakEnvelopeCurve.GetPeakEnvelopeCurve(targetPeaks, theorIntensityRatio, highestPeakIndex, peakTable, scans, precursor.ZeroBasedScanIndex, diaParam,
        //            peakFindingTolerance, maxRTRange);
        //        if (newPEC != null)
        //        {
        //            allPeakEnvelopeCurves.Add(newPEC);
        //            newPEC.MonoisotopicMass = precursor.MonoisotopicMass;
        //            newPEC.Charge = precursor.Charge;
        //            newPEC.MsLevel = 1;
        //            newPEC.IsolationRange = null;
        //            newPEC.GetFakePeakCurve();
        //            allPeakCurves.Add(newPEC.FakePeakCurve);
        //        }
        //    }
        //    //debug
        //    //for (int i = 0; i < 100; i++)
        //    //{
        //    //    Random rnd = new Random();
        //    //    int r = rnd.Next(Ms1PeakCurves.Count - 1);
        //    //    var pc = Ms1PeakCurves[r];
        //    //    pc.VisualizeRaw("point");
        //    //}
        //    return allPeakCurves;
        //}

        //public static List<PeakCurve> GetAllPeakCurves_cutPeak(MsDataScan[] scans, DIAparameters diaParam, Tolerance peakFindingTolerance, double maxRTRange
        //    , out List<Peak>[] allPeaksByScan)
        //{
        //    var allPeakCurves = new List<PeakCurve>();
        //    if (diaParam.NumScansPerCycle != 0)
        //    {
        //        allPeaksByScan = Peak.GetAllPeaksByScan(scans, diaParam.NumScansPerCycle);
        //    }
        //    else
        //    {
        //        allPeaksByScan = Peak.GetAllPeaksByScan(scans);
        //    }
        //    var allPeaks = allPeaksByScan.Where(v => v != null).SelectMany(p => p).ToList();
        //    var rankedPeaks = allPeaks.OrderByDescending(p => p.Intensity).ToList();
        //    var peakTable = Peak.GetPeakTable(allPeaks, diaParam.PeakSearchBinSize);
        //    foreach (var peak in rankedPeaks)
        //    {
        //        if (peak.PeakCurve == null)
        //        {
        //            var newPeakCurve = PeakCurve.FindPeakCurve_cutPeak(peak, peakTable, scans, scans[0].IsolationRange,
        //                diaParam.MaxNumMissedScan, peakFindingTolerance, diaParam.PeakSearchBinSize, maxRTRange);
        //            allPeakCurves.Add(newPeakCurve);
        //        }
        //    }
        //    return allPeakCurves;
        //}

        public static PrecursorFragmentsGroup PFgrouping(PeakCurve precursor, List<PeakCurve> fragments, DIAparameters diaParam)
        {
            switch (diaParam.PFGroupingType)
            {
                case PFGroupingType.RetentionTime:
                    return PrecursorFragmentsGroup.GroupPrecursorFragments(precursor, fragments, diaParam);
                //case PFGroupingType.ScanCycle:
                //    return ISDEngine.GroupPrecursorFragments_scanCycle(precursor, fragments, diaParam);
                //case PFGroupingType.OverlapAreaRatio:
                //    return ISDEngine.GroupPrecursorFragments_area(precursor, fragments, diaParam);
                //case PFGroupingType.Area_correlation:
                //    return ISDEngine.GroupPrecursorFragments_area_correlation(precursor, fragments, diaParam);
                case PFGroupingType.OverlapFirst:
                    return PrecursorFragmentsGroup.GroupPrecursorFragments_overlapFirst(precursor, fragments, diaParam);
                case PFGroupingType.Umpire:
                    return PrecursorFragmentsGroup.UmpireGrouping(precursor, fragments, diaParam);
                case PFGroupingType.SharedXIC:
                    return PrecursorFragmentsGroup.SharedXICGrouping(precursor, fragments, diaParam);
                default: return null;
            }
        }


        public static Ms2ScanWithSpecificMass ConstructNewMs2Scans(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, PseudoMs2ConstructionType pseudoMs2Type, string dataFilePath)
        {
            switch(pseudoMs2Type)
            {
                case PseudoMs2ConstructionType.mzPeak:
                    return GetPseudoMs2Scan_mzPeak(pfGroup, commonParameters, dataFilePath);
                case PseudoMs2ConstructionType.neutralMass:
                    return GetPseudoMs2Scan_neutralMass(pfGroup, commonParameters, dataFilePath);
                case PseudoMs2ConstructionType.massCurve:
                    return GetPseudoMs2Scan_massCurve(pfGroup, commonParameters, dataFilePath);

                default: return null;
            }
        }

        public static Ms2ScanWithSpecificMass GetPseudoMs2Scan_mzPeak(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, string dataFilePath)
        {
            var mzs = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedMz).ToArray();
            var intensities = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedIntensity).ToArray();
            var spectrum = new MzSpectrum(mzs, intensities, false);
            var newMs2Scan = new MsDataScan(spectrum, pfGroup.PFgroupIndex, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                        MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: pfGroup.PrecursorPeakCurve.Index);
            var neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, commonParameters);
            var charge = pfGroup.PrecursorPeakCurve.Charge;
            var monoMz = pfGroup.PrecursorPeakCurve.MonoisotopicMass.ToMz(charge);
            //should highestPeakMz used in ms2withmass???
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge, dataFilePath, 
                commonParameters, neutralExperimentalFragments);

            return scanWithprecursor;
        }

        public static Ms2ScanWithSpecificMass GetPseudoMs2Scan_neutralMass(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, string dataFilePath)
        {
            var mzs = new double[] { 1 };
            var intensities = new double[] { Double.MaxValue };
            var spectrum = new MzSpectrum(mzs, intensities, false);
            var newMs2Scan = new MsDataScan(spectrum, pfGroup.PFgroupIndex, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                        MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
            var neutralExperimentalFragments = pfGroup.PFpairs.Select(pf => new IsotopicEnvelope(1,
                            new List<(double mz, double intensity)> { (1, 1) }, pf.FragmentPeakCurve.MonoisotopicMass, pf.FragmentPeakCurve.Charge, 1, 0)).OrderBy(e => e.MonoisotopicMass).ToArray();
            var monoMz = pfGroup.PrecursorPeakCurve.MonoisotopicMass.ToMz(pfGroup.PrecursorPeakCurve.Charge);
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, pfGroup.PrecursorPeakCurve.Charge, dataFilePath, commonParameters, neutralExperimentalFragments);

            return scanWithprecursor;
        }

        public static Ms2ScanWithSpecificMass GetPseudoMs2Scan_massCurve(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, string dataFilePath)
        {
            var mzs = new double[] { 1 };
            var intensities = new double[] { Double.MaxValue };
            var spectrum = new MzSpectrum(mzs, intensities, false);
            var newMs2Scan = new MsDataScan(spectrum, pfGroup.PFgroupIndex, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                        MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
            var neutralExperimentalFragments = pfGroup.PFpairs.Select(pf => new IsotopicEnvelope(
                            new List<(double mz, double intensity)> { (1, 1) }, pf.FragmentPeakCurve.AveragedMass, pf.FragmentPeakCurve.Charge, 1, 0)).OrderBy(e => e.MonoisotopicMass).ToArray();
            var charge = pfGroup.PrecursorPeakCurve.Charge;
            var monoMz = pfGroup.PrecursorPeakCurve.AveragedMass.ToMz(charge);
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge, dataFilePath, commonParameters, neutralExperimentalFragments);

            return scanWithprecursor;
        }

        public static Dictionary<string, List<double>> TestMatchedFragments(string sequence, MsDataFile dataFile, MsDataScan[] ms1Scans, MsDataScan[] ms2Scans, CommonParameters commonParameters,
            DIAparameters diaParam, double precursorMz, int precursorCharge, double precursorApexRT)
        {
            var peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>());
            var fragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, commonParameters.DigestionParams.FragmentationTerminus, fragments);

            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, commonParameters, commonParameters.DIAparameters,
               commonParameters.DIAparameters.Ms1XICType, commonParameters.DIAparameters.Ms1PeakFindingTolerance,
               commonParameters.DIAparameters.MaxRTRangeMS1, out List<Peak>[] peaksByScan);
            var allMs2PeakCurves = ISDEngine_static.GetAllPeakCurves(ms2Scans, commonParameters, commonParameters.DIAparameters,
                commonParameters.DIAparameters.Ms2XICType, commonParameters.DIAparameters.Ms2PeakFindingTolerance,
                commonParameters.DIAparameters.MaxRTRangeMS2, out List<Peak>[] peaksByScan2);

            var preXIC = allMs1PeakCurves.Where(pc => Math.Abs(pc.AveragedMz - precursorMz) < 0.01 && precursorCharge == 2
            && Math.Round(pc.ApexRT, 2) == precursorApexRT).First();

            var pfGroup = ISDEngine_static.PFgrouping(preXIC, allMs2PeakCurves, commonParameters.DIAparameters);
            var ms2WithPre = ISDEngine_static.ConstructNewMs2Scans(pfGroup, commonParameters, commonParameters.DIAparameters.PseudoMs2ConstructionType, dataFile.FilePath);
            var matchedMs2Curves = pfGroup.PFpairs.Select(p => p.FragmentPeakCurve).ToList();
            var deconvolutedMasses = ms2WithPre.ExperimentalFragments.Select(p => Math.Round(p.MonoisotopicMass, 1)).ToList();

            var match = new Dictionary<string, List<double>>();
            match["matched"]= new List<double>();
            match["missing"] = new List<double>();
            foreach (var fragment in fragments)
            {
                if (deconvolutedMasses.Contains(Math.Round(fragment.MonoisotopicMass, 1)))
                {
                    match["matched"].Add(fragment.MonoisotopicMass);
                }
                else
                {
                    match["missing"].Add(fragment.MonoisotopicMass);
                }
            }

            return match;
        }

        public static Dictionary<double, double> GetRtMap(MsDataScan[] ms1Scans, MsDataScan[] ms2Scans)
        {
            var rtMap = new Dictionary<double, double>();
            foreach (var scan in ms2Scans)
            {
                rtMap[scan.RetentionTime] = ms1Scans.Where(s => s.OneBasedScanNumber == scan.OneBasedPrecursorScanNumber).First().RetentionTime;
            }
            return rtMap;
        }

        public static Dictionary<int, double> GetRtIndexMap(MsDataScan[] scans)
        {
            var map = new Dictionary<int, double>();
            for (int i = 0; i < scans.Length; i++)
            {
                map[i] = scans[i].RetentionTime;
            }
            return map;
        }
    }
}
    
