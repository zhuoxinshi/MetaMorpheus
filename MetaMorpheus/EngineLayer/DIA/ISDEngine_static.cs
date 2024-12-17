using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Plotly.NET.TraceObjects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http.Headers;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ISDEngine_static
    {
        public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();

            //Get ms1 XICs
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            var allMs1PeakCurves = GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance, diaParam.MaxRTRangeMS1,
                out List<Peak>[] peaksByScan);

            //calculate number of scans per cycle
            int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;
            diaParam.NumScansPerCycle = scansPerCycle;

            //Get ms2 XICs
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var isdScanVoltageMap = ConstructMs2Groups(ms2Scans);
            var allMs2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach(var ms2Group in isdScanVoltageMap)
            {
                allMs2PeakCurves[ms2Group.Key] = GetAllPeakCurves(ms2Group.Value.ToArray(), commonParameters, diaParam, diaParam.Ms2XICType,
                    diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out List<Peak>[] peaksByScan2);
            }

            //precursor fragment grouping
            var pfGroups = new List<PrecursorFragmentsGroup>();
            foreach(var precursor in allMs1PeakCurves)
            {
                foreach (var fragments in allMs2PeakCurves.Values)
                {
                    var pfGroup = PFgrouping(precursor, fragments, diaParam);
                    if (pfGroup != null && pfGroup.PFpairs.Count > 0)
                    {
                        pfGroups.Add(pfGroup);
                    }
                }
            }

            //construct new ms2Scans
            foreach(var pfGroup in pfGroups)
            {
                var newScans = ConstructNewMs2Scans(pfGroup, commonParameters, diaParam, dataFile);
                pseudoMs2Scans.Add(newScans);
            }
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
            Tolerance peakFindingTolerance, double maxRTRange, out List<Peak>[] allPeaksByScan)
        {
            switch (xicType)
            {
                case XICType.Peak:
                    var peakCurves1 = GetAllPeakCurves_Peak(scans, diaParam, peakFindingTolerance, maxRTRange, out allPeaksByScan);
                    return peakCurves1;

                case XICType.DeconHighestPeak:
                    var peakCurves2 = GetAllPeakCurves_DeconHighestPeak(scans, commonParameters, diaParam, peakFindingTolerance, maxRTRange, out allPeaksByScan);
                    return peakCurves2;

                case XICType.isoEnvelopeTotal:
                    var peakCurves3 = GetAllPeakCurves_isoEnvelopeTotal(scans, commonParameters, diaParam, peakFindingTolerance, maxRTRange, out allPeaksByScan);
                    return peakCurves3;

                default: throw new MzLibException("XICType");
            }
        }

        public static List<PeakCurve> GetAllPeakCurves_Peak(MsDataScan[] scans, DIAparameters diaParam, Tolerance peakFindingTolerance, double maxRTRange
            , out List<Peak>[] allPeaksByScan)
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
                if (peak.PeakCurve == null)
                {
                    var newPeakCurve = PeakCurve.FindPeakCurve(peak, peakTable, scans, scans[0].IsolationRange,
                        diaParam.MaxNumMissedScan, peakFindingTolerance, diaParam.PeakSearchBinSize, maxRTRange);
                    if (diaParam.SplitMS2Peak)
                    {
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            newPeakCurve.DetectPeakRegions();
                            var newPCs = newPeakCurve.SeparatePeakByRegion();
                            foreach (var pc in newPCs)
                            {
                                if (pc.Peaks.Count > 4)
                                {
                                    allPeakCurves.Add(pc);
                                }
                            }
                        }
                    }
                    else
                    {
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            allPeakCurves.Add(newPeakCurve);
                            newPeakCurve.GetScanCycleSmoothedData(diaParam.ScanCycleSplineTimeInterval);
                        }
                    }
                }
            }
            return allPeakCurves;
        }

        public static List<PeakCurve> GetAllPeakCurves_DeconHighestPeak(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam, 
            Tolerance peakFindingTolerance, double maxRTRange, out List<Peak>[] allPeaksByScan)
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
                var envelopes = Deconvoluter.Deconvolute(scans[i], commonParameters.PrecursorDeconvolutionParameters).OrderByDescending(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.Charge < 6 || envelope.Peaks.Count < 3 || envelope.MonoisotopicMass < diaParam.MinMass || envelope.MonoisotopicMass > diaParam.MaxMass)
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
                    var precursor = new DeconvolutedMass(envelope, charge, scans[i].RetentionTime, 1, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        scans[i].OneBasedScanNumber, zeroBasedScanIndex);
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
                if (peak.PeakCurve == null && peak.Intensity >= diaParam.PrecursorIntensityCutOff)
                {
                    var newPeakCurve = PeakCurve.FindPeakCurve(peak, peakTable, scans, null, diaParam.MaxNumMissedScan,
                    peakFindingTolerance, diaParam.PeakSearchBinSize, maxRTRange);
                    newPeakCurve.MonoisotopicMass = precursor.MonoisotopicMass;
                    newPeakCurve.Charge = precursor.Charge;
                    newPeakCurve.Envelope = precursor.Envelope;
                    if (diaParam.SplitMS1Peak)
                    {
                        newPeakCurve.DetectPeakRegions();
                        var newPCs = newPeakCurve.SeparatePeakByRegion();
                        foreach (var pc in newPCs)
                        {
                            if (pc.Peaks.Count > 4)
                            {
                                allPeakCurves.Add(pc);
                            }
                        }
                    }
                    else
                    {
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            allPeakCurves.Add(newPeakCurve);
                            newPeakCurve.GetScanCycleSmoothedData(diaParam.ScanCycleSplineTimeInterval);
                            newPeakCurve.Index = index;
                            index++;

                            //debug
                            //if (Math.Abs(newPeakCurve.AveragedMz - 864.60) < 0.01 && newPeakCurve.Charge == 16)
                            //{
                            //    newPeakCurve.VisualizeRaw("line");
                            //}
                        }
                    }
                }
            }
            //debug
            //for (int i = 0; i < 100; i++)
            //{
            //    Random rnd = new Random();
            //    int r = rnd.Next(Ms1PeakCurves.Count - 1);
            //    var pc = Ms1PeakCurves[r];
            //    pc.VisualizeBspline(out List<float> rtSeq).Show();
            //    pc.VisualizePeakRegions();
            //}
            //var ms1PeakCurve = Ms1PeakCurves.OrderBy(p => p.MonoisotopicMass).ToList();
            return allPeakCurves;
        }

        public static List<PeakCurve> GetAllPeakCurves_isoEnvelopeTotal(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam, 
            Tolerance peakFindingTolerance, double maxRTRange, out List<Peak>[] allPeaksByScan)
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
            var allPrecursors = new List<DeconvolutedMass>[scans.Length];
            for (int i = 0; i < scans.Length; i++)
            {
                allPrecursors[i] = new List<DeconvolutedMass>();
                var envelopes = Deconvoluter.Deconvolute(scans[i], commonParameters.PrecursorDeconvolutionParameters).OrderByDescending(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.MonoisotopicMass < diaParam.MinMass || envelope.Charge < 5)
                    {
                        continue;
                    }
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, scans[i].RetentionTime, 1, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        scans[i].OneBasedScanNumber, i);
                    allPrecursors[i].Add(precursor);
                }
                allPrecursors[i] = allPrecursors[i].OrderBy(p => p.MonoisotopicMass).ToList();
            }
            var sortedAllPrecursors = allPrecursors.SelectMany(p => p).OrderByDescending(p => p.Envelope.Peaks.Count).ThenByDescending(p => p.HighestPeakIntensity).ToList();

            int index = 1;
            var allPeakEnvelopeCurves = new List<PeakEnvelopeCurve>();
            foreach (var precursor in sortedAllPrecursors)
            {
                var targetPeaks = precursor.Envelope.Peaks.OrderBy(p => p.mz).ToArray();
                var theorIntensityRatio = targetPeaks.Select(p => p.intensity / targetPeaks.Sum(p => p.intensity)).ToArray();
                var highestPeak = targetPeaks.OrderByDescending(p => p.intensity).FirstOrDefault();
                var highestPeakIndex = Array.IndexOf(targetPeaks, highestPeak);
                var newPEC = PeakEnvelopeCurve.GetPeakEnvelopeCurve(targetPeaks, theorIntensityRatio, highestPeakIndex, peakTable, scans, precursor.ZeroBasedScanIndex, diaParam,
                    peakFindingTolerance, maxRTRange);
                if (newPEC != null && newPEC.PeakEnvelopes.Count > 4)
                {
                    allPeakEnvelopeCurves.Add(newPEC);
                    newPEC.Index = index;
                    newPEC.MonoisotopicMass = precursor.MonoisotopicMass;
                    newPEC.Charge = precursor.Charge;
                    newPEC.MsLevel = 1;
                    newPEC.IsolationRange = null;
                    newPEC.GetFakePeakCurve();
                    index++;
                    allPeakCurves.Add(newPEC.FakePeakCurve);
                    newPEC.FakePeakCurve.GetScanCycleSmoothedData(diaParam.ScanCycleSplineTimeInterval);
                }
            }
            //debug
            //for (int i = 0; i < 100; i++)
            //{
            //    Random rnd = new Random();
            //    int r = rnd.Next(Ms1PeakCurves.Count - 1);
            //    var pc = Ms1PeakCurves[r];
            //    pc.VisualizeRaw("point");
            //}
            return allPeakCurves;
        }

        public static PrecursorFragmentsGroup PFgrouping(PeakCurve precursor, List<PeakCurve> fragments, DIAparameters diaParam)
        {
            switch (diaParam.PFGroupingType)
            {
                case PFGroupingType.RetentionTime:
                    return ISDEngine.GroupPrecursorFragments(precursor, fragments, diaParam);
                case PFGroupingType.ScanCycle:
                    return ISDEngine.GroupPrecursorFragments_scanCycle(precursor, fragments, diaParam);
                case PFGroupingType.OverlapAreaRatio:
                    return ISDEngine.GroupPrecursorFragments_area(precursor, fragments, diaParam);
                case PFGroupingType.Area_correlation:
                    return ISDEngine.GroupPrecursorFragments_area_correlation(precursor, fragments, diaParam);
                case PFGroupingType.OverlapFirst:
                    return ISDEngine.GroupPrecursorFragments_overlapFirst(precursor, fragments, diaParam);
                default: return null;
            }
;
        }

        public static Ms2ScanWithSpecificMass ConstructNewMs2Scans(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, DIAparameters diaParam, MsDataFile dataFile)
        {
            switch(diaParam.PseudoMs2ConstructionType)
            {
                case PseudoMs2ConstructionType.mzPeak:
                    return GetPseudoMs2Scan_mzPeak(pfGroup, commonParameters, dataFile);
                case PseudoMs2ConstructionType.neutralMass:
                    return GetPseudoMs2Scan_neutralMass(pfGroup, commonParameters, dataFile);
                default: return null;
            }
        }

        public static Ms2ScanWithSpecificMass GetPseudoMs2Scan_mzPeak(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, MsDataFile dataFile)
        {
            var mzs = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedMz).ToArray();
            var intensities = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedIntensity).ToArray();
            var spectrum = new MzSpectrum(mzs, intensities, false);
            var newMs2Scan = new MsDataScan(spectrum, pfGroup.PrecursorPeakCurve.Index, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                        MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: pfGroup.PrecursorPeakCurve.Index);
            var neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, commonParameters);
            var charge = pfGroup.PrecursorPeakCurve.Charge;
            var monoMz = pfGroup.PrecursorPeakCurve.MonoisotopicMass.ToMz(charge);
            //should highestPeakMz used in ms2withmass???
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge
                , dataFile.FilePath, commonParameters, neutralExperimentalFragments);

            return scanWithprecursor;
        }

        public static Ms2ScanWithSpecificMass GetPseudoMs2Scan_neutralMass(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, MsDataFile dataFile)
        {
            var mzs = new double[] { 1 };
            var intensities = new double[] { Double.MaxValue };
            var spectrum = new MzSpectrum(mzs, intensities, false);
            var newMs2Scan = new MsDataScan(spectrum, pfGroup.PrecursorPeakCurve.Index, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                        MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
            var neutralExperimentalFragments = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.Envelope).ToArray();
            var charge = pfGroup.PrecursorPeakCurve.Charge;
            var monoMz = pfGroup.PrecursorPeakCurve.MonoisotopicMass.ToMz(charge);
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge
                , dataFile.FilePath, commonParameters, neutralExperimentalFragments);

            return scanWithprecursor;
        }

        
    }
}
    
