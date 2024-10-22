using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra.Storage;
using MzLibUtil;
using OpenMcdf.Extensions.OLEProperties;
using Plotly.NET.CSharp;
using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ISDEngine
    {
        public ISDEngine(MsDataFile myFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            MyMSDataFile = myFile;
            CommonParameters = commonParameters;
            DIAparameters = diaParam;
        }

        public DIAparameters DIAparameters { get; set; }
        public List<Peak>[] Ms1PeakTable { get; set; }
        public List<Peak>[] Ms2PeakTable { get; set; }
        public List<PeakCurve> Ms1PeakCurves { get; set; }
        public Dictionary<double, List<PeakCurve>> Ms2PeakCurves { get; set; }
        public MsDataFile MyMSDataFile { get; set; }
        public CommonParameters CommonParameters { get; set; }
        public List<Ms2ScanWithSpecificMass> PseudoMS2Scans { get; set; }
        public Dictionary<double, List<MsDataScan>> DIAScanWindowMap { get; set; }
        public List<PrecursorFragmentsGroup> PFgroups { get; set; }
        public double CycleTime { get; set; }
        public List<Ms2ScanWithSpecificMass> PseudoMs2WithPre { get; set; }

        public void GetPseudoMS2Scans()
        {
            Ms1PeakIndexing();
            ConstructMs2Group();
            GetMs1PeakCurves();
            GetMs2PeakCurves();
            PrecursorFragmentPairing();
            PFgroupFilter();
            //ConstructNewMs2Scans();
        }

        public void Ms1PeakIndexing()
        {
            var ms1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            var allMs1Peaks = Peak.GetAllPeaks(ms1Scans, DIAparameters.PeakSearchBinSize);
            Ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, DIAparameters.PeakSearchBinSize);
        }
        public void ConstructMs2Group()
        {
            DIAScanWindowMap = new Dictionary<double, List<MsDataScan>>();
            var ms2Scans = MyMSDataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
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
        }

        public void GetMs1PeakCurves()
        {
            var allMs1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            Ms1PeakCurves = new List<PeakCurve>();
            int index = 1;
            var allPrecursors = new List<DeconvolutedMass>();
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                //debug
                if (allMs1Scans[i].OneBasedScanNumber == 2267)
                {
                    bool stop = true;
                }
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], CommonParameters.PrecursorDeconvolutionParameters).OrderBy(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.Charge < 6  || envelope.MonoisotopicMass < 5000)
                    {
                        continue;
                    }
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
                    allPrecursors.Add(precursor);
                }
            }
            allPrecursors = allPrecursors.OrderByDescending(p => p.HighestPeakIntensity).ToList();
            var precursorGroups = allPrecursors.GroupBy(p => Math.Round(p.MonoisotopicMass, 1)).OrderByDescending(g => g.Key).Select(g => g.OrderBy(p => p.Charge)).ToList();
            foreach (var precursor in allPrecursors)
            {
                if (precursor.HighestPeakIntensity < DIAparameters.PrecursorIntensityCutOff)
                {
                    continue;
                }
                var peak = PeakCurve.GetPeakFromScan(precursor.HighestPeakMz, Ms1PeakTable, precursor.ZeroBasedScanIndex, new PpmTolerance(0),
                    DIAparameters.PeakSearchBinSize);
                if (peak.PeakCurve == null && peak.Intensity >= DIAparameters.PrecursorIntensityCutOff)
                {
                    var newPeakCurve = PeakCurve.FindPeakCurve(peak, Ms1PeakTable, allMs1Scans, null, DIAparameters.MaxNumMissedScan,
                    DIAparameters.Ms1PeakFindingTolerance, DIAparameters.PeakSearchBinSize, DIAparameters.MaxRTRangeMS1);
                    newPeakCurve.MonoisotopicMass = precursor.MonoisotopicMass;
                    newPeakCurve.Charge = precursor.Charge;
                    if (DIAparameters.SplitMS1Peak)
                    {
                        newPeakCurve.DetectPeakRegions();
                        var newPCs = newPeakCurve.SeparatePeakByRegion();
                        foreach (var pc in newPCs)
                        {
                            if (pc.Peaks.Count > 4)
                            {
                                Ms1PeakCurves.Add(pc);
                            }
                        }
                    }
                    else
                    {
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            Ms1PeakCurves.Add(newPeakCurve);
                            newPeakCurve.Index = index;
                            index++;
                        }
                    }
                }
            }
        }

        public void GetMs2PeakCurves()
        {
            var ms2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                var allMs2Peaks = Peak.GetAllPeaks(ms2scans, DIAparameters.PeakSearchBinSize);
                var rankedMs2Peaks = allMs2Peaks.OrderByDescending(p => p.Intensity).ToList();
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, DIAparameters.PeakSearchBinSize);
                ms2PeakCurves[ms2Group.Key] = new List<PeakCurve>();
                foreach (var peak in rankedMs2Peaks)
                {
                    if (peak.PeakCurve == null)
                    {
                        var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms2PeakTable, ms2scans, ms2scans[0].IsolationRange,
                            DIAparameters.MaxNumMissedScan, DIAparameters.Ms2PeakFindingTolerance, DIAparameters.PeakSearchBinSize, DIAparameters.MaxRTRangeMS2);
                        if (DIAparameters.SplitMS2Peak)
                        {
                            if (newPeakCurve.Peaks.Count > 4)
                            {
                                newPeakCurve.DetectPeakRegions();
                                var newPCs = newPeakCurve.SeparatePeakByRegion();
                                foreach (var pc in newPCs)
                                {
                                    if (pc.Peaks.Count > 4)
                                    {
                                        ms2PeakCurves[ms2Group.Key].Add(pc);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (newPeakCurve.Peaks.Count > 4)
                            {
                                ms2PeakCurves[ms2Group.Key].Add(newPeakCurve);
                            }
                        }
                    }
                }
            }
            Ms2PeakCurves = ms2PeakCurves;

            //debug
            //var testMs2PeakCurve = Ms2PeakCurves.Values.SelectMany(v => v).ToList();
            //foreach (var a in testMs2PeakCurve)
            //{
            //    Random rnd = new Random();
            //    int r = rnd.Next(testMs2PeakCurve.Count - 1);
            //    var pc = testMs2PeakCurve[r];
            //    pc.VisualizeBspline(out List<float> rtSeq).Show();
            //    pc.VisualizePeakRegions();
            //}
        }

        public void PrecursorFragmentPairing()
        {
            PFgroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2group in Ms2PeakCurves)
            {
                var precursorsInRange = Ms1PeakCurves.ToArray();

                Parallel.ForEach(Partitioner.Create(0, precursorsInRange.Length), new ParallelOptions { MaxDegreeOfParallelism = 18 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursorsInRange[i];
                        var preFragGroup = GroupPrecursorFragments(precursor, ms2group.Value, DIAparameters);

                        if (preFragGroup != null)
                        {
                            PFgroups.Add(preFragGroup);
                        }
                    }
                });
            }
            //debug
            //var rankedPFgroups = PFgroups.OrderByDescending(pf => pf.PFpairs.Count).ToList();

            //TODO?
            //Check the fragment correlations within each pfgroup for filtering
        }

        public void PFgroupFilter()
        {
            var allMatchedMs2 = Ms2PeakCurves.SelectMany(p => p.Value).Where(pc => pc.PFpairs.Count > 0).ToList();
            foreach (var ms2curve in allMatchedMs2)
            {
                ms2curve.GetPrecursorRanks();
            }
            foreach (var group in PFgroups)
            {
                group.PFpairs = group.PFpairs.Where(pf => pf.PrecursorRank <= DIAparameters.PrecursorRankCutOff).ToList();
                group.GetNumberOfHighCorrFragments(DIAparameters);
            }
            PFgroups = PFgroups.Where(pf => pf.PFpairs.Count > 0 && pf.NumHighCorrFragments >= DIAparameters.NumHighCorrFragments).ToList();
        }

        public static PrecursorFragmentsGroup GroupPrecursorFragments(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters DIAparameters)
        {
            var preFragGroup = new PrecursorFragmentsGroup(precursor);
            foreach (var ms2curve in ms2curves)
            {
                if (ms2curve.ApexRT >= precursor.StartRT && ms2curve.ApexRT <= precursor.EndRT)
                {
                    if (Math.Abs(ms2curve.ApexRT - precursor.ApexRT) <= DIAparameters.ApexRtTolerance)
                    {
                        var overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursor, ms2curve);
                        if (overlap > DIAparameters.OverlapRatioCutOff)
                        {
                            double corr = PrecursorFragmentPair.CalculatePeakCurveCorr(precursor, ms2curve);
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
            }
            if (preFragGroup.PFpairs.Count > DIAparameters.FragmentRankCutOff)
            {
                var filtered = preFragGroup.PFpairs.OrderByDescending(pair => pair.Correlation).Take(DIAparameters.FragmentRankCutOff);
                preFragGroup.PFpairs = filtered.ToList();
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
    }
    
}
