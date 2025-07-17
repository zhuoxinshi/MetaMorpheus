using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Plotly.NET.TraceObjects;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http.Headers;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DIAEngine_static
    {
        public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();

            //calculate number of scans per cycle
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;
            diaParam.NumScansPerCycle = scansPerCycle;

            //construct DIAScanWindowMap
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var DIAScanWindowMap = ConstructMs2Groups(ms2Scans);

            //Get ms1 XICs
            var allMs1PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            foreach (var ms1window in DIAScanWindowMap.Keys)
            {
                var ms1Range = new MzRange(ms1window.min, ms1window.max);
                allMs1PeakCurves[ms1window] = ISDEngine_static.GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance, 
                    diaParam.MaxRTRangeMS1, out List<Peak>[] allPeaksByScan, diaParam.CutMs1Peaks, ms1Range);
                if (allMs1PeakCurves[ms1window].Count == 0)
                {
                    DIAScanWindowMap.Remove(ms1window);
                }
            }
            ISDEngine_static.PeakCurveSpline(allMs1PeakCurves.Values.SelectMany(p => p).Where(p =>p.Peaks.Count > 4).ToList(), diaParam.Ms1SplineType, diaParam, ms1Scans, ms2Scans);
            ISDEngine_static.PeakCurveSpline(allMs1PeakCurves.Values.SelectMany(p => p).Where(p => p.Peaks.Count <= 4).ToList(), SplineType.UmpireBSpline, diaParam, ms1Scans, ms2Scans);

            //Group Ms1PeakCurves by mass and apexRT
            //var ms1PeakCurveGroups = new Dictionary<(double min, double max), List<PeakCurve>>();
            //foreach (var group in allMs1PeakCurves)
            //{
            //    var pcs = group.Value.GroupBy(pc => new { mass = Math.Round(pc.MonoisotopicMass, 1), apexRt = pc.ApexRT }).ToList();
            //}

            //Get rtMap
            var rtMap = ISDEngine_static.GetRtMap(ms1Scans, ms2Scans);

            //Get ms2 XICs and ms1SpaceSpline
            var allMs2PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                allMs2PeakCurves[ms2Group.Key] = ISDEngine_static.GetAllPeakCurves(ms2Group.Value.ToArray(), commonParameters, diaParam, diaParam.Ms2XICType,
                    diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out List<Peak>[] peaksByScan2);
                ISDEngine_static.PeakCurveSpline(allMs2PeakCurves[ms2Group.Key].Where(p => p.Peaks.Count > 4).ToList(), diaParam.Ms2SplineType, diaParam, ms1Scans, ms2Scans);
                ISDEngine_static.PeakCurveSpline(allMs2PeakCurves[ms2Group.Key].Where(p => p.Peaks.Count <= 4).ToList(), SplineType.UmpireBSpline, diaParam, ms1Scans, ms2Scans);
            }

            //precursor fragment grouping
            var pfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2group in allMs2PeakCurves)
            {
                var precursorsInRange = allMs1PeakCurves[ms2group.Key].ToArray();

                Parallel.ForEach(Partitioner.Create(0, precursorsInRange.Length), new ParallelOptions { MaxDegreeOfParallelism = 10 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursorsInRange[i];
                        var preFragGroup = ISDEngine_static.PFgrouping(precursor, allMs2PeakCurves[ms2group.Key], diaParam);

                        if (preFragGroup != null)
                        {
                            lock (pfGroups)
                            {
                                pfGroups.Add(preFragGroup);
                            }
                        }
                    }
                });
            }

            //precursor fragment grouping 2
            //var pfGroups = new List<PrecursorFragmentsGroup>();
            //var allFragments = allMs2PeakCurves.Values.SelectMany(p => p).ToList();
            //var allPrecursors = allMs1PeakCurves.Values.SelectMany(p => p).ToList();
            //Parallel.ForEach(Partitioner.Create(0, allPrecursors.Count), new ParallelOptions { MaxDegreeOfParallelism = 15 },
            //    (partitionRange, loopState) =>
            //    {
            //        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
            //        {
            //            var precursor = allPrecursors[i];
            //            var pfGroup = ISDEngine_static.PFgrouping(precursor, allFragments, diaParam);
            //            if (pfGroup != null && pfGroup.PFpairs.Count > 0)
            //            {
            //                lock (pfGroups)
            //                    pfGroups.Add(pfGroup);
            //                //pfGroup.VisualizeXYData().Show();
            //            }
            //        }
            //    });


            //pfGroups = pfGroups.Where(pf => pf.PFpairs.Count > 0 && pf.NumHighCorrFragments >= diaParam.NumHighCorrFragments).ToList();

            //combine precursor fragment groups
            //var groupedPfGroups = pfGroups
            //    .GroupBy(g => (Math.Round(g.PrecursorPeakCurve.MonoisotopicMass, 1), g.PrecursorPeakCurve.ApexRT))
            //    .ToList();
            //var combinedPFgroups = new List<PrecursorFragmentsGroup>();
            //int i = 1;
            //foreach (var group in groupedPfGroups)
            //{
            //    var highestIntensityPrecursor = group.Select(g => g.PrecursorPeakCurve).OrderByDescending(p => p.ApexIntensity).First();
            //    var combinedGroup = new PrecursorFragmentsGroup(highestIntensityPrecursor, group.SelectMany(p => p.PFpairs).ToList(), i);
            //    combinedPFgroups.Add(combinedGroup);
            //    i++;
            //}

            //precursor fragment rank filtering
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

            //write pfGroup file
            //if (diaParam.PFgroupsDictionary == null)
            //{
            //    diaParam.PFgroupsDictionary = new Dictionary<string, List<PrecursorFragmentsGroup>>();
            //}
            //diaParam.PFgroupsDictionary[dataFile.FilePath] = pfGroups;
            //if (diaParam.PeakCurveDictionary == null)
            //{
            //    diaParam.PeakCurveDictionary = new Dictionary<string, List<PeakCurve>>();
            //}
            //diaParam.PeakCurveDictionary[dataFile.FilePath] = allMs1PeakCurves.SelectMany(p => p.Value).Concat(allMs2PeakCurves.SelectMany(p => p.Value)).ToList();

            //construct new ms2Scans
            int pfGroupIndex = 1;
            foreach (var pfGroup in pfGroups)
            {
                pfGroup.PFgroupIndex = pfGroupIndex;
                var newScans = ISDEngine_static.ConstructNewMs2Scans(pfGroup, commonParameters, diaParam.PseudoMs2ConstructionType, dataFile.FilePath);
                pseudoMs2Scans.Add(newScans);
                pfGroupIndex++;
            }

            return pseudoMs2Scans;
        }

        public static Dictionary<(double min, double max), List<MsDataScan>> ConstructMs2Groups(MsDataScan[] ms2Scans)
        {
            var DIAScanWindowMap = new Dictionary<(double min, double max), List<MsDataScan>>();
            foreach (var ms2 in ms2Scans)
            {
                (double min, double max) range = new(Math.Round(ms2.IsolationRange.Minimum, 0), Math.Round(ms2.IsolationRange.Maximum, 0));
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
    }
}
