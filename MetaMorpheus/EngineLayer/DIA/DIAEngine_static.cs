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

            //Get rtMap
            var rtMap = ISDEngine_static.GetRtMap(ms1Scans, ms2Scans);

            //Get ms2 XICs and ms1SpaceSpline
            var allMs2PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                allMs2PeakCurves[ms2Group.Key] = ISDEngine_static.GetAllPeakCurves(ms2Group.Value.ToArray(), commonParameters, diaParam, diaParam.Ms2XICType,
                    diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out List<Peak>[] peaksByScan2);
                foreach(var peakCurve in allMs2PeakCurves[ms2Group.Key])
                {
                    peakCurve.GetMs1SpaceSpline(rtMap, "cubic");
                }
            }

            //precursor fragment grouping
            var pfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2group in allMs2PeakCurves)
            {
                var precursorsInRange = allMs1PeakCurves[ms2group.Key].ToArray();

                Parallel.ForEach(Partitioner.Create(0, precursorsInRange.Length), new ParallelOptions { MaxDegreeOfParallelism = 18 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursorsInRange[i];
                        var preFragGroup = ISDEngine_static.PFgrouping(precursor, allMs2PeakCurves[ms2group.Key], diaParam);

                        if (preFragGroup != null)
                        {
                            pfGroups.Add(preFragGroup);
                        }
                    }
                });
            }

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

            //construct new ms2Scans
            foreach (var pfGroup in pfGroups)
            {
                var newScans = ISDEngine_static.ConstructNewMs2Scans(pfGroup, commonParameters, diaParam.PseudoMs2ConstructionType, dataFile.FilePath);
                pseudoMs2Scans.Add(newScans);
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
