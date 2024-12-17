using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TopDownProteomics.ProForma;

namespace EngineLayer.DIA
{
    public class DIA_scanBased
    {
        public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();

            //Get ms1 peakCurves
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            var maxScanNum = ms1Scans[ms1Scans.Length - 1].OneBasedScanNumber;
            var ms1PeakList = new List<Peak>[maxScanNum + 1];
            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance,
                diaParam.MaxRTRangeMS1, out ms1PeakList);
            var allMs1Peaks = ms1PeakList.Where(v => v != null).SelectMany(p => p).ToList();
            var ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, diaParam.PeakSearchBinSize);

            //calculate number of scans per cycle
            int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;

            //get ms2withmass and ms2 peakCurves
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var ms2WithMass = ISD_scanBased.GetMs2ScansWithMass(ms1Scans, ms2Scans, dataFile.FilePath, commonParameters, diaParam);
            var diaScanWindowMap = ConstructMs2Groups(ms2Scans);
            var allMs2PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            var ms2PeakLists = new Dictionary<(double min, double max), List<Peak>[]>();
            foreach (var ms2Group in diaScanWindowMap)
            {
                var maxNum = ms2Group.Value.Last().OneBasedScanNumber;
                var ms2Peaks = new List<Peak>[maxNum + 1];
                allMs2PeakCurves[ms2Group.Key] = ISDEngine_static.GetAllPeakCurves(ms2Group.Value.ToArray(), commonParameters, diaParam, diaParam.Ms2XICType,
                    diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out ms2Peaks);
                ms2PeakLists[ms2Group.Key] = ms2Peaks;
            }
            var maxScanNumMs2 = ms2Scans[ms2Scans.Length - 1].OneBasedScanNumber;
            var allMs2Peaks = new List<Peak>[maxScanNumMs2 + 1];
            for (int i = 1; i < allMs2Peaks.Length + 1; i++)
            {
                foreach (var ms2PeakList in ms2PeakLists)
                {
                    if (i < ms2PeakList.Value.Length && ms2PeakList.Value[i] != null)
                    {
                        allMs2Peaks[i] = ms2PeakList.Value[i];
                    }
                }
            }

            //precursor fragment grouping
            var pfGroups = new List<PrecursorFragmentsGroup>();
            var allMs2WithMass = ms2WithMass.Where(v => v != null).SelectMany(s => s).ToList();
            foreach (var ms2 in allMs2WithMass)
            {
                var group = ISD_scanBased.PFgrouping_scanBased(ms2, ms1PeakTable, scansPerCycle, allMs2Peaks, diaParam);
                if (group != null)
                {
                    pfGroups.Add(group);
                }
            }

            //pfGroup filtering
            foreach(var pfGroup in pfGroups)
            {
                if (pfGroup.PFpairs.Count > diaParam.FragmentRankCutOff)
                {
                    var filtered = pfGroup.PFpairs.OrderByDescending(pair => pair.Correlation).Take(diaParam.FragmentRankCutOff);
                    pfGroup.PFpairs = filtered.ToList();
                    foreach (var ms2curve in filtered.Select(pair => pair.FragmentPeakCurve))
                    {
                        ms2curve.GetPrecursorRanks();
                    }
                    pfGroup.PFpairs = pfGroup.PFpairs.Where(pf => pf.PrecursorRank <= diaParam.PrecursorRankCutOff).ToList();
                }
            }

            //construct new ms2Scans
            foreach (var pfGroup in pfGroups)
            {
                var newScans = ISDEngine_static.ConstructNewMs2Scans(pfGroup, commonParameters, diaParam, dataFile);
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
