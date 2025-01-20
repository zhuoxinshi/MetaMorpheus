using MassSpectrometry;
using MzLibUtil;
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

            //calculate number of scans per cycle
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;
            diaParam.NumScansPerCycle = scansPerCycle;

            //construct DIAScanWindowMap
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var diaScanWindowMap = DIAEngine_static.ConstructMs2Groups(ms2Scans);

            //Get ms1 peakCurves
            var maxScanNum = ms1Scans[ms1Scans.Length - 1].OneBasedScanNumber;
            var ms1PeakList = new List<Peak>[maxScanNum + 1];
            var allMs1PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            foreach (var ms1window in diaScanWindowMap.Keys)
            {
                var ms1Range = new MzRange(ms1window.min, ms1window.max);
                allMs1PeakCurves[ms1window] = ISDEngine_static.GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance,
                    diaParam.MaxRTRangeMS1, out ms1PeakList, diaParam.CutMs1Peaks, ms1Range);
                if (allMs1PeakCurves[ms1window].Count == 0)
                {
                    diaScanWindowMap.Remove(ms1window);
                }
            }
            var allMs1Peaks = ms1PeakList.Where(v => v != null).SelectMany(p => p).ToList();
            var ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, diaParam.PeakSearchBinSize);
            
            //get ms2withmass and ms2 peakCurves
            var ms2WithMass = ISD_scanBased.GetMs2ScansWithMass(ms1Scans, ms2Scans, dataFile.FilePath, commonParameters, diaParam);
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

            //precursor fragment rank filtering
            foreach (var ms2curve in allMs2PeakCurves.Values.SelectMany(p => p))
            {
                ms2curve.GetPrecursorRanks();
            }
            foreach (var group in pfGroups)
            {
                group.PFpairs = group.PFpairs.OrderByDescending(pf => pf.Correlation).Take(diaParam.FragmentRankCutOff).ToList();
                //group.PFpairs = group.PFpairs.Where(pf => pf.PrecursorRank <= diaParam.PrecursorRankCutOff).ToList();
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

        //public static void GetPrecursorPeakCurve(List<Ms2ScanWithSpecificMass>[] scansWithPrecursor, MsDataScan[] ms1scans, List<Peak>[] allPeaks,
        //    List<Peak>[] ms1PeakTable, DIAparameters DIAparam, Dictionary<int, int> scanIndexMap)
        //{
        //    for (int i = 0; i < scansWithPrecursor.Length; i++)
        //    {
        //        foreach (var scan in scansWithPrecursor[i])
        //        {
        //            var preScan = ms1scans.Where(s => s.OneBasedScanNumber == scan.OneBasedPrecursorScanNumber).First();
        //            var precursorPeak = PeakCurve.GetPeakFromScan(scan.HighestPeakMz, ms1PeakTable, scanIndexMap[preScan.OneBasedScanNumber], new PpmTolerance(0), DIAparam.PeakSearchBinSize);
        //            if (precursorPeak.PeakCurve == null)
        //            {
        //                scan.PrecursorPeakCurve = PeakCurve.FindPeakCurve(precursorPeak, ms1PeakTable, ms1scans, null, DIAparam.MaxNumMissedScan, DIAparam.Ms1PeakFindingTolerance,
        //                    DIAparam.PeakSearchBinSize, DIAparam.MaxRTRange);
        //                scan.PrecursorPeakCurve.MonoisotopicMass = scan.PrecursorMass;
        //                scan.PrecursorPeakCurve.Charge = scan.PrecursorCharge;
        //            }
        //            else
        //            {
        //                scan.PrecursorPeakCurve = precursorPeak.PeakCurve;
        //            }
        //        }
        //    }
        //}
    }
}
