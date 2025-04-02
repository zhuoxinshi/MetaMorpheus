using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ISDEngine_massCurve
    {
        //public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        //{
        //    var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();

        //    //calculate number of scans per cycle
        //    var ms1Scans = dataFile.GetMS1Scans().ToArray();
        //    var ms2Scans = dataFile.Where(s => s.MsnOrder == 2).ToArray();
        //    int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;
        //    diaParam.NumScansPerCycle = scansPerCycle;

        //    //Get ms1 XICs
        //    var allMs1PeakCurves = GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance, diaParam.MaxRTRangeMS1,
        //        out List<Peak>[] peaksByScan, diaParam.CutMs1Peaks);
        //    PeakCurveSpline(allMs1PeakCurves, diaParam.Ms1SplineType, diaParam, ms1Scans, ms2Scans);

        //    //Get ms2 XICs
        //    var isdScanVoltageMap = ConstructMs2Groups(ms2Scans);
        //    var allMs2PeakCurves = new Dictionary<double, List<PeakCurve>>();
        //    foreach (var ms2Group in isdScanVoltageMap)
        //    {
        //        allMs2PeakCurves[ms2Group.Key] = GetAllPeakCurves(ms2Group.Value.ToArray(), commonParameters, diaParam, diaParam.Ms2XICType,
        //            diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out List<Peak>[] peaksByScan2, diaParam.CutMs2Peaks);
        //        PeakCurveSpline(allMs2PeakCurves[ms2Group.Key], diaParam.Ms2SplineType, diaParam, ms1Scans, ms2Scans);
        //    }

        //    //precursor fragment grouping
        //    var pfGroups = new List<PrecursorFragmentsGroup>();
        //    foreach (var precursor in allMs1PeakCurves)
        //    {
        //        foreach (var fragments in allMs2PeakCurves.Values)
        //        {
        //            var pfGroup = PFgrouping(precursor, fragments, diaParam);
        //            if (pfGroup != null && pfGroup.PFpairs.Count > 0)
        //            {
        //                pfGroups.Add(pfGroup);
        //            }
        //        }
        //    }

        //    //precursor fragment rank filtering
        //    //foreach (var ms2curve in allMs2PeakCurves.Values.SelectMany(p => p))
        //    //{
        //    //    ms2curve.GetPrecursorRanks();
        //    //}
        //    //foreach (var group in pfGroups)
        //    //{
        //    //    group.PFpairs = group.PFpairs.OrderByDescending(pf => pf.Correlation).Take(diaParam.FragmentRankCutOff).ToList();
        //    //    group.PFpairs = group.PFpairs.Where(pf => pf.PrecursorRank <= diaParam.PrecursorRankCutOff).ToList();
        //    //    group.GetNumberOfHighCorrFragments(diaParam);
        //    //}
        //    //pfGroups = pfGroups.Where(pf => pf.PFpairs.Count > 0 && pf.NumHighCorrFragments >= diaParam.NumHighCorrFragments).ToList();

        //    //construct new ms2Scans
        //    foreach (var pfGroup in pfGroups)
        //    {
        //        var newScans = ISDEngine_static.ConstructNewMs2Scans(pfGroup, commonParameters, diaParam.PseudoMs2ConstructionType, dataFile.FilePath);
        //        pseudoMs2Scans.Add(newScans);
        //    }
        //    return pseudoMs2Scans;
        //}

        
    }
}
