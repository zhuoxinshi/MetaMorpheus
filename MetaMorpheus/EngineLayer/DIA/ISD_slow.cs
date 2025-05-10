using MassSpectrometry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ISD_slow
    {
        public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();

            //Get ms1 XICs
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance, diaParam.MaxRTRangeMS1,
                out List<Peak>[] peaksByScan).ToArray();

            //calculate number of scans per cycle
            int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;
            diaParam.NumScansPerCycle = scansPerCycle;

            //Get ms2 scans
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var isdScanVoltageMap = ISDEngine_static.ConstructMs2Groups(ms2Scans);

            //precursor fragment grouping for each precursor
            ISDEngine_static.PeakCurveSpline(allMs1PeakCurves.ToList(), diaParam.Ms1SplineType, diaParam, ms1Scans, ms2Scans);
            var pfGroups = new List<PrecursorFragmentsGroup>();
            Parallel.ForEach(Partitioner.Create(0, allMs1PeakCurves.Length), new ParallelOptions { MaxDegreeOfParallelism = 15 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = allMs1PeakCurves[i];
                        foreach (var ms2group in isdScanVoltageMap.Values)
                        {
                            var preFragGroup = ISD_slow.FindFragments(precursor, ms1Scans, ms2group.ToArray(), commonParameters, diaParam);
                            if (preFragGroup != null)
                            {
                                lock (pfGroups)
                                {
                                    pfGroups.Add(preFragGroup);
                                }
                            }
                        }
                    }
                });

            //construct new ms2Scans
            foreach (var pfGroup in pfGroups)
            {
                var newScans = ISDEngine_static.ConstructNewMs2Scans(pfGroup, commonParameters, diaParam.PseudoMs2ConstructionType, dataFile.FilePath);
                pseudoMs2Scans.Add(newScans);
            }
            return pseudoMs2Scans;
        }

        public static PrecursorFragmentsGroup FindFragments(PeakCurve precursor, MsDataScan[] ms1scans, MsDataScan[] ms2scans, CommonParameters commonParameters, DIAparameters diaParam)
        {
            //Get all ms2 XICs in range
            double cycleTime = Math.Ceiling((ms2scans[1].RetentionTime - ms2scans[0].RetentionTime) * 100) / 100;
            double maxRTRange = precursor.EndRT - precursor.StartRT;
            var scans = ms2scans.Where(s => s.RetentionTime >= precursor.StartRT - cycleTime && s.RetentionTime <= precursor.EndRT + cycleTime).ToArray();
            var allMs2PeakCurves = ISDEngine_static.GetAllPeakCurves(scans, commonParameters, diaParam, diaParam.Ms2XICType, diaParam.Ms2PeakFindingTolerance, maxRTRange, 
                out List<Peak>[] peaksByScan);
            ISDEngine_static.PeakCurveSpline(allMs2PeakCurves, diaParam.Ms2SplineType, diaParam, ms1scans, ms2scans);

            var pfGroup = ISDEngine_static.PFgrouping(precursor, allMs2PeakCurves, diaParam);
            return pfGroup;
        }

        
    }
}
