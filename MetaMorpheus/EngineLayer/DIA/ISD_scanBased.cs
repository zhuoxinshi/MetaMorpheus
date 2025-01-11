using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using SpectralAveraging;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ISD_scanBased
    {
        public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();

            //calculate number of scans per cycle
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            int scansPerCycle = ms1Scans[1].OneBasedScanNumber - ms1Scans[0].OneBasedScanNumber;
            diaParam.NumScansPerCycle = scansPerCycle;

            //Get ms1 peakCurves
            var maxScanNum = ms1Scans[ms1Scans.Length - 1].OneBasedScanNumber;
            var ms1PeakList = new List<Peak>[maxScanNum + 1];
            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, commonParameters, diaParam, diaParam.Ms1XICType, diaParam.Ms1PeakFindingTolerance,
                diaParam.MaxRTRangeMS1, out ms1PeakList, diaParam.CutMs1Peaks);
            var allMs1Peaks = ms1PeakList.Where(v => v != null).SelectMany(p => p).ToList();
            var ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, diaParam.PeakSearchBinSize);

            //get ms2withmass and ms2 peakCurves
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var ms2WithMass = GetMs2ScansWithMass(ms1Scans, ms2Scans, dataFile.FilePath, commonParameters, diaParam);
            var isdScanVoltageMap = ISDEngine_static.ConstructMs2Groups(ms2Scans);
            var allMs2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            var ms2PeakLists = new Dictionary<double, List<Peak>[]>();
            foreach (var ms2Group in isdScanVoltageMap)
            {
                var maxNum = ms2Group.Value.Last().OneBasedScanNumber;
                var ms2Peaks = new List<Peak>[maxNum + 1];
                allMs2PeakCurves[ms2Group.Key] = ISDEngine_static.GetAllPeakCurves(ms2Group.Value.ToArray(), commonParameters, diaParam, diaParam.Ms2XICType,
                    diaParam.Ms2PeakFindingTolerance, diaParam.MaxRTRangeMS2, out ms2Peaks, diaParam.CutMs2Peaks);
                ms2PeakLists[ms2Group.Key] = ms2Peaks;
            }
            var maxScanNumMs2 = ms2Scans[ms2Scans.Length - 1].OneBasedScanNumber;
            var allMs2Peaks = new List<Peak>[maxScanNumMs2 + 1];
            for (int i = 1; i < allMs2Peaks.Length + 1; i++)
            {
                foreach(var ms2PeakList in ms2PeakLists)
                {
                    if (i < ms2PeakList.Value.Length && ms2PeakList.Value[i] != null)
                    {
                        allMs2Peaks[i] = ms2PeakList.Value[i];
                    }
                }
            }

            //precursor fragment grouping
            var pfGroups = new List<PrecursorFragmentsGroup>();
            var allMs2WithMass = ms2WithMass.Where(v => v!= null).SelectMany(s => s).ToList();
            foreach (var ms2 in allMs2WithMass)
            {
                var group = PFgrouping_scanBased(ms2, ms1PeakTable, scansPerCycle, allMs2Peaks, diaParam);
                if (group != null)
                {
                    pfGroups.Add(group);
                }
            }

            //Combine fragments for the same precursor
            //pfGroups = CombinePFGroups(pfGroups);

            //construct new ms2Scans
            foreach (var pfGroup in pfGroups)
            {
                var newScans = ISDEngine_static.ConstructNewMs2Scans(pfGroup, commonParameters, diaParam.PseudoMs2ConstructionType, dataFile.FilePath);
                pseudoMs2Scans.Add(newScans);
            }
            return pseudoMs2Scans;
        }


        public static List<Peak>[] GetPeakTable(List<Peak>[] allPeaks, int binsPerDalton, out Dictionary<int, int> scanIndexMap)
        {
            var table = new List<Peak>[(int)Math.Ceiling(allPeaks.Where(v => v != null).SelectMany(p => p).Max(p => p.Mz) * binsPerDalton) + 1];
            int zeroBasedScanIndex = 0;
            scanIndexMap = new Dictionary<int, int>();

            for (int i = 0; i < allPeaks.Length; i++)
            {
                scanIndexMap.Add(allPeaks[i].FirstOrDefault().ScanNumber, zeroBasedScanIndex);
                for (int j = 0; j < allPeaks[i].Count; j++)
                {
                    //Label the peak with zeroBasedScanIndex
                    allPeaks[i][j].ZeroBasedScanIndex = zeroBasedScanIndex;

                    int roundedMz = (int)Math.Round(allPeaks[i][j].Mz * binsPerDalton, 0);

                    if (table[roundedMz] == null)
                    {
                        table[roundedMz] = new List<Peak>();
                    }
                    table[roundedMz].Add(allPeaks[i][j]);
                }
                zeroBasedScanIndex++;
            }
            return table;
        }

        public static void GetPrecursorPeakCurve(List<Ms2ScanWithSpecificMass>[] scansWithPrecursor, MsDataScan[] ms1scans, List<Peak>[] allPeaks,
            List<Peak>[] ms1PeakTable, DIAparameters DIAparam, Dictionary<int, int> scanIndexMap)
        {
            for (int i = 0; i < scansWithPrecursor.Length; i++)
            {
                foreach (var scan in scansWithPrecursor[i])
                {
                    var preScan = ms1scans.Where(s => s.OneBasedScanNumber == scan.OneBasedPrecursorScanNumber).First();
                    var precursorPeak = PeakCurve.GetPeakFromScan(scan.HighestPrecursorPeakMz, ms1PeakTable, scanIndexMap[preScan.OneBasedScanNumber], new PpmTolerance(0), DIAparam.PeakSearchBinSize);
                    if (precursorPeak.PeakCurve == null)
                    {
                        scan.PrecursorPeakCurve = PeakCurve.FindPeakCurve(precursorPeak, ms1PeakTable, ms1scans, null, DIAparam.MaxNumMissedScan, DIAparam.Ms1PeakFindingTolerance,
                            DIAparam.PeakSearchBinSize, DIAparam.MaxRTRangeMS1);
                        scan.PrecursorPeakCurve.MonoisotopicMass = scan.PrecursorMass;
                        scan.PrecursorPeakCurve.Charge = scan.PrecursorCharge;
                    }
                    else
                    {
                        scan.PrecursorPeakCurve = precursorPeak.PeakCurve;
                    }
                }
            }
        }

        public static List<Ms2ScanWithSpecificMass>[] GetMs2ScansWithMass(MsDataScan[] ms1scans, MsDataScan[] ms2Scans, string fullFilePath, CommonParameters commonParameters,
            DIAparameters diaParam)
        {
            List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];

            if (!ms2Scans.Any())
            {
                return scansWithPrecursors;
            }

            Parallel.ForEach(Partitioner.Create(0, ms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (partitionRange, loopState) =>
                {
                    var precursors = new List<(double MonoPeakMz, int Charge, double Intensity, int PeakCount, double highestPeakMz, 
                        IsotopicEnvelope envelope, double? FractionalIntensity)>();

                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        if (GlobalVariables.StopLoops) { break; }

                        precursors.Clear();
                        MsDataScan ms2scan = ms2Scans[i];

                        if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                        {
                            MsDataScan precursorSpectrum = ms1scans.First(s => s.OneBasedScanNumber == ms2scan.OneBasedPrecursorScanNumber);

                            if (commonParameters.DoPrecursorDeconvolution)
                            {
                                foreach (IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
                                    precursorSpectrum.MassSpectrum, commonParameters.PrecursorDeconvolutionParameters))
                                {
                                    if (envelope.MonoisotopicMass < diaParam.MinMass)
                                    {
                                        continue;
                                    }
                                    double monoPeakMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
                                    int peakCount = envelope.Peaks.Count();
                                    double intensity = 1;
                                    if (commonParameters.UseMostAbundantPrecursorIntensity)
                                    {
                                        intensity = envelope.Peaks.Max(p => p.intensity);
                                    }
                                    else
                                    {
                                        intensity = envelope.Peaks.Sum(p => p.intensity);
                                    }

                                    var fractionalIntensity = envelope.TotalIntensity /
                                          (double)precursorSpectrum.MassSpectrum.YArray
                                          [
                                              precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Minimum)..precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Maximum)
                                          ].Sum();

                                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).First().mz;

                                    precursors.Add((monoPeakMz, envelope.Charge, intensity, peakCount, highestPeakMz, envelope, 
                                        fractionalIntensity));
                                }
                            }
                        }

                        scansWithPrecursors[i] = new List<Ms2ScanWithSpecificMass>();
                        IsotopicEnvelope[] neutralExperimentalFragments = null;

                        foreach (var precursor in precursors)
                        {
                            // assign precursor for this MS2 scan
                            var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.MonoPeakMz,
                                precursor.Charge, fullFilePath, commonParameters, neutralExperimentalFragments,
                                precursor.Intensity, precursor.PeakCount, precursor.FractionalIntensity, precursor.highestPeakMz);
                            scan.HighestPrecursorPeakMz = precursor.highestPeakMz;

                            scansWithPrecursors[i].Add(scan);
                        }
                    }
                });
            //debug
            var preCount = scansWithPrecursors.SelectMany(p => p).Count();

            return scansWithPrecursors;
        }

        public static PrecursorFragmentsGroup PFgrouping_scanBased(Ms2ScanWithSpecificMass ms2WithMass, List<Peak>[] ms1PeakTable, int numScansPerCycle, List<Peak>[] ms2PeaksByScan,
            DIAparameters diaParam)
        {
            int zeroBasedScanIndex = (ms2WithMass.OneBasedScanNumber - 1) / numScansPerCycle;
            var precursorPeak = PeakCurve.GetPeakFromScan(ms2WithMass.HighestPrecursorPeakMz, ms1PeakTable, zeroBasedScanIndex, diaParam.Ms1PeakFindingTolerance,
                diaParam.PeakSearchBinSize);
            if (precursorPeak.PeakCurve == null || precursorPeak.PeakCurve.Peaks.Count < 5)
            {
                return null;
            }
            var fragmentsPCs = ms2PeaksByScan[ms2WithMass.OneBasedScanNumber].Where(f => f.PeakCurve != null).
                Select(f => f.PeakCurve).Where(p => p.Peaks.Count > 4).ToList();

            var group = ISDEngine_static.PFgrouping(precursorPeak.PeakCurve, fragmentsPCs, diaParam);
            return group;
        }

        public static List<PrecursorFragmentsGroup> CombinePFGroups(List<PrecursorFragmentsGroup> pfGroups)
        {
            var newGroups = new List<PrecursorFragmentsGroup>();
            var groupByPrecursorPC = pfGroups.GroupBy(g => g.PrecursorPeakCurve).ToList();
            foreach(var group in groupByPrecursorPC)
            {
                var allFragments = group.SelectMany(g => g.PFpairs).Select(pair => pair.FragmentPeakCurve).Distinct().ToList();
                var precursor = group.Key;
                var newGroup = new PrecursorFragmentsGroup(precursor, allFragments);
                newGroups.Add(newGroup);
            }
            return newGroups;
        }
    }
}
