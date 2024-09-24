using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra.Storage;
using MzLibUtil;
using OpenMcdf.Extensions.OLEProperties;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http.Headers;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DIAEngine2
    {
        public DIAEngine2(MsDataFile myFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            MyMSDataFile = myFile;
            CommonParameters = commonParameters;
            DIAparameters = diaParam;
        }

        public DIAparameters DIAparameters { get; set; }
        public List<Peak>[] Ms1PeakTable { get; set; }
        public List<Peak>[] Ms2PeakTable { get; set; }
        public Dictionary<(double min, double max), List<PeakCurve>> Ms1PeakCurves { get; set; }
        public Dictionary<(double min, double max), List<PeakCurve>> Ms2PeakCurves { get; set; }
        public MsDataFile MyMSDataFile { get; set; }
        public CommonParameters CommonParameters { get; set; }
        public List<Ms2ScanWithSpecificMass> PseudoMS2Scans { get; set; }
        public Dictionary<(double min, double max), List<MsDataScan>> DIAScanWindowMap { get; set; }
        public List<PrecursorFragmentsGroup> PFgroups { get; set; }
        public double CycleTime { get; set; }
        public List<Ms2ScanWithSpecificMass> PseudoMs2WithPre { get; set; }

        public void GetPseudoMS2Scans()
        {
            Ms1PeakIndexing();
            ConstructMs2Group();
            GetMs1PeakCurves();
            GetMs2PeakCurves2();
            PrecursorFragmentPairing();
            PFgroupFilter();
            ConstructNewMs2Scans2();
        }

        public void Ms1PeakIndexing()
        {
            var ms1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            var allMs1Peaks = Peak.GetAllPeaks(ms1Scans, DIAparameters.PeakSearchBinSize);
            Ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, DIAparameters.PeakSearchBinSize);
        }
        public void ConstructMs2Group()
        {
            DIAScanWindowMap = new Dictionary<(double min, double max), List<MsDataScan>>();
            var ms2Scans = MyMSDataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            foreach (var ms2 in ms2Scans)
            {
                (double min, double max) range = new(ms2.IsolationRange.Minimum, ms2.IsolationRange.Maximum);
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
        }

        public void GetMs1PeakCurves()
        {
            var allMs1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            Ms1PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            int index = 1;
            foreach (var ms1window in DIAScanWindowMap.Keys)
            {
                Ms1PeakCurves[ms1window] = new List<PeakCurve>();
                var ms1Range = new MzRange(ms1window.min, ms1window.max);
                var allPrecursors = new List<DeconvolutedMass>();
                for (int i = 0; i < allMs1Scans.Length; i++)
                {
                    var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], CommonParameters.PrecursorDeconvolutionParameters, ms1Range);
                    foreach (var envelope in envelopes)
                    {
                        var charge = envelope.Charge;
                        double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                        double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                        var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                            allMs1Scans[i].OneBasedScanNumber, i);
                        allPrecursors.Add(precursor);
                    }
                }
                allPrecursors = allPrecursors.OrderByDescending(p => p.HighestPeakIntensity).ToList();
                foreach (var precursor in allPrecursors)
                {
                    if(precursor.HighestPeakIntensity < DIAparameters.PrecursorIntensityCutOff)
                    {
                        continue;
                    }
                    var peak = PeakCurve.GetPeakFromScan(precursor.HighestPeakMz, Ms1PeakTable, precursor.ZeroBasedScanIndex, new PpmTolerance(0),
                        DIAparameters.PeakSearchBinSize);
                    if (peak.PeakCurve == null)
                    {
                        var newPeakCurve = PeakCurve.FindPeakCurve(peak, Ms1PeakTable, allMs1Scans, null, DIAparameters.MaxNumMissedScan,
                        DIAparameters.Ms1PeakFindingTolerance, DIAparameters.PeakSearchBinSize, DIAparameters.MaxRTRange);
                        newPeakCurve.MonoisotopicMass = precursor.MonoisotopicMass;
                        newPeakCurve.Charge = precursor.Charge;
                        if(newPeakCurve.Peaks.Count > 4)
                        {
                            Ms1PeakCurves[ms1window].Add(newPeakCurve);
                            newPeakCurve.Index = index;
                            index++;
                        }
                    }
                }
            }
        }

        public void PrecursorFilter()
        {

        }

        public void GetMs2PeakCurves()
        {
            var ms2PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                MzRange range = ms2scans[0].IsolationRange;
                var allMs2Peaks = Peak.GetAllPeaks(ms2scans, DIAparameters.PeakSearchBinSize);
                var rankedMs2Peaks = allMs2Peaks.OrderByDescending(p => p.Intensity).ToList();
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, DIAparameters.PeakSearchBinSize);
                ms2PeakCurves[ms2Group.Key] = new List<PeakCurve>();
                foreach (var peak in rankedMs2Peaks)
                {
                    if (peak.PeakCurve == null)
                    {
                        var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms2PeakTable, ms2scans, ms2scans[0].IsolationRange,
                            DIAparameters.MaxNumMissedScan, DIAparameters.Ms2PeakFindingTolerance, DIAparameters.PeakSearchBinSize, DIAparameters.MaxRTRange);
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            ms2PeakCurves[ms2Group.Key].Add(newPeakCurve);
                        }
                    }
                }
            }
            Ms2PeakCurves = ms2PeakCurves;
        }

        public void GetMs2PeakCurves2()
        {
            var ms2PeakCurves = new Dictionary<(double min, double max), List<PeakCurve>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                var allMasses = new List<DeconvolutedMass>();
                for (int i = 0; i < ms2scans.Length; i++)
                {
                    var envelopes = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scans[i], CommonParameters);

                    foreach (var envelope in envelopes)
                    {
                        var charge = envelope.Charge;
                        double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                        double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                        var mass = new DeconvolutedMass(envelope, charge, ms2scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                            ms2scans[i].OneBasedScanNumber, i);
                        allMasses.Add(mass);
                    }
                }
                allMasses = allMasses.OrderByDescending(p => p.HighestPeakIntensity).ToList();

                MzRange range = ms2scans[0].IsolationRange;
                var allMs2Peaks = Peak.GetAllPeaks(ms2scans, DIAparameters.PeakSearchBinSize);
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, DIAparameters.PeakSearchBinSize);
                ms2PeakCurves[ms2Group.Key] = new List<PeakCurve>();
                foreach (var mass in allMasses)
                {
                    var peak = PeakCurve.GetPeakFromScan(mass.HighestPeakMz, ms2PeakTable, mass.ZeroBasedScanIndex, new PpmTolerance(0),
                        DIAparameters.PeakSearchBinSize);
                    if (peak.PeakCurve == null)
                    {
                        var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms2PeakTable, ms2scans, null, DIAparameters.MaxNumMissedScan,
                        DIAparameters.Ms1PeakFindingTolerance, DIAparameters.PeakSearchBinSize, DIAparameters.MaxRTRange);
                        newPeakCurve.MonoisotopicMass = mass.MonoisotopicMass;
                        newPeakCurve.Charge = mass.Charge;
                        newPeakCurve.Envelope = mass.Envelope;
                        if (newPeakCurve.Peaks.Count > 4)
                        {
                            ms2PeakCurves[ms2Group.Key].Add(newPeakCurve);
                        }
                    }
                }
            }
            Ms2PeakCurves = ms2PeakCurves;
        }

        public void PrecursorFragmentPairing()
        {
            PFgroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2group in Ms2PeakCurves)
            {
                var precursorsInRange = Ms1PeakCurves[ms2group.Key].ToArray();

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
            foreach(var ms2curve in allMatchedMs2)
            {
                ms2curve.GetPrecursorRanks();
            }
            foreach(var group in PFgroups)
            {
                group.PFpairs = group.PFpairs.Where(pf => pf.PrecursorRank <= DIAparameters.PrecursorRankCutOff).ToList();
                group.GetNumberOfHighCorrFragments(DIAparameters);
            }
            PFgroups = PFgroups.Where(pf => pf.PFpairs.Count > 0  && pf.NumHighCorrFragments >= DIAparameters.NumHighCorrFragments).ToList();
        }

        public void ConstructNewMs2Scans()
        {
            PseudoMs2WithPre = new List<Ms2ScanWithSpecificMass>();
            foreach (var pfGroup in PFgroups)
            {
                var mzs = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedMz).ToArray();
                var intensities = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedIntensity).ToArray();
                var spectrum = new MzSpectrum(mzs, intensities, false);
                var newMs2Scan = new MsDataScan(spectrum, pfGroup.Index, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                            MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
                var neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, CommonParameters);
                var charge = pfGroup.PrecursorPeakCurve.Charge;
                var highestPeakMz = pfGroup.PrecursorPeakCurve.AveragedMz;
                Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, highestPeakMz, charge
                    , MyMSDataFile.FilePath, CommonParameters, neutralExperimentalFragments);
                PseudoMs2WithPre.Add(scanWithprecursor);
            }
        }

        public void ConstructNewMs2Scans2()
        {
            PseudoMs2WithPre = new List<Ms2ScanWithSpecificMass>();
            foreach (var pfGroup in PFgroups)
            {
                //var peaks = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.Envelope.Peaks).SelectMany(p => p).OrderBy(p => p.mz).ToList();
                var mzs = new double[] { 1 };
                var intensities = new double[] { pfGroup.PFpairs.Sum(pf => pf.FragmentPeakCurve.Envelope.TotalIntensity) };
                var spectrum = new MzSpectrum(mzs, intensities, false);
                var newMs2Scan = new MsDataScan(spectrum, pfGroup.Index, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                            MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
                var neutralExperimentalFragments = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.Envelope).ToArray();
                var charge = pfGroup.PrecursorPeakCurve.Charge;
                var highestPeakMz = pfGroup.PrecursorPeakCurve.AveragedMz;
                Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, highestPeakMz, charge
                    , MyMSDataFile.FilePath, CommonParameters, neutralExperimentalFragments);
                PseudoMs2WithPre.Add(scanWithprecursor);
            }
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
                        var overlap = PeakCurve.CalculateRTOverlapRatio(precursor, ms2curve);
                        if (overlap > DIAparameters.OverlapRatioCutOff)
                        {
                            double corr = PeakCurve.CalculateCorr_spline(precursor, ms2curve, "cubic", 0.005);
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

