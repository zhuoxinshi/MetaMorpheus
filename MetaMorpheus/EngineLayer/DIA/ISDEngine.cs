using BayesianEstimation;
using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra.Storage;
using MzLibUtil;
using OpenMcdf.Extensions.OLEProperties;
using Plotly.NET.CSharp;
using Proteomics.AminoAcidPolymer;
using Readers;
using SpectralAveraging;
using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Reflection;
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
        public Dictionary<double, List<MsDataScan>> ISDScanVoltageMap { get; set; }
        public List<PrecursorFragmentsGroup> PFgroups { get; set; }
        public double CycleTime { get; set; }
        public List<Ms2ScanWithSpecificMass> PseudoMs2WithPre { get; set; }
        public Dictionary<int, int> Ms1ScanIndexMap { get; set; }

        public void GetPseudoMS2Scans()
        {
            //Ms1PeakIndexing();
            //ConstructMs2Group();

            //if (DIAparameters.AverageMs2Scans)
            //{
            //    AverageMs2Scans();
            //}
            ////AverageMs2Scans();
            //GetMs1PeakCurves();
            //GetMs2PeakCurves();
            ////GetMs1PeakEnvelopeCurve();
            ////GetMs2PeakEnvelopeCurve();
            //PrecursorFragmentPairing();
            ////PFgroupFilter();
            //ConstructNewMs2Scans();
            ////ConstructNewMs2Scans_peakEnvelopeCurve();

            UseExternalFeatures();
        }

        public void Ms1PeakIndexing()
        {
            var ms1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            var allMs1Peaks = Peak.GetAllPeaks(ms1Scans, DIAparameters.PeakSearchBinSize);
            //Ms1ScanIndexMap = scanIndexMap;
            Ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, DIAparameters.PeakSearchBinSize);
        }
        public void ConstructMs2Group()
        {
            ISDScanVoltageMap = new Dictionary<double, List<MsDataScan>>();
            var ms2Scans = MyMSDataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            string pattern = $@"sid=(\d+)";
            foreach (var ms2 in ms2Scans)
            {
                var match = Regex.Match(ms2.ScanFilter, pattern);
                double voltage = double.Parse(match.Groups[1].Value);
                if (!ISDScanVoltageMap.ContainsKey(voltage))
                {
                    ISDScanVoltageMap[voltage] = new List<MsDataScan>();
                    ISDScanVoltageMap[voltage].Add(ms2);
                }
                else
                {
                    ISDScanVoltageMap[voltage].Add(ms2);
                }
            }
        }

        public void AverageMs2Scans()
        {
            var averageParam = new SpectralAveragingParameters()
            {
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
                NumberOfScansToAverage = 5,
                ScanOverlap = 4,
                NormalizationType = NormalizationType.RelativeToTics,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly,
            };

            foreach (var ms2Group in ISDScanVoltageMap)
            {
                var origScans = ms2Group.Value.ToArray();
                var scansForAveraging = new List<MsDataScan>();
                for (int i = 0; i < origScans.Length; i++)
                {
                    var newScan = new MsDataScan(origScans[i].MassSpectrum, origScans[i].OneBasedScanNumber, 1, origScans[i].IsCentroid,
                        origScans[i].Polarity, origScans[i].RetentionTime, origScans[i].ScanWindowRange, origScans[i].ScanFilter, origScans[i].MzAnalyzer,
                        origScans[i].TotalIonCurrent, origScans[i].InjectionTime, origScans[i].NoiseData, origScans[i].NativeId,
                        origScans[i].SelectedIonMZ, origScans[i].SelectedIonChargeStateGuess, origScans[i].SelectedIonIntensity,
                        origScans[i].IsolationMz, origScans[i].IsolationWidth, origScans[i].DissociationType, origScans[i].OneBasedPrecursorScanNumber,
                        origScans[i].SelectedIonMonoisotopicGuessMz, origScans[i].HcdEnergy, origScans[i].ScanDescription);
                    scansForAveraging.Add(newScan);
                }
                var scans = new MsDataScan[origScans.Length];
                var averagedMS2 = SpectraFileAveraging.AverageSpectraFile(scansForAveraging, averageParam).ToList();
                for (int i = 0; i < origScans.Length; i++)
                {
                    var newScan = new MsDataScan(averagedMS2[i].MassSpectrum, origScans[i].OneBasedScanNumber, origScans[i].MsnOrder, origScans[i].IsCentroid,
                        origScans[i].Polarity, origScans[i].RetentionTime, origScans[i].ScanWindowRange, origScans[i].ScanFilter, origScans[i].MzAnalyzer,
                        origScans[i].TotalIonCurrent, origScans[i].InjectionTime, origScans[i].NoiseData, origScans[i].NativeId,
                        origScans[i].SelectedIonMZ, origScans[i].SelectedIonChargeStateGuess, origScans[i].SelectedIonIntensity,
                        origScans[i].IsolationMz, origScans[i].IsolationWidth, origScans[i].DissociationType, origScans[i].OneBasedPrecursorScanNumber,
                        origScans[i].SelectedIonMonoisotopicGuessMz, origScans[i].HcdEnergy, origScans[i].ScanDescription);
                    scans[i] = newScan;
                }
                ms2Group.Value.Clear();
                ms2Group.Value.AddRange(scans);
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
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], CommonParameters.PrecursorDeconvolutionParameters).OrderByDescending(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.Charge < 6  || envelope.Peaks.Count < 3 || envelope.MonoisotopicMass < DIAparameters.MinMass)
                    {
                        continue;
                    }
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, 1, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
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
                            newPeakCurve.GetScanCycleSmoothedData(DIAparameters.ScanCycleSplineTimeInterval);
                            newPeakCurve.Index = index;
                            index++;
                            
                            //debug
                            if(Math.Abs(newPeakCurve.AveragedMz - 864.60) < 0.01 && newPeakCurve.Charge == 16)
                            {
                                newPeakCurve.VisualizeRaw("line");
                            }
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
            var p = Ms1PeakCurves.Where(p => Math.Abs(p.AveragedMz - 864.60) < 0.01 && p.Charge == 16).First();
            p.VisualizeRaw("line");
        }

        public void GetMs2PeakCurves()
        {
            var ms2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach (var ms2Group in ISDScanVoltageMap)
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
                                newPeakCurve.GetScanCycleSmoothedData(DIAparameters.ScanCycleSplineTimeInterval);
                            }
                        }
                    }
                }
            }
            Ms2PeakCurves = ms2PeakCurves;

            //debug
            //var testMs2PeakCurve = Ms2PeakCurves.Values.SelectMany(v => v).ToList();
            //for (int i = 0; i<100; i++)
            //{
            //    Random rnd = new Random();
            //    int r = rnd.Next(testMs2PeakCurve.Count - 1);
            //    var pc = testMs2PeakCurve[r];
            //    pc.VisualizeBspline(out List<float> rtSeq).Show();
            //    pc.VisualizePeakRegions();
            //}
        }

        public void GetMs2PeakCurves_Decon()
        {
            var ms2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach (var ms2Group in ISDScanVoltageMap)
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
                        var mass = new DeconvolutedMass(envelope, charge, ms2scans[i].RetentionTime, 2, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                            ms2scans[i].OneBasedScanNumber, i);
                        allMasses.Add(mass);
                    }
                }
                allMasses = allMasses.OrderByDescending(p => p.HighestPeakIntensity).ToList();

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
                        DIAparameters.Ms1PeakFindingTolerance, DIAparameters.PeakSearchBinSize, DIAparameters.MaxRTRangeMS1);
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
            //debug
            //var testMs2PeakCurve = Ms2PeakCurves.Values.SelectMany(v => v).ToList();
            //for (int i = 0; i<100; i++)
            //{
            //    Random rnd = new Random();
            //    int r = rnd.Next(testMs2PeakCurve.Count - 1);
            //    var pc = testMs2PeakCurve[r];
            //    pc.VisualizeBspline(out List<float> rtSeq).Show();
            //    pc.VisualizePeakRegions();
            //}
        }


        public void PrecursorFragmentPairing( )
        {
            PFgroups = new List<PrecursorFragmentsGroup>();
            //foreach (var ms2group in Ms2PeakCurves)
            //{
                var precursorsInRange = Ms1PeakCurves.ToArray();
                var allMs2Curves = Ms2PeakCurves.Values.SelectMany(v => v).ToList();

            Parallel.ForEach(Partitioner.Create(0, precursorsInRange.Length), new ParallelOptions { MaxDegreeOfParallelism = 18 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursorsInRange[i];
                        var preFragGroup = GroupPrecursorFragments_scanCycle(precursor, allMs2Curves, DIAparameters);

                        if (preFragGroup != null)
                        {
                            lock (PFgroups)
                            {
                                PFgroups.Add(preFragGroup);
                            }
                        }
                    }
                });
            //}
            //debug
            //for (int i = 0; i < 100; i++)
            //{
            //    Random rnd = new Random();
            //    int r = rnd.Next(PFgroups.Count - 1);
            //    var group = PFgroups[r];
            //    group.Visualize();
            //}
            var rankedPFgroups = PFgroups.OrderByDescending(pf => pf.PFpairs.Count).ToList();
            //var groups = PFgroups.OrderBy(pf => pf.PrecursorPeakCurve.MonoisotopicMass).ToList();
            //TODO?
            //Check the fragment correlations within each pfgroup for filtering
        }

        public static PrecursorFragmentsGroup GroupPrecursorFragments(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters DIAparameters)
        {
            var preFragGroup = new PrecursorFragmentsGroup(precursor);
            foreach (var ms2curve in ms2curves)
            {
                //if (ms2curve.ApexIntensity > precursor.ApexIntensity)
                //{
                //    continue;
                //}
                if (ms2curve.ApexRT >= precursor.StartRT && ms2curve.ApexRT <= precursor.EndRT)
                {
                    if (Math.Abs(ms2curve.ApexRT - precursor.ApexRT) <= DIAparameters.ApexRtTolerance)
                    {
                        var overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursor, ms2curve);
                        if (overlap > DIAparameters.OverlapRatioCutOff)
                        {
                            //double corr = PrecursorFragmentPair.CalculateCorr_spline_scanCycle(precursor, ms2curve, DIAparameters.ScanCycleSplineTimeInterval);
                            double corr = PrecursorFragmentPair.CalculatePeakCurveCorr(precursor, ms2curve);
                            if (corr > DIAparameters.CorrelationCutOff)
                            {
                                var PFpair = new PrecursorFragmentPair(precursor, ms2curve, corr);
                                preFragGroup.PFpairs.Add(PFpair);
                            }
                        }
                    }
                }
            }
            //if (preFragGroup.PFpairs.Count > DIAparameters.FragmentRankCutOff)
            //{
            //    var filtered = preFragGroup.PFpairs.OrderByDescending(pair => pair.Correlation).Take(DIAparameters.FragmentRankCutOff);
            //    preFragGroup.PFpairs = filtered.ToList();
            //}
            if (preFragGroup.PFpairs.Count > 10)
            {
                preFragGroup.PFpairs = preFragGroup.PFpairs.OrderBy(pair => pair.FragmentPeakCurve.AveragedMz).ToList();
                return preFragGroup;
            }
            else
            {
                return null;
            }
        }

        public void ConstructNewMs2Scans()
        {
            PseudoMs2WithPre = new List<Ms2ScanWithSpecificMass>();
            int oneBasedScanNumber = 1;
            foreach (var pfGroup in PFgroups)
            {
                var mzs = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedMz).ToArray();
                var intensities = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.AveragedIntensity).ToArray();
                var spectrum = new MzSpectrum(mzs, intensities, false);
                var newMs2Scan = new MsDataScan(spectrum, oneBasedScanNumber, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                            MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: pfGroup.PrecursorPeakCurve.Index);
                oneBasedScanNumber++;
                var neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, CommonParameters);
                var charge = pfGroup.PrecursorPeakCurve.Charge;
                var highestPeakMz = pfGroup.PrecursorPeakCurve.AveragedMz;
                Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, highestPeakMz, charge
                    , MyMSDataFile.FilePath, CommonParameters, neutralExperimentalFragments);
                PseudoMs2WithPre.Add(scanWithprecursor);
            }
        }

        public void ConstructNewMs2Scans_Decon()
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

        public void GetPrecursorEnvelopeXIC()
        {
            var allMs1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            Ms1PeakCurves = new List<PeakCurve>();
            int index = 1;
            var allPrecursors = new List<DeconvolutedMass>();
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], CommonParameters.PrecursorDeconvolutionParameters).OrderByDescending(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.MonoisotopicMass < DIAparameters.MinMass || envelope.Charge < 5)
                    {
                        continue;
                    }
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, 1, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
                    allPrecursors.Add(precursor);
                }
            }
            allPrecursors = allPrecursors.OrderByDescending(p => p.Envelope.TotalIntensity).ToList();
            //var precursorGroups = allPrecursors.GroupBy(p => new { mass = Math.Round(p.MonoisotopicMass, 0), p.Charge}).ToList();
            foreach (var precursor in allPrecursors)
            {
                var targetList = precursor.Envelope.Peaks.Select(p => p.mz).ToList();
                var targetPeaks = new List<Peak>();
                foreach (var targetMz in targetList)
                {
                    var peak = PeakCurve.GetPeakFromScan(targetMz, Ms1PeakTable, precursor.ZeroBasedScanIndex, new PpmTolerance(0), DIAparameters.PeakSearchBinSize);
                    targetPeaks.Add(peak);
                }
                var isotopePeakList = PrecursorCluster.FindIsotopePeakCurve(targetPeaks, Ms1PeakTable, precursor.ZeroBasedScanIndex, allMs1Scans, DIAparameters.Ms1PeakFindingTolerance,
                    DIAparameters.PeakSearchBinSize, DIAparameters.MaxNumMissedScan);
                var fakePeaks = new List<Peak>();
                foreach (var envelopePeaks in isotopePeakList)
                {
                    var mz = envelopePeaks.OrderByDescending(p => p.Intensity).First().Mz;
                    var rt = envelopePeaks.First().RetentionTime;
                    var intensity = envelopePeaks.Sum(p => p.Intensity);
                    var fakePeak = new Peak(mz, rt, intensity, 1, envelopePeaks.First().ScanNumber, envelopePeaks.First().ZeroBasedScanIndex);
                    fakePeaks.Add(fakePeak);
                }
                var newPC = new PeakCurve(fakePeaks, 1, null, precursor.MonoisotopicMass, precursor.Charge, startMz: targetList.Min(), endMz: targetList.Max());
                Ms1PeakCurves.Add(newPC);
            }
        }

        public void GetPrecursorEnvelopeXIC2()
        {
            var allMs1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            Ms1PeakCurves = new List<PeakCurve>();
            var allPrecursors = new List<DeconvolutedMass>[allMs1Scans.Length];
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                allPrecursors[i] = new List<DeconvolutedMass>();
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], CommonParameters.PrecursorDeconvolutionParameters).OrderByDescending(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.MonoisotopicMass < DIAparameters.MinMass || envelope.Charge < 5)
                    {
                        continue;
                    }
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, 1,highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
                    allPrecursors[i].Add(precursor);
                }
                allPrecursors[i] = allPrecursors[i].OrderBy(p => p.MonoisotopicMass).ToList();
            }
            var sortedAllPrecursors = allPrecursors.SelectMany(p => p).OrderByDescending(p => p.Envelope.TotalIntensity).ToList();

            int index = 1;
            foreach (var precursor in sortedAllPrecursors)
            {
                if (precursor.EnvelopeCurve == null)
                {
                    var newEC = EnvelopeCurve.GetEnvelopeCurve(precursor, allPrecursors, Ms1PeakTable ,DIAparameters);
                    newEC.GetFakePeakCurve();
                    if(newEC.FakePeakCurve.Peaks.Count > 4)
                    {
                        Ms1PeakCurves.Add(newEC.FakePeakCurve);
                        newEC.Index = index;
                        index++;
                    }
                }
            }
           
        }

        
        public void GetEnvelopeXICForMs2()
        {
            var ms2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach (var ms2Group in ISDScanVoltageMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                var allMasses = new List<DeconvolutedMass>[ms2scans.Length];
                for (int i = 0; i < ms2scans.Length; i++)
                {
                    allMasses[i] = new List<DeconvolutedMass>();
                    var envelopes = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scans[i], CommonParameters);

                    foreach (var envelope in envelopes)
                    {
                        var charge = envelope.Charge;
                        double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                        double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                        var mass = new DeconvolutedMass(envelope, charge, ms2scans[i].RetentionTime, 2, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                            ms2scans[i].OneBasedScanNumber, i);
                        allMasses[i].Add(mass);
                    }
                }
                var sortedMasses = allMasses.SelectMany(p => p).OrderByDescending(p => p.TotalIntensity).ToList();

                var allMs2Peaks = Peak.GetAllPeaks(ms2scans, DIAparameters.PeakSearchBinSize);
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, DIAparameters.PeakSearchBinSize);
                ms2PeakCurves[ms2Group.Key] = new List<PeakCurve>();

                foreach (var mass in sortedMasses)
                {
                    //var targetList = mass.Envelope.Peaks.Select(p => p.mz).ToList();
                    //var targetPeaks = new List<Peak>();
                    //foreach (var targetMz in targetList)
                    //{
                    //    var peak = PeakCurve.GetPeakFromScan(targetMz, ms2PeakTable, mass.ZeroBasedScanIndex, new PpmTolerance(0), DIAparameters.PeakSearchBinSize);
                    //    targetPeaks.Add(peak);
                    //}
                    //var isotopePeakList = PrecursorCluster.FindIsotopePeakCurve(targetPeaks, ms2PeakTable, mass.ZeroBasedScanIndex, ms2scans, DIAparameters.Ms1PeakFindingTolerance,
                    //    DIAparameters.PeakSearchBinSize, DIAparameters.MaxNumMissedScan);
                    //var fakePeaks = new List<Peak>();
                    //foreach (var envelopePeaks in isotopePeakList)
                    //{
                    //    var mz = envelopePeaks.OrderByDescending(p => p.Intensity).First().Mz;
                    //    var rt = envelopePeaks.First().RetentionTime;
                    //    var intensity = envelopePeaks.Sum(p => p.Intensity);
                    //    var fakePeak = new Peak(mz, rt, intensity, 1, envelopePeaks.First().ScanNumber, envelopePeaks.First().ZeroBasedScanIndex);
                    //    fakePeaks.Add(fakePeak);
                    //}
                    //var newPC = new PeakCurve(fakePeaks, 1, null, mass.MonoisotopicMass, mass.Charge, startMz: targetList.Min(), endMz: targetList.Max());
                    //ms2PeakCurves[ms2Group.Key].Add(newPC);

                    if (mass.EnvelopeCurve == null)
                    {
                        var newEC = EnvelopeCurve.GetEnvelopeCurve(mass, allMasses, ms2PeakTable, DIAparameters);
                        newEC.GetFakePeakCurve();
                        if (newEC.FakePeakCurve.Peaks.Count > 4)
                        {
                            ms2PeakCurves[ms2Group.Key].Add(newEC.FakePeakCurve);
                        }
                    }
                }
            }
            Ms2PeakCurves = ms2PeakCurves;
        }

        public void GetMs1PeakEnvelopeCurve()
        {
            var allMs1Scans = MyMSDataFile.GetMS1Scans().ToArray();
            Ms1PeakCurves = new List<PeakCurve>();
            var allPrecursors = new List<DeconvolutedMass>[allMs1Scans.Length];
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                allPrecursors[i] = new List<DeconvolutedMass>();
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], CommonParameters.PrecursorDeconvolutionParameters).OrderByDescending(E => E.MonoisotopicMass);
                foreach (var envelope in envelopes)
                {
                    if (envelope.MonoisotopicMass < DIAparameters.MinMass || envelope.Charge < 5)
                    {
                        continue;
                    }
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, 1, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
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
                var newPEC = PeakEnvelopeCurve.GetPeakEnvelopeCurve(targetPeaks, theorIntensityRatio, highestPeakIndex, Ms1PeakTable, allMs1Scans, precursor.ZeroBasedScanIndex, DIAparameters);
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
                    Ms1PeakCurves.Add(newPEC.FakePeakCurve);
                    newPEC.FakePeakCurve.GetScanCycleSmoothedData(DIAparameters.ScanCycleSplineTimeInterval);
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
        }

        public void GetMs2PeakEnvelopeCurve()
        {
            var ms2PeakCurves = new Dictionary<double, List<PeakCurve>>();
            foreach (var ms2Group in ISDScanVoltageMap)
            {
                var ms2scans = ms2Group.Value.ToArray();
                var allMasses = new List<DeconvolutedMass>[ms2scans.Length];
                for (int i = 0; i < ms2scans.Length; i++)
                {
                    allMasses[i] = new List<DeconvolutedMass>();
                    var envelopes = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scans[i], CommonParameters);

                    foreach (var envelope in envelopes)
                    {
                        var charge = envelope.Charge;
                        double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                        double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                        var mass = new DeconvolutedMass(envelope, charge, ms2scans[i].RetentionTime, 2, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                            ms2scans[i].OneBasedScanNumber, i);
                        allMasses[i].Add(mass);
                    }
                }
                var sortedMasses = allMasses.SelectMany(p => p).OrderByDescending(p => p.TotalIntensity).ToList();

                var allMs2Peaks = Peak.GetAllPeaks(ms2scans, DIAparameters.PeakSearchBinSize);
                var ms2PeakTable = Peak.GetPeakTable(allMs2Peaks, DIAparameters.PeakSearchBinSize);
                ms2PeakCurves[ms2Group.Key] = new List<PeakCurve>();

                foreach (var mass in sortedMasses)
                {
                    var targetPeaks = mass.Envelope.Peaks.OrderBy(p => p.mz).ToArray();
                    var theorIntensityRatio = targetPeaks.Select(p => p.intensity / targetPeaks.Sum(p => p.intensity)).ToArray();
                    var highestPeak = targetPeaks.OrderByDescending(p => p.intensity).FirstOrDefault();
                    var highestPeakIndex = Array.IndexOf(targetPeaks, highestPeak);
                    var newPEC = PeakEnvelopeCurve.GetPeakEnvelopeCurve(targetPeaks, theorIntensityRatio, highestPeakIndex, ms2PeakTable, ms2scans, mass.ZeroBasedScanIndex, DIAparameters);
                    if (newPEC != null && newPEC.PeakEnvelopes.Count > 4)
                    {
                        newPEC.MonoisotopicMass = mass.MonoisotopicMass;
                        newPEC.Charge = mass.Charge;
                        newPEC.MsLevel = 2;
                        newPEC.IsolationRange = null;
                        newPEC.GetFakePeakCurve();
                        ms2PeakCurves[ms2Group.Key].Add(newPEC.FakePeakCurve);
                        newPEC.FakePeakCurve.GetScanCycleSmoothedData(DIAparameters.ScanCycleSplineTimeInterval);
                    }
                }
            }
            Ms2PeakCurves = ms2PeakCurves;
            var ms2PC = Ms2PeakCurves.Values.SelectMany(v => v).ToList();
            //debug
            //for (int i = 0; i <100; i++)
            //{
            //    Random rnd = new Random();
            //    int r = rnd.Next(ms2PC.Count - 1);
            //    var pc = ms2PC[r];
            //    pc.VisualizeRaw("line");
            //}
        }

        public void ConstructNewMs2Scans_peakEnvelopeCurve()
        {
            PseudoMs2WithPre = new List<Ms2ScanWithSpecificMass>();
            int oneBasedScanNum = 1;
            foreach (var pfGroup in PFgroups)
            {
                //var peaks = pfGroup.PFpairs.Select(pf => pf.FragmentPeakCurve.Envelope.Peaks).SelectMany(p => p).OrderBy(p => p.mz).ToList();
                var mzs = new double[] { 1 };
                var intensities = new double[] { pfGroup.PFpairs.Sum(pf => pf.FragmentPeakCurve.TotalIntensity) };
                var spectrum = new MzSpectrum(mzs, intensities, false);
                var newMs2Scan = new MsDataScan(spectrum, oneBasedScanNum, 2, true, Polarity.Positive, pfGroup.PrecursorPeakCurve.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null,
                            MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: pfGroup.Index);
                oneBasedScanNum++;
                var neutralExperimentalFragments = new List<IsotopicEnvelope>();
                foreach (var pf in pfGroup.PFpairs)
                {
                    var envelope = new IsotopicEnvelope(new List<(double mz, double intensity)> { (pf.FragmentPeakCurve.AveragedMz, pf.FragmentPeakCurve.AveragedIntensity)}, 
                        pf.FragmentPeakCurve.MonoisotopicMass, pf.FragmentPeakCurve.Charge, pf.FragmentPeakCurve.TotalIntensity, 0, 0);
                    neutralExperimentalFragments.Add(envelope);
                }
                Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, pfGroup.PrecursorPeakCurve.AveragedMz, 
                    pfGroup.PrecursorPeakCurve.Charge, MyMSDataFile.FilePath, CommonParameters, neutralExperimentalFragments.ToArray());
                PseudoMs2WithPre.Add(scanWithprecursor);
            }
        }

        public static PrecursorFragmentsGroup GroupPrecursorFragments_scanCycle(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters DIAparameters)
        {
            var preFragGroup = new PrecursorFragmentsGroup(precursor);
            foreach (var ms2curve in ms2curves)
            {
                if (ms2curve.ApexScanCycle >= precursor.StartCycle && ms2curve.ApexScanCycle <= precursor.EndCycle)
                {
                    if (Math.Abs(ms2curve.ApexScanCycle - precursor.ApexScanCycle) <= DIAparameters.ApexCycleTolerance)
                    {
                        var overlap = PrecursorFragmentPair.CalculateRTOverlapRatio_scanCycle(precursor, ms2curve);
                        if (overlap > DIAparameters.OverlapRatioCutOff)
                        {
                            double corr = PrecursorFragmentPair.CalculateCorr_spline_scanCycle(precursor, ms2curve);
                            if (corr > DIAparameters.CorrelationCutOff)
                            {
                                var overlapAreaRatio = PrecursorFragmentPair.CalculateOverlapAreaRatio(precursor, ms2curve);
                                var PFpair = new PrecursorFragmentPair(precursor, ms2curve, corr);
                                //lock (ms2curve.PFpairs)
                                //{
                                //    ms2curve.PFpairs.Add(PFpair);
                                //}
                                preFragGroup.PFpairs.Add(PFpair);
                            }
                        }
                    }
                }
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

        public static PrecursorFragmentsGroup GroupPrecursorFragments_area(PeakCurve precursor, List<PeakCurve> ms2curves, DIAparameters DIAparameters)
        {
            var preFragGroup = new PrecursorFragmentsGroup(precursor);
            foreach (var ms2curve in ms2curves)
            {
                if (ms2curve.ApexScanCycle >= precursor.StartCycle && ms2curve.ApexScanCycle <= precursor.EndCycle)
                {
                    if (Math.Abs(ms2curve.ApexScanCycle - precursor.ApexScanCycle) <= DIAparameters.ApexCycleTolerance)
                    {
                        var overlap = PrecursorFragmentPair.CalculateOverlapAreaRatio(precursor, ms2curve);
                        if (overlap > DIAparameters.OverlapRatioCutOff)
                        {
                            var PFpair = new PrecursorFragmentPair(precursor, ms2curve, 0);
                                //lock (ms2curve.PFpairs)
                                //{
                                //    ms2curve.PFpairs.Add(PFpair);
                                //}
                            preFragGroup.PFpairs.Add(PFpair);
                        }
                    }
                }
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

        public void UseExternalFeatures()
        {
            var ms1FeatureFile = @"E:\toppic-windows-1.6.5\ISD\test1\id_08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_MS1_ms1.feature";
            var featureFile60 = @"E:\toppic-windows-1.6.5\ISD\test1\id_08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_60_ms1.feature";
            var featureFile80 = @"E:\toppic-windows-1.6.5\ISD\test1\id_08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_80_ms1.feature";
            var featureFile100 = @"E:\toppic-windows-1.6.5\ISD\test1\id_08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_100_ms1.feature";
            var allMs1Features = new Ms1FeatureFile(ms1FeatureFile);
            var all60Features = new Ms1FeatureFile(featureFile60);
            var all80Features = new Ms1FeatureFile(featureFile80);
            var all100Features = new Ms1FeatureFile(featureFile100);

            var ms1features = allMs1Features.Where(f => f.RetentionTimeBegin <= 32.69*60 && f.RetentionTimeEnd >= 28.01*60).ToList();
            var features60 = all60Features.Where(f => f.RetentionTimeBegin <= 32.69*60 && f.RetentionTimeEnd >= 28.01*60).ToList();
            var features80 = all80Features.Where(f => f.RetentionTimeBegin <= 32.69*60 && f.RetentionTimeEnd >= 28.01*60).ToList();
            var features100 = all100Features.Where(f => f.RetentionTimeBegin <= 32.69*60 && f.RetentionTimeEnd >= 28.01*60).ToList();
            var ms2features = new List<Ms1Feature> { features60, features80, features100 };

            PseudoMs2WithPre = new List<Ms2ScanWithSpecificMass>();
            foreach(var ms1 in ms1features)
            {
                var fragments = new List<Ms1Feature>();
                foreach(var ms2 in ms2features)
                {
                    if (OverlapRatio(ms1, ms2) > 0.3 && Math.Abs(ms1.RetentionTimeApex - ms2.RetentionTimeApex) < 0.3*60)
                    {
                        double overlap = PrecursorFragmentPair.CalculateOverlapRatioGaussian(ms1, ms2);
                        if (overlap > 0.2)
                        {
                            fragments.Add(ms2);
                        }
                    }
                }
                fragments = fragments.OrderBy(f => f.Mass).ToList();
                var mzs = new double[] { 1 };
                var intensities = new double[] { fragments.Sum(f => f.Intensity) };
                var spectrum = new MzSpectrum(mzs, intensities, false);
                var newMs2Scan = new MsDataScan(spectrum, ms1.Id, 2, true, Polarity.Positive, ms1.RetentionTimeApex, new MzRange(mzs.Min(), mzs.Max()), null,
                MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null);
                var neutralExperimentalFragments = new List<IsotopicEnvelope>();
                if (fragments.Count > 0)
                {
                    foreach (var ms2 in fragments)
                    {
                        var envelope = new IsotopicEnvelope(new List<(double, double)> { (ms2.Mass.ToMz(ms2.ChargeStateMax), ms2.Intensity)}, ms2.Mass, ms2.ChargeStateMax, 
                            ms2.Intensity, 0.1, 1);
                        neutralExperimentalFragments.Add(envelope);
                    }
                    var newMs2WithMass = new Ms2ScanWithSpecificMass(newMs2Scan, ms1.Mass.ToMz(ms1.ChargeStateMax), ms1.ChargeStateMax, MyMSDataFile.FilePath,
                    CommonParameters, neutralExperimentalFragments.ToArray());
                    PseudoMs2WithPre.Add(newMs2WithMass);
                }
            }
        }

        public static double OverlapRatio(Ms1Feature feature1, Ms1Feature feature2)
        {
            var start = Math.Min(feature1.RetentionTimeBegin, feature2.RetentionTimeBegin);
            var end = Math.Max(feature1.RetentionTimeEnd, feature2.RetentionTimeEnd);
            var overlapStart = Math.Max(feature1.RetentionTimeBegin, feature2.RetentionTimeBegin);
            var overlapEnd = Math.Min(feature1.RetentionTimeEnd, feature2.RetentionTimeEnd);
            double ratio = (overlapEnd - overlapStart) / (end - start);

            return ratio;
        }
    }

}
