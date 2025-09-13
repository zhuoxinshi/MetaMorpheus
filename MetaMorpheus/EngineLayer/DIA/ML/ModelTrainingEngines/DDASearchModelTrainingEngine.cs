using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FlashLFQ;
using Chemistry;
using System.Runtime.InteropServices.ObjectiveC;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.DIA
{
    public class DDASearchModelTrainingEngine : ModelTrainingEngine
    {
        public MsDataScan[] Ms1Scans { get; set; }
        public MsDataScan[] Ms2Scans { get; set; }
        public Dictionary<object, object> AllMs1PeakIndexingEngines { get; set; }
        public Dictionary<object, object> AllMs2PeakIndexingEngines { get; set; }
        public Dictionary<IIndexedPeak, ExtractedIonChromatogram> Ms1PeakXicDictionary { get; set; }
        public Dictionary<IIndexedPeak, ExtractedIonChromatogram> Ms2PeakXicDictionary { get; set; }
        public Dictionary<int, object> ScanNumberWindowMap { get; set; }
        public DDASearchModelTrainingEngine(MLbasedDIAparameters mlDIAparams, CommonParameters commonParameters, MsDataScan[] ms1Scans, MsDataScan[] ms2Scans, Dictionary<object, object> allMs1PeakIndexingEngines, Dictionary<object, object> allMs2PeakIndexingEngines, Dictionary<IIndexedPeak, ExtractedIonChromatogram> ms1PeakXicDictionary, Dictionary<IIndexedPeak, ExtractedIonChromatogram> ms2PeakXicDictionary, Dictionary<(double min, double max), List<MsDataScan>> diaScanWindowMap = null, Dictionary<double, List<MsDataScan>> isdScanVoltageMap = null) : base(mlDIAparams, commonParameters)
        {
            Ms1Scans = ms1Scans;
            Ms2Scans = ms2Scans;
            AllMs1PeakIndexingEngines = allMs1PeakIndexingEngines;
            AllMs2PeakIndexingEngines = allMs2PeakIndexingEngines;
            Ms1PeakXicDictionary = ms1PeakXicDictionary;
            Ms2PeakXicDictionary = ms2PeakXicDictionary;
            ScanNumberWindowMap = new Dictionary<int, object>();
            if (diaScanWindowMap != null)
            {
                foreach (var kvp in diaScanWindowMap)
                {
                    foreach (var scan in kvp.Value)
                    {
                        ScanNumberWindowMap[scan.OneBasedScanNumber] = kvp.Key;
                    }
                }
            }
            else if (isdScanVoltageMap != null  )
            {
                foreach (var kvp in isdScanVoltageMap)
                {
                    foreach (var scan in kvp.Value)
                    {
                        ScanNumberWindowMap[scan.OneBasedScanNumber] = kvp.Key;
                    }
                }
            }
        }

        public override void GeneratePseudoSearchScans()
        {
            PseudoSearchMs2Scans = GetMs2Scans(Ms1Scans, Ms2Scans, CommonParameters).ToArray();
        }

        public override IEnumerable<PfPairTrainingSample> GetTrainingSamples()
        {
            GeneratePseudoSearchScans();
            Psms = RunClassicSearch(PseudoSearchMs2Scans).Where(p => p.Score >= MlDIAparams.PsmScoreCutOff).ToArray();
            var allSamples = new List<PfPairTrainingSample>();
            var groupedPsms = Psms.GroupBy(p => new { p.FullSequence, p.ScanPrecursorCharge }).Select(g => g.OrderByDescending(p => p.PrecursorScanIntensity));

            if (MlDIAparams.AnalysisType == AnalysisType.MLbased_bottomUp)
            {
                foreach (var psm in groupedPsms)
                {
                    //var samples = GetTrainingSamplesFromDDASearchPsm_bottomUp(psm, PseudoSearchMs2Scans.FirstOrDefault(s => s.OneBasedScanNumber == psm.ScanNumber && psm.ScanPrecursorCharge == s.PrecursorCharge && psm.ScanPrecursorMass == s.PrecursorMass));
                    //if (samples.Count() > 0) allSamples.AddRange(samples);
                }
            } 
            else if (MlDIAparams.AnalysisType == AnalysisType.MLbased_topDown)
            {
                foreach (var psms in groupedPsms)
                {
                    //var samples = GetTrainingSamplesFromDDASearchPsm_topDown(psm, PseudoSearchMs2Scans.FirstOrDefault(s => s.OneBasedScanNumber == psm.ScanNumber && psm.ScanPrecursorCharge == s.PrecursorCharge && psm.ScanPrecursorMass == s.PrecursorMass));
                    var samples = GetTrainingSamplesFromDDASearchPsm_topDown(psms);
                    if (samples.Count() > 0) allSamples.AddRange(samples);
                }
            }
            return allSamples;
        }

        public IEnumerable<PfPairTrainingSample> GetTrainingSamplesFromDDASearchPsm_bottomUp(SpectralMatch psm, Ms2ScanWithSpecificMass ms2WithPrecursor)
        {
            var allSamples = new List<PfPairTrainingSample>();
            var key = ScanNumberWindowMap[psm.ScanNumber];
            var ms1PeakIndexingEngine = AllMs1PeakIndexingEngines[key] as PeakIndexingEngine;
            var ms2PeakIndexingEngine = AllMs2PeakIndexingEngines[key] as PeakIndexingEngine;

            //find precursor xic
            int numberOfScansPerCycle = 3;
            int zeroBasedPrecursorScanIndex = (psm.PrecursorScanNumber.Value - 1) / numberOfScansPerCycle;
            int zeroBasedMs2ScanIndex = (psm.ScanNumber - 1) / numberOfScansPerCycle;
            var precursorPeak = ms1PeakIndexingEngine.GetIndexedPeak(psm.ScanPrecursorHighestIsotopeMz, zeroBasedPrecursorScanIndex, MlDIAparams.Ms1XicConstructor.PeakFindingTolerance);
            if (precursorPeak != null && Ms1PeakXicDictionary.ContainsKey(precursorPeak))
            {
                var precursorXic = Ms1PeakXicDictionary[precursorPeak];
                if (precursorXic == null) return allSamples;

                //find all fragment xics
                //create training samples with positive and negative labels
                var visitedXics = new HashSet<ExtractedIonChromatogram>();
                foreach (var ion in psm.MatchedFragmentIons)
                {
                    double minMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMinimumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                    double maxMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMaximumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                    var masses = ms2WithPrecursor.GetClosestExperimentalIsotopicEnvelopeList(minMass, maxMass);
                    var targetMass = masses.FirstOrDefault(m => m.Charge == ion.Charge);
                    if (targetMass != null)
                    {
                        foreach (var peak in targetMass.Peaks)
                        {
                            var fragmentPeak = ms2PeakIndexingEngine.GetIndexedPeak(peak.mz, zeroBasedMs2ScanIndex, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                            if (fragmentPeak != null && Ms2PeakXicDictionary.ContainsKey(fragmentPeak))
                            {
                                var fragmentXic = Ms2PeakXicDictionary[fragmentPeak];
                                if (fragmentXic == null) continue;
                                var newSample = new PfPairTrainingSample(precursorXic, fragmentXic, true, psm);
                                allSamples.Add(newSample);
                                visitedXics.Add(fragmentXic);
                            }
                        }
                    }
                }

                foreach (var mz in ms2WithPrecursor.TheScan.MassSpectrum.XArray)
                {
                    var peak = ms2PeakIndexingEngine.GetIndexedPeak(mz, zeroBasedMs2ScanIndex, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                    if (peak != null && Ms2PeakXicDictionary.ContainsKey(peak))
                    {
                        var fragmentXic = Ms2PeakXicDictionary[peak];
                        if (!visitedXics.Contains(fragmentXic) && fragmentXic != null)
                        {
                            var newSample = new PfPairTrainingSample(precursorXic, fragmentXic, false, psm);
                            allSamples.Add(newSample);
                            visitedXics.Add(fragmentXic);
                        }
                    }
                }
            }
            return allSamples;
        }

        public IEnumerable<PfPairTrainingSample> GetTrainingSamplesFromDDASearchPsm_topDown(SpectralMatch psm, Ms2ScanWithSpecificMass ms2WithPrecursor)
        {
            var allSamples = new List<PfPairTrainingSample>();
            var key = ScanNumberWindowMap[psm.ScanNumber];
            var ms1PeakIndexingEngine = AllMs1PeakIndexingEngines[key] as MassIndexingEngine;
            var ms2PeakIndexingEngine = AllMs2PeakIndexingEngines[key] as MassIndexingEngine;

            //find precursor xic
            int numberOfScansPerCycle = ScanNumberWindowMap.Values.Distinct().Count() + 1;
            int zeroBasedPrecursorScanIndex = (psm.PrecursorScanNumber.Value - 1) / numberOfScansPerCycle;
            int zeroBasedMs2ScanIndex = (psm.ScanNumber - 1) / numberOfScansPerCycle;
            var precursorPeak = ms1PeakIndexingEngine.GetIndexedPeak(psm.ScanPrecursorMass, zeroBasedPrecursorScanIndex, MlDIAparams.Ms1XicConstructor.PeakFindingTolerance);
            if (precursorPeak != null && Ms1PeakXicDictionary.ContainsKey(precursorPeak))
            {
                var precursorXic = Ms1PeakXicDictionary[precursorPeak];
                if (precursorXic == null) return allSamples;

                //find all fragment xics
                //create training samples with positive and negative labels
                var visitedXics = new HashSet<ExtractedIonChromatogram>();
                foreach (var ion in psm.MatchedFragmentIons)
                {
                    double minMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMinimumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                    double maxMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMaximumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                    var masses = ms2WithPrecursor.GetClosestExperimentalIsotopicEnvelopeList(minMass, maxMass);
                    var targetMass = masses.FirstOrDefault(m => m.Charge == ion.Charge);
                    if (targetMass != null)
                    {
                        var fragmentPeak = ms2PeakIndexingEngine.GetIndexedPeak(targetMass.MonoisotopicMass, zeroBasedMs2ScanIndex, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                        if (fragmentPeak != null && Ms2PeakXicDictionary.ContainsKey(fragmentPeak))
                        {
                            var fragmentXic = Ms2PeakXicDictionary[fragmentPeak];
                            if (fragmentXic == null) continue;
                            var newSample = new PfPairTrainingSample(precursorXic, fragmentXic, true, psm);
                            allSamples.Add(newSample);
                            visitedXics.Add(fragmentXic);
                        }
                    }
                }

                foreach (var mass in ms2WithPrecursor.ExperimentalFragments)
                {
                    var peak = ms2PeakIndexingEngine.GetIndexedPeak(mass.MonoisotopicMass, zeroBasedMs2ScanIndex, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                    if (peak != null && Ms2PeakXicDictionary.ContainsKey(peak))
                    {
                        var fragmentXic = Ms2PeakXicDictionary[peak];
                        if (!visitedXics.Contains(fragmentXic) && fragmentXic != null)
                        {
                            var newSample = new PfPairTrainingSample(precursorXic, fragmentXic, false, psm);
                            allSamples.Add(newSample);
                            visitedXics.Add(fragmentXic);
                        }
                    }
                }
            }
            return allSamples;
        }

        public IEnumerable<PfPairTrainingSample> GetTrainingSamplesFromDDASearchPsm_topDown(IEnumerable<SpectralMatch> psms)
        {
            var allSamples = new List<PfPairTrainingSample>();
            var key = ScanNumberWindowMap[psms.First().ScanNumber];
            var ms1PeakIndexingEngine = AllMs1PeakIndexingEngines[key] as MassIndexingEngine;
            var ms2PeakIndexingEngine = AllMs2PeakIndexingEngines[key] as MassIndexingEngine;

            //find precursor xic
            int numberOfScansPerCycle = ScanNumberWindowMap.Values.Distinct().Count() + 1;
            int zeroBasedPrecursorScanIndex = (psms.First().PrecursorScanNumber.Value - 1) / numberOfScansPerCycle;
            int zeroBasedMs2ScanIndex = (psms.First().ScanNumber - 1) / numberOfScansPerCycle;
            var precursorPeak = ms1PeakIndexingEngine.GetIndexedPeak(psms.First().ScanPrecursorMass, zeroBasedPrecursorScanIndex, MlDIAparams.Ms1XicConstructor.PeakFindingTolerance);

            //find all fragment xics
            //create training samples with positive and negative labels
            var visitedXics = new HashSet<ExtractedIonChromatogram>();
            if (precursorPeak != null && Ms1PeakXicDictionary.ContainsKey(precursorPeak))
            {
                var precursorXic = Ms1PeakXicDictionary[precursorPeak];
                if (precursorXic == null) return allSamples;

                foreach(var psm in psms)
                {
                    var ms2WithPrecursor = PseudoSearchMs2Scans.FirstOrDefault(s => s.OneBasedScanNumber == psm.ScanNumber && psm.ScanPrecursorCharge == s.PrecursorCharge && psm.ScanPrecursorMass == s.PrecursorMass);
                    if (ms2WithPrecursor == null) continue;
                    foreach (var ion in psm.MatchedFragmentIons)
                    {
                        double minMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMinimumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                        double maxMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMaximumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                        var masses = ms2WithPrecursor.GetClosestExperimentalIsotopicEnvelopeList(minMass, maxMass);
                        var targetMass = masses.FirstOrDefault(m => m.Charge == ion.Charge);
                        if (targetMass != null)
                        {
                            var fragmentPeak = ms2PeakIndexingEngine.GetIndexedPeak(targetMass.MonoisotopicMass, zeroBasedMs2ScanIndex, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                            if (fragmentPeak != null && Ms2PeakXicDictionary.ContainsKey(fragmentPeak))
                            {
                                var fragmentXic = Ms2PeakXicDictionary[fragmentPeak];
                                if (fragmentXic == null || visitedXics.Contains(fragmentXic)) continue;
                                var newSample = new PfPairTrainingSample(precursorXic, fragmentXic, true, psm);
                                allSamples.Add(newSample);
                                visitedXics.Add(fragmentXic);
                            }
                        }
                    }

                    foreach (var mass in ms2WithPrecursor.ExperimentalFragments)
                    {
                        var peak = ms2PeakIndexingEngine.GetIndexedPeak(mass.MonoisotopicMass, zeroBasedMs2ScanIndex, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                        if (peak != null && Ms2PeakXicDictionary.ContainsKey(peak))
                        {
                            var fragmentXic = Ms2PeakXicDictionary[peak];
                            if (!visitedXics.Contains(fragmentXic) && fragmentXic != null)
                            {
                                var newSample = new PfPairTrainingSample(precursorXic, fragmentXic, false, psm);
                                allSamples.Add(newSample);
                                visitedXics.Add(fragmentXic);
                            }
                        }
                    }
                }
            }
            return allSamples;
        }

        private List<Ms2ScanWithSpecificMass> GetMs2Scans(MsDataScan[] ms1Scans, MsDataScan[] ms2Scans, CommonParameters commonParameters)
        {
            List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];
            Parallel.ForEach(Partitioner.Create(0, ms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (partitionRange, loopState) =>
                {
                    var precursors = new List<(double MonoPeakMz, int Charge, double HighestIsotopeMz)>();

                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        precursors.Clear();
                        MsDataScan ms2scan = ms2Scans[i];

                        if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                        {
                            MsDataScan precursorScan = ms1Scans.First(s => s.OneBasedScanNumber == ms2scan.OneBasedPrecursorScanNumber);

                            foreach (MassSpectrometry.IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
                                    precursorScan.MassSpectrum, commonParameters.PrecursorDeconvolutionParameters))
                            {
                                double monoPeakMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
                                var highestMz = envelope.Peaks.MaxBy(p => p.intensity).mz;
                                precursors.Add((monoPeakMz, envelope.Charge, highestMz));
                            }
                        }
                        scansWithPrecursors[i] = new List<Ms2ScanWithSpecificMass>();
                        MassSpectrometry.IsotopicEnvelope[] neutralExperimentalFragments = null;

                        if (commonParameters.DissociationType != DissociationType.LowCID)
                        {
                            neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scan, commonParameters);
                        }

                        foreach (var precursor in precursors)
                        {
                            // assign precursor for this MS2 scan
                            var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.MonoPeakMz,
                                precursor.Charge, null, commonParameters, neutralExperimentalFragments, precursorHighestIsotopeMz: precursor.HighestIsotopeMz);
                            scansWithPrecursors[i].Add(scan);
                        }
                    }
                });

            var parentScans = scansWithPrecursors.Where(p => p.Any()).SelectMany(v => v).OrderBy(p => p.PrecursorMass).ToList();
            return parentScans;
        }
    }
}
