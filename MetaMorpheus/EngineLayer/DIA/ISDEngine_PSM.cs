using Easy.Common.Extensions;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using ThermoFisher.CommonCore.Data.Business;
using UsefulProteomicsDatabases;
using static Microsoft.FSharp.Core.ByRefKinds;
namespace EngineLayer.DIA
{
    public class ISDEngine_PSM : ISDEngine
    {
        private readonly MsDataFile DataFile;
        public MLbasedDIAparameters MlDIAparams => (MLbasedDIAparameters)DIAparams;
        public ISDEngine_PSM(DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) : base(DIAparameters, dataFile, commonParameters, fileSpecificParameters, nestedIds)
        {
            DIAparams = DIAparameters;
            DataFile = dataFile;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //read in scans and isd scan pre-process
            var allScans = DataFile.GetAllScansList().ToArray();
            var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
            ReLabelIsdScans(isdVoltageMap, allScans);

            //Get all MS1 and MS2 XICs
            var allMs1Xics = new Dictionary<object, List<ExtractedIonChromatogram>>();
            var allMs2Xics = new Dictionary<object, List<ExtractedIonChromatogram>>();
            var ms2PeakEngines = new List<MassIndexingEngine>();
            var ms1PeakXicDictionary = new Dictionary<IIndexedPeak, ExtractedIonChromatogram>();
            var ms2PeakXicDictionary = new Dictionary<IIndexedPeak, ExtractedIonChromatogram>();
            var ms1Xics = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks1, out var indexingEngine1);
            foreach (var kvp in matchedPeaks1)
            {
                if (!ms1PeakXicDictionary.ContainsKey(kvp.Key))
                    ms1PeakXicDictionary.Add(kvp.Key, kvp.Value);
            }
            foreach (var ms2Group in isdVoltageMap)
            {
                allMs1Xics[ms2Group.Key] = ms1Xics;
                allMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(ms2Group.Value.ToArray(), out var matchedPeaks2, out var indexingEngine2);
                ms2PeakEngines.Add(indexingEngine2 as MassIndexingEngine);
                foreach (var kvp in matchedPeaks2)
                {
                    if (!ms2PeakXicDictionary.ContainsKey(kvp.Key))
                        ms2PeakXicDictionary.Add(kvp.Key, kvp.Value);
                }
            }
            int numberOfScansPerCycle = isdVoltageMap.Count() + 1;

            //Get pseudo MS2 scans for classic search
            var pseudoSearchScans = DDASearchModelTrainingEngine.GetMs2Scans(ms1Scans, allScans.Where(s => s.MsnOrder == 2).ToArray(), CommonParameters).ToArray();

            //Run classic search
            SpectralMatch[] allPsms = new SpectralMatch[pseudoSearchScans.Length];
            var fixedModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var variableModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var proteins = ModelTrainingEngine.LoadProteinDb(MlDIAparams.ProteinDb, true, DecoyType.Reverse, null, false, out var um, out int count, CommonParameters);
            var newClassicSearchEngine = new ClassicSearchEngine(allPsms, pseudoSearchScans, variableModifications, fixedModifications, null, null, null, proteins, new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, CommonParameters.PrecursorMassTolerance), CommonParameters, null, null, null, true);
            var result = newClassicSearchEngine.Run();

            //Precursor-fragment Grouping
            var filteredPsms = new List<SpectralMatch>();
            var groupedPsms = allPsms.GroupBy(p => p.FullSequence).ToArray();
            var ms1PeakIndexingEngine = indexingEngine1 as MassIndexingEngine;
            foreach (var group in groupedPsms)
            {
                var allMatchedIons = group.SelectMany(p => p.MatchedFragmentIons).GroupBy(i => i.NeutralTheoreticalProduct.FragmentNumber).ToArray();
                if (allMatchedIons.Length < 5) continue;
                var highestIntensityPrecursorPsm = group.MaxBy(p => p.PrecursorScanIntensity);
                int zeroBasedPrecursorScanIndex = (highestIntensityPrecursorPsm.PrecursorScanNumber.Value - 1) / numberOfScansPerCycle;
                var precursorPeak = ms1PeakIndexingEngine.GetIndexedPeak(highestIntensityPrecursorPsm.ScanPrecursorMass, zeroBasedPrecursorScanIndex, MlDIAparams.Ms1XicConstructor.PeakFindingTolerance);

                if (precursorPeak == null && !matchedPeaks1.ContainsKey(precursorPeak)) continue;
                var precursorXic = matchedPeaks1[precursorPeak];
                if (precursorXic == null) continue;

                //find all fragment xics
                var filteredFragments = new List<MatchedFragmentIon>();
                var uniqueFragments = new Dictionary<(double mass, int charge), (int, MatchedFragmentIon)>();
                foreach(var psm in group)
                {
                    foreach(var ion in psm.MatchedFragmentIons)
                    {
                        if (!uniqueFragments.ContainsKey((ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge)))
                        {
                            uniqueFragments.Add((ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge), (psm.ScanNumber, ion));
                        }
                        else
                        {
                            if (ion.Intensity > uniqueFragments[(ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge)].Item2.Intensity)
                            {
                                uniqueFragments[(ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge)] = (psm.ScanNumber, ion);
                            }
                        }
                    }
                }

                foreach (var kvp in uniqueFragments)
                {
                    double minMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMinimumValue(kvp.Value.Item2.NeutralTheoreticalProduct.MonoisotopicMass);
                    double maxMass = MlDIAparams.Ms2XicConstructor.PeakFindingTolerance.GetMaximumValue(kvp.Value.Item2.NeutralTheoreticalProduct.MonoisotopicMass);
                    int zeroBasedScanIndex = (kvp.Value.Item1 - 1) / numberOfScansPerCycle;
                    int key = kvp.Value.Item1 % numberOfScansPerCycle;
                    var ms2PeakIndexingEngine = ms2PeakEngines[key];
                    var fragmentPeak = ms2PeakIndexingEngine.GetIndexedPeak(kvp.Key.mass, zeroBasedScanIndex, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance, kvp.Key.charge);
                    if (fragmentPeak == null || !ms2PeakXicDictionary.ContainsKey(fragmentPeak)) continue;
                    var fragmentXic = ms2PeakXicDictionary[fragmentPeak];
                    if (fragmentXic == null) continue;
                    if (Math.Abs(precursorXic.ApexRT - fragmentXic.ApexRT) < MlDIAparams.ApexRtTolerance)
                    {
                        var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(new List<ExtractedIonChromatogram> { precursorXic}, new List<ExtractedIonChromatogram> { fragmentXic });
                        if (pfGroups.Count > 0)
                        {
                            filteredFragments.Add(kvp.Value.Item2);
                        }
                    }
                }

                var ms2WithMass = new Ms2ScanWithSpecificMass(null, highestIntensityPrecursorPsm.ScanPrecursorMonoisotopicPeakMz, highestIntensityPrecursorPsm.ScanPrecursorCharge, highestIntensityPrecursorPsm.FullFilePath, CommonParameters);
                var newPsm = new PeptideSpectralMatch(highestIntensityPrecursorPsm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer, highestIntensityPrecursorPsm.Notch.Value, filteredFragments.Count, highestIntensityPrecursorPsm.ScanIndex, ms2WithMass, CommonParameters, filteredFragments);
            }

            return new MetaMorpheusEngineResults(this);
        }
    }
}
