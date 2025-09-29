using Easy.Common.Extensions;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using Readers;
using Proteomics.ProteolyticDigestion;
using static Microsoft.FSharp.Core.ByRefKinds;
using Chemistry;

namespace EngineLayer.DIA
{
    public class ISDEngine_PSM : ISDEngine
    {
        private readonly MsDataFile DataFile;
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
            var ms2PeakEngines = new Dictionary<int, object>();
            var ms1PeakXicDictionary = new Dictionary<IIndexedPeak, ExtractedIonChromatogram>();
            var ms2PeakXicDictionary = new Dictionary<IIndexedPeak, ExtractedIonChromatogram>();
            var ms1Xics = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks1, out var indexingEngine1);
            foreach (var kvp in matchedPeaks1)
            {
                if (!ms1PeakXicDictionary.ContainsKey(kvp.Key))
                    ms1PeakXicDictionary.Add(kvp.Key, kvp.Value);
            }
            int keyIndex = 0;
            foreach (var ms2Group in isdVoltageMap)
            {
                allMs1Xics[ms2Group.Key] = ms1Xics;
                allMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(ms2Group.Value.ToArray(), out var matchedPeaks2, out var indexingEngine2);
                ms2PeakEngines[keyIndex] = indexingEngine2;
                foreach (var kvp in matchedPeaks2)
                {
                    if (!ms2PeakXicDictionary.ContainsKey(kvp.Key))
                        ms2PeakXicDictionary.Add(kvp.Key, kvp.Value);
                }
                keyIndex++;
            }
            int numberOfScansPerCycle = isdVoltageMap.Count() + 1;

            //Run classic search or load in PSMs
            PsmFromTsv[] allPsms = null;
            if (DIAparams.PsmPath == null)
            {
                var pseudoSearchScans = DDASearchModelTrainingEngine.GetMs2Scans(ms1Scans, allScans.Where(s => s.MsnOrder == 2).ToArray(), CommonParameters).ToArray();
                var allPsms1 = new SpectralMatch[pseudoSearchScans.Length];
                var fixedModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
                var variableModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
                var proteins = ModelTrainingEngine.LoadProteinDb(DIAparams.DbPath, true, DecoyType.Reverse, null, false, out var um, out int count, CommonParameters);
                var newClassicSearchEngine = new ClassicSearchEngine(allPsms1, pseudoSearchScans, variableModifications, fixedModifications, null, null, null, proteins, new DotMassDiffAcceptor("PlusOrMinus3Da", new List<double> { -3 * Chemistry.Constants.C13MinusC12, -2 * Chemistry.Constants.C13MinusC12, -1 * Chemistry.Constants.C13MinusC12, 0, 1 * Chemistry.Constants.C13MinusC12, 2 * Chemistry.Constants.C13MinusC12, 3 * Chemistry.Constants.C13MinusC12 },
                            CommonParameters.PrecursorMassTolerance), CommonParameters, null, null, null, true);
                var result = newClassicSearchEngine.Run();
            }
            else
            {
                var minMass = matchedPeaks1.Keys.Min(p => p.M);
                allPsms = SpectrumMatchTsvReader.ReadPsmTsv(DIAparams.PsmPath, out var warning).Where(p => p.PrecursorMass >= minMass).ToArray();
            }

            //Precursor-fragment Grouping
            var pairedXics = new Dictionary<ExtractedIonChromatogram, string>();
            var filteredPsms = new List<(string, int, int, double)>();
            var groupedPsms = allPsms.GroupBy(p => p.FullSequence).OrderByDescending(g => g.Max(p => p.PrecursorIntensity)).ToArray();
            var ms1PeakIndexingEngine = indexingEngine1 as MassIndexingEngine;
            foreach (var group in groupedPsms)
            {
                var allMatchedIons = group.SelectMany(p => p.MatchedIons.Where(i => i.IsInternalFragment == false)).GroupBy(i => (i.NeutralTheoreticalProduct.MonoisotopicMass, i.Charge)).ToArray();
                if (allMatchedIons.Length < 5) continue;
                var highestIntensityPrecursorPsm = group.MaxBy(p => p.PrecursorIntensity);
                int zeroBasedPrecursorScanIndex = (highestIntensityPrecursorPsm.PrecursorScanNum - 1) / numberOfScansPerCycle;
                var precursorPeak = ms1PeakIndexingEngine.GetIndexedPeak(highestIntensityPrecursorPsm.PrecursorMass, zeroBasedPrecursorScanIndex, DIAparams.Ms1XicConstructor.PeakFindingTolerance);

                if (precursorPeak == null || !matchedPeaks1.ContainsKey(precursorPeak)) continue;
                var precursorXic = matchedPeaks1[precursorPeak];
                if (precursorXic == null) continue;

                //find all fragment xics
                var filteredFragments = new List<MatchedFragmentIon>();
                var uniqueFragments = new Dictionary<(double mass, int charge), (int, MatchedFragmentIon)>();
                foreach(var psm in group)
                {
                    foreach(var ion in psm.MatchedIons.Where(i => i.NeutralTheoreticalProduct.IsInternalFragment == false))
                    {
                        if (!uniqueFragments.ContainsKey((ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge)))
                        {
                            uniqueFragments.Add((ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge), (psm.Ms2ScanNumber, ion));
                        }
                        else
                        {
                            if (ion.Intensity > uniqueFragments[(ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge)].Item2.Intensity)
                            {
                                uniqueFragments[(ion.NeutralTheoreticalProduct.MonoisotopicMass, ion.Charge)] = (psm.Ms2ScanNumber, ion);
                            }
                        }
                    }
                }

                foreach (var kvp in uniqueFragments)
                {
                    int zeroBasedScanIndex = (kvp.Value.Item1 - 1) / numberOfScansPerCycle;
                    int key = (kvp.Value.Item1 - 2) % numberOfScansPerCycle;
                    var ms2PeakIndexingEngine = ms2PeakEngines[key] as MassIndexingEngine;
                    var fragmentPeak = ms2PeakIndexingEngine.GetIndexedPeak(kvp.Value.Item2.Mz.ToMass(kvp.Value.Item2.Charge), zeroBasedScanIndex, DIAparams.Ms2XicConstructor.PeakFindingTolerance, kvp.Key.charge);
                    if (fragmentPeak == null || !ms2PeakXicDictionary.ContainsKey(fragmentPeak)) continue;
                    var fragmentXic = ms2PeakXicDictionary[fragmentPeak];
                    if (fragmentXic == null) continue;
                    //if (pairedXics.ContainsKey(fragmentXic) && pairedXics[fragmentXic] != highestIntensityPrecursorPsm.BaseSeq) continue;
                    if (Math.Abs(precursorXic.ApexRT - fragmentXic.ApexRT) < 0.3f)
                    {
                        var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(new List<ExtractedIonChromatogram> { precursorXic}, new List<ExtractedIonChromatogram> { fragmentXic });
                        if (pfGroups.Count > 0)
                        {
                            filteredFragments.Add(kvp.Value.Item2);
                            //if (!pairedXics.ContainsKey(fragmentXic)) pairedXics[fragmentXic] = highestIntensityPrecursorPsm.BaseSeq;
                        }
                    }
                }

                //var ms2WithMass = new Ms2ScanWithSpecificMass(ms1Scans.First(), highestIntensityPrecursorPsm.PrecursorMz, highestIntensityPrecursorPsm.PrecursorCharge, highestIntensityPrecursorPsm.DecoyContamTarget, CommonParameters);
                //var newPsm = new PeptideSpectralMatch(highestIntensityPrecursorPsm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer, highestIntensityPrecursorPsm.Notch.Value, filteredFragments.Count, highestIntensityPrecursorPsm.ScanIndex, ms2WithMass, CommonParameters, filteredFragments);
                if (filteredFragments.Count < 5) continue;
                filteredPsms.Add((highestIntensityPrecursorPsm.DecoyContamTarget, filteredFragments.Count, filteredFragments.GroupBy(f => f.NeutralTheoreticalProduct.Annotation).Count(), highestIntensityPrecursorPsm.PrecursorIntensity.Value));
            }

            var orderedPsms = filteredPsms.OrderByDescending(p => p.Item3).ToList();
            //new FdrAnalysisEngine(filteredPsms, 3, CommonParameters, this.FileSpecificParameters, null).Run();
            return new MetaMorpheusEngineResults(this);
        }
    }
}
