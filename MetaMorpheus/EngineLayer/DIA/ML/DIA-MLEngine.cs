using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;
using Chemistry;
using System.Collections.Concurrent;
using EngineLayer.ClassicSearch;
using Omics;
using Microsoft.ML;
using Proteomics;
using Omics.Modifications;
using System.IO;
using UsefulProteomicsDatabases;
using CsvHelper.Configuration.Attributes;

namespace EngineLayer.DIA
{
    public class DIA_MLEngine : DIAEngine
    {
        private readonly MsDataFile DataFile;
        public MLbasedDIAparameters MlDIAparams => (MLbasedDIAparameters)DIAparams;
        public Dictionary<object, List<ExtractedIonChromatogram>> AllMs1Xics { get; private set; }
        public Dictionary<object, List<ExtractedIonChromatogram>> AllMs2Xics { get; private set; }
        public DIA_MLEngine(DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) : base(DIAparameters, dataFile, commonParameters, fileSpecificParameters, nestedIds)
        {
            DIAparams = DIAparameters;
            DataFile = dataFile;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //read in scans
            var ms1Scans = DataFile.GetMS1Scans().ToArray();
            var ms2Scans = DataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var DIAScanWindowMap = ConstructMs2Groups(ms2Scans);
            AllMs1Xics = new Dictionary<object, List<ExtractedIonChromatogram>>();
            AllMs2Xics = new Dictionary<object, List<ExtractedIonChromatogram>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                AllMs1Xics[ms2Group.Key] = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks, new MzRange(ms2Group.Key.min, ms2Group.Key.max));
                AllMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(ms2Group.Value.ToArray(), out matchedPeaks);
            }
            
            //construct pseudo scans and run a pseudo search
            var pseudoSearchScans = GeneratePseudoSearchScans(MlDIAparams.PseudoSearchType, ms1Scans, ms2Scans, out List <PrecursorFragmentsGroup> pseudoPfGroups).OrderBy(p => p.PrecursorMass).ToArray();
            var pseudoPsms = RunClassicSearch(pseudoSearchScans);

            //get normalized points if needed
            if (MlDIAparams.Features.Contains("SharedXIC"))
            {
                var allXics = AllMs1Xics.Values.SelectMany(p => p).Concat(AllMs2Xics.Values.SelectMany(p => p));
                foreach(var xic in allXics) xic.SetNormalizedPeakIntensities();
            }

            //train model
            var modelTrainingEngine = new ModelTrainingEngine(MlDIAparams, pseudoPsms, pseudoPfGroups.ToArray(), pseudoSearchScans);
            var model = modelTrainingEngine.TrainModel();

            //ml based pf grouping
            var groupingEngine = new MLgroupingEngine(model, MlDIAparams.PredictionScoreThreshold);
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                var pfGroups = groupingEngine.PrecursorFragmentGrouping(AllMs1Xics[ms2Group.Key], AllMs2Xics[ms2Group.Key]);
                allPfGroups.AddRange(pfGroups);
            }

            //make pseudo ms2 scans
            int oneBasedNumber = 1;
            PseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();
            foreach (var group in allPfGroups)
            {
                group.PFgroupIndex = oneBasedNumber;
                var scan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(group, MlDIAparams.PseudoMs2ConstructionType, CommonParameters, DataFile.FilePath);
                PseudoMs2Scans.Add(scan);
                oneBasedNumber++;
            }

            return new MetaMorpheusEngineResults(this);
        }

        public IEnumerable<Ms2ScanWithSpecificMass> GeneratePseudoSearchScans(PseudoSearchScanType pseudoSearchScanType, MsDataScan[] ms1Scans, MsDataScan[] ms2Scans, out List<PrecursorFragmentsGroup> pseudoPfGroups)
        {
            IEnumerable<Ms2ScanWithSpecificMass> scansWithPrecursors = null;
            pseudoPfGroups = new List<PrecursorFragmentsGroup>();
            switch (pseudoSearchScanType)
            {
                case (PseudoSearchScanType.DirectSearch):
                    scansWithPrecursors = GetMs2Scans(ms1Scans,ms2Scans, CommonParameters);
                    return scansWithPrecursors;

                case (PseudoSearchScanType.AllOverlap):
                    scansWithPrecursors = GeneratePseudoSearchScans_AllOverlap(out pseudoPfGroups);
                    return scansWithPrecursors;
                default:
                    throw new NotImplementedException();
            }
        }

        private IEnumerable<Ms2ScanWithSpecificMass> GeneratePseudoSearchScans_AllOverlap(out List<PrecursorFragmentsGroup> pseudoPfGroups)
        {
            var xicGroupingEngine = new XicGroupingEngine(0.5f, 0, -1, CommonParameters.MaxThreadsToUsePerFile, 1, 100, 500);
            pseudoPfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2Group in AllMs2Xics)
            {
                var groups = xicGroupingEngine.PrecursorFragmentGrouping(AllMs1Xics[ms2Group.Key], AllMs2Xics[ms2Group.Key]);
                pseudoPfGroups.AddRange(groups);
            }

            var pseudoSearchScans = new List<Ms2ScanWithSpecificMass>();
            int oneBasedNumber = 1;
            foreach (var group in pseudoPfGroups)
            {
                group.PFgroupIndex = oneBasedNumber;
                var pseudoSearchScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(group, MlDIAparams.PseudoMs2ConstructionType, CommonParameters, null);
                oneBasedNumber++;
                pseudoSearchScans.Add(pseudoSearchScan);
            }
            return pseudoSearchScans;
        }

        public SpectralMatch[] RunClassicSearch(Ms2ScanWithSpecificMass[] pseudoSearchScans)
        {
            SpectralMatch[] fileSpecificPsms = new SpectralMatch[pseudoSearchScans.Length];
            var fixedModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var variableModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            var decoyType = MlDIAparams.SearchDecoys ? DecoyType.Reverse : DecoyType.None;
            var proteins = LoadProteinDb(MlDIAparams.ProteinDb, true, decoyType, null, false, out var um, out int count, CommonParameters);
            var newClassicSearchEngine = new ClassicSearchEngine(fileSpecificPsms, pseudoSearchScans, variableModifications, fixedModifications, null,
                       null, null, proteins, new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, CommonParameters.PrecursorMassTolerance), CommonParameters, this.FileSpecificParameters, null, NestedIds, true);
            var result = newClassicSearchEngine.Run();
            return fileSpecificPsms.Where(p => p != null).ToArray();
        }

        public ITransformer TrainModel_DirectSearch(SpectralMatch[] psms)
        {
            var trainingSamples = new List<PfPairTrainingSample>();
            foreach (var psm in psms)
            {

            }
            var mlContext = new MLContext();
            var trainingData = mlContext.Data.LoadFromEnumerable(psms);
            var featureColumns = MlDIAparams.Features.ToArray();
            var pipeline = mlContext.Transforms.Concatenate("Features", featureColumns)
                .Append(mlContext.BinaryClassification.Trainers.LbfgsLogisticRegression(labelColumnName: "Label", featureColumnName: "Features"));
            var model = pipeline.Fit(trainingData);
            return model;
        }

        //Helper methods inaccessible to this class
        private List<Ms2ScanWithSpecificMass> GetMs2Scans(MsDataScan[] ms1Scans, MsDataScan[] ms2Scans, CommonParameters commonParameters)
        {
            List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];
            Parallel.ForEach(Partitioner.Create(0, ms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (partitionRange, loopState) =>
                {
                    var precursors = new List<(double MonoPeakMz, int Charge)>();

                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        precursors.Clear();
                        MsDataScan ms2scan = ms2Scans[i];

                        if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                        {
                            MsDataScan precursorScan = ms1Scans.First(s => s.OneBasedScanNumber == ms2scan.OneBasedPrecursorScanNumber);

                            foreach (IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
                                    precursorScan.MassSpectrum, commonParameters.PrecursorDeconvolutionParameters))
                            {
                                double monoPeakMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
                                precursors.Add((monoPeakMz, envelope.Charge));
                            }
                        }
                        scansWithPrecursors[i] = new List<Ms2ScanWithSpecificMass>();
                        IsotopicEnvelope[] neutralExperimentalFragments = null;

                        if (commonParameters.DissociationType != DissociationType.LowCID)
                        {
                            neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scan, commonParameters);
                        }

                        foreach (var precursor in precursors)
                        {
                            // assign precursor for this MS2 scan
                            var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.MonoPeakMz,
                                precursor.Charge, null, commonParameters, neutralExperimentalFragments);
                            scansWithPrecursors[i].Add(scan);
                        }
                    }
                });

            var parentScans = scansWithPrecursors.Where(p => p.Any()).SelectMany(v => v).OrderBy(p => p.OneBasedScanNumber);
            return parentScans.ToList();
        }

        protected static List<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, out Dictionary<string, Modification> um,
            out int emptyEntriesCount, CommonParameters commonParameters)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                um = null;
                proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, out dbErrors,
                    ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, commonParameters.MaxThreadsToUsePerFile, addTruncations: commonParameters.AddTruncations);
            }
            else
            {
                List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
                proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, commonParameters.MaxHeterozygousVariants, commonParameters.MinVariantDepth, addTruncations: commonParameters.AddTruncations);
            }

            emptyEntriesCount = proteinList.Count(p => p.BaseSequence.Length == 0);
            return proteinList.Where(p => p.BaseSequence.Length > 0).ToList();
        }
    }
}
