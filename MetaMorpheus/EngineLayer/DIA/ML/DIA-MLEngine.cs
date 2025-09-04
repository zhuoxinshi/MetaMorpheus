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
            var allMs1Xics = new Dictionary<object, List<ExtractedIonChromatogram>>();
            var allMs2Xics = new Dictionary<object, List<ExtractedIonChromatogram>>();
            var ms1PeakEngines = new Dictionary<object, object>();
            var ms2PeakEngines = new Dictionary<object, object>();
            var peakXicDictionary = new Dictionary<IIndexedPeak, ExtractedIonChromatogram>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                allMs1Xics[ms2Group.Key] = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks1, out var indexingEngine1, new MzRange(ms2Group.Key.min, ms2Group.Key.max));
                ms1PeakEngines[ms2Group.Key] = indexingEngine1;
                allMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(ms2Group.Value.ToArray(), out var matchedPeaks2, out var indexingEngine2);
                ms2PeakEngines[ms2Group.Key] = indexingEngine2;

                var allKeys = matchedPeaks1.Select(p => p.Key).Concat(matchedPeaks2.Select(p => p.Key));
                var duplicateKeys = allKeys
                    .GroupBy(k => k)
                    .Where(g => g.Count() > 1)
                    .Select(g => g.Key).ToList();
            }

            //train model
            ITransformer model = null;
            ModelTrainingEngine modelTrainingEngine = null;
            switch (MlDIAparams.PseudoSearchType)
            {
                case (PseudoSearchScanType.DirectSearch):
                    modelTrainingEngine = new DDASearchModelTrainingEngine(MlDIAparams, CommonParameters, ms1Scans, ms2Scans, ms1PeakEngines, ms2PeakEngines, peakXicDictionary, DIAScanWindowMap);
                    model = modelTrainingEngine.TrainModel();
                    break;
                case (PseudoSearchScanType.AllOverlap):
                    modelTrainingEngine = new PfPairModelTrainingEngine(MlDIAparams, CommonParameters, allMs1Xics, allMs2Xics);
                    model = modelTrainingEngine.TrainModel();
                    break;
                default:
                    throw new NotImplementedException();
            }

            //ml based pf grouping
            var groupingEngine = new MLgroupingEngine(model, MlDIAparams.PredictionScoreThreshold);
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                var pfGroups = groupingEngine.PrecursorFragmentGrouping(allMs1Xics[ms2Group.Key], allMs2Xics[ms2Group.Key]);
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
    }
}
