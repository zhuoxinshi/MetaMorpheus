using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;
using EngineLayer.DIA.ML;
using Chemistry;
using System.Collections.Concurrent;
using EngineLayer.ClassicSearch;
using Readers.SpectralLibrary;
using Omics;
using Microsoft.ML;

namespace EngineLayer.DIA
{
    public class DIA_MLEngine : DIAEngine
    {
        private readonly MsDataFile DataFile;
        public PseudoSearchScanType PseudoSearchScanType { get; set; }
        private readonly List<IBioPolymer> Proteins;
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
            var allMs1Xics = new Dictionary<(double min, double max), List<ExtractedIonChromatogram>>();
            var allMs2Xics = new Dictionary<(double min, double max), List<ExtractedIonChromatogram>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                allMs1Xics[ms2Group.Key] = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, new MzRange(ms2Group.Key.min, ms2Group.Key.max));
                allMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(ms2Group.Value.ToArray());
            }

            var pseudoSearchScans = GeneratePseudoSearchScans(PseudoSearchScanType, ms1Scans, ms2Scans);
            SpectralMatch[] fileSpecificPsms = new SpectralMatch[pseudoSearchScans.Length];

            var fixedModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var variableModifications = GlobalVariables.AllModsKnown.Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            var newClassicSearchEngine = new ClassicSearchEngine(fileSpecificPsms, pseudoSearchScans, variableModifications, fixedModifications, null,
                       null, null, Proteins, new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, CommonParameters.PrecursorMassTolerance), CommonParameters, this.FileSpecificParameters, null, NestedIds, false);
            var result = newClassicSearchEngine.Run();

            var model = TrainModel(fileSpecificPsms);



            throw new NotImplementedException();
        }

        public Ms2ScanWithSpecificMass[] GeneratePseudoSearchScans(PseudoSearchScanType pseudoSearchScanType, MsDataScan[] ms1Scans, MsDataScan[] ms2Scans)
        {
            switch (pseudoSearchScanType)
            {
                case (PseudoSearchScanType.DirectSearch):
                    var scansWithPrecursors = GetMs2Scans(ms1Scans,ms2Scans, CommonParameters).OrderBy(p => p.PrecursorMass).ToArray();
                    return scansWithPrecursors;
                default:
                    throw new NotImplementedException();
            }
        }

        public ITransformer TrainModel(SpectralMatch[] psms)
        {

            return null;
        }

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
    }
}
