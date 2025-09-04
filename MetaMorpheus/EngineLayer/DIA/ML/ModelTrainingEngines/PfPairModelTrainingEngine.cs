using FlashLFQ;
using MassSpectrometry;
using Microsoft.ML;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PfPairModelTrainingEngine : ModelTrainingEngine
    {
        public PrecursorFragmentsGroup[] PfGroups { get; set; }
        public Dictionary<object, List<ExtractedIonChromatogram>> AllMs1Xics { get; set; }
        public Dictionary<object, List<ExtractedIonChromatogram>> AllMs2Xics { get; set; }

        public PfPairModelTrainingEngine(MLbasedDIAparameters mlDIAparams, CommonParameters commonParameters, Dictionary<object, List<ExtractedIonChromatogram>> allMs1Xics, Dictionary<object, List<ExtractedIonChromatogram>> allMs2Xics) :base(mlDIAparams, commonParameters)
        {
            AllMs1Xics = allMs1Xics;
            AllMs2Xics = allMs2Xics;
        }

        public override void GeneratePseudoSearchScans()
        {
            var xicGroupingEngine = new XicGroupingEngine(0.5f, 0, -1, CommonParameters.MaxThreadsToUsePerFile, 1, 100, 500);
            var pseudoPfGroups = new List<PrecursorFragmentsGroup>();
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
            PseudoSearchMs2Scans = pseudoSearchScans.ToArray();
        }

        public override IEnumerable<PfPairTrainingSample> GetTrainingSamples()
        {
            var filteredPsms = Psms.Where(p => p.Score >= MlDIAparams.PsmScoreCutOff).ToArray();
            var peptides = filteredPsms.GroupBy(p => p.FullSequence).ToList();
            var allSamples = new List<PfPairTrainingSample>();
            foreach (var psm in filteredPsms)
            {
                var correspondingPfGroup = PfGroups.FirstOrDefault(g => g.PFgroupIndex == psm.ScanNumber);
                var correspondingMs2Scan = PseudoSearchMs2Scans.FirstOrDefault(s => s.OneBasedScanNumber == psm.ScanNumber);
                var samplesFromThisPsm = GetTrainingSamplesFromPfGroup(psm, correspondingPfGroup, correspondingMs2Scan, MlDIAparams.Ms2XicConstructor.PeakFindingTolerance);
                allSamples.AddRange(samplesFromThisPsm);
            }
            return allSamples;
        }

        public static List<PfPairTrainingSample> GetTrainingSamplesFromPfGroup(SpectralMatch psm, PrecursorFragmentsGroup pfGroup, Ms2ScanWithSpecificMass ms2WithPrecursor, Tolerance tol)
        {
            var samples = new List<PfPairTrainingSample>();
            var sortedPfPairs = pfGroup.PFpairs.OrderByDescending(pf => pf.FragmentXic.AveragedMassOrMz).ToArray();
            var sortedPfPairMzs = sortedPfPairs.Select(pf => (double)pf.FragmentXic.ApexPeak.M).ToArray();

            var positiveIndices = new List<int>();
            foreach (var ion in psm.MatchedFragmentIons)
            {
                double minMass = tol.GetMinimumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                double maxMass = tol.GetMaximumValue(ion.NeutralTheoreticalProduct.MonoisotopicMass);
                var masses = ms2WithPrecursor.GetClosestExperimentalIsotopicEnvelopeList(minMass, maxMass);
                var targetMass = masses.FirstOrDefault(m => m.Charge == ion.Charge);
                if (targetMass != null)
                {
                    foreach (var peak in targetMass.Peaks)
                    {
                        var closestIndex = FindClosestIndexOfPfPair(sortedPfPairMzs, peak.mz);
                        positiveIndices.Add(closestIndex);
                    }
                }
            }

            for (int i = 0; i < sortedPfPairs.Length; i++)
            {
                var newSample = new PfPairTrainingSample(sortedPfPairs[i]);
                if (positiveIndices.Contains(i))
                {
                    newSample.Label = true;
                }
                else
                {
                    newSample.Label = false;
                }
                samples.Add(newSample);
            }
            return samples;
        }


        private static int FindClosestIndexOfPfPair(double[] sortedPfPairMzs, double targetMz)
        {
            int closestIndex = -1;
            var index = Array.BinarySearch(sortedPfPairMzs, targetMz);
            if (index >= 0)
            {
                return index;
            }
            else
            {
                index = ~index;
                if (index == 0)
                {
                    closestIndex = 0;
                }
                else if (index >= sortedPfPairMzs.Length)
                {
                    closestIndex = sortedPfPairMzs.Length - 1;
                }
                else
                {
                    if (Math.Abs(sortedPfPairMzs[index] - targetMz) < Math.Abs(sortedPfPairMzs[index - 1] - targetMz))
                    {
                        closestIndex = index;
                    }
                    else
                    {
                        closestIndex = index - 1;
                    }
                }
            }
            return closestIndex;
        }
    }
}
