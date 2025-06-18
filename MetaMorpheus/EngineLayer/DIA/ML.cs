using Chemistry;
using Microsoft.ML;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer.DIA;
using System.Reflection.Emit;
using MathNet.Numerics.Statistics;
using MathNet.Numerics.Random;
using MzLibUtil;
using System.Linq.Expressions;

namespace EngineLayer.DIA
{
    public class PFPair
    {
        public float Correlation { get; set; }   // e.g., XIC correlation
        public float ApexRtDelta { get; set; }   // e.g., RT diff
        public float Overlap { get; set; }      // e.g., XIC overlap
        public float FragmentIntensity { get; set; } // e.g., fragment intensity ratio
        public bool Label { get; set; }       // true = target, false = decoy
        public float Weight { get; set; }
    }

    public class PredictionResult
    {
        public bool PredictedLabel { get; set; }
        public float Probability { get; set; }
        public float Score { get; set; }
    }

    public class TrainModel
    {
        public static void Train(IEnumerable<PFPair> yourDataList, bool weight = false)
        {
            var mlContext = new MLContext();

            if (weight)
            {
                AssignWeightsToPFPairs(yourDataList, ComputeWeights(yourDataList.ToList()));
            }
            IDataView data = mlContext.Data.LoadFromEnumerable(yourDataList);

            // Split into train/test
            var split = mlContext.Data.TrainTestSplit(data, testFraction: 0.2);
            var trainData = split.TrainSet;
            var testData = split.TestSet;

            var pipeline = mlContext.Transforms.Concatenate("Features", nameof(PFPair.Correlation), nameof(PFPair.ApexRtDelta))
                .Append(mlContext.BinaryClassification.Trainers.SdcaLogisticRegression(labelColumnName: "Label", featureColumnName: "Features"));
            var model = pipeline.Fit(trainData);
            var predictions = model.Transform(testData);
            var metrics = mlContext.BinaryClassification.Evaluate(predictions);

            var predictor = mlContext.Model.CreatePredictionEngine<PrecursorFragmentPair, PredictionResult>(model);
            Console.WriteLine($"AUC: {metrics.AreaUnderRocCurve:F4}");
            Console.WriteLine($"F1 Score: {metrics.F1Score:F4}");
            Console.WriteLine($"Accuracy: {metrics.Accuracy:P2}");
        }

        public static PFPair GetPFPairsFromPfPairMetrics(PFpairMetrics pfPairMetrics)
        {
            var pfPair = new PFPair
            {
                Correlation = (float)pfPairMetrics.Correlation,
                ApexRtDelta = (float)pfPairMetrics.ApexRtDelta,
                Overlap = (float)pfPairMetrics.Overlap,
                FragmentIntensity = (float)pfPairMetrics.FragmentIntensity,
                Label = pfPairMetrics.MatchedIonType == "NA" ? false : true
            };
            return pfPair;
        }

        public static List<PFPair> GetPFPairsFromPFMetricsExcludingInternal(IEnumerable<PFpairMetrics> pfPairMetricsList)
        {
            var list = new List<PFPair>();
            foreach (var pfPairMetrics in pfPairMetricsList)
            {
                if (pfPairMetrics.MatchedIonType == "Internal") continue;
                var pfPair = GetPFPairsFromPfPairMetrics(pfPairMetrics);
                list.Add(pfPair);
            }
            return list;
        }

        public static (float, float) ComputeWeights(List<PFPair> pfPairs)
        {
            float w0 = pfPairs.Count / (2 * pfPairs.Where(p => p.Label == false).Count());
            float w1 = pfPairs.Count / (2 * pfPairs.Where(p => p.Label == true).Count());
            return (w0, w1);
        }

        public static void AssignWeightsToPFPairs(IEnumerable<PFPair> pfPairs, (float, float) weights)
        {
            foreach (var pfPair in pfPairs)
            {
                pfPair.Weight = pfPair.Label ? weights.Item2 : weights.Item1;
            }
        }

        //this method assumes that Xic construction is already done
        public static PrecursorFragmentsGroup FindPfGroupBasedOnPsm(PsmFromTsv psmTsv, List<DeconvolutedMass>[] ms1PeakTable, Dictionary<double, List<DeconvolutedMass>[]> allMs2MassTables,
            DIAparameters diaParam)
        {
            var cycle = (psmTsv.Ms2ScanNumber - 1) % 4;
            List<DeconvolutedMass>[] ms2PeakTable = null;
            switch (cycle)
            {
                case 1:
                    ms2PeakTable = allMs2MassTables[60];
                    break;
                case 2:
                    ms2PeakTable = allMs2MassTables[80];
                    break;
                case 3:
                    ms2PeakTable = allMs2MassTables[100];
                    break;
            }
            var precursorMass = MassCurve.GetMassFromScan(psmTsv.PrecursorMass, psmTsv.PrecursorCharge, ms1PeakTable, (psmTsv.PrecursorScanNum - 1) / 4,
                diaParam.Ms1PeakFindingTolerance, diaParam.PeakSearchBinSize);
            if (precursorMass == null || precursorMass.PeakCurve == null)
            {
                return null;
            }
            var pfGroup = new PrecursorFragmentsGroup(precursorMass.PeakCurve);
            precursorMass.PeakCurve.GetRawXYData();

            foreach (var matchedIon in psmTsv.MatchedIons)
            {
                var fragmentMass = MassCurve.GetMassFromScan(matchedIon.Mz.ToMass(matchedIon.Charge), matchedIon.Charge, ms2PeakTable, (psmTsv.Ms2ScanNumber - 1) / 4,
                                       diaParam.Ms2PeakFindingTolerance, diaParam.PeakSearchBinSize);
                if (fragmentMass.PeakCurve != null)//fragmentMass == null ||
                {
                    double overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    fragmentMass.PeakCurve.GetRawXYData();
                    double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    var pfPair = new PrecursorFragmentPair(precursorMass.PeakCurve, fragmentMass.PeakCurve, overlap, corr, psmTsv.DecoyContamTarget);
                    pfGroup.PFpairs.Add(pfPair);
                }
            }
            return pfGroup;
        }

        public static List<PFpairMetrics> GetPFPairMetricsFromPsm(PsmFromTsv psmTsv, List<Peak>[] ms1PeaksByScan, List<Peak>[] ms2PeaksByScan, DIAparameters diaParam)
        {
            var pfPairMetricsList = new List<PFpairMetrics>();
            var precursorMass = SearchPeakByScan(psmTsv.PrecursorMass, psmTsv.PrecursorCharge, ms1PeaksByScan[psmTsv.PrecursorScanNum], diaParam.Ms1PeakFindingTolerance);
            if (precursorMass == null || precursorMass.PeakCurve == null)
            {
                return null;
            }
            var pfGroup = new PrecursorFragmentsGroup(precursorMass.PeakCurve);
            precursorMass.PeakCurve.GetRawXYData();

            var visited = new HashSet<Peak>();
            foreach (var matchedIon in psmTsv.MatchedIons.Where(i => i.IsInternalFragment == false))
            {
                var fragmentMass = SearchPeakByScan(matchedIon.Mz.ToMass(matchedIon.Charge), matchedIon.Charge, ms2PeaksByScan[psmTsv.Ms2ScanNumber],
                                       diaParam.Ms2PeakFindingTolerance);
                if (fragmentMass != null && fragmentMass.PeakCurve != null)//fragmentMass == null ||
                {
                    double overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    fragmentMass.PeakCurve.GetRawXYData();
                    double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorMass.PeakCurve, fragmentMass.PeakCurve);
                    var pfPair = new PrecursorFragmentPair(precursorMass.PeakCurve, fragmentMass.PeakCurve, overlap, corr, psmTsv.DecoyContamTarget);
                    var pfPairMetrics = new PFpairMetrics(pfPair);
                    pfPairMetrics.MatchedIonType = "Terminal";
                    visited.Add(fragmentMass);
                    pfPairMetricsList.Add(pfPairMetrics);
                }
            }

            var umatchedMasses = ms2PeaksByScan[psmTsv.Ms2ScanNumber].Where(p => !visited.Contains(p)).ToList();
            foreach(var mass in umatchedMasses)
            {
                double overlap = PrecursorFragmentPair.CalculateRTOverlapRatio(precursorMass.PeakCurve, mass.PeakCurve);
                mass.PeakCurve.GetRawXYData();
                double corr = PrecursorFragmentPair.CalculatePeakCurveCorrXYData(precursorMass.PeakCurve, mass.PeakCurve);
                var pfPair = new PrecursorFragmentPair(precursorMass.PeakCurve, mass.PeakCurve, overlap, corr, psmTsv.DecoyContamTarget);
                var pfPairMetrics = new PFpairMetrics(pfPair);
                pfPairMetrics.MatchedIonType = "NA";
                pfPairMetricsList.Add(pfPairMetrics);
            }
            return pfPairMetricsList;
        }

        public static Peak SearchPeakByScan(double targetMass, int targetCharge, List<Peak> peaks, Tolerance tolerance)
        {
            var allPossiblePeaks = peaks.Where(p => p.Charge == targetCharge).ToArray();
            if (allPossiblePeaks == null || allPossiblePeaks.Length == 0)
            {
                return null; // No peaks found for the specified charge
            }
            int index = Array.BinarySearch(allPossiblePeaks.Select(p => p.MonoisotopicMass).ToArray(), targetMass);
            if (index < 0)
            {
                index = ~index; 
            }
            if (index == 0)
            {
                if (tolerance.Within(targetMass, allPossiblePeaks[0].MonoisotopicMass))
                    return allPossiblePeaks[0];
                else
                    return null;
            }
            if (index == allPossiblePeaks.Length)
            {
                if (tolerance.Within(targetMass, allPossiblePeaks[allPossiblePeaks.Length - 1].MonoisotopicMass))
                    return allPossiblePeaks[allPossiblePeaks.Length - 1];
                else
                    return null;
            }

            Peak bestMass = null;
            if (tolerance.Within(targetMass, allPossiblePeaks[index - 1].MonoisotopicMass))
            {
                bestMass = allPossiblePeaks[index - 1];
            }
            if (tolerance.Within(targetMass, allPossiblePeaks[index].MonoisotopicMass)
                && (bestMass == null || Math.Abs(bestMass.MonoisotopicMass - targetMass) > Math.Abs(allPossiblePeaks[index].MonoisotopicMass - targetMass)))
            {
                bestMass = allPossiblePeaks[index];
            }
            return bestMass;
        }
    }
    
}
