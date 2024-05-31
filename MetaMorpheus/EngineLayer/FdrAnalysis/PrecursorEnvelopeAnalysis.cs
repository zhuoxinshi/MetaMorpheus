
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;
using Chemistry;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using ThermoFisher.CommonCore.Data.Business;
using Microsoft.ML;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearRegression;
using SharpLearning.Containers.Matrices;
using Newtonsoft.Json.Linq;
using Easy.Common.Extensions;
using System.Reflection;
using MassSpectrometry.MzSpectra;
using static MassSpectrometry.MzSpectra.SpectralSimilarity;

namespace EngineLayer.FdrAnalysis
{
    public class PrecursorEnvelopeAnalysis
    {
        public MsDataFile RawScans;
        public List<SpectralMatch> AllPsms;
        public double MzTolerance;

        public PrecursorEnvelopeAnalysis(MsDataFile rawScans, List<SpectralMatch> psms, double tolerance)
        {
            RawScans = rawScans;
            AllPsms = psms;
            MzTolerance = tolerance;
        }

        public static List<(double mz, double intensity)> FindTheoreticalIsotopicPeaks(SpectralMatch psm, int charge)
        {
            ChemicalFormula formula = new Peptide(psm.FullSequence).GetChemicalFormula();
            IsotopicDistribution theoreticalDistribution = IsotopicDistribution.GetDistribution(formula);

            List<(double mz, double intensity)> theoreticalPeaks = new List<(double mz, double intensity)>();
            foreach (double mass in theoreticalDistribution.Masses)
            {
                double mz = mass.ToMz(charge);
                if (mz >= psm.MsDataScan.IsolationRange.Minimum && mz <= psm.MsDataScan.IsolationRange.Maximum)
                {
                    double intensity = theoreticalDistribution.Intensities[Array.IndexOf(theoreticalDistribution.Masses, mass)];
                    theoreticalPeaks.Add(new(mz, intensity));
                }
            }
            return theoreticalPeaks;
        }

        public static List<List<(double mz, double intensity)>> FindTheoreticalMs1Peaks(List<SpectralMatch> psms, int maxCharge)
        {
            //find all theoretical peaks for the set of PSMs and store them in a list of list of peaks
            //outside list is for the sequence from each PSM, inside list is the peaks for each charge state of a specific sequence
            var allTheoreticalPeaks = new List<List<(double mz, double intensity)>>();
            foreach (var psm in psms)
            {
                for (int charge = 1; charge <= maxCharge; charge++)
                {
                    var theoreticalPeaks = FindTheoreticalIsotopicPeaks(psm, charge);
                    allTheoreticalPeaks.Add(theoreticalPeaks);
                }
            }
            return allTheoreticalPeaks;
        }

        public static bool WithinTolerance(double mz1, double mz2, double tolerance)
        {
            return Math.Abs(mz1 - mz2) < tolerance;
        }

        public static List<(double experimentalMz, double theoreticalMz)> MatchedMzs(double[] experimentalMzs, double[] experimentalIntensities, double[] theoreticalMzs, double tolerance)
        {
            List<(double, double)> mzPairs = new List<(double, double)>();
            List<(double mz, double intensity)> experimental = new();
            for (int i = 0; i < experimentalMzs.Length; i++)
            {
                experimental.Add((experimentalMzs[i], experimentalIntensities[i]));
            }
            experimental = experimental.OrderByDescending(i => i.intensity).ToList();

            foreach(var theoreticalMz in theoreticalMzs)
            {
                int index = 0;
                while (experimental.Count > 0 && index < experimental.Count)
                {
                    if (WithinTolerance(experimental[index].mz, theoreticalMz, tolerance))
                    {
                        mzPairs.Add((experimental[index].mz, theoreticalMz));
                        experimental.RemoveAt(index);
                        index = -1;
                        break;
                    }
                    index++;
                }
                if (experimental.Count == 0)
                {
                    index++;
                }
                if (index > 0)
                {
                    //didn't find a experimental mz in range
                    mzPairs.Add((-1, theoreticalMz));
                }
            }

            return mzPairs;

        }
        public static List<double> FindMatchedIntensities (List<SpectralMatch> psms, int maxCharge, MsDataFile dataFile, double tolerance)
        {
            var allTheoreticalPeaks = FindTheoreticalMs1Peaks(psms, maxCharge);
            double[] allTheoreticalMzs = allTheoreticalPeaks.SelectMany(peaks => peaks.Select(p => p.mz)).Distinct().OrderBy(mz => mz).ToArray();
            List<double> matchedIndex = new List<double>();
            var precursorSpectrum = dataFile.Scans.Where(scan => scan.OneBasedScanNumber == psms.First().PrecursorScanNumber).First().MassSpectrum;

            var mzPairs = MatchedMzs(precursorSpectrum.XArray, precursorSpectrum.YArray, allTheoreticalMzs, tolerance);
            var matchedMzsExperimental = mzPairs.Select(mz => mz.Item1).OrderBy(mz => mz);
            List<double> matchedIntensities = new List<double>();

            foreach(double mz in matchedMzsExperimental)
            {
                if (mz >= 0)
                {
                    int index = precursorSpectrum.XArray.IndexOf(mz);
                    matchedIntensities.Add(precursorSpectrum.YArray[index]);
                }
                else
                {
                    matchedIntensities.Add(0);
                }
            }

            return matchedIntensities;
        }

        public static double FindFractionOfMatchedIntensities(List<SpectralMatch> psms, int maxCharge, MsDataFile dataFile, double tolerance)
        {
            var precursorSpectrum = dataFile.Scans.Where(scan => scan.OneBasedScanNumber == psms.First().PrecursorScanNumber).First().MassSpectrum;
            var matchedIntensities = FindMatchedIntensities(psms, maxCharge, dataFile, tolerance);
            double fractionMatched = matchedIntensities.Sum() / precursorSpectrum.SumOfAllY;

            return fractionMatched;
        }

        public static MzSpectrum GetTheoreticalMs1Spectrum(List<SpectralMatch> psms, int maxCharge, MsDataFile dataFile, double tolerance)
        {
            var allTheoreticalPeaks = FindTheoreticalMs1Peaks(psms, maxCharge);
            double[] allTheoreticalMzs = allTheoreticalPeaks.SelectMany(peaks => peaks.Select(p => p.mz)).Distinct().OrderBy(mz => mz).ToArray();

            //make a feature matrix with each row one mz and each column as a peptide ion
            double[][] features = new double[allTheoreticalMzs.Length][];
            for (int i = 0; i < allTheoreticalMzs.Length; i++)
            {
                for (int j = 0; j < allTheoreticalPeaks.Count; j++)
                {
                    if (allTheoreticalPeaks[j].Select(p => p.mz).Contains(allTheoreticalMzs[i]))
                    {
                        int indexOfMz = allTheoreticalPeaks[j].Select(p => p.mz).ToList().IndexOf(allTheoreticalMzs[i]);
                        features[i][j] = allTheoreticalPeaks[j][indexOfMz].intensity;
                    }
                    else
                    {
                        features[i][j] = 0;
                    }
                }
            }

            double[] matchedIntensities = FindMatchedIntensities(psms, maxCharge, dataFile, tolerance).ToArray();

            //make a multiple linear regression model to find the relative intensities for each peptide ion
            var coefficients = MultipleRegression.QR(features, matchedIntensities, false);

            //Find the normalized theoretical spectrum by matrix multiplication
            var featureMatrix = new double[allTheoreticalMzs.Length, allTheoreticalPeaks.Count()];
            for (int i = 0; i < allTheoreticalMzs.Length; i++)
            {
                for (int j = 0; j < allTheoreticalPeaks.Count(); j++)
                {
                    featureMatrix[i, j] = features[i][j];
                }
            }
            var matrix = Matrix<double>.Build.DenseOfArray(featureMatrix);
            var vector = Vector<double>.Build.Dense(coefficients);
            var resultIntensities = matrix.Multiply(vector).ToArray();
            MzSpectrum theoreticalSpectrum = new MzSpectrum(allTheoreticalMzs, resultIntensities, true);

            return theoreticalSpectrum;
        }

        public static double? CalculateSimilarityScore (List<SpectralMatch> psms, int maxCharge, MsDataFile dataFile, double tolerance, 
            SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks) 
        {
            MzSpectrum theoreticalSpectrum = GetTheoreticalMs1Spectrum(psms, maxCharge, dataFile, tolerance);
            var precursorSpectrum = dataFile.Scans.Where(scan => scan.OneBasedScanNumber == psms.First().PrecursorScanNumber).First().MassSpectrum;
            var allTheoreticalMzs = theoreticalSpectrum.XArray;
            var mzPairs = MatchedMzs(precursorSpectrum.XArray, precursorSpectrum.YArray, allTheoreticalMzs, tolerance);
            var matchedMzsExperimental = mzPairs.Where(pair => pair.Item1 >= 0).Select(pair =>pair.Item1).OrderBy(mz => mz).ToArray();
            var matchedIntensitiesExperimental = FindMatchedIntensities(psms, maxCharge, dataFile, tolerance).Where(i => i > 0).ToArray();
            MzSpectrum experimentalSpectrum = new MzSpectrum(matchedMzsExperimental, matchedIntensitiesExperimental, true);

            SpectralSimilarity similarity = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, scheme, toleranceInPpm, allPeaks);

            return similarity.CosineSimilarity();
        }


        public void RunAnalysis()
        {
            var psmSets = AllPsms.GroupBy(p => p.ScanNumber);
            foreach (var psmSet in psmSets)
            {
                var psmList = psmSet.ToList();
                //double? similarityScore = CalculateSimilarityScore(RawScans, psmList);
                //pmList.ForEach(psm => psm.PrecursorScanSmilarityScore = similarityScore);
            }
        }
    }
}
