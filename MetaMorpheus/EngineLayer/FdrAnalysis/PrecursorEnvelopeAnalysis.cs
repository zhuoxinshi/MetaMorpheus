
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
using TopDownProteomics.MassSpectrometry;
using MathNet.Numerics.LinearAlgebra.Complex;
using Proteomics.ProteolyticDigestion;
using System.Runtime.InteropServices;
using 

namespace EngineLayer.FdrAnalysis
{
    public class PrecursorEnvelopeAnalysis
    {
        public MsDataFile RawScans;
        public List<SpectralMatch> AllPsms;
        public double MzTolerance;
        public int[] ChargeStates;
        MsDataFile AllDataScans;

        public PrecursorEnvelopeAnalysis(MsDataFile rawScans, List<SpectralMatch> psms, double tolerance)
        {
            RawScans = rawScans;
            AllPsms = psms;
            MzTolerance = tolerance;
        }

        public static List<(double mz, double intensity)> FindTheoreticalIsotopicPeaks(string sequence, int charge, MzRange range)
        {
            if(sequence.Contains("|"))
            {
                return null;
            }
            var peptide = new PeptideWithSetModifications(sequence, GlobalVariables.AllModsKnownDictionary);
            ChemicalFormula formula = peptide.FullChemicalFormula;
            var theoreticalDistribution = Chemistry.IsotopicDistribution.GetDistribution(formula);

            List<(double mz, double intensity)> theoreticalPeaks = new List<(double mz, double intensity)>();
            foreach (double mass in theoreticalDistribution.Masses)
            {
                double mz = mass.ToMz(charge);
                if (mz >= range.Minimum && mz <= range.Maximum)
                {
                    double intensity = theoreticalDistribution.Intensities[Array.IndexOf(theoreticalDistribution.Masses, mass)];
                    theoreticalPeaks.Add(new(mz, intensity));
                }
            }
            return theoreticalPeaks.Where(p => p.intensity >= 0.01).ToList();
        }

        public static List<List<(double mz, double intensity)>> FindTheoreticalMs1Peaks(List<string> sequences, int[] charges, MzRange range)
        {
            //find all theoretical peaks for the set of PSMs and store them in a list of list of peaks
            //outside list is for the sequence from each PSM, inside list is the peaks for each charge state of a specific sequence
            var allTheoreticalPeaks = new List<List<(double mz, double intensity)>>();
            foreach (var sequence in sequences)
            {
                foreach (int charge in charges)
                {
                    var theoreticalPeaks = FindTheoreticalIsotopicPeaks(sequence, charge, range);
                    if (theoreticalPeaks != null)
                    {
                        allTheoreticalPeaks.Add(theoreticalPeaks);
                    }
                }
            }
            return allTheoreticalPeaks;
        }

        public static List<(double experimentalMz, double theoreticalMz)> MatchedMzs(double[] experimentalMzs, double[] experimentalIntensities, double[] theoreticalMzs, Tolerance tolerance)
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
                    if (tolerance.Within(experimental[index].mz, theoreticalMz))
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
        public static List<double> FindMatchedIntensities (List<string> sequences, int[] charges, MzRange range, MzSpectrum experimentalSpectrum, Tolerance tolerance, out List<string> matchedSequences)
        {
            var allTheoreticalPeaks = FindTheoreticalMs1Peaks(sequences, charges , range);
            matchedSequences = new List<string>();
            var allseq = sequences.SelectMany(s => Enumerable.Repeat(s, charges.Length)).ToList();  
            List<int> indicesToRemove = new List<int>();
            for (int i = allTheoreticalPeaks.Count - 1; i >= 0; i--) 
            {
                var match = MatchedMzs(experimentalSpectrum.XArray, experimentalSpectrum.YArray, allTheoreticalPeaks[i].Select(p => p.mz).ToArray(), tolerance);
                if (match.All(m => m.experimentalMz == -1))
                {
                    allTheoreticalPeaks.RemoveAt(i);
                }
                else
                {
                    matchedSequences.Add(allseq[i]);
                }
            }
            double[] allTheoreticalMzs = allTheoreticalPeaks.SelectMany(peaks => peaks.Select(p => p.mz)).Distinct().OrderBy(mz => mz).ToArray();

            var mzPairs = MatchedMzs(experimentalSpectrum.XArray, experimentalSpectrum.YArray, allTheoreticalMzs, tolerance);
            var matchedMzs = mzPairs.OrderBy(pair => pair.theoreticalMz).Select(pair => pair.experimentalMz);

            List<double> matchedIntensities = new List<double>();

            foreach (double mz in matchedMzs)
            {
                if (mz >= 0)
                {
                    int index = experimentalSpectrum.XArray.IndexOf(mz);
                    matchedIntensities.Add(experimentalSpectrum.YArray[index]);
                }
                else
                {
                    matchedIntensities.Add(0);
                }
            }

            return matchedIntensities;
        }

        public static double FindFractionOfMatchedIntensities(List<string> sequences, int[] charges, MzRange range, MzSpectrum experimentalSpectrum, Tolerance tolerance, out List<string> matchedSequences)
        {
            var matchedIntensities = FindMatchedIntensities(sequences, charges, range, experimentalSpectrum, tolerance, out List<string> PSMsequences);
            double fractionMatched = matchedIntensities.Sum() / experimentalSpectrum.SumOfAllY;
            matchedSequences = PSMsequences;

            return fractionMatched;
        }

        public static MzSpectrum GetTheoreticalMs1Spectrum(List<string> sequences, int[] charges, MzRange range, MzSpectrum experimentalSpectrum, Tolerance tolerance, out double[] coefficients)
        {
            var allTheoreticalPeaks = FindTheoreticalMs1Peaks(sequences, charges, range);
            List<int> indicesToRemove = new List<int>();
            for (int i = allTheoreticalPeaks.Count - 1; i >= 0; i--)
            {
                var match = MatchedMzs(experimentalSpectrum.XArray, experimentalSpectrum.YArray, allTheoreticalPeaks[i].Select(p => p.mz).ToArray(), tolerance);
                if (match.All(m => m.experimentalMz == -1))
                {
                    allTheoreticalPeaks.RemoveAt(i);
                }
            }
            double[] allTheoreticalMzs = allTheoreticalPeaks.SelectMany(peaks => peaks.Select(p => p.mz)).Distinct().OrderBy(mz => mz).ToArray();

            //make a feature matrix with each row one mz and each column as a peptide ion
            double[][] features = new double[allTheoreticalMzs.Length][];
            for (int i = 0; i < allTheoreticalMzs.Length; i++)
            {
                features[i] = new double[allTheoreticalPeaks.Count];
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

            double[] matchedIntensities = FindMatchedIntensities(sequences, charges, range, experimentalSpectrum, tolerance, out List<string> matchedSequences).ToArray();
            double[] matchedIntensities_normalized = matchedIntensities.Select(m => m/matchedIntensities.Sum()).ToArray();
            //make a multiple linear regression model to find the relative intensities for each peptide ion
            coefficients = MultipleRegression.QR(features, matchedIntensities, false).Select(c => Math.Round(c, 2)).ToArray();
            //var coef_normalized = MultipleRegression.QR(features, matchedIntensities_normalized, false).Select(c => Math.Round(c, 2)).ToArray();

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

        public static double? CalculateSimilarityScore (List<string> sequences, int[] charges, MzRange range, MzSpectrum precursorSpectrum, Tolerance tolerance, 
            SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks) 
        {
            double[] coefficients;
            MzSpectrum theoreticalSpectrum = GetTheoreticalMs1Spectrum(sequences, charges, range, precursorSpectrum, tolerance, out coefficients);
            var allTheoreticalMzs = theoreticalSpectrum.XArray;
            var mzPairs = MatchedMzs(precursorSpectrum.XArray, precursorSpectrum.YArray, allTheoreticalMzs, tolerance);
            var matchedMzsExperimental = mzPairs.Where(pair => pair.experimentalMz >= 0).Select(pair => pair.experimentalMz).OrderBy(mz => mz).ToArray();
            var matchedIntensitiesExperimental = FindMatchedIntensities(sequences, charges, range, precursorSpectrum, tolerance, out List<string> matchedSequences).Where(i => i > 0).ToArray();
            MzSpectrum experimentalSpectrum = new MzSpectrum(matchedMzsExperimental, matchedIntensitiesExperimental, true);

            SpectralSimilarity similarity = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, scheme, toleranceInPpm, allPeaks);

            return similarity.CosineSimilarity();
        }

        public static List<double> FindMzsForPrecursors(double precursorMz, int precursorCharge, int[] charges, MzRange range)
        {
            var theoreticalMzs = new List<double>();
            double precursorMass = precursorMz.ToMass(precursorCharge);
            foreach (int charge in charges)
            {
                double mz = precursorMass.ToMz(charge);
                if (mz >= range.Minimum && mz <= range.Maximum)
                {
                    theoreticalMzs.Add(mz);
                }
            }
            return theoreticalMzs;
            ;
        }

        public static List<double> FindMatchedIntensities2(List<string> sequences, int[] charges, MzRange range, MzSpectrum experimentalSpectrum, Tolerance tolerance, out List<string> matchedSequences, 
            out List<(double experimentalMz, double theoreticalMz)> mzPairs)
        {
            var allTheoreticalPeaks = FindTheoreticalMs1Peaks(sequences, charges, range);
            matchedSequences = new List<string>();
            var allseq = sequences.SelectMany(s => Enumerable.Repeat(s, charges.Length)).ToList();
            List<int> indicesToRemove = new List<int>();
            for (int i = allTheoreticalPeaks.Count - 1; i >= 0; i--)
            {
                var match = MatchedMzs(experimentalSpectrum.XArray, experimentalSpectrum.YArray, allTheoreticalPeaks[i].Select(p => p.mz).ToArray(), tolerance);
                if (match.All(m => m.experimentalMz == -1))
                {
                    allTheoreticalPeaks.RemoveAt(i);
                }
                else
                {
                    matchedSequences.Add(allseq[i]);
                }
            }
            double[] allTheoreticalMzs = allTheoreticalPeaks.SelectMany(peaks => peaks.Select(p => p.mz)).Distinct().OrderBy(mz => mz).ToArray();

            mzPairs = MatchedMzs(experimentalSpectrum.XArray, experimentalSpectrum.YArray, allTheoreticalMzs, tolerance);
            var matchedMzs = mzPairs.OrderBy(pair => pair.theoreticalMz).Select(pair => pair.experimentalMz);

            List<double> matchedIntensities = new List<double>();

            foreach (double mz in matchedMzs)
            {
                if (mz >= 0)
                {
                    int index = experimentalSpectrum.XArray.IndexOf(mz);
                    matchedIntensities.Add(experimentalSpectrum.YArray[index]);
                }
                else
                {
                    matchedIntensities.Add(0);
                }
            }

            return matchedIntensities;
        }


    }
}
