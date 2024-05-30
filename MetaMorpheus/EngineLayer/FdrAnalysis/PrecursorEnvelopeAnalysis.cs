
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
        
        public static double[] FindRelativeIntensities(List<SpectralMatch> psms, int maxCharge, MsDataFile dataFile, double tolerance)
        {
            //find all theoretical peaks for the set of PSMs and store them in a list of list of peaks
            var allTheoreticalPeaks = new List<List<(double mz, double intensity)>>();
            foreach (var psm in psms)
            {
                for (int charge = 1; charge <= maxCharge; charge++)
                {
                    var theoreticalPeaks = FindTheoreticalIsotopicPeaks(psm, charge);
                    allTheoreticalPeaks.Add(theoreticalPeaks);
                }
            }
            double[] allTheoreticalMzs = allTheoreticalPeaks.SelectMany(peaks => peaks.Select(p => p.mz)).Distinct().OrderBy(mz => mz).ToArray();
            double[][] features = new double[allTheoreticalMzs.Length][];
            
            //make a matrix for theoretical mz and all peptide ions
            for (int i = 0; i < allTheoreticalMzs.Length; i++)
            {
                for (int j = 0; j < allTheoreticalPeaks.Count; j++)
                {
                    if (allTheoreticalPeaks[j].Select(p => p.mz).Contains(allTheoreticalMzs[i]))
                    {
                        int indexOfMz = allTheoreticalPeaks[j].Select(p => p.mz).ToList().IndexOf(allTheoreticalMzs[i]);
                        features[i][j] = allTheoreticalPeaks[j][indexOfMz].intensity;
                    } else
                    {
                        features[i][j] = 0;
                    }
                }
            }

            //match experimental mz with theoretical mz
            //what if a theoretical peak is matched to an experimental peak but there is another theoretical peak closer to that experimental peak?
            double[] matchedIntensities = new double[allTheoreticalMzs.Length];
            List<double> matchedIndex = new List<double>();
            var precursorSpectrum = dataFile.Scans.Where(scan => scan.OneBasedScanNumber == psms.First().PrecursorScanNumber).First().MassSpectrum;
            for (int i = 0; i < allTheoreticalMzs.Length; i++)
            {
                int closestTheoreticalPeakIndex = precursorSpectrum.GetClosestPeakIndex(allTheoreticalMzs[i]);
                if (Math.Abs(allTheoreticalMzs[i] - precursorSpectrum.XArray[closestTheoreticalPeakIndex]) <= tolerance && matchedIndex.Contains(closestTheoreticalPeakIndex))
                {
                    matchedIntensities[i] = precursorSpectrum.YArray[closestTheoreticalPeakIndex];
                    matchedIndex.Add(closestTheoreticalPeakIndex);
                } else
                {
                    matchedIntensities[i] = 0;
                }
            }

            var ratios = MultipleRegression.QR(features, matchedIntensities, false);
            return ratios;
        }

        /*
        public static List<(double mz, double intensity)> FindAllMs1TheoreticalPeaks(List<SpectralMatch> psms)
        {
            List<double> allIntensities = psms.Select(p => p.PrecursorScanIntensity).ToList();
            double[] relativeIntensities = allIntensities.Select(i => i / allIntensities.Sum()).ToArray();
            List<(double mz, double intensity)> allPeaks = new List<(double mz, double intensity)>();
            SpectralMatch[] psmsArray = psms.ToArray();

            for (int psmIndex = 0; psmIndex < psms.Count; psmIndex++)
            {
                List<(double mz, double intensity)> theoreticalPeaks = FindTheoreticalIsotopicPeaks(psmsArray[psmIndex]);
                List<(double mz, double intensity)> normalizedPeaks = theoreticalPeaks.Select(p => (p.mz, p.intensity * relativeIntensities[psmIndex])).ToList();
                allPeaks.AddRange(normalizedPeaks);
            }

            //Can use matrices to calculate the theoretical intensities
            Dictionary<double, double> dict = new Dictionary<double, double>();

            foreach(var peak in allPeaks)
            {
                if (dict.ContainsKey(peak.mz))
                {
                    dict[peak.mz] += peak.intensity;
                }
                else
                {
                    dict[peak.mz] = peak.intensity;
                }
            }

            List <(double mz, double intensity)> normalizedAllPeaks = dict.Select(p => (p.Key, p.Value)).ToList();

            return normalizedAllPeaks;
        }
        */

        public static double CalculatePrecursorScanSimilarityScore(MsDataFile dataScans, List<SpectralMatch> psms)
        {
            return 1;
        }

        public void RunAnalysis()
        {
            var psmSets = AllPsms.GroupBy(p => p.ScanNumber);
            foreach (var psmSet in psmSets)
            {
                var psmList = psmSet.ToList();
                double similarityScore = CalculatePrecursorScanSimilarityScore(RawScans, psmList);
                //pmList.ForEach(psm => psm.PrecursorScanSmilarityScore = similarityScore);
            }
        }
    }
}
