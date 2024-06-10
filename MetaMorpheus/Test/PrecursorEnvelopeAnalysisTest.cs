using EngineLayer.ClassicSearch;
using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;
using Chemistry;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using System.IO;
using Omics;
using Proteomics.AminoAcidPolymer;
using EngineLayer.FdrAnalysis;
using MassSpectrometry.MzSpectra;

namespace Test
{
    public static class PrecursorEnvelopeAnalysisTest
    {
        [Test]
        public static void RandomTest()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters();
            var myMsDataFile = myFileManager.LoadFile(filePath, commonParameters);
            SearchTask task = new SearchTask();
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(myDatabase, true, DecoyType.Reverse, false, out List<string> errors);
            var fsp = new List<(string, CommonParameters)>();
            fsp.Add(("SmallCalibratible_Yeast.mzML", commonParameters));
            var arrayOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, filePath, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            var variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => commonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var fixedModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => commonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            bool writeSpectralLibrary = false;
            MassDiffAcceptor massDiffAcceptor = new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, new PpmTolerance(5));
            SpectralMatch[] allPsmsArray = new SpectralMatch[arrayOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, arrayOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, massDiffAcceptor, commonParameters, fsp, null, new List<string>(), writeSpectralLibrary).Run();

            List<SpectralMatch> psms = new List<SpectralMatch>();
            foreach (SpectralMatch psm in allPsmsArray)
            {
                if (psm != null)
                {
                    psms.Add(psm);
                }
            }
            SpectralMatch psmScan23 = psms.ToArray()[33];
            ChemicalFormula formula = new Peptide(psmScan23.FullSequence).GetChemicalFormula();
            IsotopicDistribution theoreticalDistribution = IsotopicDistribution.GetDistribution(formula);
            var sum = theoreticalDistribution.Intensities.Sum();
        }

        [Test]
        public static void TestOnSimpleSpectrum()
        {
            List<string> sequences = new List<string>() 
            {"KYDNSLKIISNASCTTNCLAPLA", "QYDNSLKIISNASCTTNCLAPLA", "KYDNSLKIISNASCTTNCLAPNA" };
            MzRange range = new MzRange(0, 1500);
            double[] myRatio = new double[] { 1, 2, 3 };

            int[] charges = new int[] { 3 };
            var peaks = PrecursorEnvelopeAnalysis.FindTheoreticalMs1Peaks(sequences, charges, range);
            var allPeaks = new List<(double mz, double intensity)>();

            for (int i = 0; i < sequences.Count; i++)
            {
                var modifiedPeaks = peaks[i].Where(p => p.intensity >= 0.0001).Select(p => (p.mz, p.intensity * myRatio[i])).ToList();
                allPeaks.AddRange(modifiedPeaks);
            }

            var sortedPeaks = allPeaks.OrderBy(p => p.mz);
            double[] allMzs = sortedPeaks.Select(p => p.mz).ToArray();
            double[] allIntensities = sortedPeaks.Select(p => p.intensity).ToArray();
            MzSpectrum mySpectrum = new MzSpectrum(allMzs, allIntensities, true);
            double tolerance = 0.0001;

            MzSpectrum resultSpectrum = PrecursorEnvelopeAnalysis.GetTheoreticalMs1Spectrum
                (sequences, charges, range, mySpectrum, tolerance, out double[] coefficients);
            Assert.AreEqual(mySpectrum.XArray, resultSpectrum.XArray);
            for (int i = 0; i < mySpectrum.YArray.Length; i++)
            {
                Assert.IsTrue(Math.Abs(mySpectrum.YArray[i] - resultSpectrum.YArray[i]) < 0.0001);
            }

            //Test for FindFractionOfMatchedIntensities
            for (int i = 0; i < myRatio.Length; i++)
                Assert.That(myRatio[i], Is.EqualTo(coefficients[i]).Within(tolerance));
            
            double fraction = PrecursorEnvelopeAnalysis.FindFractionOfMatchedIntensities(sequences, charges, range, mySpectrum, tolerance);
            Assert.AreEqual(1, fraction);
            var similarityScore = PrecursorEnvelopeAnalysis.CalculateSimilarityScore(sequences, charges, range, mySpectrum, tolerance, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, 0.01, false);
            Assert.AreEqual(1, similarityScore);
        }

        [Test]
        public static void TestOnRealData()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters();
            var myMsDataFile = myFileManager.LoadFile(filePath, commonParameters);
            SearchTask task = new SearchTask();
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(myDatabase, true, DecoyType.Reverse, false, out List<string> errors);
            var fsp = new List<(string, CommonParameters)>();
            fsp.Add(("SmallCalibratible_Yeast.mzML", commonParameters));
        }
    }
}
