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
using Nett;
using NUnit.Framework.Constraints;
using System.Runtime.InteropServices;
using Microsoft.VisualStudio.TestPlatform.ObjectModel;

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
            var tolerance = new PpmTolerance(5);

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

            //Test for charge states that do not exist
            MzSpectrum resultSpectrum2 = PrecursorEnvelopeAnalysis.GetTheoreticalMs1Spectrum
                (sequences, new int[] {2,3}, range, mySpectrum, tolerance, out double[] coefficients2);
            var similarityScore2 = PrecursorEnvelopeAnalysis.CalculateSimilarityScore(sequences, new int[] { 2, 3 }, range, mySpectrum, tolerance, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, 0.01, false);
            Assert.AreEqual(1, similarityScore);
        }

        [Test]
        public static void TestOnRealData()
        {
            string tomlFilePath = @"\\bison.chem.wisc.edu\share\Users\Nic\Chimeras\Mann_11cell_analysis\Jurkat\SearchResults\MetaMorpheusWithLibrary\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFilePath, MetaMorpheusTask.tomlConfig);
            string outputPath = @"E:\Test Data";
            task.SearchParameters.DoLabelFreeQuantification = true;
            string myDatabase = @"\\bison.chem.wisc.edu\share\Users\Nic\Chimeras\Mann_11cell_analysis\Jurkat\human_UP000005640_reviewedGPTMD.xml";
            //List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(myDatabase, true, DecoyType.Reverse, false, out List<string> errors);
            List<string> dataFilePaths = new List<string>
                {
                @"E:\Test Data\20100614_Velos1_TaGe_SA_Jurkat_1-calib-averaged.mzML"
                };
            task.RunTask(outputPath, new List<DbForTask> { new DbForTask(myDatabase, false) }, dataFilePaths, "testOnRealData");
            var allPsms = PsmTsvReader.ReadTsv(Path.Combine(outputPath, @"AllPSMs.psmtsv"), out var warnings);

            //27534
            var psms = allPsms.GroupBy(p => new { p.PrecursorScanNum, p.FileNameWithoutExtension }).Where(g => g.ToList().First().PrecursorScanNum == 27534).SelectMany(p => p).ToList();
            MyFileManager myFileManager = new MyFileManager(true);
            MsDataFile scans = myFileManager.LoadFile(dataFilePaths[0], task.CommonParameters);
            var experimentalMs1 = scans.Scans.Where(scan => scan.OneBasedPrecursorScanNumber == 27534).First().MassSpectrum;
            List<string> sequences = psms.Select(p => p.FullSequence).ToList();
            int[] charges = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            var range = scans.Scans.Where(scan => scan.OneBasedPrecursorScanNumber == 27534).First().ScanWindowRange;
            var tolerance = task.CommonParameters.PrecursorMassTolerance;
            MzSpectrum resultSpectrum = PrecursorEnvelopeAnalysis.GetTheoreticalMs1Spectrum
            (sequences, charges, range, experimentalMs1, tolerance, out double[] coefficients);
            double fraction = PrecursorEnvelopeAnalysis.FindFractionOfMatchedIntensities(sequences, charges, range, experimentalMs1, tolerance);
            var similarityScore = PrecursorEnvelopeAnalysis.CalculateSimilarityScore(sequences, charges, range, experimentalMs1, tolerance, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, tolerance.Value, false);
        }

        [Test]
        public static void TestOnSimpleSpectrum2()
        {

        }
    }
}
