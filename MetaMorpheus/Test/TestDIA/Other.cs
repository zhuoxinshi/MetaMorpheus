using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using EngineLayer;
using NUnit.Framework;
using System.IO;
using TaskLayer;
using EngineLayer.DIA;
using MzLibUtil;
using Plotly.NET;
using Nett;
using Omics.SpectrumMatch;
using MassSpectrometry.MzSpectra;
using Readers;

namespace Test.TestDIA
{
    public class Other
    {
        [Test]
        public static void SimilarityCalculation()
        {
            var libraryPath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\1614_HCDch2_decoys_predictions.msp";
            var neutralLossLibPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\spectralLibraryNeutralLossTest.msp");
            var pathList = new List<string> { libraryPath };
            var library = new SpectralLibrary(pathList);
            var librarySpectra = library.GetAllLibrarySpectra().ToList();

            var rawPath = @"E:\Aneuploidy\DDA\071525\07-15-25_1614-R1-Q_E1+5-calib.mzML";
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(rawPath, new CommonParameters());
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();

            var psmTsvPath = @"E:\Aneuploidy\DDA\071525\1614_E1-8_calied-generalGPTMD+1NAsub_noTrunc\Task2-SearchTask\Individual File Results\07-15-25_1614-R1-Q_E1+5-calib_PSMs.psmtsv";
            var allPsmTsv = SpectrumMatchTsvReader.ReadTsv(psmTsvPath, out List<string> warnings).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var allPsmTsv_decoy = SpectrumMatchTsvReader.ReadTsv(psmTsvPath, out List<string> warnings2).Where(p => p.DecoyContamTarget == "D").ToList();
            var allSequences = librarySpectra.Select(s => s.Sequence).ToList();
            var psmToLook = allPsmTsv_decoy.Where(p => allSequences.Contains(p.FullSequence)).ToList();
            var cosineSimilarity = new List<double>();
            foreach (var psmTsv in psmToLook)
            {
                if (library.TryGetSpectrum(psmTsv.FullSequence, psmTsv.PrecursorCharge, out LibrarySpectrum libSpectrum))
                {
                    var rawScan = ms2Scans.FirstOrDefault(s => s.OneBasedScanNumber == psmTsv.Ms2ScanNumber);
                    var similarity = new SpectralSimilarity(rawScan.MassSpectrum, libSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, 20, false);
                    cosineSimilarity.Add(similarity.CosineSimilarity().Value);
                }
            }
            var densityPlot = Chart2D.Chart.Histogram<double, string>(
                    cosineSimilarity.ToArray(), orientation: StyleParam.Orientation.Vertical,
                    HistNorm: StyleParam.HistNorm.ProbabilityDensity,Opacity: 0.6);
            densityPlot.Show();
        }
    }
}
