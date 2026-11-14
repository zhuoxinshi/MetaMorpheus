using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.DIATests
{
    /// <summary>
    /// Generates ISD pseudo MS2 spectra (msAlign) from instrument .raw files and TopFD feature files. All inputs and
    /// settings are set in <see cref="WriteHcdYbcIsdPseudoMs2"/>; as written they reproduce the pseudo MS2 spectra of
    /// the HCD_YBC benchmark (the six HCD-YBC ISD60-80-100 runs).
    ///
    /// Inputs, for each run (stem = .raw file name without extension):
    ///   {raw folder}\{stem}.raw
    ///   {TopFD folder}\{stem}_ms1_feature.xml                     precursor channel (sid=15)
    ///   {TopFD folder}\{stem}_isd60/_isd80/_isd100_feature.xml    fragment channels (sid=60, 80, 100)
    /// The feature files are from TopFD 1.8.1, run on per-channel mzML files that keep the .raw scan numbers
    /// (mzLib branch ISD_pub, LabelCorrection.SeparateFilesForTopFD).
    ///
    /// Output: {output folder}\{stem}_pseudo_ms2.msalign for each run (with an empty {stem}_pseudo_ms1.msalign beside
    /// it, which TopPIC's combined-file mode expects) and DIAparameters.txt, a record of the settings.
    /// Explicit: needs the local data files.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class IsdPseudoMs2Reproduction
    {
        [Test]
        [Explicit("needs the local HCD_YBC .raw files and TopFD feature files")]
        public static void WriteHcdYbcIsdPseudoMs2()
        {
            // ---- Inputs and output ----
            string rawFolder = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
            string topFdFolder = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\Reproduce_test\SeparateISD";
            string outputFolder = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\Reproduce_test\PseudoMs2";
            var rawFiles = Directory.GetFiles(rawFolder, "*.raw")
                .Where(f => Path.GetFileName(f).Contains("ISD60-80-100") && Path.GetFileName(f).Contains("HCD-YBC"))
                .OrderBy(f => f)
                .ToList();
            Assert.That(rawFiles, Has.Count.EqualTo(6));

            // ---- Settings shared by all runs ----
            // XIC spline: cubic, on the scan-index axis, sampled every 0.05 scan, with one zero-intensity point added
            // one scan beyond each end. Pearson correlations are computed on the splined profiles.
            var xicSpline = new XicCubicSpline(splineRtInterval: 0.05, numberOfPeaksToAdd: 1, gap: 1, scanIndexBased: true);

            // Grouping: a fragment XIC is paired with a precursor XIC when their apexes are at most 2 scans
            // (acquisition cycles) apart and their Pearson correlation is at least 0.5; no overlap requirement.
            // Each precursor keeps its 100 best-correlated fragments and needs more than 10 to give a pseudo MS2 scan.
            var groupingEngine = new XicGroupingEngine(apexRTTolerance: 2, overlapThreshold: 0, correlationThreshold: 0.5,
                maxThreadsForGrouping: 10, minFragmentCountForGrouping: 10, precursorRankThreshold: null, fragmentRankThreshold: 100,
                apexToleranceUnit: ApexToleranceUnit.ScanIndex);

            foreach (var rawFile in rawFiles)
            {
                string stem = Path.GetFileNameWithoutExtension(rawFile);

                // ---- Settings of each run (the XIC constructors read that run's TopFD feature files) ----
                // Precursor XICs: TopFD features of the MS1 channel, 4,000 to 70,000 Da, observed in at least 3 scans
                // and spanning at most 2 min.
                var ms1XicConstructor = new FeatureXmlXicConstructor(Path.Combine(topFdFolder, stem + "_ms1_feature.xml"),
                    minNumberOfPeaks: 3, minMass: 4000, maxMass: 70000, minEcScore: 0, maxXicWidth: 2, xicSpline: xicSpline);

                // Fragment XICs: TopFD features of each fragment channel (60, 80, 100 V), keyed by source voltage, below
                // 70,000 Da, with the same limits.
                var fragmentFeatureXmlByVoltage = new[] { 60, 80, 100 }
                    .ToDictionary(voltage => voltage, voltage => Path.Combine(topFdFolder, $"{stem}_isd{voltage}_feature.xml"));
                var ms2XicConstructor = FeatureXmlXicConstructor.BySourceVoltage(fragmentFeatureXmlByVoltage,
                    minNumberOfPeaks: 3, minMass: 0, maxMass: 70000, minEcScore: 0, maxXicWidth: 2, xicSpline: xicSpline);

                // Fragment XICs of the three voltages are pooled before grouping. Each pseudo MS2 peak is a fragment's
                // mass with its summed XIC intensity.
                var diaParameters = new DIAparameters(DIAanalysisType.ISD, ms1XicConstructor, ms2XicConstructor, groupingEngine,
                    PseudoMs2ConstructionType.Mass, combineFragments: true, pseudoMs2IntensityType: PseudoMs2IntensityType.XicSum);

                WritePseudoMs2MsAlign(rawFile, diaParameters, outputFolder);
            }
        }

        /// <summary>
        /// Builds the pseudo MS2 scans of one ISD run with <paramref name="diaParameters"/> and writes
        /// {outputFolder}\{stem}_pseudo_ms2.msalign, an empty {stem}_pseudo_ms1.msalign and DIAparameters.txt.
        /// </summary>
        public static void WritePseudoMs2MsAlign(string rawFile, DIAparameters diaParameters, string outputFolder)
        {
            string stem = Path.GetFileNameWithoutExtension(rawFile);
            var dataFile = MsDataFileReader.GetDataFile(rawFile).LoadAllStaticData();
            var commonParameters = new CommonParameters { DIAparameters = diaParameters };
            var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParameters).ToArray();

            Directory.CreateDirectory(outputFolder);
            ISDEngine.WriteMsAlignFile(Path.Combine(outputFolder, stem + "_pseudo_ms2.msalign"), pseudoScans);
            File.WriteAllText(Path.Combine(outputFolder, stem + "_pseudo_ms1.msalign"), string.Empty);
            diaParameters.WriteSettingsToFile(Path.Combine(outputFolder, "DIAparameters.txt"));
        }
    }
}
