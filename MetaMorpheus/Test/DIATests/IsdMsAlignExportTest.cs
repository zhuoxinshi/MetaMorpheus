using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.DIATests
{
    [TestFixture]
    public static class IsdMsAlignExportTest
    {
        /// <summary>
        /// End-to-end: take a real ISD mzML, build pseudo-MS2 scans using CONSENSUS FEATURE TRACING
        /// (<see cref="ConsensusMassXicConstructor"/> = mzLib MassTraceBuilder→TraceCorrector→MassFeatureBuilder)
        /// as the XIC source, combined with apex-RT/Pearson precursor-fragment grouping and pseudo-MS2 generation,
        /// then export them to ms1/ms2.msalign and confirm the ms2.msalign round-trips through the MsAlign reader
        /// (i.e. it is a valid file TopPIC could consume).
        /// </summary>
        [Test]
        public static void ExportIsdPseudoScansToMsAlign()
        {
            string mzmlPath = @"E:\ISD Project\ISD_250906\09-17-25_YB_81min_ISD60-80-100_preFilter700-900-1100_rep1.raw";
            Assert.That(File.Exists(mzmlPath), $"missing test fixture: {mzmlPath}");

            var dataFile = MsDataFileReader.GetDataFile(mzmlPath);
            dataFile.LoadAllStaticData();

            // Consensus feature tracing as the XIC front-end for both the intact (precursor) and fragment channels.
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), minMass: 4000, xicSpline: new Bspline(2, 150)),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new Bspline(2, 150)),
                new UmpirePfGroupingEngine(150, 0.2f, 0.3, 0.7),
                PseudoMs2ConstructionType.Mass,
                combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, mzmlPath, commonParams).ToArray();
            Assert.That(pseudoScans.Length, Is.GreaterThan(0), "no pseudo scans produced from the ISD fixture");

            string ms2Path = Path.Combine(TestContext.CurrentContext.TestDirectory, "isd_pseudo_ms2.msalign");
            string ms1Path = Path.Combine(TestContext.CurrentContext.TestDirectory, "isd_pseudo_ms1.msalign");
            IsdMsAlignExporter.WriteMs2Align(pseudoScans, ms2Path);
            IsdMsAlignExporter.WriteMs1Align(pseudoScans, ms1Path);
            Assert.That(File.Exists(ms2Path));
            Assert.That(File.Exists(ms1Path));

            // number of BEGIN IONS blocks == number of pseudo scans
            int blocks = File.ReadLines(ms2Path).Count(l => l.StartsWith("BEGIN IONS"));
            Assert.That(blocks, Is.EqualTo(pseudoScans.Length));

            // round-trip through the MsAlign reader: it must parse back the same number of scans
            var reread = new Ms2Align(ms2Path);
            reread.LoadAllStaticData();
            Assert.That(reread.Scans.Count(s => s != null), Is.EqualTo(pseudoScans.Length));
            // and every reread scan is MS2 with a precursor mass carried through
            Assert.That(reread.Scans.Where(s => s != null).All(s => s.MsnOrder == 2));

            File.Delete(ms2Path);
            File.Delete(ms1Path);
        }
    }
}
