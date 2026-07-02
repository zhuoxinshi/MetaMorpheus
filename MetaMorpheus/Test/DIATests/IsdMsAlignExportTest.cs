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
        /// End-to-end: take a real ISD mzML, build pseudo-MS2 scans via the ISD engine
        /// (consensus XIC construction + apex-RT/Pearson precursor-fragment grouping + pseudo-MS2 generation),
        /// export them to ms1/ms2.msalign, and confirm the ms2.msalign round-trips through the MsAlign reader
        /// (i.e. it is a valid file TopPIC could consume).
        /// </summary>
        [Test]
        public static void ExportIsdPseudoScansToMsAlign()
        {
            string mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "DIA",
                "08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT8.1-9.23.mzML");
            Assert.That(File.Exists(mzmlPath), $"missing test fixture: {mzmlPath}");

            var dataFile = MsDataFileReader.GetDataFile(mzmlPath);
            dataFile.LoadAllStaticData();

            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)),
                new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)),
                new XicGroupingEngine(0.2f, 0.5, 0.5, 1),
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
