using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
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
            string mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "DIA",
                "08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT8.1-9.23.mzML");
            Assert.That(File.Exists(mzmlPath), $"missing test fixture: {mzmlPath}");

            var dataFile = MsDataFileReader.GetDataFile(mzmlPath);
            dataFile.LoadAllStaticData();

            // Consensus feature tracing as the XIC front-end for both the intact (precursor) and fragment channels.
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)),
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

        /// <summary>
        /// Environment-driven harness: generate consensus-path pseudo-MS2 MGFs at several precursor-fragment
        /// correlation thresholds, reusing one consensus-XIC build. Set ISD_MZML (input), ISD_OUTDIR (output),
        /// and optionally ISD_CORRS (comma list, default "0,0.4,0.6,0.7"). Writes consensus_corr&lt;t&gt;.mgf per
        /// threshold; each MGF is then searched externally with the top-down td_pseudoMS2 config to count IDs.
        /// Ignored unless ISD_MZML is set, so it never runs in normal CI.
        /// </summary>
        [Test]
        public static void ConsensusIdSweep_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("ISD_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set ISD_MZML to run this harness"); return; }
            string outDir = Environment.GetEnvironmentVariable("ISD_OUTDIR") ?? TestContext.CurrentContext.TestDirectory;
            Directory.CreateDirectory(outDir);
            var corrs = (Environment.GetEnvironmentVariable("ISD_CORRS") ?? "0,0.4,0.6,0.7")
                .Split(',').Select(s => double.Parse(s, System.Globalization.CultureInfo.InvariantCulture)).ToArray();
            float apex = 0.2f;
            double overlap = 0.5;
            var decon = new ClassicDeconvolutionParameters(1, 30, 4, 3);

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // consensus feature tracing -> XICs, built ONCE and reused across thresholds
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon);
            var ms2Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon);
            var ms1Xics = ms1Ctor.GetAllXics(ms1Scans);
            var ms2Xics = isdMap.Values.SelectMany(v => ms2Ctor.GetAllXics(v.ToArray())).ToList();
            TestContext.WriteLine($"consensus XICs: {ms1Xics.Count} precursor, {ms2Xics.Count} fragment");

            foreach (var corr in corrs)
            {
                var grouper = new XicGroupingEngine(apex, overlap, corr, Environment.ProcessorCount);
                var groups = grouper.PrecursorFragmentGrouping(ms1Xics, ms2Xics).ToList();
                int idx = 0;
                var scans = new List<Ms2ScanWithSpecificMass>();
                foreach (var g in groups)
                {
                    idx++;
                    g.PFgroupIndex = idx;
                    scans.Add(g.GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.Mass, new CommonParameters(), mzml));
                }
                string outMgf = Path.Combine(outDir, $"consensus_corr{corr.ToString(System.Globalization.CultureInfo.InvariantCulture)}.mgf");
                IsdMsAlignExporter.WriteMgf(scans, outMgf);
                TestContext.WriteLine($"corr={corr}: {scans.Count} pseudo scans -> {outMgf}");
            }
        }
    }
}
