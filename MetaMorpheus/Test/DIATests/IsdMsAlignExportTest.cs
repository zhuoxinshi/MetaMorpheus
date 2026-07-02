using Chemistry;
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

            // consensus feature tracing -> XICs, built ONCE and reused across thresholds.
            // PRECURSOR channel: keep only intact proteoforms (>=3 kDa, >=3 charge states) so small noise
            // features do not become precursors. FRAGMENT channel: unfiltered (fragments are low mass/charge).
            double precMinMass = double.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINMASS") ?? "3000",
                System.Globalization.CultureInfo.InvariantCulture);
            int precMinCharge = int.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINCHARGE") ?? "3");
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
                minMass: precMinMass, minChargeCount: precMinCharge);
            var ms2Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon);
            // ISD_FRAG_AGG=false -> fragments use mass tracing only (one XIC per charge-locked trace, no
            // cross-charge aggregation); precursors always use full feature tracing.
            bool fragAgg = (Environment.GetEnvironmentVariable("ISD_FRAG_AGG") ?? "true").ToLowerInvariant() != "false";
            ms2Ctor.AggregateCharges = fragAgg;
            var ms1Xics = ms1Ctor.GetAllXics(ms1Scans);
            var ms2Xics = isdMap.Values.SelectMany(v => ms2Ctor.GetAllXics(v.ToArray())).ToList();
            TestContext.WriteLine($"consensus XICs: {ms1Xics.Count} precursor (>= {precMinMass} Da, >= {precMinCharge} charges), {ms2Xics.Count} fragment (chargeAggregation={fragAgg})");

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

        /// <summary>
        /// Search DDA runs using CONSENSUS MASS TRACING for the precursors: trace the DDA MS1 scans with
        /// <see cref="ConsensusMassXicConstructor"/> to get clean intact precursor features, pair each real DDA
        /// MS2 scan to the consensus feature matching its isolation m/z + retention time, and build a pseudo scan
        /// from the consensus precursor + the scan's REAL (deconvoluted) fragments. Writes an MGF that is then
        /// searched with the same td_pseudoMS2 config as the ISD path (apples-to-apples ISD-consensus vs
        /// DDA-consensus). Env-driven: DDA_MZML (input), DDA_OUT (output MGF). Ignored unless DDA_MZML is set.
        /// </summary>
        [Test]
        public static void SearchDdaWithConsensusPrecursors_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("DDA_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set DDA_MZML to run this harness"); return; }
            string outMgf = Environment.GetEnvironmentVariable("DDA_OUT")
                ?? Path.Combine(TestContext.CurrentContext.TestDirectory, "dda_consensus.mgf");
            const double proton = 1.007276466879;
            var decon = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            // real DDA MS2 fragments are multiply charged -> deconvolute them to neutral masses
            var fragParams = new CommonParameters(productDeconParams: new ClassicDeconvolutionParameters(1, 20, 4, 3));

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2Scans = allScans.Where(s => s.MsnOrder == 2).ToArray();

            // consensus mass tracing on the MS1 scans -> intact precursor features (>=3 kDa, >=3 charges)
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
                minMass: 3000, minChargeCount: 3);
            var precursorXics = ms1Ctor.GetAllXics(ms1Scans);

            // pair each real DDA MS2 scan with the consensus feature matching its isolation m/z + RT
            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            int assigned = 0;
            foreach (var ms2 in ms2Scans)
            {
                if (ms2.IsolationMz == null) continue;
                double isoMz = ms2.IsolationMz.Value;
                double rt = ms2.RetentionTime;
                ExtractedIonChromatogram best = null; double bestErr = double.MaxValue; int bestZ = 1;
                foreach (var xic in precursorXics)
                {
                    if (rt < xic.StartRT - 0.5 || rt > xic.EndRT + 0.5) continue;
                    double mass = xic.ApexPeak.M;
                    int z = (int)Math.Round(mass / (isoMz - proton)); // charge that places this mass at the iso m/z
                    if (z < 1 || z > 60) continue;
                    double err = Math.Abs(mass.ToMz(z) - isoMz);
                    if (err < bestErr && err < 2.5) { bestErr = err; best = xic; bestZ = z; }
                }
                if (best == null) continue;
                assigned++;
                double precMz = ((double)best.ApexPeak.M).ToMz(bestZ);
                var frags = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2, fragParams);
                pseudoScans.Add(new Ms2ScanWithSpecificMass(ms2, precMz, bestZ, mzml, fragParams, frags));
            }
            IsdMsAlignExporter.WriteMgf(pseudoScans, outMgf);
            TestContext.WriteLine($"DDA consensus: {precursorXics.Count} MS1 features, {ms2Scans.Length} MS2, {assigned} assigned -> {outMgf}");
            Assert.That(assigned, Is.GreaterThan(0), "no MS2 scans matched a consensus precursor feature");
        }
    }
}
