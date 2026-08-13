using EngineLayer.DIA;
using EngineLayer.DIA.FeatureAlign;
using iText.Forms.Form.Element;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    internal class FeatureQuant
    {
        /// <summary>
        /// Manual runner (not an assertion test): aligns the feature tables listed below with
        /// <see cref="FeatureAligner"/> and writes the consensus tsv. Marked [Explicit] so it is skipped by the
        /// normal test run and only executes when you select it directly (right-click -> Run in the IDE, or
        /// dotnet test --filter). Edit the input files, output path, and tolerances below.
        /// </summary>
        [Test, Explicit("manual runner: aligns feature tables and writes a tsv")]
        public void RunFeatureAligner()
        {
            // ---------- EDIT: input feature tables (FLASHQuant or FLASHDeconv .tsv) ----------
            var dir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HC_YC";
            // filter on the FILE NAME (not the full path, which contains the folder "HC_YC" and would match YC too)
            var files = Directory.GetFiles(dir, "*ms1.feature").Where(f => f.Contains("DDA") & !f.Contains("_YC_")).ToList();

            double massTol    = 20;
            var    massUnit   = MassToleranceUnit.Ppm;   // or MassToleranceUnit.Dalton
            double rtTol      = 120;                      // seconds, +/-

            var outFolder = @"E:\ISD Project\071626_ISD\Quant\HC_YC\TopFD";
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);
            string outputPath = Path.Combine(outFolder, $"FeatureAlign_HC_DDA_11files_{massTol}_{massUnit}_{rtTol}.tsv");

            string written = FeatureAligner.Align(files, outputPath, massTol, massUnit, rtTol, grouping: GroupingMethod.ReferenceAnchor, positiveNotch:2, negativeNotch:2);

            TestContext.Out.WriteLine("FeatureAligner wrote: " + written);
            TestContext.Out.WriteLine("Rows: " + (File.ReadAllLines(written).Length - 1));
        }

        /// <summary>
        /// Manual runner (not an assertion test): cross-run feature alignment of the FLASHDeconv-on-averaged-mzML
        /// outputs for the 11 files whose name contains "DDA" and "HC" (9 HC-YC mix + 2 pure-HC; excludes the 2 YC).
        /// FLASHDeconv's -out tsv mixes MS1/MS2/MS3 features (ini max_MS_level=3); intact-proteoform alignment uses
        /// MS1 only, so this writes MS1-filtered copies first, then aligns them with <see cref="FeatureAligner"/>.
        /// [Explicit] so it runs only when selected.
        /// </summary>
        [Test, Explicit("manual runner: FLASHDeconv averaged MS1 feature alignment, 11 DDA+HC files")]
        public void RunAveragedHCFeatureAlignment()
        {
            var dir = @"E:\ISD Project\071626_ISD\Quant\HC_YC\FlashDeconv\Averaged";
            double massTol = 20; var massUnit = MassToleranceUnit.Ppm; double rtTol = 120;   // ppm / seconds
            int posNotch = 2, negNotch = 2;                                                  // +/-2 monoisotope notches

            // 11 files: file name contains "DDA" and "HC" (HC-YC mixes + pure-HC; the pure-YC files lack "HC").
            var selected = Directory.GetFiles(dir, "*.tsv")
                .Where(f => { var n = Path.GetFileName(f); return n.Contains("_DDA_") && n.Contains("HC"); })
                .OrderBy(f => f).ToList();
            TestContext.Out.WriteLine($"selected {selected.Count} files:");
            foreach (var f in selected) TestContext.Out.WriteLine("  " + Path.GetFileName(f));

            // MS1-only filtered copies (keep header + rows where MSLevel==1).
            var ms1dir = Path.Combine(dir, "ms1_only");
            Directory.CreateDirectory(ms1dir);
            var ms1files = new List<string>();
            foreach (var f in selected)
            {
                var lines = File.ReadAllLines(f);
                if (lines.Length == 0) continue;
                var header = lines[0].Split('\t');
                int ms = System.Array.FindIndex(header, h => h.Trim().Equals("MSLevel", System.StringComparison.OrdinalIgnoreCase));
                var outF = Path.Combine(ms1dir, Path.GetFileName(f));
                using (var sw = new StreamWriter(outF, false))
                {
                    sw.WriteLine(lines[0]);
                    for (int i = 1; i < lines.Length; i++)
                    {
                        if (lines[i].Length == 0) continue;
                        var t = lines[i].Split('\t');
                        if (ms >= 0 && ms < t.Length && t[ms].Trim() == "1") sw.WriteLine(lines[i]);
                    }
                }
                ms1files.Add(outF);
            }

            string outputPath = Path.Combine(dir, "FeatureAlign_averaged_11_DDA_HC_MS1_20ppm_120s_notch2.tsv");
            FeatureAligner.Align(ms1files, outputPath, massTol, massUnit, rtTol,
                positiveNotch: posNotch, negativeNotch: negNotch, grouping: GroupingMethod.ReferenceAnchor);

            var l = File.ReadAllLines(outputPath);
            int nfCol = System.Array.IndexOf(l[0].Split('\t'), "n_files");
            var dist = l.Skip(1).Select(x => int.Parse(x.Split('\t')[nfCol]))
                        .GroupBy(x => x).OrderBy(g => g.Key).Select(g => $"{g.Key}:{g.Count()}");
            int shared = l.Skip(1).Count(x => int.Parse(x.Split('\t')[nfCol]) >= 2);
            int inAll  = l.Skip(1).Count(x => int.Parse(x.Split('\t')[nfCol]) == ms1files.Count);
            TestContext.Out.WriteLine($"consensus rows: {l.Length - 1}; shared(>=2): {shared}; in all {ms1files.Count}: {inAll}");
            TestContext.Out.WriteLine("n_files dist -> " + string.Join("  ", dist));
        }

        /// <summary>
        /// Manual runner (not an assertion test): aligns one feature table to a MetaMorpheus proteoform/PSM
        /// .psmtsv by mass + RT and writes one row per feature with its best-matching identification.
        /// [Explicit] so it runs only when you select it. Edit the paths / tolerances below.
        /// </summary>
        [Test, Explicit("manual runner: aligns features to identifications and writes a tsv")]
        public void RunFeatureToIdAlignment()
        {
            // ---------- EDIT: one feature table + one identification .psmtsv ----------
            // A FLASHQuant/FLASHDeconv .tsv OR a TopFD .feature file both work (columns auto-detected).
            string featureFile        = @"E:\ISD Project\071626_ISD\Quant\HC_YC\FlashDeconv\Ratios\07-23-26_HC-YC_1-1_81min_DDA_3iso_mscan4_rep1.tsv";
            // TopFD example: @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HC_YC\07-23-26_HC-YC_1-1_81min_DDA_3iso_mscan4_rep1_ms1.feature";
            string identificationFile = @"E:\ISD Project\071626_ISD\Quant\HC_YC\HC_YC_DDA_xml-cali-xml\Task1-SearchTask\Individual File Results\07-23-26_HC-YC_1-1_81min_DDA_3iso_mscan4_rep1_Proteoforms.psmtsv";
            string outputPath         = @"E:\ISD Project\071626_ISD\Quant\HC_YC\FlashDeconv\Ratios\07-23-26_HC-YC_1-1_rep1_features_to_ids.tsv";

            // ---------- EDIT: tolerances ----------
            double massTol      = 20;
            var    massUnit     = MassToleranceUnit.Ppm;   // or MassToleranceUnit.Dalton
            double rtTol        = 150;                      // seconds, +/-
            double qValueThresh = 0.01;                     // ignore IDs with QValue above this
            // -------------------------------------

            string written = FeatureAligner.AlignToIdentifications(
                featureFile, identificationFile, outputPath, massTol, massUnit, rtTol, qValueThresh);

            TestContext.Out.WriteLine("Feature->ID alignment wrote: " + written);
            TestContext.Out.WriteLine("Rows: " + (File.ReadAllLines(written).Length - 1));
        }

        /// <summary>
        /// Manual runner (not an assertion test): aligns TopPIC proteoform identifications across the DDA runs and
        /// reports how many proteoforms are shared. Produces FOUR outputs: two file groupings
        /// (the 9 HC-YC mix files + the 2 HC-only files; and the 9 mix files + the 2 YC-only files) × two grouping
        /// methods (accession+mass ppm, and exact proteoform string). Reads *_ms2_toppic_proteoform_single.tsv.
        /// [Explicit] so it runs only when selected.
        /// </summary>
        [Test, Explicit("manual runner: aligns TopPIC proteoforms across DDA runs (4 outputs)")]
        public void RunToppicProteoformAligner()
        {
            // ---------- EDIT: folder with TopPIC per-run outputs + where to write ----------
            var dir    = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HC_YC";
            var outDir = Path.Combine(dir, "proteoform_alignments");
            double massPpm = 10;   // MassPpm tolerance (matches the search mass tol)
            // -------------------------------------------------------------------------------

            var single = Directory.GetFiles(dir, "*_ms2_toppic_proteoform_single.tsv").ToList();
            // 9 HC-YC mix DDA files (ratios 1-1/1-2/1-3, rep1/2/3); 2 HC-only; 2 YC-only. Distinguish by name:
            // mix = "HC-YC_1-" (hyphen), HC-only = "_HC_1-3_", YC-only = "_YC_1-8_"; all DDA.
            bool IsDda(string f) => f.Contains("_DDA_");
            var mix = single.Where(f => IsDda(f) && Path.GetFileName(f).Contains("HC-YC_1-")).OrderBy(f => f).ToList();
            var hc  = single.Where(f => IsDda(f) && Path.GetFileName(f).Contains("_HC_1-3_")).OrderBy(f => f).ToList();
            var yc  = single.Where(f => IsDda(f) && Path.GetFileName(f).Contains("_YC_1-8_")).OrderBy(f => f).ToList();
            TestContext.Out.WriteLine($"mix={mix.Count}  HC-only={hc.Count}  YC-only={yc.Count}");

            var groupA = mix.Concat(hc).ToList();   // 9 + 2 HC  = 11
            var groupB = mix.Concat(yc).ToList();   // 9 + 2 YC  = 11

            void RunOne(string tag, List<string> files, ToppicMatchKey key)
            {
                if (files.Count < 2) { TestContext.Out.WriteLine($"[skip] {tag}: <2 files"); return; }
                string keyTag = key == ToppicMatchKey.MassPpm ? "massppm" : "exact";
                string outPath = Path.Combine(outDir, $"toppic_align_{tag}_{keyTag}.tsv");
                ToppicProteoformAligner.Align(files, outPath, key, massPpm);

                var lines = File.ReadAllLines(outPath);
                int nRunsCol = System.Array.IndexOf(lines[0].Split('\t'), "n_runs");
                var dist = lines.Skip(1).Select(l => int.Parse(l.Split('\t')[nRunsCol]))
                                .GroupBy(x => x).OrderBy(g => g.Key)
                                .Select(g => $"{g.Key}:{g.Count()}");
                int total = lines.Length - 1;
                int shared = lines.Skip(1).Count(l => int.Parse(l.Split('\t')[nRunsCol]) >= 2);
                int allRuns = files.Count;
                int inAll = lines.Skip(1).Count(l => int.Parse(l.Split('\t')[nRunsCol]) == allRuns);
                TestContext.Out.WriteLine($"{tag}/{keyTag}: {total} proteoforms, {shared} shared(>=2 runs), {inAll} in all {allRuns}. n_runs dist -> {string.Join("  ", dist)}");
            }

            RunOne("HC", groupA, ToppicMatchKey.MassPpm);
            RunOne("HC", groupA, ToppicMatchKey.ExactProteoform);
            RunOne("YC", groupB, ToppicMatchKey.MassPpm);
            RunOne("YC", groupB, ToppicMatchKey.ExactProteoform);
        }

        /// <summary>
        /// Manual runner (not an assertion test): aligns MetaMorpheus proteoform IDENTIFICATIONS across runs by
        /// identity (default: Accession + Base Sequence + mods ignoring position) and writes a proteoform × run
        /// matrix (per-run intensity, RT and q-value). Input is one proteoform *_Proteoforms.psmtsv per run, OR a
        /// single pooled AllPSMs.psmtsv (auto-split by its File Name column). [Explicit] so it runs only when
        /// selected. Edit the paths / options below.
        /// </summary>
        [Test, Explicit("manual runner: aligns proteoform IDs across runs and writes a tsv")]
        public void RunProteoformAligner()
        {
            // ---------- EDIT: per-run proteoform .psmtsv files (Individual File Results) ----------
            var dir = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\HC_YC_DDA_xml-cali-xml\Task3-SearchTask\Individual File Results";
            var files = Directory.GetFiles(dir, "*_Proteoforms.psmtsv").ToList();
            // Single pooled file instead? e.g. AllPSMs.psmtsv (auto-splits by File Name):
            // var files = new List<string> { @"...\Task3-SearchTask\AllPSMs.psmtsv" };

            // ---------- EDIT: output + options ----------
            string outputPath   = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\ProteoformAlign_DDA_modsIgnorePos.tsv";
            var    matchKey     = ProteoformMatchKey.ModsIgnorePosition;  // or FullSequence / BaseSequence
            double qValueThresh = 0.01;                                   // drop IDs above this QValue
            // --------------------------------------------

            string written = ProteoformAligner.Align(files, outputPath, matchKey, qValueThresh);

            TestContext.Out.WriteLine("ProteoformAligner wrote: " + written);
            TestContext.Out.WriteLine("Proteoform rows: " + (File.ReadAllLines(written).Length - 1));
        }

        /// <summary>
        /// Manual runner (not an assertion test): runs OpenMS's ConsensusFeatureGroupDetector.exe directly on the
        /// FLASHQuant output tsv files below and reports the consensus row count, for side-by-side comparison with
        /// <see cref="FeatureAligner"/>. [Explicit] so it runs only when selected; Ignored if the exe isn't installed.
        /// NOTE: this tool takes FLASHQuant *_out tsv files (not FLASHDeconv/TopFD). Edit paths/params below.
        /// </summary>
        [Test, Explicit("manual runner: runs OpenMS ConsensusFeatureGroupDetector.exe")]
        public void RunOpenMSConsensusDetector()
        {
            // ---------- EDIT: exe + FLASHQuant output tsv files (one per run) ----------
            string exe = @"C:\Program Files\OpenMS-3.5.0\bin\ConsensusFeatureGroupDetector.exe";
            var dir1 = @"E:\ISD Project\071626_ISD\Quant\HC_YC\FlashDeconv\HC_DDA_11_files";
            var dir2 = @"E:\ISD Project\071626_ISD\Quant\HC_YC\FlashDeconv\YC_DDA_11_files";
            // filter on the FILE NAME (not the full path, which contains the folder "HC_YC" and would match YC too)
            var files1 = Directory.GetFiles(dir1, "*.tsv").ToList();
            var files2 = Directory.GetFiles(dir2, "*.tsv").ToList();

            // ---------- EDIT: output + parameters (defaults match the tool's defaults) ----------
            string outputPath  = @"E:\ISD Project\071626_ISD\Quant\HC_YC\FlashDeconv\HC_11files_DDA_20ppm_120s.tsv";
            string massTol     = "20";             // default 100
            string massTolUnit = "ppm";             // ppm or dalton
            string rtTol       = "120";             // seconds
            string whenDup     = "max_abundance";   // max_abundance or nearest_mass
            // -----------------------------------------------------------------------------------

            var files = files1.ToList();
            if (!File.Exists(exe)) Assert.Ignore("ConsensusFeatureGroupDetector.exe not found: " + exe);
            if (files.Count < 2) Assert.Ignore("Fewer than 2 input files matched.");

            // OpenMS reads FLASHQuant's fixed 17-column layout and CRASHES on FLASHDeconv/TopFD files, so convert
            // each input to FLASHQuant format first (a no-op copy if it is already FLASHQuant).
            var convDir = Path.Combine(Path.GetTempPath(), "fq_converted");
            Directory.CreateDirectory(convDir);
            var inputFiles = files
                .Select(f => FeatureAligner.ConvertFlashDeconvToFlashQuant(f, Path.Combine(convDir, Path.GetFileName(f))))
                .ToList();

            var psi = new ProcessStartInfo
            {
                FileName = exe,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false,
                CreateNoWindow = true
            };
            psi.ArgumentList.Add("-in");
            foreach (var f in inputFiles) psi.ArgumentList.Add(f);   // ArgumentList quotes paths with spaces for us
            psi.ArgumentList.Add("-out");            psi.ArgumentList.Add(outputPath);
            psi.ArgumentList.Add("-mass_tol");       psi.ArgumentList.Add(massTol);
            psi.ArgumentList.Add("-mass_tol_unit");  psi.ArgumentList.Add(massTolUnit);
            psi.ArgumentList.Add("-rt_tol");         psi.ArgumentList.Add(rtTol);
            psi.ArgumentList.Add("-when_duplicate"); psi.ArgumentList.Add(whenDup);

            using var proc = Process.Start(psi);
            var outTask = proc.StandardOutput.ReadToEndAsync();   // read both streams concurrently to avoid deadlock
            var errTask = proc.StandardError.ReadToEndAsync();
            proc.WaitForExit();

            TestContext.Out.WriteLine(outTask.Result);
            TestContext.Out.WriteLine(errTask.Result);            // OpenMS logs progress to stderr
            TestContext.Out.WriteLine("exit code: " + proc.ExitCode);
            if (File.Exists(outputPath))
                TestContext.Out.WriteLine($"Consensus rows: {File.ReadAllLines(outputPath).Length - 1} -> {outputPath}");
        }

        /// <summary>
        /// Assertion test (runs in the normal suite): builds a synthetic multi-charge <see cref="MassFeature"/>
        /// with known per-scan intensities and verifies the two top-3 quant metrics computed in
        /// <see cref="MassFeature.Finalise"/>:
        ///   * <see cref="MassFeature.SummedIntensityTop3"/> — sum of the 3 most abundant traces' sum-over-scans.
        ///   * <see cref="MassFeature.ApexIntensityTop3"/>   — the 3 most abundant traces' intensities read at the
        ///     feature's shared apex scan (the scan where intensity summed across ALL traces is highest).
        /// The scenario uses 4 charge traces so the weakest is correctly dropped from the top-3, and the apex scan
        /// is chosen using all 4 traces while the reported apex value sums only the top-3 (the excluded trace's
        /// intensity counts toward CHOOSING the apex scan but not toward the value).
        /// </summary>
        [Test]
        public void Top3Quant_ComputesSharedApexAndSum()
        {
            // Build one charge-locked trace from (scanIndex, intensity) points, added in scan order.
            static CorrectedTrace Trace(int charge, double mass, params (int scan, double intensity)[] pts)
            {
                var t = new CorrectedTrace { Charge = charge, ConsensusMass = mass };
                foreach (var (scan, intensity) in pts)
                    t.Envelopes.Add(new CorrectedEnvelope
                    {
                        ScanIndex = scan,
                        ScanNumber = scan + 1,
                        RT = scan * 0.1,
                        OriginalMass = mass,
                        CorrectedMass = mass,
                        Charge = charge,
                        Intensity = intensity,
                    });
                return t;
            }

            // Four charge traces. Per-scan intensities chosen so the arithmetic is easy to verify by hand;
            // trailing comment is each trace's TotalIntensity (its sum over scans).
            var A = Trace(10, 10000.0, (0, 100), (1, 500), (2, 1000), (3, 400)); // total 2000
            var B = Trace(11, 10000.1, (0,  50), (1, 300), (2,  800), (3, 200)); // total 1350
            var C = Trace(12, 10000.2,           (1, 100), (2,  600), (3, 100)); // total  800
            var D = Trace( 9, 10000.3,                     (2,   50), (3,  50)); // total  100  (weakest)

            var feature = new MassFeature { Traces = new List<CorrectedTrace> { A, B, C, D } };
            feature.Finalise();

            // Full summed intensity is every trace: 2000 + 1350 + 800 + 100.
            Assert.That(feature.SummedIntensity, Is.EqualTo(4250.0).Within(1e-6));

            // (1) Summed intensity of the top-3 (A, B, C): 2000 + 1350 + 800.
            Assert.That(feature.SummedIntensityTop3, Is.EqualTo(4150.0).Within(1e-6));

            // (2) Shared apex scan = scan 2. Per-scan totals over ALL traces: s0=150, s1=900, s2=2450, s3=750.
            //     Apex value = top-3 (A, B, C) intensities at scan 2 = 1000 + 800 + 600 = 2400.
            //     D's 50 at scan 2 helps pick scan 2 as the apex but is excluded from the reported value.
            Assert.That(feature.ApexIntensityTop3, Is.EqualTo(2400.0).Within(1e-6));

            // Fallback: a feature with fewer than 3 traces uses all it has (here A + C).
            var twoTrace = new MassFeature { Traces = new List<CorrectedTrace> { A, C } };
            twoTrace.Finalise();
            Assert.That(twoTrace.SummedIntensityTop3, Is.EqualTo(2800.0).Within(1e-6)); // 2000 + 800
            // apex scan over {A, C}: s0=100, s1=600, s2=1600, s3=500 -> scan 2; value = 1000 + 600 = 1600.
            Assert.That(twoTrace.ApexIntensityTop3, Is.EqualTo(1600.0).Within(1e-6));
        }

        /// <summary>
        /// Manual runner (not an assertion test): runs consensus feature finding on ONE real MS data file and
        /// writes the `_ms1.feature` table — now including the Apex_intensity_top3 / Sum_intensity_top3 columns.
        /// Point <c>inputFile</c> at a .raw/.mzML, set <c>outputPath</c> + the parameters, then run it directly
        /// (right-click -> Run in the IDE, or dotnet test --filter). [Explicit] so the normal suite skips it;
        /// Ignored if the input file is missing.
        /// </summary>
        [Test, Explicit("manual runner: consensus feature finding on a real file -> _ms1.feature tsv")]
        public void RunConsensusFeatureFinding()
        {
            // ---------- EDIT: input spectra file (.raw or .mzML) ----------
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HC_YC";
            var fileList = Directory.GetFiles(dataDir, "*DDA*").Where(f => f.Contains(".raw") & !f.Contains("_YC_")).ToList();

            var outDir = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\Consensus";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);
            var outFolder = Path.Combine(outDir, "HC_DDA_11files");
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);
            //string outputPath = Path.Combine(outDir, "FeatureAlign_HC_DDA_11files_1.5Da.tsv");

            var    decon          = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            double traceTolDa     = 1.5;    // Da tolerance for stitching envelopes into a charge-locked trace
            int    maxGap         = 1;      // allowed missed scans within a trace
            double featureMassPpm = 15.0;   // ppm tolerance for stitching per-charge traces into a cross-charge feature

            double minMass        = 4000;   // drop features below this neutral mass
            int    minChargeCount = 3;      // drop features seen at fewer than this many charge states
            
            foreach(var file in fileList)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();
                var ms1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToArray();

                var features = ConsensusMassXicConstructor.TraceFeatures(ms1Scans, decon, traceTolDa, maxGap, featureMassPpm);
                var outPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(file) + "_consensus_ms1.feature");
                IsdMsAlignExporter.WriteMs1FeatureFile(features, outPath, minMass, minChargeCount);
            }
        }

        /// <summary>
        /// Manual runner (not an assertion test): aligns TWO OR MORE consensus `_ms1.feature` tables (written by
        /// <see cref="IsdMsAlignExporter.WriteMs1FeatureFile"/> / the RunConsensusFeatureFinding runner) into one
        /// cross-run consensus tsv — per-run SummedIntensity / SummedIntensityTop3 / ApexIntensityTop3, plus the
        /// feature id each run contributed. Point <c>dir</c> at a folder of `_ms1.feature` files, set the output
        /// path + tolerances, then run it directly. [Explicit] so the normal suite skips it; Ignored if fewer than
        /// 2 files match.
        /// </summary>
        [Test, Explicit("manual runner: aligns consensus _ms1.feature files across runs -> consensus tsv")]
        public void RunConsensusFeatureAligner()
        {
            // ---------- EDIT: folder of consensus feature tables to align (_ms1.feature), >= 2 ----------
            var dir = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\Consensus\HC_DDA_11files";
            var files = Directory.GetFiles(dir, "*ms1.feature").ToList();

            // ---------- EDIT: alignment tolerances ----------
            double massTol  = 1.5;
            var    massUnit = MassToleranceUnit.Dalton;
            double rtTol    = 120;          // SECONDS

            var outFolder = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\Consensus";
            string outputPath = Path.Combine(outFolder, $"FeatureAlign_HC_DDA_11files_{massTol}_{massUnit}_{rtTol}.tsv");

            string written = ConsensusFeatureAligner.AlignFeatureFiles(
                files, outputPath, massTol, massUnit, rtTol, positiveNotch: 2, negativeNotch: 2,
                grouping: GroupingMethod.ReferenceAnchor, whenDuplicate: DuplicateResolution.MaxAbundance,
                quantMeasure: ConsensusQuantMeasure.SummedIntensityTop3);

            TestContext.Out.WriteLine("ConsensusFeatureAligner wrote: " + written);
            TestContext.Out.WriteLine("Rows: " + (File.ReadAllLines(written).Length - 1));
        }

        /// <summary>
        /// Manual runner (not an assertion test): answers "of the fully-shared consensus proteoforms
        /// (n_files == <c>requiredNFiles</c>) in a <see cref="ConsensusFeatureAligner"/> table, how many are
        /// identified at 1% FDR (target, QValue &lt;= threshold) in a MetaMorpheus proteoform .psmtsv?".
        ///
        /// Process:
        ///   1. Read the consensus tsv; keep rows whose n_files == requiredNFiles. Each is a "shared" proteoform
        ///      with a consensus mass (mass_mean) and apex RT (rt_apex_mean, already in SECONDS).
        ///   2. Read the psmtsv; keep confident IDs = Decoy/Contaminant/Target == "T" AND QValue &lt;= threshold.
        ///      Each carries the experimental Precursor Mass and Scan Retention Time (MINUTES -> *60 to seconds).
        ///   3. A shared feature is "identified" if some confident ID agrees within the mass tolerance (applied
        ///      around +/- isotope notches, since deconvolutors can disagree on the monoisotope by 1-2 Da) AND
        ///      within the RT tolerance. Reuses <see cref="FeatureAligner"/>'s exact notch/tolerance matching.
        ///   4. Report the count (and a mass-only count, to show how much the RT gate removes).
        ///
        /// Observed 2026-08-12 on FeatureAlign_HC_DDA_11files_20_Ppm_120.tsv vs AllProteoforms.psmtsv
        /// (20 ppm, 120 s, notch +/-2): 15 of 97 shared features identified; 16 with mass only. The notch made no
        /// difference (0/+-1/+-2/+-3 all give 15) -> every match is an exact monoisotope, no off-by-one.
        /// [Explicit] so the normal suite skips it; Ignored if either file is missing.
        /// </summary>
        [Test, Explicit("manual runner: counts identified fully-shared consensus proteoforms vs a psmtsv")]
        public void CountIdentifiedSharedProteoforms()
        {
            var ci = System.Globalization.CultureInfo.InvariantCulture;
            var num = System.Globalization.NumberStyles.Any;

            // ---------- EDIT: consensus alignment tsv + identification psmtsv ----------
            string consensusFile = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\Consensus\FeatureAlign_HC_DDA_11files_20_Ppm_120.tsv";
            string psmFile       = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\HC_YC_DDA_xml-cali-xml\Task1-SearchTask\AllProteoforms.psmtsv";

            // ---------- EDIT: definition of "shared", ID confidence, and match tolerances ----------
            int    requiredNFiles = 11;    // "shared" = present in this many runs (the n_files column)
            double qValueThresh   = 0.01;  // keep only target IDs at/below this QValue
            double massTolPpm     = 20;    // ppm
            double rtTolSeconds   = 120;   // +/- seconds
            int    positiveNotch  = 2, negativeNotch = 2;
            // ---------------------------------------------------------------------------------------

            if (!File.Exists(consensusFile)) Assert.Ignore("Consensus file not found: " + consensusFile);
            if (!File.Exists(psmFile)) Assert.Ignore("psmtsv not found: " + psmFile);

            // ---- (1) shared consensus features: (mass_mean, rt_apex_mean in seconds) ----
            var cLines = FeatureAligner.ReadAllLinesShared(consensusFile);
            var cHead  = cLines[0].Split('\t');
            int cMass  = FeatureAligner.FindColumn(cHead, new[] { "mass_mean" });
            int cRt    = FeatureAligner.FindColumn(cHead, new[] { "rt_apex_mean" });
            int cN     = FeatureAligner.FindColumn(cHead, new[] { "n_files" });
            int cNeed  = System.Math.Max(cMass, System.Math.Max(cRt, cN));
            var shared = new List<(double mass, double rt)>();
            for (int i = 1; i < cLines.Length; i++)
            {
                if (cLines[i].Length == 0) continue;
                var t = cLines[i].Split('\t');
                if (t.Length <= cNeed) continue;
                if (!int.TryParse(t[cN], out int nf) || nf != requiredNFiles) continue;
                if (double.TryParse(t[cMass], num, ci, out double m) && double.TryParse(t[cRt], num, ci, out double r))
                    shared.Add((m, r));
            }

            // ---- (2) confident target IDs: (Precursor Mass, Scan Retention Time -> seconds) ----
            var pLines = FeatureAligner.ReadAllLinesShared(psmFile);
            var pHead  = pLines[0].Split('\t');
            int pMass  = FeatureAligner.FindColumn(pHead, new[] { "Precursor Mass" });
            int pRt    = FeatureAligner.FindColumn(pHead, new[] { "Scan Retention Time" });
            int pQ     = FeatureAligner.FindColumn(pHead, new[] { "QValue" });
            int pTgt   = FeatureAligner.FindColumn(pHead, new[] { "Decoy/Contaminant/Target" });
            int pNeed  = new[] { pMass, pRt, pQ, pTgt }.Max();
            var ids = new List<(double mass, double rt)>();
            for (int i = 1; i < pLines.Length; i++)
            {
                if (pLines[i].Length == 0) continue;
                var t = pLines[i].Split('\t');
                if (t.Length <= pNeed) continue;
                if (t[pTgt].Trim() != "T") continue;                                      // target only
                if (!double.TryParse(t[pQ], num, ci, out double q) || q > qValueThresh) continue;
                if (double.TryParse(t[pMass], num, ci, out double m) && double.TryParse(t[pRt], num, ci, out double r))
                    ids.Add((m, r * 60.0));                                                // psmtsv RT is in minutes
            }

            // ---- (3) match each shared feature to a confident ID by mass(+notch) [and RT] ----
            var tol    = FeatureAligner.MakeTolerance(massTolPpm, MassToleranceUnit.Ppm);
            var shifts = FeatureAligner.NotchShifts(positiveNotch, negativeNotch);
            int idMatch = 0, massOnly = 0;
            foreach (var f in shared)
            {
                bool massHit = false, rtHit = false;
                foreach (var id in ids)
                {
                    if (!FeatureAligner.WithinNotch(tol, f.mass, id.mass, shifts)) continue;
                    massHit = true;
                    if (System.Math.Abs(f.rt - id.rt) <= rtTolSeconds) { rtHit = true; break; }
                }
                if (massHit) massOnly++;
                if (rtHit)  idMatch++;
            }

            // ---- (4) report ----
            TestContext.Out.WriteLine($"shared consensus features (n_files=={requiredNFiles}): {shared.Count}");
            TestContext.Out.WriteLine($"confident target IDs (target, QValue<={qValueThresh}): {ids.Count}");
            TestContext.Out.WriteLine($"IDENTIFIED (mass {massTolPpm}ppm +/-{positiveNotch} notch AND RT +/-{rtTolSeconds}s): {idMatch} / {shared.Count}");
            TestContext.Out.WriteLine($"identified (mass only, RT ignored): {massOnly} / {shared.Count}");

            Assert.That(idMatch, Is.LessThanOrEqualTo(shared.Count));   // sanity: can't identify more than exist
        }
    }
}
