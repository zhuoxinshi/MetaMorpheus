using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using EngineLayer.DIA.FeatureAlign;
using MassSpectrometry;
using Microsoft.Win32;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;
using UsefulProteomicsDatabases;
using static iText.StyledXmlParser.Jsoup.Select.Evaluator;

namespace Test.DIATests.Manuscript
{
    internal class Draft
    {
        [Test]
        public static void SearchAllPseudoScansMgf()
        {
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
            var fileList = Directory.GetFiles(dataDir, "*ISD60-80-100*").Where(f => f.Contains(".raw") & f.Contains("HCD-YBC")).ToList();
            //var fileList = new List<string>
            //{
            //    @"E:\ISD Project\071626_ISD\Quant\4proYB\08-02-26_4proYB_1-10_ISD_65min_rep1.raw"
            //};

            var outDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\MM\ISD_search";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);
            var outFolder = Path.Combine(outDir, "neutralMass10");
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

            bool group = true;
            bool write = true;
            bool search = true;

            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), traceToleranceDa: 1.5, minMass: 4000, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new XicGroupingEngine(0.2f, 0.3, 0.5, maxThreadsForGrouping: 10),
                new XicGroupingEngine(0.1f, 0, 0.5, maxThreadsForGrouping: 10, fragmentRankThreshold: 200), PseudoMs2ConstructionType.Mass, combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            if (group)
            {
                foreach (var file in fileList)
                {
                    var dataFile = MsDataFileReader.GetDataFile(file);
                    dataFile.LoadAllStaticData();
                    var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();
                    var msAlignOutPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.msalign");
                    var mgfOutPath = Path.Combine(outFolder, Path.GetFileNameWithoutExtension(dataFile.FilePath) + "_pseudo_ms2.mgf");
                    if (write)
                    {
                        diaParams.WriteSettingsToFile(Path.Combine(outFolder, "DIAparameters.txt"));
                        ISDEngine.WriteMsAlignFile(msAlignOutPath, pseudoScans);
                        IsdMsAlignExporter.WriteMgf(pseudoScans, mgfOutPath);
                    }
                }
            }

            if (search)
            {
                string calibrationToml = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\HC_YC_DDA_xml-cali-xml\Task Settings\Task2-CalibrateTaskconfig.toml";
                var calibrationTask = Toml.ReadFile<CalibrationTask>(calibrationToml, MetaMorpheusTask.tomlConfig);
                string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
                var searchTomlNoFixed = @"E:\ISD Project\071626_ISD\Quant\PierceYB\MM\DDA_6files_6protDb\Task Settings\Task1-SearchTaskconfig.toml";
                var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DoPrecursorDeconvolution = false; // already deconvoluted via FromFile
                searchTask.CommonParameters.UseProvidedPrecursorInfo = true; // FromFile assembly provides the precursor info
                searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;
                //searchTask.SearchParameters.MinAllowedInternalFragmentLength = 0;

                var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters;

                string yeast_xml = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.xml";
                var human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
                var yeast_fasta = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.fasta";
                string pierceDb = @"E:\ISD Project\071626_ISD\Quant\Pierce\Pierce_6prot_new.fasta";
                string human_fasta = @"E:\Databases\Human_9606.fasta";
                var std4pro_xml = @"E:\ISD Project\Std_4pro.xml";
                var databases = new List<DbForTask> { new DbForTask(yeast_xml, false), new DbForTask(human_xml, false) };//, new DbForTask(human_xml, false)

                var searchFileList = Directory.GetFiles(outFolder, "*_pseudo_ms2.mgf").ToList();
                var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) }; // ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) 
                var engine = new EverythingRunnerEngine(taskList, searchFileList, databases, outFolder);
                engine.Run();
            }
        }

        [Test]
        public static void SearchDDA()
        {
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
            var fileList = Directory.GetFiles(dataDir, "*raw").Where(f => !f.Contains("ISD60-80-100") & f.Contains("HCD-YBC_1-2")).ToList();
            //var fileList = new List<string>
            //{
            //    @"E:\ISD Project\071626_ISD\Quant\4proYB\08-02-26_4proYB_1-10_ISD_65min_rep1.raw"
            //};
            var outFolder = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\MM\DDA_custom";
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

            string calibrationToml = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\HC_YC_DDA_xml-cali-xml\Task Settings\Task2-CalibrateTaskconfig.toml";
            var calibrationTask = Toml.ReadFile<CalibrationTask>(calibrationToml, MetaMorpheusTask.tomlConfig);
            string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
            var searchTomlNoFixed = @"E:\ISD Project\071626_ISD\Quant\PierceYB\MM\DDA_6files_6protDb\Task Settings\Task1-SearchTaskconfig.toml";
            var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);

            var lessGPTMD_noFilter_toml = @"E:\ISD Project\ISD_250428\new_topDown_lessGPTMD_noFilter_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
            var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
            var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);

            string yeast_xml = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.xml";
            var human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
            var yeast_fasta = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.fasta";
            string pierceDb = @"E:\ISD Project\071626_ISD\Quant\Pierce\Pierce_6prot_new.fasta";
            string human_fasta = @"E:\Databases\Human_9606.fasta";
            var std4pro_xml = @"E:\ISD Project\Std_4pro.xml";
            var databases = new List<DbForTask> { new DbForTask(yeast_xml, false), new DbForTask(human_xml, false) };//, new DbForTask(human_xml, false)

            var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) }; // ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) 
            var engine = new EverythingRunnerEngine(taskList, fileList, databases, outFolder);
            engine.Run();
        }

        // Voltage-contribution analysis (manuscript point 2): rebuild the nm5 pseudo-MS2 exactly as the search
        // did, but keep which ISD source voltage (60/80/100) each fragment XIC came from, and dump one row per
        // (pseudo-scan, fragment) with that voltage. Precursor mass + apex RT are written so the table can be
        // joined to the existing nm5 search identifications without re-searching.
        [Test, Explicit("manuscript: per-fragment source-voltage table for the nm5 ISD pseudo-MS2")]
        public static void VoltageContributionFragmentTable()
        {
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
            var files = Directory.GetFiles(dataDir, "*ISD60-80-100*")
                .Where(f => f.EndsWith(".raw") && f.Contains("HCD-YBC")).OrderBy(f => f).ToList();
            string outPath = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\MM\ISD_search\neutralMass_venn\voltage_fragments_nm5.tsv";

            using var w = new StreamWriter(outPath, false);
            w.WriteLine("file	scan_id	precursor_mass	precursor_rt	precursor_charge	n_fragments	frag_mass	frag_intensity	voltage");

            foreach (var file in files)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();
                var allScans = dataFile.GetAllScansList().ToArray();
                var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
                ISDEngine.ReLabelIsdScans(isdMap, allScans);

                // nm5 constructors: MS1 (precursor) and MS2 (fragment), NeutralMass, exactly as SearchAllPseudoScansMgf
                var ms1Ctor = new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                    new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3,
                    xicSpline: new XicCubicSpline(0.05, 1, 1, true));
                var ms2Ctor = new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                    new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new XicCubicSpline(0.05, 1, 1, true));

                var ms1Xics = ms1Ctor.GetAllXicsWithXicSpline(ms1Scans);

                // build fragment XICs one voltage channel at a time, remembering each XIC's voltage (reference key)
                var fragToVoltage = new Dictionary<ExtractedIonChromatogram, int>();
                var allFrag = new List<ExtractedIonChromatogram>();
                foreach (var kvp in isdMap.OrderBy(k => k.Key))
                {
                    foreach (var x in ms2Ctor.GetAllXicsWithXicSpline(kvp.Value.ToArray()))
                    {
                        fragToVoltage[x] = kvp.Key;
                        allFrag.Add(x);
                    }
                }

                // nm5 grouping: apexRT 0.1, overlap 0, correlation 0, fragment-rank cap 200
                var grouper = new XicGroupingEngine(0.1f, 0, 0, maxThreadsForGrouping: 10, fragmentRankThreshold: 200);
                var pfGroups = grouper.PrecursorFragmentGrouping(ms1Xics, allFrag).ToList();

                string stem = Path.GetFileNameWithoutExtension(file);
                int scanId = 0;
                foreach (var g in pfGroups)
                {
                    scanId++;
                    int pcharge = g.PrecursorXic.ApexPeak is IndexedMass im ? im.Charge : 1;
                    double pmass = g.PrecursorXic.ApexPeak.M;
                    double prt = g.PrecursorXic.ApexRT;
                    int nfrag = g.PFpairs.Count;
                    foreach (var pf in g.PFpairs)
                    {
                        int volt = fragToVoltage.TryGetValue(pf.FragmentXic, out var vv) ? vv : 0;
                        double fmass = pf.FragmentXic.ApexPeak.M;
                        double fint = pf.FragmentXic.ApexPeak.Intensity;
                        w.WriteLine($"{stem}	{scanId}	{pmass:R}	{prt:R}	{pcharge}	{nfrag}	{fmass:R}	{fint:R}	{volt}");
                    }
                }
                TestContext.Out.WriteLine($"{stem}: {pfGroups.Count} pseudo-scans, " +
                    $"{pfGroups.Sum(g => g.PFpairs.Count)} fragment-pairs, voltages {string.Join("/", isdMap.Keys.OrderBy(v => v))}");
            }
            TestContext.Out.WriteLine("wrote " + outPath);
        }

        // Export the MS1 precursor XICs (mass, apex RT, and the full elution profile as scanIndex:intensity)
        // for the nm5 setting, so precursor-to-precursor co-elution can be tested in Python: how many
        // UNIDENTIFIED precursors co-elute tightly (apex RT + elution-profile correlation) with an
        // IDENTIFIED precursor. Aligning two XICs by scan index makes the correlation a simple vector op.
        [Test, Explicit("manuscript: export MS1 precursor XICs (elution profiles) for precursor co-elution grouping")]
        public static void PrecursorCoelutionExport()
        {
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
            var files = Directory.GetFiles(dataDir, "*ISD60-80-100*")
                .Where(f => f.EndsWith(".raw") && f.Contains("HCD-YBC")).OrderBy(f => f).ToList();
            string outPath = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\MM\ISD_search\neutralMass_venn\precursor_xics_nm5.tsv";

            using var w = new StreamWriter(outPath, false);
            w.WriteLine("file\txic_id\tmass\tapex_rt\tapex_scan\tapex_intensity\tcharge\tprofile");

            foreach (var file in files)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                dataFile.LoadAllStaticData();
                var allScans = dataFile.GetAllScansList().ToArray();
                var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
                ISDEngine.ReLabelIsdScans(isdMap, allScans);

                // nm5 MS1 (precursor) constructor, exactly as SearchAllPseudoScansMgf
                var ms1Ctor = new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                    new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3,
                    xicSpline: new XicCubicSpline(0.05, 1, 1, true));
                var xics = ms1Ctor.GetAllXicsWithXicSpline(ms1Scans);

                string stem = Path.GetFileNameWithoutExtension(file);
                int id = 0;
                foreach (var x in xics)
                {
                    id++;
                    int charge = x.ApexPeak is IndexedMass im ? im.Charge : 1;
                    var prof = string.Join(",", x.Peaks.Select(p => p.ZeroBasedScanIndex + ":" + ((float)p.Intensity)));
                    w.WriteLine($"{stem}\t{id}\t{x.ApexPeak.M:R}\t{x.ApexRT:R}\t{x.ApexPeak.ZeroBasedScanIndex}\t" +
                                $"{x.ApexPeak.Intensity:R}\t{charge}\t{prof}");
                }
                TestContext.Out.WriteLine($"{stem}: {xics.Count} precursor XICs");
            }
            TestContext.Out.WriteLine("wrote " + outPath);
        }

        [Test]
        public void RunFeatureAligner()
        {
            // ---------- EDIT: input feature tables (FLASHQuant or FLASHDeconv .tsv) ----------
            var dir = @"E:\ISD Project\071626_ISD\Quant\Manuscript\TopFD\ms1_feature";
            // filter on the FILE NAME (not the full path, which contains the folder "HC_YC" and would match YC too)
            var files = Directory.GetFiles(dir, "*ms1.feature").Where(f => f.Contains("HCD-YBC_1-2")).ToList();// & !f.Contains("ISD60-80-100")

            double massTol = 20;
            var massUnit = MassToleranceUnit.Ppm;   // or MassToleranceUnit.Dalton
            double rtTol = 120;                      // seconds, +/-

            // ---------- EDIT: OUTPUT FILTERS (write out fewer features) ----------
            // Keep only consensus rows whose n_files and mean mass fall inside these INCLUSIVE ranges. Any bound left
            // null is not enforced, so all four null => write every row (the original behaviour).
            int?    minNFiles = 1;      // e.g. 8 keeps only features present in all 8 runs; null = no lower bound
            int?    maxNFiles = null;   // null = no upper bound
            double? minMass   = 4000;   // Da; null = no lower bound
            double? maxMass   = null;   // Da; null = no upper bound
            // --------------------------------------------------------------------

            var outFolder = @"E:\ISD Project\071626_ISD\Quant\Manuscript\FeatureAlign";
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);
            string outputPath = Path.Combine(outFolder, $"FeatureAlign_HCD-YBC_1-2_DDA&ISD_1-6files_{massTol}_{massUnit}_{rtTol}.tsv");

            string written = FeatureAligner.Align(files, outputPath, massTol, massUnit, rtTol,
                grouping: GroupingMethod.ReferenceAnchor, positiveNotch: 2, negativeNotch: 2,
                minNFiles: minNFiles, maxNFiles: maxNFiles, minMass: minMass, maxMass: maxMass);

            TestContext.Out.WriteLine("FeatureAligner wrote: " + written);
            TestContext.Out.WriteLine($"Filters -> n_files [{minNFiles?.ToString() ?? "-"}, {maxNFiles?.ToString() ?? "-"}], " +
                $"mass [{minMass?.ToString() ?? "-"}, {maxMass?.ToString() ?? "-"}]");
            TestContext.Out.WriteLine("Rows: " + (File.ReadAllLines(written).Length - 1));
        }

        /// <summary>
        /// Runs OUR ReferenceAnchor aligner on the same 3 FLASHDeconv mix files that OpenMS
        /// ConsensusFeatureGroupDetector was run on, with OpenMS-matched tolerances (100 ppm, 180 s, notch 0),
        /// so the two consensus tables can be diffed. [Explicit].
        /// </summary>
        [Test, Explicit("head-to-head vs OpenMS ConsensusFeatureGroupDetector on 3 FLASHDeconv mix files")]
        public void OursVsOpenMS_ThreeMixFlashDeconv()
        {
            var fd = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\FeatureFiles";
            var files = new[]
            {
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep1.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep2.tsv",
                "07-27-26_HCD-YBC_1-1_81min_DDA_3iso_mscan4_rep3.tsv",
            }.Select(f => Path.Combine(fd, f)).ToList();

            string outPath = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\FlashDeconv\_compare\ours_consensus.tsv";
            FeatureAligner.Align(files, outPath,
                massTolerance: 100, massUnit: MassToleranceUnit.Ppm, rtTolerance: 90,
                preferredQuantColumn: "FeatureQuantity",
                positiveNotch: 0, negativeNotch: 0,
                grouping: GroupingMethod.ReferenceAnchor, whenDuplicate: DuplicateResolution.MaxAbundance);

            var lines = File.ReadAllLines(outPath);
            int nf = System.Array.IndexOf(lines[0].Split('\t'), "n_files");
            var dist = lines.Skip(1).Where(l => l.Length > 0).Select(l => int.Parse(l.Split('\t')[nf]))
                            .GroupBy(x => x).OrderBy(g => g.Key).Select(g => $"{g.Key}:{g.Count()}");
            TestContext.Out.WriteLine($"OURS consensus rows: {lines.Length - 1}");
            TestContext.Out.WriteLine("n_files: " + string.Join("  ", dist));
        }

    }
}
