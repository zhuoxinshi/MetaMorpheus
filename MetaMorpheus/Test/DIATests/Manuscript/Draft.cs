using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using EngineLayer.DIA.FeatureAlign;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
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

            bool group = false;
            bool write = false;
            bool search = true;

            var apexTolerance = 0.1f; var corrThreshold = 0.75; var fragmentRank = 200;
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), traceToleranceDa: 1.5, minMass: 4000, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 5, new ClassicDeconvolutionParameters(1, 60, 4, 3), minMass: 4000, minCharge: 3, xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                new NeutralMassXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), xicSpline: new XicCubicSpline(0.05, 1, 1, true)),
                //new XicGroupingEngine(0.2f, 0.3, 0.5, maxThreadsForGrouping: 10),
                new XicGroupingEngine(apexTolerance, 0, corrThreshold, maxThreadsForGrouping: 10, fragmentRankThreshold: fragmentRank), PseudoMs2ConstructionType.Mass, combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            var outDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\MM\ISD_search";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);
            //var outFolder = Path.Combine(outDir, $"{apexTolerance}_{corrThreshold}_{fragmentRank}");
            var outFolder = Path.Combine(outDir, "neutralMass6");
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

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
                searchTask.SearchParameters.WritePrunedDatabase = true;

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

                var searchFileList = Directory.GetFiles(outFolder, "*_pseudo_ms2.msalign").ToList();
                var taskList = new List<(string, MetaMorpheusTask)> { ("search_msalign", searchTask), ("GPTMD_msalign", gptmdTask), ("postGPTMD-search_msalign", searchTask) }; // ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) 
                var engine = new EverythingRunnerEngine(taskList, searchFileList, databases, outFolder);
                engine.Run();
            }
        }

        [Test]
        public static void SearchAllPseudoScansMsalign_ExternalDecon()
        {
            var topFdRoot = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\SeparateISD\TopFD";
            var dataDir = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
            var rawFiles = Directory.GetFiles(dataDir, "*ISD60-80-100*").Where(f => f.Contains(".raw") & f.Contains("HCD-YBC")).ToList();

            var apexTolerance = 2f; var corrThreshold = 0.5; int? fragmentRank = 100; int? precursorRank = null;int minFrags = 10;
            var outDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);
            var outFolder = Path.Combine(outDir, $"{apexTolerance}_{corrThreshold}_{fragmentRank}_{precursorRank}_{minFrags}");
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);

            bool group = true;
            bool write = true; 
            bool search = false;

            if (group)
            {
                foreach (var file in rawFiles)
                {
                    string stem = Path.GetFileNameWithoutExtension(file);

                    // TopFD was run SEPARATELY on each ISD channel, so every run has four *_feature.xml under
                    // topFdRoot: one for the MS1 (precursor) channel and one per source-fragmentation voltage.
                    // Build all of them from stem (Path.Combine, not string concat: topFdRoot is a directory).
                    string precursorFeatureXml = Path.Combine(topFdRoot, stem + "_ms1_feature.xml");
                    var fragmentFeatureXmlByVoltage = new Dictionary<int, string>
                    {
                        [60]  = Path.Combine(topFdRoot, stem + "_isd60_feature.xml"),
                        [80]  = Path.Combine(topFdRoot, stem + "_isd80_feature.xml"),
                        [100] = Path.Combine(topFdRoot, stem + "_isd100_feature.xml"),
                    };

                    var dataFile = MsDataFileReader.GetDataFile(file);
                    dataFile.LoadAllStaticData();

                    // Precursor XICs read straight from the MS1 channel's feature.xml (intact proteoforms, >= 4 kDa).
                    var ms1Constructor = new FeatureXmlXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                        precursorFeatureXml, minMass: 4000, maxMass: 70000, minEcScore: 0,
                        xicSpline: new XicCubicSpline(0.05, 1, 1, true));

                    // Fragment XICs built the SAME way, but keyed by source voltage so each channel is read from its
                    // own feature.xml. ISDEngine hands one voltage's scans at a time and BySourceVoltage dispatches
                    // on the sid= in the scan filter; the feature.xml scan numbers are the original raw numbers
                    // (each channel steps by 4), so they join to those scans by scan number. No mass floor on
                    // fragments.
                    var ms2Constructor = FeatureXmlXicConstructor.BySourceVoltage(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                        fragmentFeatureXmlByVoltage, minMass: 0, maxMass: 70000, minEcScore: 0,
                        xicSpline: new XicCubicSpline(0.05, 1, 1, true));

                    var diaParams = new DIAparameters(
                        DIAanalysisType.ISD, ms1Constructor, ms2Constructor,
                        new XicGroupingEngine(apexTolerance, 0, corrThreshold, maxThreadsForGrouping: 10, fragmentRankThreshold: fragmentRank, precursorRankThreshold: precursorRank, minFragmentCountForGrouping: minFrags), PseudoMs2ConstructionType.Mass, combineFragments: true);
                    var commonParams = new CommonParameters { DIAparameters = diaParams };

                    var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();

                    var msAlignOutPath = Path.Combine(outFolder, stem + "_pseudo_ms2.msalign");
                    if (write)
                    {
                        diaParams.WriteSettingsToFile(Path.Combine(outFolder, "DIAparameters.txt"));
                        // Fragments no longer come from a TopFD msAlign, so the per-peak EnvCNN 4th column is not
                        // available here; write the standard 3-column msAlign. (To restore EnvCNN, pass the per-
                        // voltage *_ms2.msalign as the score source.)
                        ISDEngine.WriteMsAlignFile(msAlignOutPath, pseudoScans);
                    }
                }
            }

            if (search)
            {
                string searchTomlCommonFixedVariable = @"E:\ISD Project\ISD_250428\new_topDown_search_toml_commonFixedVariable\Task Settings\Task1-SearchTaskconfig.toml";
                var searchTask = Toml.ReadFile<SearchTask>(searchTomlCommonFixedVariable, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DoPrecursorDeconvolution = false; // already deconvoluted via FromFile
                searchTask.CommonParameters.UseProvidedPrecursorInfo = true; // FromFile assembly provides the precursor info
                searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;

                var moreGPTMD_toml = @"E:\ISD Project\ISD_250428\moreGPTMD_toml\Task Settings\Task1-GPTMDTaskconfig.toml";
                var gptmdTask = Toml.ReadFile<GptmdTask>(moreGPTMD_toml, MetaMorpheusTask.tomlConfig);
                gptmdTask.CommonParameters = searchTask.CommonParameters;

                // YB is yeast; add the human xml back if the data files are switched to the HC_YC mixes
                string yeast_xml = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2026_08_05.xml";
                var human_xml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
                var databases = new List<DbForTask> { new DbForTask(yeast_xml, false) };

                var searchFileList = Directory.GetFiles(outFolder, "*_pseudo_ms2.mgf").ToList();
                var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask), ("GPTMD", gptmdTask), ("postGPTMD-search", searchTask) };
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
            var outFolder = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC\MM\DDA_custom\Normal";
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
            var dir = @"E:\ISD Project\071626_ISD\Quant\Manuscript\TopFD\ms1_feature";
            var files = Directory.GetFiles(dir, "*ms1.feature").Where(f => f.Contains("HCD-YBC_1-1") & f.Contains("ISD60-80-100")).ToList();// f.Contains("HCD-YBC_1-") & !f.Contains("ISD60-80-100")

            //var files = new List<string> { @"E:\ISD Project\071626_ISD\Quant\Manuscript\FeatureAlign\FeatureAlign_YBC_DDA_8files_20_Ppm_120.tsv", @"E:\ISD Project\071626_ISD\Quant\Manuscript\FeatureAlign\FeatureAlign_YBC_ISD_8files_20_Ppm_120.tsv" };

            double massTol = 20;
            var massUnit = MassToleranceUnit.Ppm;   // or MassToleranceUnit.Dalton
            double rtTol = 60;                      // seconds, +/-

            // ---------- EDIT: OUTPUT FILTERS (write out fewer features) ----------
            // Keep only consensus rows whose n_files and mean mass fall inside these INCLUSIVE ranges. Any bound left
            // null is not enforced, so all four null => write every row (the original behaviour).
            int?    minNFiles = 3;      // e.g. 8 keeps only features present in all 8 runs; null = no lower bound
            int?    maxNFiles = null;   // null = no upper bound
            double? minMass   = 4000;   // Da; null = no lower bound
            double? maxMass   = null;   // Da; null = no upper bound
            int notches = 2;
            // --------------------------------------------------------------------

            var outFolder = @"E:\ISD Project\071626_ISD\Quant\Manuscript\FeatureAlign_new";
            if (!Directory.Exists(outFolder)) Directory.CreateDirectory(outFolder);
            string outputPath = Path.Combine(outFolder, $"FeatureAlign_HCD-YBC_1-1_ISD_3files_{massTol}_{massUnit}_{rtTol}_{notches}.tsv");

            string written = FeatureAligner.Align(files, outputPath, massTol, massUnit, rtTol,
                grouping: GroupingMethod.ReferenceAnchor, positiveNotch: notches, negativeNotch: notches,
                minNFiles: minNFiles, maxNFiles: maxNFiles, minMass: minMass, maxMass: maxMass);

            TestContext.Out.WriteLine("FeatureAligner wrote: " + written);
            TestContext.Out.WriteLine($"Filters -> n_files [{minNFiles?.ToString() ?? "-"}, {maxNFiles?.ToString() ?? "-"}], " +
                $"mass [{minMass?.ToString() ?? "-"}, {maxMass?.ToString() ?? "-"}]");
            TestContext.Out.WriteLine("Rows: " + (File.ReadAllLines(written).Length - 1));
        }
    }
}
