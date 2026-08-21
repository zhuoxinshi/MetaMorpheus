using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Nett;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;

namespace Test.DIATests
{
    /// <summary>
    /// Builds ISD pseudo-MS2 scans with TopFD feature XICs on BOTH sides.
    ///
    /// Earlier variants traced the fragment channels out of an ms2.msalign. Here each ISD source voltage was
    /// split into its own mzML and run through TopFD as if it were an MS1 acquisition, so every channel has its
    /// own <c>*_feature.xml</c> with per-charge envelope lists. That means fragments get the same treatment as
    /// precursors: features and elution profiles straight from TopFD, no peak tracing and no tolerance matching
    /// anywhere in this test.
    ///
    /// The channel splits keep the parent file's scan numbering - isd60 holds scans 2, 6, 10 ..., isd80 holds
    /// 3, 7, 11 ..., isd100 holds 4, 8, 12 ... - so each channel's feature xml joins only to its own scans and
    /// the voltages can never contaminate one another.
    ///
    /// Every knob is a constant at the top of the test so thresholds can be changed without hunting through the
    /// body.
    /// </summary>
    [TestFixture]
    public static class FeatureXmlPseudoMs2Test
    {
        /// <summary>
        /// The same feature-xml-on-both-sides workflow as
        /// <see cref="PseudoMs2FromFeatureXmlPrecursorsAndFragments"/>, but driven through the NORMAL ISD
        /// pipeline: the constructors go onto <see cref="DIAparameters"/> and the pseudo scans come out of
        /// <see cref="MetaMorpheusTask.GetMs2Scans"/>, exactly as a real search would obtain them.
        ///
        /// The one thing that made this awkward is that <see cref="ISDEngine"/> has a single
        /// <see cref="DIAparameters.Ms2XicConstructor"/> slot while the fragment channels live in separate files.
        /// <see cref="FeatureXmlXicConstructor.BySourceVoltage"/> solves it inside the constructor: ISDEngine
        /// hands over one channel at a time and the constructor picks that channel's features from the sid= scan
        /// filter, so no change to ISDEngine or DIAparameters was needed.
        ///
        /// Set <c>search</c> to run MetaMorpheus on the result instead of only writing the msAlign/mgf.
        /// </summary>
        [Test, Explicit]
        public static void PseudoMs2FromFeatureXmlThroughIsdEngine()
        {
            // ------------------------------------------------------------------ knobs
            const double precursorMinMass = 4000;
            const double precursorMaxMass = 70000;
            const double fragmentMaxMass = 70000;
            const double minEcScore = 0;
            const int minNumberOfPeaks = 3;
            const bool aggregateFragmentCharges = true;

            const float apexRtTolerance = 0.2f;
            const double overlapThreshold = 0;
            const double correlationThreshold = 0.0;

            bool write = true;
            bool search = false;

            // ------------------------------------------------------------------ files
            string toppicDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string stem = @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected";
            string parentMzml = Path.Combine(toppicDir, stem + ".mzML");
            string precursorFeatureXml = Path.Combine(toppicDir, stem + "_feature.xml");
            string channelDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\HCD_YBC_ISDseparate";
            string searchToml = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\DDA_1-1\Task Settings\Task1-SearchTaskconfig.toml";
            string databaseXml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
            string outDir = Path.Combine(Path.GetDirectoryName(channelDir), "PseudoMs2_isdEngine");

            var fragmentFeatureXmlByVoltage = Directory.GetFiles(channelDir, "*_feature.xml")
                .Select(f => (File: f, Match: Regex.Match(Path.GetFileName(f), @"_isd(\d+)_feature\.xml$", RegexOptions.IgnoreCase)))
                .Where(t => t.Match.Success)
                .ToDictionary(t => int.Parse(t.Match.Groups[1].Value), t => t.File);
            Assert.That(fragmentFeatureXmlByVoltage, Is.Not.Empty, $"no *_isd<NN>_feature.xml under {channelDir}");

            // ------------------------------------------------------------------ DIA parameters
            // MS1: the whole-file feature xml, one XIC per proteoform feature.
            // MS2: one feature xml per source voltage, dispatched on the sid= of whichever channel ISDEngine
            //      passes in. Fragment charges are aggregated so each fragment mass appears once per pseudo scan.
            var ms1Constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                precursorFeatureXml, minMass: precursorMinMass, maxMass: precursorMaxMass, minEcScore: minEcScore);

            var ms2Constructor = FeatureXmlXicConstructor.BySourceVoltage(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                fragmentFeatureXmlByVoltage, maxMass: fragmentMaxMass, minEcScore: minEcScore,
                aggregateCharges: aggregateFragmentCharges);

            var diaParams = new DIAparameters(DIAanalysisType.ISD, ms1Constructor, ms2Constructor,
                new XicGroupingEngine(apexRtTolerance, overlapThreshold, correlationThreshold,
                    maxThreadsForGrouping: Math.Max(1, Environment.ProcessorCount - 2)),
                PseudoMs2ConstructionType.Mass, combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            // ------------------------------------------------------------------ run it the normal way
            var dataFile = MsDataFileReader.GetDataFile(parentMzml);
            dataFile.LoadAllStaticData();
            var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, dataFile.FilePath, commonParams).ToArray();

            var fragmentCounts = pseudoScans.Select(s => s.ExperimentalFragments.Length).OrderBy(n => n).ToArray();
            TestContext.Progress.WriteLine($"{stem}: {pseudoScans.Length} pseudo-MS2 scans via ISDEngine; " +
                $"MS1 {ms1Constructor.FeatureCount} features -> {ms1Constructor.LastXicCount} XICs; " +
                $"fragment channels {string.Join(", ", fragmentFeatureXmlByVoltage.Keys.OrderBy(v => v))}");
            TestContext.Progress.WriteLine($"   fragments/precursor: mean={fragmentCounts.Average():F1} " +
                $"median={fragmentCounts[fragmentCounts.Length / 2]} max={fragmentCounts[^1]}; " +
                $"total={fragmentCounts.Sum():N0}");
            Assert.That(pseudoScans.Length, Is.GreaterThan(0));

            Directory.CreateDirectory(outDir);
            string mgfPath = Path.Combine(outDir, stem + "_pseudo_ms2.mgf");
            if (write)
            {
                diaParams.WriteSettingsToFile(Path.Combine(outDir, "DIAparameters.txt"));
                ISDEngine.WriteMsAlignFile(Path.Combine(outDir, stem + "_pseudo_ms2.msalign"), pseudoScans);
                IsdMsAlignExporter.WriteMgf(pseudoScans, mgfPath);
                TestContext.Progress.WriteLine($"wrote msAlign + mgf to {outDir}");
            }

            if (search)
            {
                var searchTask = Toml.ReadFile<SearchTask>(searchToml, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DoPrecursorDeconvolution = false;
                searchTask.CommonParameters.UseProvidedPrecursorInfo = true;
                searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;

                new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("search", searchTask) },
                    new List<string> { mgfPath },
                    new List<DbForTask> { new DbForTask(databaseXml, false) },
                    outDir).Run();
                TestContext.Progress.WriteLine("search finished; results under " + outDir);
            }
        }

        [Test, Explicit]
        public static void PseudoMs2FromFeatureXmlPrecursorsAndFragments()
        {
            // ------------------------------------------------------------------ knobs
            // precursor (MS1) feature filters
            const double precursorMinMass = 4000;
            const double precursorMaxMass = 70000;
            const int precursorMinChargeCount = 1;
            const double precursorMinEcScore = 0;

            // fragment feature filters
            const double fragmentMinMass = 0;
            const double fragmentMaxMass = 70000;
            const int fragmentMinChargeCount = 1;
            const double fragmentMinEcScore = 0;

            // A fragment feature is one neutral mass; aggregating sums its charge states into a single XIC, which
            // is what the pseudo scan wants (one line per fragment mass). Set false to keep each charge state as
            // its own XIC, the way the msAlign-tracing route behaves.
            const bool aggregateFragmentCharges = true;

            // minimum points an XIC must have, applied to precursors and fragments alike
            const int minNumberOfPeaks = 3;

            // precursor-fragment grouping
            const float apexRtTolerance = 0.1f;
            const double overlapThreshold = 0;
            const double correlationThreshold = 0.75;
            int? fragmentRankThreshold = null; // e.g. 200 to cap fragments per pseudo scan

            // Short output stem: TopPIC rejects any spectrum path over 200 characters
            // (toppic_argument.cpp), and the full run stem plus this directory exceeds it.
            const string outputStem = "fxml_a0p1_c0p75";

            // ------------------------------------------------------------------ files
            // the parent file, the one with all four channels interleaved and the sid= scan filters
            string toppicDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string stem = @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected";
            string parentMzml = Path.Combine(toppicDir, stem + ".mzML");
            // precursors keep coming from the whole-file TopFD run, whose feature xml covers the MS1 scans
            string precursorFeatureXml = Path.Combine(toppicDir, stem + "_feature.xml");

            // fragments: one TopFD run per source voltage, discovered by the _isd<NN>_ in the file name
            string channelDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\HCD_YBC_ISDseparate";
            var fragmentFeatureXmlByVoltage = Directory.GetFiles(channelDir, "*_feature.xml")
                .Select(f => (File: f, Match: Regex.Match(Path.GetFileName(f), @"_isd(\d+)_feature\.xml$", RegexOptions.IgnoreCase)))
                .Where(t => t.Match.Success)
                .ToDictionary(t => int.Parse(t.Match.Groups[1].Value), t => t.File);

            string outDir = Path.Combine(Path.GetDirectoryName(channelDir), "PseudoMs2_featureXml");

            // ------------------------------------------------------------------ scans
            var dataFile = MsDataFileReader.GetDataFile(parentMzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // ------------------------------------------------------------------ precursor XICs
            var ms1Constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                precursorFeatureXml, minMass: precursorMinMass, maxMass: precursorMaxMass,
                minChargeCount: precursorMinChargeCount, minEcScore: precursorMinEcScore, aggregateCharges: true);
            var precursorXics = ms1Constructor.GetAllXics(ms1Scans);

            // ------------------------------------------------------------------ fragment XICs, one channel at a time
            var fragmentXics = new List<ExtractedIonChromatogram>();
            foreach (var channel in isdMap.OrderBy(k => k.Key))
            {
                var constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                    fragmentFeatureXmlByVoltage[channel.Key], minMass: fragmentMinMass, maxMass: fragmentMaxMass,
                    minChargeCount: fragmentMinChargeCount, minEcScore: fragmentMinEcScore,
                    aggregateCharges: aggregateFragmentCharges);

                var channelXics = constructor.GetAllXics(channel.Value.ToArray());
                fragmentXics.AddRange(channelXics);
            }

            // ------------------------------------------------------------------ grouping -> pseudo scans -> msAlign
            var grouper = new XicGroupingEngine(apexRtTolerance, overlapThreshold, correlationThreshold,
                maxThreadsForGrouping: Math.Max(1, Environment.ProcessorCount - 2),
                fragmentRankThreshold: fragmentRankThreshold);
            var pfGroups = grouper.PrecursorFragmentGrouping(precursorXics, fragmentXics).ToList();

            // stamp each group with the TopFD feature id behind its precursor, so the exported msAlign can be
            // joined straight back to the precursor feature xml
            PrecursorFragmentsGroup.AssignPrecursorFeatureIds(pfGroups, ms1Constructor.LastFeatureByXic);

            // PFgroupIndex doubles as the pseudo scan's one-based scan number, so assign it before building
            var commonParams = new CommonParameters();
            var pseudoScans = new List<Ms2ScanWithSpecificMass>(pfGroups.Count);
            for (int i = 0; i < pfGroups.Count; i++)
            {
                pfGroups[i].PFgroupIndex = i + 1;
                pseudoScans.Add(pfGroups[i].GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.Mass, commonParams, dataFile.FilePath));
            }

            Directory.CreateDirectory(outDir);
            string msAlignOut = Path.Combine(outDir, outputStem + "_ms2.msalign");
            ISDEngine.WriteMsAlignFile(msAlignOut, pseudoScans);

            var fragmentCounts = pseudoScans.Select(s => s.ExperimentalFragments.Length).OrderBy(n => n).ToArray();
            TestContext.Progress.WriteLine($"apex {apexRtTolerance} / corr {correlationThreshold}: " +
                $"{pseudoScans.Count} pseudo scans, fragments/precursor mean={fragmentCounts.Average():F1} " +
                $"median={fragmentCounts[fragmentCounts.Length / 2]}, total={fragmentCounts.Sum():N0}");
            TestContext.Progress.WriteLine("wrote " + msAlignOut);
        }

        /// <summary>
        /// Sweeps the precursor-fragment grouping thresholds and reports what each combination produces. The XICs
        /// are built ONCE and only the grouping is repeated, so a whole grid costs little more than a single run.
        ///
        /// IMPORTANT - what this can and cannot tell you. Nothing here measures identifications, so it cannot
        /// name an optimum on its own; the real criterion is confident proteoform IDs out of a search, which
        /// needs a search per condition. What it does give is the trade-off surface plus one internal quality
        /// signal, FRAGMENT PROMISCUITY: the mean number of distinct precursors a fragment gets assigned to.
        /// A fragment ion physically belongs to one proteoform, so a setting that hands the same fragment to
        /// dozens of precursors is buying fragment count with co-elution coincidence. Reading the table, prefer
        /// the loosest setting that still keeps promiscuity low and leaves most precursors with enough fragments
        /// to score - then confirm the shortlist with real searches.
        /// </summary>
        [Test, Explicit]
        public static void ScreenXicGroupingParameters()
        {
            // ------------------------------------------------------------------ grid
            float[] apexRtTolerances = { 0.05f, 0.1f, 0.2f };
            double[] correlationThresholds = { 0, 0.3, 0.5, 0.7, 0.9 };
            const double overlapThreshold = 0;
            const int identifiableFragmentCount = 10; // "enough fragments to be worth searching"

            const double precursorMinMass = 4000;
            const double precursorMaxMass = 70000;
            const int minNumberOfPeaks = 3;
            const bool aggregateFragmentCharges = true;

            // ------------------------------------------------------------------ files (same set as the test above)
            string toppicDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string stem = @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected";
            string parentMzml = Path.Combine(toppicDir, stem + ".mzML");
            string precursorFeatureXml = Path.Combine(toppicDir, stem + "_feature.xml");
            string channelDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\HCD_YBC_ISDseparate";
            var fragmentFeatureXmlByVoltage = Directory.GetFiles(channelDir, "*_feature.xml")
                .Select(f => (File: f, Match: Regex.Match(Path.GetFileName(f), @"_isd(\d+)_feature\.xml$", RegexOptions.IgnoreCase)))
                .Where(t => t.Match.Success)
                .ToDictionary(t => int.Parse(t.Match.Groups[1].Value), t => t.File);

            // ------------------------------------------------------------------ build the XICs once
            var dataFile = MsDataFileReader.GetDataFile(parentMzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            var ms1Constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                precursorFeatureXml, minMass: precursorMinMass, maxMass: precursorMaxMass);
            var precursorXics = ms1Constructor.GetAllXics(ms1Scans);

            var fragmentXics = new List<ExtractedIonChromatogram>();
            foreach (var channel in isdMap.OrderBy(k => k.Key))
            {
                fragmentXics.AddRange(new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                        fragmentFeatureXmlByVoltage[channel.Key], maxMass: 70000,
                        aggregateCharges: aggregateFragmentCharges)
                    .GetAllXics(channel.Value.ToArray()));
            }
            TestContext.Progress.WriteLine($"built once: {precursorXics.Count} precursor and {fragmentXics.Count} fragment XICs\n");

            // ------------------------------------------------------------------ sweep
            TestContext.Progress.WriteLine(
                "  apex  corr | groups  medFrag  meanFrag    total | %>=10frag | fragUsed  promiscuity");
            TestContext.Progress.WriteLine(
                "  -----------|--------------------------------------|-----------|----------------------");

            foreach (var apex in apexRtTolerances)
            {
                foreach (var correlation in correlationThresholds)
                {
                    var grouper = new XicGroupingEngine(apex, overlapThreshold, correlation,
                        maxThreadsForGrouping: Math.Max(1, Environment.ProcessorCount - 2));
                    var groups = grouper.PrecursorFragmentGrouping(precursorXics, fragmentXics).ToList();
                    if (groups.Count == 0)
                    {
                        TestContext.Progress.WriteLine($"  {apex,5:F2}  {correlation,4:F1} | no groups");
                        continue;
                    }

                    var counts = groups.Select(g => g.PFpairs.Count).OrderBy(n => n).ToArray();
                    long total = counts.Sum(n => (long)n);

                    // how many distinct precursors does each fragment end up in? reference identity is enough,
                    // since every group holds references into the one fragment XIC list built above
                    var precursorsPerFragment = new Dictionary<ExtractedIonChromatogram, int>();
                    foreach (var group in groups)
                        foreach (var pair in group.PFpairs)
                            precursorsPerFragment[pair.FragmentXic] =
                                precursorsPerFragment.TryGetValue(pair.FragmentXic, out int c) ? c + 1 : 1;

                    double identifiable = 100.0 * counts.Count(n => n >= identifiableFragmentCount) / counts.Length;
                    double promiscuity = precursorsPerFragment.Count == 0 ? 0 : precursorsPerFragment.Values.Average();

                    TestContext.Progress.WriteLine(
                        $"  {apex,5:F2}  {correlation,4:F1} | {groups.Count,6} {counts[counts.Length / 2],8} " +
                        $"{counts.Average(),9:F1} {total,8:N0} | {identifiable,8:F1}% | {precursorsPerFragment.Count,8} " +
                        $"{promiscuity,12:F1}");
                }
            }
        }

        /// <summary>
        /// The screen that actually answers the question: for each grouping setting, build the pseudo-MS2 scans,
        /// search them in MetaMorpheus, and count confident proteoform identifications. Unlike
        /// <see cref="ScreenXicGroupingParameters"/>, this measures the thing being optimised rather than a proxy.
        ///
        /// One search per condition, so keep <see cref="Conditions"/> short - this is hours, not minutes. The XICs
        /// are still built only once and shared across conditions.
        /// </summary>
        [Test, Explicit]
        public static void ScreenXicGroupingParametersBySearch()
        {
            // ------------------------------------------------------------------ conditions to compare
            // (apex RT tolerance, correlation threshold). First entry is the incumbent setting.
            var conditions = new (float Apex, double Correlation)[]
            {
                // Rounds 1-2: IDs tracked mean fragment count almost perfectly (r ~ 0.99) regardless of which
                // knob produced it, best so far (0.2, 0.0) -> 110 proteoforms. Round 3 fills in the apex axis at
                // correlation 0, which is also a test of that relationship: if fragment count is all that
                // matters, these should order 0.05 < 0.1 < 0.2 < 0.3 < 0.5 by both fragments AND proteoforms.
                (0.1f, 0.0),
                (0.1f, 0.3),
                (0.3f, 0.0),
                (0.5f, 0.0),
            };
            const double overlapThreshold = 0;
            const double qValueCutoff = 0.01;

            const double precursorMinMass = 4000;
            const double precursorMaxMass = 70000;
            const int minNumberOfPeaks = 3;
            const bool aggregateFragmentCharges = true;

            // ------------------------------------------------------------------ files
            string toppicDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string stem = @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected";
            string parentMzml = Path.Combine(toppicDir, stem + ".mzML");
            string precursorFeatureXml = Path.Combine(toppicDir, stem + "_feature.xml");
            string channelDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\HCD_YBC_ISDseparate";
            string searchToml = @"E:\ISD Project\071626_ISD\Quant\HC_YC\MM\DDA_1-1\Task Settings\Task1-SearchTaskconfig.toml";
            string databaseXml = @"E:\Databases\uniprotkb_taxonomy_id_9606_AND_reviewed_2026_07_23.xml";
            string screenRoot = Path.Combine(Path.GetDirectoryName(channelDir), "GroupingScreen");

            foreach (var f in new[] { parentMzml, precursorFeatureXml, searchToml, databaseXml })
                Assert.That(File.Exists(f), $"missing {f}");

            var fragmentFeatureXmlByVoltage = Directory.GetFiles(channelDir, "*_feature.xml")
                .Select(f => (File: f, Match: Regex.Match(Path.GetFileName(f), @"_isd(\d+)_feature\.xml$", RegexOptions.IgnoreCase)))
                .Where(t => t.Match.Success)
                .ToDictionary(t => int.Parse(t.Match.Groups[1].Value), t => t.File);

            // ------------------------------------------------------------------ XICs, built once
            var dataFile = MsDataFileReader.GetDataFile(parentMzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // held in a variable so its LastFeatureByXic map is available for the feature-id stamping below
            var ms1Constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                precursorFeatureXml, minMass: precursorMinMass, maxMass: precursorMaxMass);
            var precursorXics = ms1Constructor.GetAllXics(ms1Scans);
            var fragmentXics = new List<ExtractedIonChromatogram>();
            foreach (var channel in isdMap.OrderBy(k => k.Key))
                fragmentXics.AddRange(new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                        fragmentFeatureXmlByVoltage[channel.Key], maxMass: 70000,
                        aggregateCharges: aggregateFragmentCharges)
                    .GetAllXics(channel.Value.ToArray()));
            TestContext.Progress.WriteLine($"built once: {precursorXics.Count} precursor and {fragmentXics.Count} fragment XICs");

            // ------------------------------------------------------------------ one search per condition
            var results = new List<string>();
            foreach (var (apex, correlation) in conditions)
            {
                string label = $"apex{apex:0.00}_corr{correlation:0.0}".Replace('.', 'p');
                string conditionDir = Path.Combine(screenRoot, label);
                Directory.CreateDirectory(conditionDir);

                var grouper = new XicGroupingEngine(apex, overlapThreshold, correlation,
                    maxThreadsForGrouping: Math.Max(1, Environment.ProcessorCount - 2));
                var pfGroups = grouper.PrecursorFragmentGrouping(precursorXics, fragmentXics).ToList();
                PrecursorFragmentsGroup.AssignPrecursorFeatureIds(pfGroups, ms1Constructor.LastFeatureByXic);

                var commonParams = new CommonParameters();
                var pseudoScans = new List<Ms2ScanWithSpecificMass>(pfGroups.Count);
                for (int i = 0; i < pfGroups.Count; i++)
                {
                    pfGroups[i].PFgroupIndex = i + 1;
                    pseudoScans.Add(pfGroups[i].GetPseudoMs2ScanFromPfGroup(
                        PseudoMs2ConstructionType.Mass, commonParams, dataFile.FilePath));
                }

                string mgfPath = Path.Combine(conditionDir, stem + "_pseudo_ms2.mgf");
                ISDEngine.WriteMsAlignFile(Path.Combine(conditionDir, stem + "_pseudo_ms2.msalign"), pseudoScans);
                IsdMsAlignExporter.WriteMgf(pseudoScans, mgfPath);

                double meanFragments = pfGroups.Count == 0 ? 0 : pfGroups.Average(g => g.PFpairs.Count);
                TestContext.Progress.WriteLine($"\n[{label}] {pfGroups.Count} groups, mean {meanFragments:F1} fragments " +
                    "-> searching ...");

                // the pseudo scans are already deconvoluted neutral masses, so the search must not redo it
                var searchTask = Toml.ReadFile<SearchTask>(searchToml, MetaMorpheusTask.tomlConfig);
                searchTask.CommonParameters.DoPrecursorDeconvolution = false;
                searchTask.CommonParameters.UseProvidedPrecursorInfo = true;
                searchTask.CommonParameters.ProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 1, 4, 3);
                searchTask.CommonParameters.AssumeOrphanPeaksAreZ1Fragments = true;

                var engine = new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("search", searchTask) },
                    new List<string> { mgfPath },
                    new List<DbForTask> { new DbForTask(databaseXml, false) },
                    conditionDir);
                engine.Run();

                int proteoforms = CountTargetsBelowQValue(conditionDir, "AllProteoforms.psmtsv", qValueCutoff);
                int psms = CountTargetsBelowQValue(conditionDir, "AllPSMs.psmtsv", qValueCutoff);
                string row = $"  {apex,5:F2}  {correlation,4:F1} | {pfGroups.Count,6} groups | {meanFragments,7:F1} frag " +
                             $"| {psms,7:N0} PSMs | {proteoforms,7:N0} proteoforms";
                results.Add(row);
                TestContext.Progress.WriteLine($"[{label}] -> {psms:N0} PSMs, {proteoforms:N0} proteoforms at q<={qValueCutoff}");
            }

            TestContext.Progress.WriteLine($"\n=== proteoform IDs by grouping setting (q <= {qValueCutoff}) ===");
            TestContext.Progress.WriteLine("  apex  corr | groups | mean frag |    PSMs | proteoforms");
            foreach (var row in results) TestContext.Progress.WriteLine(row);
        }

        /// <summary>
        /// Counts target (non-decoy, non-contaminant) rows at or below a q-value in a MetaMorpheus result file,
        /// searching the whole task output tree because the task folder name is generated by the runner.
        /// </summary>
        private static int CountTargetsBelowQValue(string searchRoot, string fileName, double qValueCutoff)
        {
            var file = Directory.GetFiles(searchRoot, fileName, SearchOption.AllDirectories).FirstOrDefault();
            if (file == null)
                return -1;

            var lines = File.ReadAllLines(file);
            if (lines.Length < 2) return 0;
            var header = lines[0].Split('\t');
            int qIndex = Array.FindIndex(header, h => h.Equals("QValue", StringComparison.OrdinalIgnoreCase));
            int targetIndex = Array.FindIndex(header, h => h.Equals("Decoy/Contaminant/Target", StringComparison.OrdinalIgnoreCase));
            if (qIndex < 0 || targetIndex < 0) return -1;

            int count = 0;
            foreach (var line in lines.Skip(1))
            {
                if (string.IsNullOrWhiteSpace(line)) continue;
                var splits = line.Split('\t');
                if (splits.Length <= Math.Max(qIndex, targetIndex)) continue;
                if (!splits[targetIndex].Trim().Equals("T", StringComparison.OrdinalIgnoreCase)) continue;
                if (double.TryParse(splits[qIndex], System.Globalization.NumberStyles.Any,
                        System.Globalization.CultureInfo.InvariantCulture, out double q) && q <= qValueCutoff)
                    count++;
            }
            return count;
        }
    }
}
