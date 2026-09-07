using EngineLayer.DIA.FeatureAlign;
using NUnit.Framework;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Self-contained proofs that <see cref="FeatureAligner"/> aligns features correctly.
    ///
    /// Unlike the sibling *_FeatureAlignTest fixtures — which are [Explicit] runners over local TopFD
    /// output and only report row counts — every test here builds tiny feature tables whose correct
    /// alignment is known by construction, runs the aligner, and asserts the consensus table exactly.
    /// Each test isolates ONE property of the matching rule (a match = mass within tolerance around a
    /// monoisotope notch AND apex RT within tolerance), so a regression names the property it broke.
    /// </summary>
    [TestFixture]
    public static class FeatureAlignerUnitTests
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

        // ================================================================================================
        //  POSITIVE CONTROL — features that ARE the same proteoform collapse into one consensus row
        // ================================================================================================

        [Test]
        public static void SameMassAndRt_AlignAcrossRuns_IntoOneConsensusRow()
        {
            // Two runs, each with one feature at (essentially) the same mass and RT: 5 ppm and 10 s apart,
            // both inside 20 ppm / 30 s. They must become ONE row spanning both files, and the consensus
            // mass/RT must be the mean of the two.
            string dir = FreshDir();
            string a = WriteFlashQuant(dir, "runA", (idx: "1", mass: 10000.00, rtSec: 1000.0, quant: 100, zMin: 10, zMax: 12));
            string b = WriteFlashQuant(dir, "runB", (idx: "7", mass: 10000.05, rtSec: 1010.0, quant: 200, zMin: 11, zMax: 13));
            string outp = Path.Combine(dir, "consensus.tsv");

            FeatureAligner.Align(new[] { a, b }, outp, massTolerance: 20, rtTolerance: 30);

            var rows = ReadConsensus(outp);
            Assert.That(rows, Has.Count.EqualTo(1), "the two identical features should collapse to a single consensus row");

            var row = rows[0];
            Assert.That(Int(row, "n_files"), Is.EqualTo(2), "row should span both runs");
            Assert.That(Dbl(row, "mass_mean"), Is.EqualTo(10000.025).Within(1e-6), "consensus mass is the mean of the two feature masses");
            Assert.That(Dbl(row, "rt_mean"), Is.EqualTo(1005.0).Within(1e-6), "consensus RT is the mean of the two apex RTs");
            Assert.That(Int(row, "max_features_per_run"), Is.EqualTo(1), "one feature per run = a clean group");

            // Each run's abundance and original feature index land in that run's own columns.
            Assert.That(Dbl(row, "quant_runA"), Is.EqualTo(100).Within(1e-9));
            Assert.That(Dbl(row, "quant_runB"), Is.EqualTo(200).Within(1e-9));
            Assert.That(row["idx_runA"], Is.EqualTo("1"));
            Assert.That(row["idx_runB"], Is.EqualTo("7"));

            // Charge range spans the union of the two runs' observed charges.
            Assert.That(row["charge_range"], Is.EqualTo("10-13"));
        }

        // ================================================================================================
        //  RETENTION-TIME GATE — RT must be within tolerance, even when the mass matches perfectly
        // ================================================================================================

        [Test]
        public static void RetentionTimeGate_SeparatesFeaturesBeyondRtTolerance()
        {
            // Identical mass in both runs; only the RT differs. Inside the RT window -> align; outside -> not.
            // Proves RT is a required, independent gate and that the boundary is inclusive.
            AssertPairAlignsOrNot(rtGapSeconds: 29, rtTol: 30, expectAligned: true,  "29 s < 30 s tolerance must align");
            AssertPairAlignsOrNot(rtGapSeconds: 30, rtTol: 30, expectAligned: true,  "30 s == tolerance is inclusive, must align");
            AssertPairAlignsOrNot(rtGapSeconds: 31, rtTol: 30, expectAligned: false, "31 s > 30 s tolerance must NOT align");
        }

        // ================================================================================================
        //  MASS GATE — mass must be within tolerance, even when the RT matches perfectly
        // ================================================================================================

        [Test]
        public static void MassGate_SeparatesFeaturesBeyondMassTolerance()
        {
            // Same RT in both runs; only the mass differs. 20 ppm at 10000 Da is 0.2 Da.
            string dir = FreshDir();
            string a = WriteFlashQuant(dir, "runA", (idx: "1", mass: 10000.000, rtSec: 1000, quant: 100, zMin: 0, zMax: 0));
            // 10000.10 -> 10 ppm (inside); 10000.30 -> 30 ppm (outside).
            string bIn  = WriteFlashQuant(dir, "runIn",  (idx: "1", mass: 10000.100, rtSec: 1000, quant: 100, zMin: 0, zMax: 0));
            string bOut = WriteFlashQuant(dir, "runOut", (idx: "1", mass: 10000.300, rtSec: 1000, quant: 100, zMin: 0, zMax: 0));

            string inOut = Path.Combine(dir, "in.tsv");
            FeatureAligner.Align(new[] { a, bIn }, inOut, massTolerance: 20, rtTolerance: 30);
            Assert.That(ReadConsensus(inOut).Single()["n_files"], Is.EqualTo("2"), "10 ppm apart must align");

            string outOut = Path.Combine(dir, "out.tsv");
            FeatureAligner.Align(new[] { a, bOut }, outOut, massTolerance: 20, rtTolerance: 30);
            var rows = ReadConsensus(outOut);
            Assert.That(rows, Has.Count.EqualTo(2), "30 ppm apart exceeds 20 ppm, must stay two rows");
            Assert.That(rows.All(r => r["n_files"] == "1"), Is.True, "each unmatched feature spans exactly one run");
        }

        // ================================================================================================
        //  MONOISOTOPE NOTCH — bridges ±1 Da monoisotope-assignment disagreements only when enabled
        // ================================================================================================

        [Test]
        public static void Notch_BridgesOneDaltonOffset_OnlyWhenNotchAllowed()
        {
            // Two features one C13 step apart (~1.0034 Da). At 10000 Da that is ~100 ppm, far outside 20 ppm,
            // so with notch 0 they are two proteoforms; with +1 notch the aligner recognises the +1 monoisotope
            // and merges them. This is the behaviour that lets a deconvolutor's off-by-one calls still align.
            string dir = FreshDir();
            string a = WriteFlashQuant(dir, "runA", (idx: "1", mass: 10000.0000000, rtSec: 1000, quant: 100, zMin: 0, zMax: 0));
            string b = WriteFlashQuant(dir, "runB", (idx: "1", mass: 10000.0000000 + FeatureAligner.NotchStep, rtSec: 1005, quant: 100, zMin: 0, zMax: 0));

            string noNotch = Path.Combine(dir, "notch0.tsv");
            FeatureAligner.Align(new[] { a, b }, noNotch, massTolerance: 20, rtTolerance: 30, positiveNotch: 0, negativeNotch: 0);
            Assert.That(ReadConsensus(noNotch), Has.Count.EqualTo(2), "1 Da apart, notch 0: must stay separate (~100 ppm)");

            string withNotch = Path.Combine(dir, "notch1.tsv");
            FeatureAligner.Align(new[] { a, b }, withNotch, massTolerance: 20, rtTolerance: 30, positiveNotch: 1, negativeNotch: 0);
            var rows = ReadConsensus(withNotch);
            Assert.That(rows, Has.Count.EqualTo(1), "1 Da apart, +1 notch: must merge as the same proteoform");
            Assert.That(rows[0]["n_files"], Is.EqualTo("2"));
        }

        // ================================================================================================
        //  MULTI-RUN — shared and run-specific proteoforms get the right n_files and column placement
        // ================================================================================================

        [Test]
        public static void ThreeRuns_SharedAndUniqueProteoforms_ReportCorrectOccupancy()
        {
            // One proteoform present in all three runs, plus one feature unique to each run.
            // Expect: the shared one has n_files == 3 with all three quant columns filled; each unique one
            // has n_files == 1 with only its own run's quant filled and the other runs blank.
            string dir = FreshDir();
            string a = WriteFlashQuant(dir, "r1",
                (idx: "S", mass: 10000.00, rtSec: 1000, quant: 111, zMin: 0, zMax: 0),
                (idx: "U1", mass: 12000.00, rtSec: 1500, quant: 10, zMin: 0, zMax: 0));
            string b = WriteFlashQuant(dir, "r2",
                (idx: "S", mass: 10000.02, rtSec: 1005, quant: 222, zMin: 0, zMax: 0),
                (idx: "U2", mass: 14000.00, rtSec: 2000, quant: 20, zMin: 0, zMax: 0));
            string c = WriteFlashQuant(dir, "r3",
                (idx: "S", mass: 9999.98, rtSec: 995, quant: 333, zMin: 0, zMax: 0),
                (idx: "U3", mass: 16000.00, rtSec: 2500, quant: 30, zMin: 0, zMax: 0));

            string outp = Path.Combine(dir, "consensus.tsv");
            FeatureAligner.Align(new[] { a, b, c }, outp, massTolerance: 20, rtTolerance: 30);

            var rows = ReadConsensus(outp);
            Assert.That(rows, Has.Count.EqualTo(4), "one shared + three unique proteoforms = 4 rows");

            var shared = RowNear(rows, 10000.0);
            Assert.That(Int(shared, "n_files"), Is.EqualTo(3), "the shared proteoform spans all three runs");
            Assert.That(Dbl(shared, "quant_r1"), Is.EqualTo(111).Within(1e-9));
            Assert.That(Dbl(shared, "quant_r2"), Is.EqualTo(222).Within(1e-9));
            Assert.That(Dbl(shared, "quant_r3"), Is.EqualTo(333).Within(1e-9));

            var unique = RowNear(rows, 12000.0);
            Assert.That(Int(unique, "n_files"), Is.EqualTo(1), "a run-specific feature spans one run");
            Assert.That(Dbl(unique, "quant_r1"), Is.EqualTo(10).Within(1e-9));
            Assert.That(row_blank(unique, "quant_r2"), Is.True, "runs without this feature leave the quant cell blank");
            Assert.That(row_blank(unique, "quant_r3"), Is.True);
        }

        // ================================================================================================
        //  SAME-RUN DUPLICATE — a run reporting one species twice is flagged, not double-counted
        // ================================================================================================

        [Test]
        public static void SameRunDuplicate_FlagsMaxPerRun_AndKeepsMostAbundant()
        {
            // runA reports the SAME proteoform twice (a deconvolution artifact): two features within tolerance
            // of each other in ONE run. runB reports it once. The consensus must:
            //   - keep runA's MOST ABUNDANT copy as runA's quantity (not sum, not the weaker one),
            //   - set max_features_per_run = 2 so the artifact is auditable,
            //   - still count n_files = 2 (two distinct runs), not 3.
            string dir = FreshDir();
            string a = WriteFlashQuant(dir, "runA",
                (idx: "weak",   mass: 10000.00, rtSec: 1000, quant: 100, zMin: 0, zMax: 0),
                (idx: "strong", mass: 10000.05, rtSec: 1002, quant: 500, zMin: 0, zMax: 0));
            string b = WriteFlashQuant(dir, "runB",
                (idx: "1", mass: 10000.02, rtSec: 1001, quant: 300, zMin: 0, zMax: 0));

            string outp = Path.Combine(dir, "consensus.tsv");
            FeatureAligner.Align(new[] { a, b }, outp, massTolerance: 20, rtTolerance: 30);

            var row = ReadConsensus(outp).Single();
            Assert.That(Int(row, "n_files"), Is.EqualTo(2), "two distinct runs, despite runA's duplicate");
            Assert.That(Int(row, "max_features_per_run"), Is.EqualTo(2), "runA's same-run duplicate must be flagged");
            Assert.That(Dbl(row, "quant_runA"), Is.EqualTo(500).Within(1e-9), "the more abundant of runA's two copies is kept");
            Assert.That(row["idx_runA"], Is.EqualTo("strong"), "and its index, not the weaker copy's");
            Assert.That(Dbl(row, "quant_runB"), Is.EqualTo(300).Within(1e-9));
        }

        // ================================================================================================
        //  RT UNIT CONVERSION — a TopFD table (Apex_time in MINUTES) is read in seconds like everything else
        // ================================================================================================

        [Test]
        public static void TopFdMinutes_AreConvertedToSeconds_BeforeMatching()
        {
            // TopFD writes Apex_time in MINUTES. Two features 1.0 minute apart are 60 s apart after conversion,
            // which exceeds a 30 s RT tolerance -> they must NOT align. If the aligner forgot to convert, it
            // would see "1.0" vs "0.0" ~ 1 unit apart, well inside 30, and wrongly merge them. So a NON-match
            // here proves the minutes->seconds conversion happened.
            string dir = FreshDir();
            string a = WriteTopFd(dir, "runA", (idx: "1", mass: 10000.00, rtMin: 16.0000, apexInt: 100, zMin: 10, zMax: 12));
            string bFar = WriteTopFd(dir, "runFar", (idx: "1", mass: 10000.05, rtMin: 17.0000, apexInt: 100, zMin: 10, zMax: 12)); // +1.0 min = +60 s

            string far = Path.Combine(dir, "far.tsv");
            FeatureAligner.Align(new[] { a, bFar }, far, massTolerance: 20, rtTolerance: 30);
            Assert.That(ReadConsensus(far), Has.Count.EqualTo(2), "60 s apart (1 min) exceeds 30 s -> conversion applied, no merge");

            // And a 0.25 min (= 15 s) gap DOES align, confirming the converted RT is what is compared.
            string bNear = WriteTopFd(dir, "runNear", (idx: "1", mass: 10000.05, rtMin: 16.2500, apexInt: 100, zMin: 10, zMax: 12));
            string near = Path.Combine(dir, "near.tsv");
            FeatureAligner.Align(new[] { a, bNear }, near, massTolerance: 20, rtTolerance: 30);
            var row = ReadConsensus(near).Single();
            Assert.That(Int(row, "n_files"), Is.EqualTo(2), "15 s apart (0.25 min) is within 30 s -> merge");
            // rt_mean should be the mean IN SECONDS: (16.0*60 + 16.25*60)/2 = (960 + 975)/2 = 967.5 s.
            Assert.That(Dbl(row, "rt_mean"), Is.EqualTo(967.5).Within(1e-6), "consensus RT is in seconds, i.e. minutes were scaled x60");
        }

        // ================================================================================================
        //  GROUPING CONTRACT — ReferenceAnchor caps a group at the tolerance; UnionFind chains a mass ladder
        // ================================================================================================

        [Test]
        public static void GroupingMethods_DifferOnAMassLadder_AsDocumented()
        {
            // Four features in four runs at the same RT, each 15 ppm above the previous (so consecutive pairs
            // match at 20 ppm, but F1..F4 span 45 ppm). This is the ladder the class comments describe.
            //   UnionFind: transitive links F1-F2-F3-F4 -> ONE cluster of all four runs.
            //   ReferenceAnchor: the anchor's ±20 ppm window cannot reach both ends -> the ladder is split.
            // Abundances make F2 the anchor deterministically (highest), so ReferenceAnchor yields {F1,F2,F3}
            // then {F4}.
            string dir = FreshDir();
            const double m = 10000.0;
            double step = m * 15e-6; // 15 ppm in Da
            string a = WriteFlashQuant(dir, "rA", (idx: "1", mass: m + 0 * step, rtSec: 1000, quant: 100, zMin: 0, zMax: 0));
            string b = WriteFlashQuant(dir, "rB", (idx: "1", mass: m + 1 * step, rtSec: 1000, quant: 1000, zMin: 0, zMax: 0)); // anchor
            string c = WriteFlashQuant(dir, "rC", (idx: "1", mass: m + 2 * step, rtSec: 1000, quant: 100, zMin: 0, zMax: 0));
            string d = WriteFlashQuant(dir, "rD", (idx: "1", mass: m + 3 * step, rtSec: 1000, quant: 500, zMin: 0, zMax: 0));
            var files = new[] { a, b, c, d };

            string uf = Path.Combine(dir, "unionfind.tsv");
            FeatureAligner.Align(files, uf, massTolerance: 20, rtTolerance: 30, grouping: GroupingMethod.UnionFind);
            var ufRows = ReadConsensus(uf);
            Assert.That(ufRows, Has.Count.EqualTo(1), "UnionFind chains the 45 ppm ladder into one cluster");
            Assert.That(Int(ufRows[0], "n_files"), Is.EqualTo(4), "and it spans all four runs");

            string ra = Path.Combine(dir, "refanchor.tsv");
            FeatureAligner.Align(files, ra, massTolerance: 20, rtTolerance: 30, grouping: GroupingMethod.ReferenceAnchor);
            var raRows = ReadConsensus(ra);
            Assert.That(raRows, Has.Count.EqualTo(2), "ReferenceAnchor caps the group at the tolerance width, splitting the ladder");
            Assert.That(raRows.Max(r => Int(r, "n_files")), Is.EqualTo(3), "the anchor gathers three consecutive runs");
            Assert.That(raRows.Min(r => Int(r, "n_files")), Is.EqualTo(1), "the far end (F4) falls into its own group");
        }

        // ================================================================================================
        //  FEATURE -> IDENTIFICATION — the second entry point matches a feature to a psmtsv by mass + RT
        // ================================================================================================

        [Test]
        public static void AlignToIdentifications_MatchesByMassAndRt_ConvertsIdMinutes_AndFiltersByQValue()
        {
            // One feature; a psmtsv with two identifications:
            //   - a real match (0.5 ppm, RT equal after the ID's minutes are converted to seconds, QValue 0.001)
            //   - a decoy-ish high-QValue row that must be filtered out before matching.
            // Proves: mass+RT matching, ID minutes->seconds conversion, the QValue threshold, and that the ID's
            // metadata (accession/sequence) is written onto the feature row.
            string dir = FreshDir();
            string feat = WriteFlashQuant(dir, "feat", (idx: "1", mass: 10000.000, rtSec: 1000.0, quant: 999, zMin: 10, zMax: 15));

            // Feature RT is 1000 s = 16.66667 min. The good ID is given at 16.66667 min so it matches only if
            // the aligner converts minutes -> seconds. The bad ID is at the same mass/RT but QValue 0.5.
            string psm = Path.Combine(dir, "ids.psmtsv");
            File.WriteAllLines(psm, new[]
            {
                string.Join("\t", "Precursor Mass", "Scan Retention Time", "Full Sequence", "Accession", "Gene Name", "Name", "QValue", "Precursor Charge"),
                string.Join("\t", (10000.0 + 0.005).ToString(Ci), (1000.0 / 60.0).ToString(Ci), "PEPTIDEK", "P12345", "GENE1", "prot1", "0.001", "12"),
                string.Join("\t", (10000.0 + 0.005).ToString(Ci), (1000.0 / 60.0).ToString(Ci), "DECOYSEQ", "XXXXX",  "GENEX", "protX", "0.500", "12"),
            });

            string outp = Path.Combine(dir, "feat_to_id.tsv");
            FeatureAligner.AlignToIdentifications(feat, psm, outp, massTolerance: 20, rtTolerance: 30, qValueThreshold: 0.01);

            var rows = ReadTsv(outp);
            Assert.That(rows, Has.Count.EqualTo(1), "one row per feature");
            var r = rows[0];
            Assert.That(r["num_id_matches"], Is.EqualTo("1"), "only the QValue<=0.01 identification is a candidate");
            Assert.That(r["accession"], Is.EqualTo("P12345"), "the surviving ID's accession is written onto the feature");
            Assert.That(r["full_sequence"], Is.EqualTo("PEPTIDEK"));
            Assert.That(double.Parse(r["qvalue"], Ci), Is.EqualTo(0.001).Within(1e-9));
            // The ID is 0.005 Da HEAVIER than the feature, so the signed residual (feature - id)/id is
            // -0.5 ppm. Asserting the sign as well as the magnitude proves the reported error is directional.
            Assert.That(double.Parse(r["mass_error_ppm"], Ci), Is.EqualTo(-0.5).Within(0.05), "0.005 Da at 10000 = -0.5 ppm residual (feature lighter than ID)");
        }

        [Test]
        public static void AlignToIdentifications_LeavesUnmatchedFeatureBlank()
        {
            // A feature whose mass is nowhere near any ID: the row is emitted with num_id_matches 0 and empty
            // identification columns (the feature is not silently dropped).
            string dir = FreshDir();
            string feat = WriteFlashQuant(dir, "feat", (idx: "1", mass: 5000.0, rtSec: 1000.0, quant: 999, zMin: 0, zMax: 0));
            string psm = Path.Combine(dir, "ids.psmtsv");
            File.WriteAllLines(psm, new[]
            {
                string.Join("\t", "Precursor Mass", "Scan Retention Time", "Full Sequence", "Accession", "QValue"),
                string.Join("\t", "20000.0", (1000.0 / 60.0).ToString(Ci), "SEQ", "P1", "0.001"),
            });

            string outp = Path.Combine(dir, "feat_to_id.tsv");
            FeatureAligner.AlignToIdentifications(feat, psm, outp, massTolerance: 20, rtTolerance: 30);

            var r = ReadTsv(outp).Single();
            Assert.That(r["num_id_matches"], Is.EqualTo("0"));
            Assert.That(r["accession"], Is.Empty, "no ID -> the accession column is blank");
            Assert.That(r["full_sequence"], Is.Empty);
        }

        // ================================================================================================
        //  ALIGN-OF-ALIGN — a consensus result file is itself readable, so aligned features can be re-aligned
        // ================================================================================================

        [Test]
        public static void ConsensusResultFiles_CanBeReAligned_AveragingMassRtAndSummingQuant()
        {
            // Build two independent consensus tables (each from two runs), then align the two consensus tables to
            // each other. This exercises reading the aligner's OWN output: mass_mean/rt_mean become the feature's
            // mass/RT, and the per-run quant_<run> columns are summed into its abundance.
            string dir = FreshDir();

            // condition A: two runs -> consensusA.tsv (shared feature S + unique UA)
            string a1 = WriteFlashQuant(dir, "a1",
                (idx: "1", mass: 10000.00, rtSec: 1000, quant: 100, zMin: 0, zMax: 0),
                (idx: "2", mass: 12000.00, rtSec: 1500, quant: 10, zMin: 0, zMax: 0));
            string a2 = WriteFlashQuant(dir, "a2",
                (idx: "1", mass: 10000.02, rtSec: 1005, quant: 200, zMin: 0, zMax: 0),
                (idx: "2", mass: 12000.01, rtSec: 1500, quant: 20, zMin: 0, zMax: 0));
            string consensusA = Path.Combine(dir, "consensusA.tsv");
            FeatureAligner.Align(new[] { a1, a2 }, consensusA, massTolerance: 20, rtTolerance: 30);

            // condition B: two runs -> consensusB.tsv (same shared feature S + unique UB)
            string b1 = WriteFlashQuant(dir, "b1",
                (idx: "1", mass: 10000.05, rtSec: 1010, quant: 1000, zMin: 0, zMax: 0),
                (idx: "2", mass: 14000.00, rtSec: 2000, quant: 5, zMin: 0, zMax: 0));
            string b2 = WriteFlashQuant(dir, "b2",
                (idx: "1", mass: 10000.03, rtSec: 1008, quant: 2000, zMin: 0, zMax: 0),
                (idx: "2", mass: 14000.01, rtSec: 2000, quant: 6, zMin: 0, zMax: 0));
            string consensusB = Path.Combine(dir, "consensusB.tsv");
            FeatureAligner.Align(new[] { b1, b2 }, consensusB, massTolerance: 20, rtTolerance: 30);

            // align the two consensus tables (align aligned features)
            string outp = Path.Combine(dir, "consensus_of_consensus.tsv");
            FeatureAligner.Align(new[] { consensusA, consensusB }, outp, massTolerance: 20, rtTolerance: 30);

            var rows = ReadConsensus(outp);
            Assert.That(rows, Has.Count.EqualTo(3), "shared S + unique UA + unique UB");

            var shared = RowNear(rows, 10000.0);
            Assert.That(Int(shared, "n_files"), Is.EqualTo(2), "S is present in both consensus tables");
            Assert.That(Dbl(shared, "mass_mean"), Is.EqualTo(10000.025).Within(0.02), "mass is the average of the two consensus masses");
            // Each side's abundance is the SUM of that consensus's per-run quant columns: A = 100+200 = 300, B = 1000+2000 = 3000.
            Assert.That(Dbl(shared, "quant_consensusA"), Is.EqualTo(300).Within(1e-6));
            Assert.That(Dbl(shared, "quant_consensusB"), Is.EqualTo(3000).Within(1e-6));

            var ua = RowNear(rows, 12000.0);
            Assert.That(Int(ua, "n_files"), Is.EqualTo(1), "UA is only in consensusA");
            Assert.That(Dbl(ua, "quant_consensusA"), Is.EqualTo(30).Within(1e-6), "10 + 20 summed");
            Assert.That(row_blank(ua, "quant_consensusB"), Is.True);

            var ub = RowNear(rows, 14000.0);
            Assert.That(Int(ub, "n_files"), Is.EqualTo(1), "UB is only in consensusB");
            Assert.That(Dbl(ub, "quant_consensusB"), Is.EqualTo(11).Within(1e-6), "5 + 6 summed");
            Assert.That(row_blank(ub, "quant_consensusA"), Is.True);
        }

        // ================================================================================================
        //  GUARD RAILS — the aligner refuses under-specified input rather than producing a wrong table
        // ================================================================================================

        [Test]
        public static void Align_RequiresAtLeastTwoExistingInputFiles()
        {
            string dir = FreshDir();
            string a = WriteFlashQuant(dir, "runA", (idx: "1", mass: 10000.0, rtSec: 1000, quant: 100, zMin: 0, zMax: 0));
            string outp = Path.Combine(dir, "consensus.tsv");

            // One real file (+ one that does not exist) is not enough to form a consensus.
            Assert.That(() => FeatureAligner.Align(new[] { a, Path.Combine(dir, "missing.tsv") }, outp),
                Throws.ArgumentException, "a consensus needs at least two existing runs");
        }

        // ================================================================================================
        //  helpers
        // ================================================================================================

        /// <summary>A fresh, empty scratch directory unique to the calling test.</summary>
        private static string FreshDir()
        {
            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "FeatureAlignerUnitTests", TestContext.CurrentContext.Test.Name);
            if (Directory.Exists(dir)) Directory.Delete(dir, true);
            Directory.CreateDirectory(dir);
            return dir;
        }

        /// <summary>Write a FLASHQuant-style feature table (RT already in SECONDS). Returns its path.</summary>
        private static string WriteFlashQuant(string dir, string label,
            params (string idx, double mass, double rtSec, double quant, int zMin, int zMax)[] rows)
        {
            string path = Path.Combine(dir, label + ".tsv");
            var lines = new List<string>
            {
                string.Join("\t", "FeatureGroupIndex", "MonoisotopicMass", "ApexRetentionTime", "FeatureGroupQuantity", "MinCharge", "MaxCharge")
            };
            foreach (var r in rows)
                lines.Add(string.Join("\t", r.idx, r.mass.ToString("R", Ci), r.rtSec.ToString("R", Ci),
                    r.quant.ToString("R", Ci), r.zMin.ToString(Ci), r.zMax.ToString(Ci)));
            File.WriteAllLines(path, lines);
            return path;
        }

        /// <summary>Write a TopFD-style feature table (Apex_time in MINUTES). Returns its path.</summary>
        private static string WriteTopFd(string dir, string label,
            params (string idx, double mass, double rtMin, double apexInt, int zMin, int zMax)[] rows)
        {
            string path = Path.Combine(dir, label + ".tsv");
            var lines = new List<string>
            {
                string.Join("\t", "Feature_ID", "Mass", "Apex_time", "Apex_intensity", "Min_charge", "Max_charge")
            };
            foreach (var r in rows)
                lines.Add(string.Join("\t", r.idx, r.mass.ToString("R", Ci), r.rtMin.ToString("R", Ci),
                    r.apexInt.ToString("R", Ci), r.zMin.ToString(Ci), r.zMax.ToString(Ci)));
            File.WriteAllLines(path, lines);
            return path;
        }

        /// <summary>Read a tab file into a list of column-name -> value dictionaries.</summary>
        private static List<Dictionary<string, string>> ReadTsv(string path)
        {
            var lines = File.ReadAllLines(path);
            var header = lines[0].Split('\t');
            var rows = new List<Dictionary<string, string>>();
            for (int i = 1; i < lines.Length; i++)
            {
                if (lines[i].Length == 0) continue;
                var t = lines[i].Split('\t');
                var d = new Dictionary<string, string>();
                for (int c = 0; c < header.Length; c++) d[header[c]] = c < t.Length ? t[c] : "";
                rows.Add(d);
            }
            return rows;
        }

        private static List<Dictionary<string, string>> ReadConsensus(string path) => ReadTsv(path);

        private static int Int(Dictionary<string, string> row, string col) => int.Parse(row[col], Ci);
        private static double Dbl(Dictionary<string, string> row, string col) => double.Parse(row[col], Ci);
        private static bool row_blank(Dictionary<string, string> row, string col) => string.IsNullOrEmpty(row[col]);

        /// <summary>The single consensus row whose mass_mean is closest to <paramref name="mass"/>.</summary>
        private static Dictionary<string, string> RowNear(List<Dictionary<string, string>> rows, double mass) =>
            rows.OrderBy(r => System.Math.Abs(Dbl(r, "mass_mean") - mass)).First();

        /// <summary>Align two same-mass features whose RTs differ by <paramref name="rtGapSeconds"/>; assert merged or not.</summary>
        private static void AssertPairAlignsOrNot(double rtGapSeconds, double rtTol, bool expectAligned, string because)
        {
            string dir = Path.Combine(FreshDir(), rtGapSeconds.ToString(Ci));
            Directory.CreateDirectory(dir);
            string a = WriteFlashQuant(dir, "runA", (idx: "1", mass: 10000.0, rtSec: 1000.0, quant: 100, zMin: 0, zMax: 0));
            string b = WriteFlashQuant(dir, "runB", (idx: "1", mass: 10000.0, rtSec: 1000.0 + rtGapSeconds, quant: 100, zMin: 0, zMax: 0));
            string outp = Path.Combine(dir, "consensus.tsv");
            FeatureAligner.Align(new[] { a, b }, outp, massTolerance: 20, rtTolerance: rtTol);

            var rows = ReadConsensus(outp);
            if (expectAligned)
                Assert.That(rows.Single()["n_files"], Is.EqualTo("2"), because);
            else
                Assert.That(rows, Has.Count.EqualTo(2), because);
        }
    }
}
