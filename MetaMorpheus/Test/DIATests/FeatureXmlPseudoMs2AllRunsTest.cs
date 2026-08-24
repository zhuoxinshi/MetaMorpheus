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

namespace Test.DIATests
{
    /// <summary>
    /// Batch pseudo-MS2 generation over ALL 10 ISD runs (2 pure HCD, 2 pure YBC, 6 HCD-YBC mixtures),
    /// using the per-voltage TopFD runs in SeparateISD\TopFD.
    ///
    /// Differences from the OneFileTest generators:
    ///   - the PRECURSOR feature xml is the split MS1-channel TopFD run (<stem>_ms1_feature.xml), not a
    ///     whole-parent run. The split preserves parent scan numbers (all envelopes are ≡1 mod 4), so the
    ///     scan-number join against the parent's MS1 scans holds - verified before this was written.
    ///   - grouping is the STRICT operating point, XicGroupingEngine(0.1, 0, 0.75), per instruction.
    ///
    /// The <see cref="FeatureChargeMode"/> knob chooses how a feature's charge states become one XIC: summing them
    /// (the original behaviour) or taking only the most abundant charge as the feature's representative. Each mode
    /// writes to its own directory, so the two can be generated and compared without overwriting each other.
    ///
    /// Output: SeparateISD\PseudoMs2[_maxZ]\<stem>_pseudo_ms2.msalign, written with the fixed writer
    /// (real XicSum intensities, PRECURSOR_FEATURE_ID from the ms1-channel TopFD feat_id).
    /// </summary>
    [TestFixture]
    public static class FeatureXmlPseudoMs2AllRunsTest
    {
        [Test]
        [Explicit("Long batch over local data; run by name.")]
        public static void WritePseudoMs2ForAllIsdRuns()
        {
            // ------------------------------------------------------------------ knobs (match OneFileTest)
            const double precursorMinMass = 4000;
            const double precursorMaxMass = 70000;
            const double fragmentMaxMass = 70000;
            const double minEcScore = 0;
            const int minNumberOfPeaks = 3;

            const float apexRtTolerance = 0.1f;
            const double overlapThreshold = 0;
            const double correlationThreshold = 0.75;

            // How each feature's charge states collapse into one XIC, applied to precursors and fragments alike.
            //   SumAllCharges      - the original behaviour: intensities summed across charges within each scan
            //   MostAbundantCharge - the feature is represented by its strongest charge state only, so the XIC is
            //                        one charge's measured elution instead of a blend of several
            const FeatureChargeMode chargeMode = FeatureChargeMode.MostAbundantCharge;

            // ------------------------------------------------------------------ files
            string parentDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string topfdDir = Path.Combine(parentDir, @"SeparateISD\TopFD");
            // each mode gets its own directory; the unsuffixed one stays the summed output already on disk. The
            // tag is kept short because TopPIC rejects a spectrum path over 200 characters and these stems are long
            string modeTag = chargeMode == FeatureChargeMode.SumAllCharges ? ""
                : chargeMode == FeatureChargeMode.MostAbundantCharge ? "_maxZ" : "_sepZ";
            string outDir = Path.Combine(parentDir, @"SeparateISD\PseudoMs2" + modeTag);
            var voltages = new[] { 60, 80, 100 };

            // run stems: derived from the ms1 feature xmls so the list always matches what TopFD produced
            var stems = Directory.GetFiles(topfdDir, "*_ms1_feature.xml")
                .Select(f => Path.GetFileName(f).Replace("_ms1_feature.xml", ""))
                .OrderBy(s => s)
                .ToList();
            TestContext.Progress.WriteLine($"{stems.Count} runs found; charge mode {chargeMode} -> {outDir}");
            Assert.That(stems.Count, Is.EqualTo(10));

            Directory.CreateDirectory(outDir);
            foreach (var stem in stems)
            {
                // the parent mzML keeps the _labelCorrected suffix the split files drop
                string parentMzml = Path.Combine(parentDir, stem + "_labelCorrected.mzML");
                Assert.That(File.Exists(parentMzml), $"missing {parentMzml}");

                var dataFile = MsDataFileReader.GetDataFile(parentMzml);
                dataFile.LoadAllStaticData();
                var allScans = dataFile.GetAllScansList().ToArray();
                var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
                ISDEngine.ReLabelIsdScans(isdMap, allScans);

                // ---------------------------------------------------------- precursor XICs from the ms1 split
                var ms1Constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                    Path.Combine(topfdDir, stem + "_ms1_feature.xml"),
                    minMass: precursorMinMass, maxMass: precursorMaxMass,
                    minEcScore: minEcScore, chargeMode: chargeMode);
                var precursorXics = ms1Constructor.GetAllXics(ms1Scans);
                // how often the charge picked on the MS1 channel is the one TopFD called representative over the
                // whole run - a sanity check that the pick is not being driven by the channel subsetting
                int ms1ChargesPicked = ms1Constructor.LastMostAbundantChargeByFeature.Count;
                int ms1ChargesMatchingTopFd = ms1Constructor.LastMostAbundantChargeMatchingTopFdCount;
                // snapshot before the fragment constructors run - the map is per-call
                var featureByPrecursorXic = new Dictionary<ExtractedIonChromatogram, TopFdFeature>(ms1Constructor.LastFeatureByXic);

                // ---------------------------------------------------------- fragment XICs, one voltage at a time
                var fragmentXics = new List<ExtractedIonChromatogram>();
                foreach (var channel in isdMap.OrderBy(k => k.Key))
                {
                    Assert.That(voltages.Contains(channel.Key), $"unexpected voltage {channel.Key} in {stem}");
                    fragmentXics.AddRange(new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                            Path.Combine(topfdDir, stem + "_isd" + channel.Key + "_feature.xml"),
                            maxMass: fragmentMaxMass, minEcScore: minEcScore, chargeMode: chargeMode)
                        .GetAllXics(channel.Value.ToArray()));
                }

                // ---------------------------------------------------------- grouping -> pseudo scans -> msalign
                var grouper = new XicGroupingEngine(apexRtTolerance, overlapThreshold, correlationThreshold,
                    maxThreadsForGrouping: Math.Max(1, Environment.ProcessorCount - 2));
                var pfGroups = grouper.PrecursorFragmentGrouping(precursorXics, fragmentXics).ToList();
                PrecursorFragmentsGroup.AssignPrecursorFeatureIds(pfGroups, featureByPrecursorXic);

                var commonParams = new CommonParameters();
                var pseudoScans = new List<Ms2ScanWithSpecificMass>(pfGroups.Count);
                for (int i = 0; i < pfGroups.Count; i++)
                {
                    pfGroups[i].PFgroupIndex = i + 1;
                    pseudoScans.Add(pfGroups[i].GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.Mass, commonParams, dataFile.FilePath));
                }

                string outPath = Path.Combine(outDir, stem + "_pseudo_ms2.msalign");
                ISDEngine.WriteMsAlignFile(outPath, pseudoScans);

                var counts = pseudoScans.Select(s => s.ExperimentalFragments.Length).OrderBy(n => n).ToArray();
                TestContext.Progress.WriteLine(
                    $"{stem}: precursors {precursorXics.Count:N0}, fragments {fragmentXics.Count:N0} -> " +
                    $"{pseudoScans.Count:N0} pseudo scans, fragments/scan median {(counts.Length > 0 ? counts[counts.Length / 2] : 0)}, " +
                    $"total peaks {counts.Sum():N0}");
                if (chargeMode == FeatureChargeMode.MostAbundantCharge && ms1ChargesPicked > 0)
                    TestContext.Progress.WriteLine(
                        $"   representative charge chosen for {ms1ChargesPicked:N0} precursor features, " +
                        $"{ms1ChargesMatchingTopFd:N0} ({100.0 * ms1ChargesMatchingTopFd / ms1ChargesPicked:F1}%) " +
                        "agreeing with TopFD rep_charge");
                TestContext.Progress.WriteLine($"   wrote {outPath}");
            }
        }
    }
}
