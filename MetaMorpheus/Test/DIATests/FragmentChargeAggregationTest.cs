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
using System.Text.RegularExpressions;

namespace Test.DIATests
{
    /// <summary>
    /// Generates two pseudo-MS2 msAlign files that differ in ONE setting - whether fragment charge states are
    /// aggregated onto their feature - so a TopPIC search of each isolates the effect of that choice.
    ///
    /// aggregateCharges TRUE  : one fragment XIC per feature, envelope intensities summed across charge states.
    ///                          One peak line per fragment feature.
    /// aggregateCharges FALSE : one fragment XIC per charge state. The same fragment mass can appear several
    ///                          times in a spectrum, once per charge it was observed at.
    ///
    /// Everything else is shared: the SAME precursor XICs, the same grouping engine and thresholds, the same
    /// intensity convention. Only the fragment XIC list is rebuilt, so any difference in identifications is
    /// attributable to the aggregation choice alone.
    /// </summary>
    [TestFixture]
    public static class FragmentChargeAggregationTest
    {
        [Test]
        [Explicit("Writes msAlign files over local TopFD output; run by name, then search both with TopPIC.")]
        public static void WriteAggregatedAndUnaggregatedMsAlign()
        {
            // ------------------------------------------------------------------ knobs (the loose operating point)
            const double precursorMinMass = 4000;
            const double precursorMaxMass = 70000;
            const double fragmentMaxMass = 70000;
            const double minEcScore = 0;
            const int minNumberOfPeaks = 3;

            const float apexRtTolerance = 0.2f;
            const double overlapThreshold = 0;
            const double correlationThreshold = 0.50;

            // ------------------------------------------------------------------ files
            string toppicDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string stem = @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected";
            string parentMzml = Path.Combine(toppicDir, stem + ".mzML");
            string precursorFeatureXml = Path.Combine(toppicDir, stem + "_feature.xml");
            string channelDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\HCD_YBC_ISDseparate";
            string outDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\ChargeAggregation";

            var fragmentFeatureXmlByVoltage = Directory.GetFiles(channelDir, "*_feature.xml")
                .Select(f => (File: f, Match: Regex.Match(Path.GetFileName(f), @"_isd(\d+)_feature\.xml$", RegexOptions.IgnoreCase)))
                .Where(t => t.Match.Success)
                .ToDictionary(t => int.Parse(t.Match.Groups[1].Value), t => t.File);

            // ------------------------------------------------------------------ scans
            var dataFile = MsDataFileReader.GetDataFile(parentMzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // ------------------------------------------------------------------ precursor XICs, built ONCE and shared
            var ms1Constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                precursorFeatureXml, minMass: precursorMinMass, maxMass: precursorMaxMass,
                minEcScore: minEcScore, aggregateCharges: true);
            var precursorXics = ms1Constructor.GetAllXics(ms1Scans);
            // snapshot before the fragment constructors run - the map is per-call
            var featureByPrecursorXic = new Dictionary<ExtractedIonChromatogram, TopFdFeature>(ms1Constructor.LastFeatureByXic);
            TestContext.Progress.WriteLine($"precursors: {precursorXics.Count:N0} XICs (shared by both arms)");

            Directory.CreateDirectory(outDir);
            foreach (bool aggregate in new[] { true, false })
            {
                // ---------------------------------------------------------- fragment XICs, the ONLY thing that varies
                var fragmentXics = new List<ExtractedIonChromatogram>();
                foreach (var channel in isdMap.OrderBy(k => k.Key))
                    fragmentXics.AddRange(new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                            fragmentFeatureXmlByVoltage[channel.Key], maxMass: fragmentMaxMass,
                            minEcScore: minEcScore, aggregateCharges: aggregate)
                        .GetAllXics(channel.Value.ToArray()));

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

                string label = aggregate ? "aggregated" : "perCharge";
                string outPath = Path.Combine(outDir, label + "_ms2.msalign");
                ISDEngine.WriteMsAlignFile(outPath, pseudoScans);

                var counts = pseudoScans.Select(s => s.ExperimentalFragments.Length).OrderBy(n => n).ToArray();
                // how much of the extra is the same mass repeated at another charge?
                double dupFraction = pseudoScans.Average(s =>
                {
                    int n = s.ExperimentalFragments.Length;
                    if (n == 0) return 0.0;
                    int distinct = s.ExperimentalFragments.Select(e => Math.Round(e.MonoisotopicMass, 3)).Distinct().Count();
                    return 1.0 - (double)distinct / n;
                });

                TestContext.Progress.WriteLine(
                    $"{label,-11}: fragment XICs {fragmentXics.Count,7:N0} | pseudo scans {pseudoScans.Count,6:N0} | " +
                    $"fragments/scan mean {counts.Average(),6:F1} median {counts[counts.Length / 2],4} | " +
                    $"total peaks {counts.Sum(),9:N0} | repeated-mass fraction {dupFraction:P1}");
                TestContext.Progress.WriteLine($"   wrote {outPath}");
            }
        }
    }
}
