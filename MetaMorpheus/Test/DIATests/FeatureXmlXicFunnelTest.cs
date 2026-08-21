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
    /// Diagnostic, not an assertion suite: accounts for the gap between how many features TopFD wrote into a
    /// <c>*_feature.xml</c> and how many XICs <see cref="FeatureXmlXicConstructor"/> actually emits from it.
    ///
    /// The constructor drops features at four separate places and only reports the total, so a large shortfall is
    /// easy to misread as a bug. This walks the funnel stage by stage, for the precursor file and for each ISD
    /// fragment channel, using the same settings the pseudo-MS2 generators use.
    /// </summary>
    [TestFixture]
    public static class FeatureXmlXicFunnelTest
    {
        [Test]
        [Explicit("Diagnostic over local TopFD output; run by name.")]
        public static void WhereDoTheFeaturesGo()
        {
            // same knobs as PseudoMs2FromFeatureXmlThroughIsdEngine
            const double precursorMinMass = 4000;
            const double precursorMaxMass = 70000;
            const double fragmentMaxMass = 70000;
            const double minEcScore = 0;
            const int minNumberOfPeaks = 3;

            string toppicDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string stem = @"07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1_labelCorrected";
            string parentMzml = Path.Combine(toppicDir, stem + ".mzML");
            string precursorFeatureXml = Path.Combine(toppicDir, stem + "_feature.xml");
            string channelDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\OneFileTest\HCD_YBC_ISDseparate";

            var fragmentFeatureXmlByVoltage = Directory.GetFiles(channelDir, "*_feature.xml")
                .Select(f => (File: f, Match: Regex.Match(Path.GetFileName(f), @"_isd(\d+)_feature\.xml$", RegexOptions.IgnoreCase)))
                .Where(t => t.Match.Success)
                .ToDictionary(t => int.Parse(t.Match.Groups[1].Value), t => t.File);

            var dataFile = MsDataFileReader.GetDataFile(parentMzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            TestContext.Progress.WriteLine($"scan array sizes: MS1 {ms1Scans.Length}, " +
                string.Join(", ", isdMap.OrderBy(k => k.Key).Select(k => $"isd{k.Key} {k.Value.Count}")));
            TestContext.Progress.WriteLine("");

            Report("PRECURSOR (MS1)", precursorFeatureXml, ms1Scans, precursorMinMass, precursorMaxMass, true);
            foreach (var channel in isdMap.OrderBy(k => k.Key))
                Report($"FRAGMENT isd{channel.Key}", fragmentFeatureXmlByVoltage[channel.Key],
                    channel.Value.ToArray(), 0, fragmentMaxMass, true);

            // how much of the precursor loss is the mass floor alone?
            TestContext.Progress.WriteLine("---- precursor with the mass floor removed, for contrast ----");
            Report("PRECURSOR (no mass filter)", precursorFeatureXml, ms1Scans, 0, double.MaxValue, false);

            void Report(string label, string featureXml, MsDataScan[] scans, double minMass, double maxMass, bool detail)
            {
                var constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                    featureXml, minMass: minMass, maxMass: maxMass, minEcScore: minEcScore, aggregateCharges: true);

                int inFile = constructor.FeatureCountBeforeFilters;
                int afterFilters = constructor.FeatureCount;
                var xics = constructor.GetAllXics(scans);

                int noScans = constructor.LastFeaturesWithNoScansInRangeCount;
                int tooShort = constructor.LastFeaturesBelowMinPeaksCount;

                TestContext.Progress.WriteLine($"== {label} ==  {Path.GetFileName(featureXml)}");
                TestContext.Progress.WriteLine($"   features in file          {inFile,8:N0}");
                TestContext.Progress.WriteLine($"   after mass/charge/EC      {afterFilters,8:N0}   (-{inFile - afterFilters:N0})");
                TestContext.Progress.WriteLine($"   no envelope in this array {noScans,8:N0}");
                TestContext.Progress.WriteLine($"   fewer than {minNumberOfPeaks} points      {tooShort,8:N0}");
                TestContext.Progress.WriteLine($"   XICs emitted              {xics.Count,8:N0}   ({100.0 * xics.Count / Math.Max(1, inFile):F1}% of file)");

                if (detail && xics.Count > 0)
                {
                    var lens = xics.Select(x => x.Peaks.Count).OrderBy(n => n).ToArray();
                    TestContext.Progress.WriteLine($"   XIC length: median {lens[lens.Length / 2]}, max {lens[^1]}");
                }
                TestContext.Progress.WriteLine("");
            }
        }
    }
}
