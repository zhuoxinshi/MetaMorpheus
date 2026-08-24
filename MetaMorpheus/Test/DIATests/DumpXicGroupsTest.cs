using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Dumps the EXACT contents of each precursor-fragment group, so downstream analysis never has to
    /// guess which feature a fragment came from.
    ///
    /// The msalign records only fragment MASSES. Because the same fragment mass legitimately appears in
    /// several ISD voltage channels (measured: ~61% of fragment masses match a feature in more than one
    /// channel), reconstructing group membership by mass lookup picks the wrong trace most of the time.
    /// This writes the source voltage, feature id, correlation and overlap for every fragment as the
    /// grouping actually assigned them.
    ///
    /// Two outputs:
    ///   *_group_membership.tsv  one row per (group, fragment): provenance + correlation + overlap
    ///   *_group_xics.tsv        XIC traces (scan, rt, intensity) for a sample of groups, for plotting
    /// </summary>
    [TestFixture]
    public static class DumpXicGroupsTest
    {
        [Test]
        [Explicit("manual dump over local TopFD output")]
        public static void DumpGroups()
        {
            const double precursorMinMass = 4000, precursorMaxMass = 70000, fragmentMaxMass = 70000;
            const int minNumberOfPeaks = 3;
            const float apexRtTolerance = 0.1f;
            const double overlapThreshold = 0, correlationThreshold = 0.75;
            const int xicSampleGroups = 300;   // groups to dump full traces for

            string parentDir = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD";
            string topfdDir = Path.Combine(parentDir, @"SeparateISD\TopFD");
            string outDir = Path.Combine(parentDir, @"SeparateISD\GroupDumps");
            string stem = "07-26-26_HCD_2-3_81min_ISD60-80-100_preFilter800-1000-1200_rep1";
            Directory.CreateDirectory(outDir);

            var dataFile = MsDataFileReader.GetDataFile(Path.Combine(parentDir, @"TopFD_labelCorrected", stem + "_labelCorrected.mzML"));
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // ---------------- precursors
            var ms1Constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                Path.Combine(topfdDir, stem + "_ms1_feature.xml"),
                minMass: precursorMinMass, maxMass: precursorMaxMass, aggregateCharges: true);
            var precursorXics = ms1Constructor.GetAllXics(ms1Scans);
            var precFeature = new Dictionary<ExtractedIonChromatogram, TopFdFeature>(ms1Constructor.LastFeatureByXic);

            // ---------------- fragments, remembering VOLTAGE + FEATURE ID per XIC
            var fragmentXics = new List<ExtractedIonChromatogram>();
            var fragSource = new Dictionary<ExtractedIonChromatogram, (int Voltage, int FeatureId, double Mass)>();
            foreach (var channel in isdMap.OrderBy(k => k.Key))
            {
                var c = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks,
                    Path.Combine(topfdDir, $"{stem}_isd{channel.Key}_feature.xml"),
                    maxMass: fragmentMaxMass, aggregateCharges: true);
                var xics = c.GetAllXics(channel.Value.ToArray());
                foreach (var x in xics)
                    if (c.LastFeatureByXic.TryGetValue(x, out var f))
                        fragSource[x] = (channel.Key, f.Id, f.Mass);
                fragmentXics.AddRange(xics);
            }

            // ---------------- group
            var grouper = new XicGroupingEngine(apexRtTolerance, overlapThreshold, correlationThreshold,
                maxThreadsForGrouping: Math.Max(1, Environment.ProcessorCount - 2));
            var pfGroups = grouper.PrecursorFragmentGrouping(precursorXics, fragmentXics).ToList();
            for (int i = 0; i < pfGroups.Count; i++) pfGroups[i].PFgroupIndex = i + 1;
            PrecursorFragmentsGroup.AssignPrecursorFeatureIds(pfGroups, precFeature);
            TestContext.Progress.WriteLine($"{pfGroups.Count:N0} groups");

            var ci = CultureInfo.InvariantCulture;

            // ---------------- membership: every fragment of every group
            string memPath = Path.Combine(outDir, stem + "_group_membership.tsv");
            using (var w = new StreamWriter(memPath))
            {
                w.WriteLine("scan\tprecursor_feature_id\tprecursor_mass\tprecursor_apex_rt\t" +
                            "frag_voltage\tfrag_feature_id\tfrag_mass\tfrag_apex_rt\tapex_rt_diff\tcorrelation\toverlap\tn_frag_points");
                foreach (var g in pfGroups)
                {
                    double pApex = g.PrecursorXic.ApexRT;
                    double pMass = g.PrecursorXic.ApexPeak.M;
                    foreach (var pf in g.PFpairs)
                    {
                        var src = fragSource.TryGetValue(pf.FragmentXic, out var s) ? s : (-1, -1, double.NaN);
                        w.WriteLine(string.Join("\t",
                            g.PFgroupIndex,
                            g.PrecursorFeatureId?.ToString() ?? "",
                            pMass.ToString("F4", ci),
                            pApex.ToString("F4", ci),
                            src.Item1, src.Item2,
                            pf.FragmentXic.ApexPeak.M.ToString("F4", ci),
                            pf.FragmentXic.ApexRT.ToString("F4", ci),
                            Math.Abs(pf.FragmentXic.ApexRT - pApex).ToString("F4", ci),
                            pf.Correlation?.ToString("F4", ci) ?? "",
                            pf.Overlap?.ToString("F4", ci) ?? "",
                            pf.FragmentXic.Peaks.Count));
                    }
                }
            }
            TestContext.Progress.WriteLine("wrote " + memPath);

            // ---------------- XIC traces for a sample of groups
            var rng = new Random(11);
            var sample = pfGroups.OrderBy(_ => rng.Next()).Take(xicSampleGroups).ToList();
            string xicPath = Path.Combine(outDir, stem + "_group_xics.tsv");
            using (var w = new StreamWriter(xicPath))
            {
                w.WriteLine("scan\trole\tfrag_voltage\tfeature_id\tmass\tcorrelation\tpoint_rt\tpoint_intensity");
                foreach (var g in sample)
                {
                    foreach (var pk in g.PrecursorXic.Peaks)
                        w.WriteLine(string.Join("\t", g.PFgroupIndex, "precursor", "", g.PrecursorFeatureId,
                            g.PrecursorXic.ApexPeak.M.ToString("F4", ci), "",
                            pk.RetentionTime.ToString("F4", ci), pk.Intensity.ToString("F1", ci)));
                    foreach (var pf in g.PFpairs)
                    {
                        var src = fragSource.TryGetValue(pf.FragmentXic, out var s) ? s : (-1, -1, double.NaN);
                        foreach (var pk in pf.FragmentXic.Peaks)
                            w.WriteLine(string.Join("\t", g.PFgroupIndex, "fragment", src.Item1, src.Item2,
                                pf.FragmentXic.ApexPeak.M.ToString("F4", ci),
                                pf.Correlation?.ToString("F4", ci) ?? "",
                                pk.RetentionTime.ToString("F4", ci), pk.Intensity.ToString("F1", ci)));
                    }
                }
            }
            TestContext.Progress.WriteLine($"wrote {xicPath} ({sample.Count} groups)");
        }
    }
}
