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
    /// Throwaway extractor (not CI): runs the ISD consensus tracer on a real raw file, groups
    /// precursor/fragment XICs (fragments with AggregateCharges = false), picks one strong
    /// co-eluting precursor+fragment pair, and dumps their RAW traced peaks to CSV for plotting.
    /// </summary>
    [TestFixture]
    public static class DumpRealXicTest
    {
        [Test, Explicit]
        public static void DumpRealIsdXicPair()
        {
            string mzml = @"E:\ISD Project\ISD_250428\05-04-25_PEPPI-YB_81min_ISD60-80-100_preFilter700-900-1100_rep1.raw";
            string outDir = @"C:\Users\Zhuoxin Shi\AppData\Local\Temp\claude\E--Islets-Brian-data-PTM\b851b496-68df-489a-95bd-116185bab683\scratchpad\realxic";
            Directory.CreateDirectory(outDir);

            var decon = new ClassicDeconvolutionParameters(1, 30, 4, 3);
            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // precursor channel: full consensus features (>=3 kDa, >=3 charges); fragment channel: mass-tracing only
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
                minMass: 3000, minChargeCount: 3, traceToleranceDa: 1.5);
            var ms2Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon, traceToleranceDa: 1.5)
            { AggregateCharges = false };

            var ms1Xics = ms1Ctor.GetAllXics(ms1Scans);
            var ms2Xics = isdMap.Values.SelectMany(v => ms2Ctor.GetAllXics(v.ToArray())).ToList();
            TestContext.Progress.WriteLine($"consensus XICs: {ms1Xics.Count} precursor, {ms2Xics.Count} fragment");

            var grouper = new XicGroupingEngine(0.2f, 0.3, 0.7, Environment.ProcessorCount);
            var groups = grouper.PrecursorFragmentGrouping(ms1Xics, ms2Xics).ToList();
            TestContext.Progress.WriteLine($"groups: {groups.Count}");

            // pick a pair where BOTH the precursor and fragment have plenty of points (>= minPts),
            // so Akima is valid on both and the spline/padding effect is meaningful.
            int minPts = int.TryParse(Environment.GetEnvironmentVariable("MIN_PTS"), out var mp) ? mp : 8;
            var pairs = groups
                .SelectMany(x => x.PFpairs.Select(p => (grp: x, pair: p)))
                .Where(t => t.grp.PrecursorXic.Peaks.Count >= minPts && t.pair.FragmentXic.Peaks.Count >= minPts)
                .ToList();
            TestContext.Progress.WriteLine($"pairs with both >= {minPts} pts: {pairs.Count}");

            var pick = pairs
                .Where(t => (t.pair.Correlation ?? 0) >= 0.9)
                .OrderByDescending(t => Math.Min(t.grp.PrecursorXic.Peaks.Count, t.pair.FragmentXic.Peaks.Count))
                .ThenByDescending(t => t.pair.Correlation ?? 0)
                .DefaultIfEmpty(pairs.OrderByDescending(t => t.pair.Correlation ?? 0).First())
                .First();

            var g = pick.grp;
            var prec = pick.grp.PrecursorXic;
            var bestPair = pick.pair;
            var frag = pick.pair.FragmentXic;

            Dump(Path.Combine(outDir, "precursor.csv"), prec);
            Dump(Path.Combine(outDir, "fragment.csv"), frag);

            using (var meta = new StreamWriter(Path.Combine(outDir, "meta.txt")))
            {
                meta.WriteLine($"file={Path.GetFileName(mzml)}");
                meta.WriteLine($"groups={groups.Count}");
                meta.WriteLine($"precursor mass={prec.ApexPeak.M:F3} apexRT={prec.ApexRT:F3} pts={prec.Peaks.Count} scanIdx[{prec.StartScanIndex}..{prec.EndScanIndex}]");
                meta.WriteLine($"fragment  mass={frag.ApexPeak.M:F3} apexRT={frag.ApexRT:F3} pts={frag.Peaks.Count} scanIdx[{frag.StartScanIndex}..{frag.EndScanIndex}]");
                meta.WriteLine($"pair correlation={(bestPair.Correlation ?? 0):F4} overlap={(bestPair.Overlap ?? 0):F4} fragmentsInGroup={g.PFpairs.Count}");
            }
            TestContext.Progress.WriteLine("dumped precursor/fragment CSV to " + outDir);
        }

        private static void Dump(string path, ExtractedIonChromatogram xic)
        {
            using var w = new StreamWriter(path);
            w.WriteLine("scanIndex,rt,intensity,mass");
            foreach (var p in xic.Peaks.OrderBy(p => p.ZeroBasedScanIndex))
                w.WriteLine(string.Join(",",
                    p.ZeroBasedScanIndex,
                    p.RetentionTime.ToString(CultureInfo.InvariantCulture),
                    p.Intensity.ToString(CultureInfo.InvariantCulture),
                    p.M.ToString(CultureInfo.InvariantCulture)));
        }
    }
}
