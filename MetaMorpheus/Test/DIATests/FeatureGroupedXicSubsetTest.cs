using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Runs the MS1 path - trace XICs from an ms1.msalign, group them onto ms1.feature entries - over a SMALL
    /// REAL SUBSET so the whole funnel can be inspected in a second instead of a minute.
    ///
    /// The fixture is 38 consecutive MS1 scans (RT 30.05-31.95 min) cut out of a TopFD 1.8.1 ISD run, plus the
    /// 1139 features whose elution windows overlap that slice. No raw file is needed: the scan numbers and
    /// retention times come from the msAlign itself, and the constructor only ever reads those two fields off
    /// each <see cref="MsDataScan"/> - the spectra are placeholders because every mass comes from the msAlign.
    /// </summary>
    [TestFixture]
    public static class FeatureGroupedXicSubsetTest
    {
        private static string DataDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "DIATests", "ExternalDeconTestData");

        private static string MsAlignPath => Path.Combine(DataDirectory, "HcdYbc_subset_ms1.msalign");
        private static string FeaturePath => Path.Combine(DataDirectory, "HcdYbc_subset_ms1.feature");

        /// <summary>
        /// Rebuilds the scan array the constructor would normally get from the raw file, using the scan numbers
        /// and retention times recorded in the msAlign.
        /// </summary>
        private static MsDataScan[] SynthesiseScans(ExternalDeconResults results)
        {
            return results.PeaksByOneBasedScanNumber.Keys
                .OrderBy(scanNumber => scanNumber)
                .Select(scanNumber => new MsDataScan(
                    new MzSpectrum(new double[] { 700 }, new double[] { 1 }, false),
                    oneBasedScanNumber: scanNumber, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                    retentionTime: results.RetentionTimeByOneBasedScanNumber[scanNumber],
                    scanWindowRange: new MzRange(200, 2000), scanFilter: "FTMS + p NSI sid=15.00 Full ms",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1, injectionTime: 1.0, noiseData: null,
                    nativeId: "scan=" + scanNumber))
                .ToArray();
        }

        /// <summary>
        /// Prints the whole funnel - masses in, masses indexed, XICs traced, XICs assigned, features hit - with
        /// the reasons XICs failed to land on a feature. This is the test to run when the feature-XIC count looks
        /// too low: it says whether the loss is in tracing or in assignment.
        /// </summary>
        [Test]
        public static void FeatureGroupingFunnelOnRealSubset()
        {
            Assert.That(File.Exists(MsAlignPath), $"missing fixture {MsAlignPath}");

            var results = ExternalDeconResults.ReadMsAlign(MsAlignPath);
            var scans = SynthesiseScans(results);
            var allFeatures = Ms1FeatureRecord.ReadFile(FeaturePath);

            TestContext.Progress.WriteLine($"fixture: {results.ScanCount} scans, {results.PeakCount} masses, " +
                $"RT {scans.First().RetentionTime:F2}-{scans.Last().RetentionTime:F2} min " +
                $"({(scans.Last().RetentionTime - scans.First().RetentionTime) / (scans.Length - 1) * 60:F1} s/scan); " +
                $"{allFeatures.Count} features");

            var constructor = new FeatureGroupedXicConstructor(new PpmTolerance(20), maxMissedScansAllowed: 2,
                maxPeakHalfWidth: 1, minNumberOfPeaks: 3, results, allFeatures,
                featureMassTolerance: new PpmTolerance(20), minMass: 3000, maxMass: 50000);

            TestContext.Progress.WriteLine($"features after [{constructor.MinMass}, {constructor.MaxMass}) mass filter: " +
                $"{constructor.FeatureCount} of {constructor.FeatureCountBeforeMassFilter}");

            var featureXics = constructor.GetAllXics(scans);

            TestContext.Progress.WriteLine(
                $"  indexed masses      : {constructor.LastIndexedPeakCount} " +
                $"(of {results.PeakCount}; {constructor.LastDuplicatePeakCount} dropped as duplicates)\n" +
                $"  traced XICs         : {constructor.LastTracedXicCount}\n" +
                $"  assigned to features: {constructor.LastAssignedXicCount}\n" +
                $"    unassigned, no feature within mass tolerance : {constructor.LastUnassignedNoMassMatchCount}\n" +
                $"    unassigned, mass matched but RT/charge failed: {constructor.LastUnassignedRejectedOnRtOrChargeCount}\n" +
                $"  feature XICs        : {featureXics.Count} " +
                $"({constructor.LastFeaturesWithXicsCount}/{constructor.FeatureCount} features hit)\n" +
                $"  features too short to ever reach {constructor.MinNumberOfPeaks} points on this scan axis: " +
                $"{constructor.LastFeaturesTooShortToTraceCount}");

            if (featureXics.Count > 0)
            {
                var lengths = featureXics.Select(x => x.Peaks.Count).OrderBy(n => n).ToArray();
                TestContext.Progress.WriteLine($"  feature XIC length: min {lengths.First()}, " +
                    $"median {lengths[lengths.Length / 2]}, max {lengths.Last()}");
                var chargesPerFeature = constructor.LastFeatureByXic
                    .Where(kvp => kvp.Value != null)
                    .Select(kvp => (double)constructor.LastAssignedXicCount / constructor.LastFeaturesWithXicsCount)
                    .FirstOrDefault();
                TestContext.Progress.WriteLine($"  mean traces merged per feature: {chargesPerFeature:F1}");
            }

            Assert.That(constructor.LastMatchedScanCount, Is.EqualTo(scans.Length), "every scan should match itself");
            Assert.That(constructor.FeatureCount, Is.LessThanOrEqualTo(constructor.FeatureCountBeforeMassFilter));
            Assert.That(featureXics.Count, Is.GreaterThan(0));
        }

        /// <summary>
        /// The mass filter on the feature list is what makes "features hit" mean anything: without it the
        /// denominator includes features the traced XICs could never reach, because those XICs are themselves
        /// filtered to [minMass, maxMass) before indexing.
        /// </summary>
        [Test]
        public static void FeaturesAreFilteredToTheSameMassRangeAsTheXics()
        {
            var results = ExternalDeconResults.ReadMsAlign(MsAlignPath);
            var allFeatures = Ms1FeatureRecord.ReadFile(FeaturePath);

            var constructor = new FeatureGroupedXicConstructor(new PpmTolerance(20), 2, 1, 3, results, allFeatures,
                minMass: 3000, maxMass: 20000);

            Assert.That(constructor.FeatureCountBeforeMassFilter, Is.EqualTo(allFeatures.Count));
            Assert.That(constructor.Features, Is.Not.Empty);
            Assert.That(constructor.Features.Select(f => f.Mass), Has.All.InRange(3000.0, 20000.0));
            Assert.That(constructor.FeatureCount, Is.EqualTo(allFeatures.Count(f => f.Mass >= 3000 && f.Mass < 20000)));

            // and the features outside the window really were present to begin with, i.e. the filter did work
            Assert.That(allFeatures.Any(f => f.Mass < 3000 || f.Mass >= 20000),
                "fixture should contain features outside the mass window for this to be a meaningful check");
            TestContext.Progress.WriteLine($"features {allFeatures.Count} -> {constructor.FeatureCount} " +
                $"after the [3000, 20000) filter");
        }

        /// <summary>
        /// TopFD flips the monoisotope assignment of a species between adjacent scans - the same charge is
        /// reported at M and at M+1.003 Da - so a plain ppm window breaks the trace at every flip. This sweeps a
        /// notched tolerance (which admits +/-1 isotope) against a plain one, and a wider missed-scan allowance,
        /// because each individual charge state is also reported only intermittently.
        /// </summary>
        [Test]
        public static void NotchedToleranceAndMissedScansRecoverTraces()
        {
            var results = ExternalDeconResults.ReadMsAlign(MsAlignPath);
            var scans = SynthesiseScans(results);
            var allFeatures = Ms1FeatureRecord.ReadFile(FeaturePath);

            foreach (var (label, tolerance) in new (string, Tolerance)[]
                     {
                         ("plain 20 ppm     ", new PpmTolerance(20)),
                         ("notched 20 ppm +/-1", new PpmToleranceWithNotch(20, 1, 1)),
                     })
            {
                foreach (int maxMissed in new[] { 2, 4, 6 })
                {
                    var constructor = new FeatureGroupedXicConstructor(tolerance, maxMissed, 1, 3,
                        results, allFeatures, featureMassTolerance: tolerance, minMass: 3000, maxMass: 50000);
                    var xics = constructor.GetAllXics(scans);
                    TestContext.Progress.WriteLine($"{label} maxMissed={maxMissed}: " +
                        $"{constructor.LastTracedXicCount,5} traced -> {constructor.LastAssignedXicCount,5} assigned -> " +
                        $"{xics.Count,4} feature XICs ({constructor.LastFeaturesWithXicsCount}/{constructor.FeatureCount} features)");
                }
            }
        }

        /// <summary>
        /// How the feature-XIC yield responds to minNumberOfPeaks. On a run whose features elute over only a few
        /// MS1 scans this threshold, not the mass or RT matching, is what governs the yield.
        /// </summary>
        [Test]
        public static void MinNumberOfPeaksDrivesTheYield()
        {
            var results = ExternalDeconResults.ReadMsAlign(MsAlignPath);
            var scans = SynthesiseScans(results);
            var allFeatures = Ms1FeatureRecord.ReadFile(FeaturePath);

            foreach (int minPeaks in new[] { 2, 3, 5, 8 })
            {
                var constructor = new FeatureGroupedXicConstructor(new PpmTolerance(20), 2, 1, minPeaks,
                    results, allFeatures, featureMassTolerance: new PpmTolerance(20), minMass: 3000, maxMass: 50000);
                var xics = constructor.GetAllXics(scans);
                TestContext.Progress.WriteLine($"minNumberOfPeaks={minPeaks}: {constructor.LastTracedXicCount} traced -> " +
                    $"{xics.Count} feature XICs ({constructor.LastFeaturesWithXicsCount}/{constructor.FeatureCount} features), " +
                    $"{constructor.LastFeaturesTooShortToTraceCount} features too short to qualify");
            }
        }
    }
}
