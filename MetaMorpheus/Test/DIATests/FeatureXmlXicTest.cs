using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace Test.DIATests
{
    /// <summary>
    /// Tests for reading feature XICs straight out of TopFD's <c>*_feature.xml</c>, which stores a per-charge
    /// <c>envelope_list</c> (scan, retention time, intensity) for every proteoform feature - i.e. the elution
    /// profiles the msAlign path has to reconstruct by tracing.
    ///
    /// The fixture is 230 features (every 4th of those lying entirely within scans 2305-2453 of a TopFD 1.8.1 ISD
    /// run), so it covers exactly the same slice as the msAlign fixture and the two approaches can be compared
    /// head to head on identical scans.
    /// </summary>
    [TestFixture]
    public static class FeatureXmlXicTest
    {
        private static string DataDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "DIATests", "ExternalDeconTestData");

        private static string FeatureXmlPath => Path.Combine(DataDirectory, "HcdYbc_subset_feature.xml");
        private static string MsAlignPath => Path.Combine(DataDirectory, "HcdYbc_subset_ms1.msalign");
        private static string FeatureTsvPath => Path.Combine(DataDirectory, "HcdYbc_subset_ms1.feature");

        /// <summary>Scan array rebuilt from the msAlign fixture, so both paths share one time axis.</summary>
        private static MsDataScan[] SynthesiseScans()
        {
            var results = ExternalDeconResults.ReadMsAlign(MsAlignPath);
            return results.PeaksByOneBasedScanNumber.Keys.OrderBy(s => s)
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
        /// The reader reproduces the file: counts of features, single-charge features and envelopes match a plain
        /// text count of the tags, and the nested values land on the right objects.
        /// </summary>
        [Test]
        public static void ReaderReproducesTheFeatureXml()
        {
            Assert.That(File.Exists(FeatureXmlPath), $"missing fixture {FeatureXmlPath}");
            string text = File.ReadAllText(FeatureXmlPath);
            int expectedFeatures = Regex.Matches(text, "<frac_feature>").Count;
            int expectedChargeFeatures = Regex.Matches(text, "<single_charge_feature>").Count;
            int expectedEnvelopes = Regex.Matches(text, "<envelope>").Count;

            var features = TopFdFeatureXmlFile.ReadFile(FeatureXmlPath);

            TestContext.Progress.WriteLine($"text tags : {expectedFeatures} features, " +
                $"{expectedChargeFeatures} single-charge, {expectedEnvelopes} envelopes");
            TestContext.Progress.WriteLine($"reader    : {features.Count} features, " +
                $"{features.Sum(f => f.SingleChargeFeatures.Count)} single-charge, " +
                $"{features.Sum(f => f.AllEnvelopes.Count())} envelopes");

            Assert.That(features.Count, Is.EqualTo(expectedFeatures));
            Assert.That(features.Sum(f => f.SingleChargeFeatures.Count), Is.EqualTo(expectedChargeFeatures));
            Assert.That(features.Sum(f => f.AllEnvelopes.Count()), Is.EqualTo(expectedEnvelopes));

            // every feature must be internally consistent
            foreach (var feature in features)
            {
                Assert.That(feature.Mass, Is.GreaterThan(0), $"feature {feature.Id} has no mass");
                Assert.That(feature.SingleChargeFeatures, Is.Not.Empty, $"feature {feature.Id} has no charge states");
                Assert.That(feature.ChargeStateMin, Is.LessThanOrEqualTo(feature.ChargeStateMax));
                foreach (var chargeFeature in feature.SingleChargeFeatures)
                {
                    Assert.That(chargeFeature.Charge, Is.InRange(feature.ChargeStateMin, feature.ChargeStateMax),
                        $"feature {feature.Id} has a charge outside its own declared range");
                    Assert.That(chargeFeature.Envelopes.Count, Is.EqualTo(chargeFeature.EnvelopeNumber),
                        $"feature {feature.Id} z{chargeFeature.Charge}: envelope_num disagrees with the list length");
                    Assert.That(chargeFeature.Envelopes.Select(e => e.OneBasedScanNumber),
                        Is.Ordered.Ascending, $"feature {feature.Id} z{chargeFeature.Charge}: envelopes out of order");
                }
                // the top-level envelope_num is the total over charge states
                Assert.That(feature.AllEnvelopes.Count(), Is.EqualTo(feature.EnvelopeNumber),
                    $"feature {feature.Id}: total envelope_num disagrees with the flattened list");
            }
        }

        /// <summary>
        /// The xml stores SECONDS while the companion tsv stores minutes. The reader must hand back minutes, so
        /// retention times line up with MsDataScan and with everything built from the msAlign.
        /// </summary>
        [Test]
        public static void RetentionTimesAreConvertedToMinutes()
        {
            var features = TopFdFeatureXmlFile.ReadFile(FeatureXmlPath);
            var scans = SynthesiseScans();
            double firstScanRt = scans.First().RetentionTime, lastScanRt = scans.Last().RetentionTime;

            TestContext.Progress.WriteLine($"scan RT range from msAlign : {firstScanRt:F2}-{lastScanRt:F2} min");
            TestContext.Progress.WriteLine($"feature apex RT range      : " +
                $"{features.Min(f => f.RetentionTimeApex):F2}-{features.Max(f => f.RetentionTimeApex):F2} min");

            // the fixture's features lie inside the fixture's scan slice, so their apexes must too
            Assert.That(features.Min(f => f.RetentionTimeApex), Is.GreaterThan(firstScanRt - 0.5));
            Assert.That(features.Max(f => f.RetentionTimeApex), Is.LessThan(lastScanRt + 0.5));

            // and each envelope's own RT must agree with the scan it names
            var rtByScanNumber = scans.ToDictionary(s => s.OneBasedScanNumber, s => s.RetentionTime);
            int compared = 0;
            foreach (var envelope in features.SelectMany(f => f.AllEnvelopes))
            {
                if (!rtByScanNumber.TryGetValue(envelope.OneBasedScanNumber, out double scanRt)) continue;
                Assert.That(envelope.RetentionTime, Is.EqualTo(scanRt).Within(0.01),
                    $"envelope RT disagrees with scan {envelope.OneBasedScanNumber}");
                compared++;
            }
            TestContext.Progress.WriteLine($"{compared} envelope retention times agree with their scans");
            Assert.That(compared, Is.GreaterThan(0));

            // reading without conversion should give seconds, i.e. ~60x larger
            var raw = TopFdFeatureXmlFile.ReadFile(FeatureXmlPath, convertSecondsToMinutes: false);
            Assert.That(raw.First().RetentionTimeApex,
                Is.EqualTo(features.First().RetentionTimeApex * 60).Within(1e-6));
        }

        /// <summary>
        /// The xml and the tsv describe the same features, so masses and apex times must agree between them.
        /// </summary>
        [Test]
        public static void XmlAgreesWithTheFeatureTsv()
        {
            var xmlFeatures = TopFdFeatureXmlFile.ReadFile(FeatureXmlPath);
            var tsvFeatures = Ms1FeatureRecord.ReadFile(FeatureTsvPath);
            var tsvById = tsvFeatures.GroupBy(f => f.Id).ToDictionary(g => g.Key, g => g.First());

            int matched = 0;
            foreach (var xmlFeature in xmlFeatures)
            {
                if (!tsvById.TryGetValue(xmlFeature.Id, out var tsvFeature)) continue;
                Assert.That(xmlFeature.Mass, Is.EqualTo(tsvFeature.Mass).Within(0.001),
                    $"feature {xmlFeature.Id}: mass differs between xml and tsv");
                Assert.That(xmlFeature.RetentionTimeApex, Is.EqualTo(tsvFeature.RetentionTimeApex).Within(0.01),
                    $"feature {xmlFeature.Id}: apex time differs between xml and tsv (unit mismatch?)");
                Assert.That(xmlFeature.ChargeStateMin, Is.EqualTo(tsvFeature.ChargeStateMin));
                Assert.That(xmlFeature.ChargeStateMax, Is.EqualTo(tsvFeature.ChargeStateMax));
                matched++;
            }
            TestContext.Progress.WriteLine($"{matched} features cross-checked against the tsv: mass, apex time and " +
                "charge range all agree");
            Assert.That(matched, Is.GreaterThan(0), "no feature ids were shared between xml and tsv");
        }

        /// <summary>
        /// The whole point: feature XICs come out with no tracing. Prints the yield and length distribution, and
        /// the same numbers for the msAlign-tracing path over the identical scans, for comparison.
        /// </summary>
        [Test]
        public static void FeatureXmlXicsVersusMsAlignTracing()
        {
            var scans = SynthesiseScans();
            var features = TopFdFeatureXmlFile.ReadFile(FeatureXmlPath);

            var constructor = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, minNumberOfPeaks: 3,
                features, minMass: 3000, maxMass: 50000);
            var xics = constructor.GetAllXics(scans);

            TestContext.Progress.WriteLine($"feature.xml : {constructor.FeatureCountBeforeFilters} features -> " +
                $"{constructor.FeatureCount} after mass filter -> {xics.Count} XICs " +
                $"({constructor.LastFeaturesBelowMinPeaksCount} below minNumberOfPeaks, " +
                $"{constructor.LastFeaturesWithNoScansInRangeCount} with no scans in range)");
            if (xics.Count > 0)
            {
                var lengths = xics.Select(x => x.Peaks.Count).OrderBy(n => n).ToArray();
                TestContext.Progress.WriteLine($"              XIC length min {lengths.First()}, " +
                    $"median {lengths[lengths.Length / 2]}, max {lengths.Last()}");
            }

            // same scans, same mass window, but reconstructed by tracing the msAlign
            var traced = new FeatureGroupedXicConstructor(new PpmToleranceWithNotch(20, 1, 1), 2, 1, 3,
                ExternalDeconResults.ReadMsAlign(MsAlignPath), Ms1FeatureRecord.ReadFile(FeatureTsvPath),
                minMass: 3000, maxMass: 50000);
            var tracedXics = traced.GetAllXics(scans);
            TestContext.Progress.WriteLine($"msAlign     : {traced.LastTracedXicCount} traced -> {tracedXics.Count} " +
                $"feature XICs ({traced.LastFeaturesWithXicsCount}/{traced.FeatureCount} features)");

            Assert.That(xics.Count, Is.GreaterThan(0));
            // every emitted XIC must be traceable back to the feature it came from
            foreach (var xic in xics)
                Assert.That(constructor.LastFeatureByXic.ContainsKey(xic));
        }

        /// <summary>
        /// Charge handling: aggregated gives one XIC per feature with intensities summed per scan; unaggregated
        /// gives one XIC per charge state. The summed intensity of the aggregate must equal the sum over charges.
        /// </summary>
        [Test]
        public static void AggregateChargesSumsTheChargeStates()
        {
            var scans = SynthesiseScans();
            var features = TopFdFeatureXmlFile.ReadFile(FeatureXmlPath);

            var aggregated = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, 3, features,
                minMass: 3000, maxMass: 50000, aggregateCharges: true);
            var perCharge = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, 3, features,
                minMass: 3000, maxMass: 50000, aggregateCharges: false);

            var aggregatedXics = aggregated.GetAllXics(scans);
            var perChargeXics = perCharge.GetAllXics(scans);
            TestContext.Progress.WriteLine($"aggregated: {aggregatedXics.Count} XICs (one per feature); " +
                $"per charge: {perChargeXics.Count} XICs");
            Assert.That(perChargeXics.Count, Is.GreaterThanOrEqualTo(aggregatedXics.Count));

            // pick a feature that produced both, and check the intensities reconcile
            var feature = aggregated.LastFeatureByXic.Values.FirstOrDefault(f =>
                perCharge.LastFeatureByXic.Values.Count(v => v == f) > 1);
            if (feature == null) Assert.Ignore("no feature produced more than one charge XIC in this fixture");

            double aggregateTotal = aggregated.LastFeatureByXic.Where(kvp => kvp.Value == feature)
                .Sum(kvp => kvp.Key.Peaks.Sum(p => (double)p.Intensity));
            double perChargeTotal = perCharge.LastFeatureByXic.Where(kvp => kvp.Value == feature)
                .Sum(kvp => kvp.Key.Peaks.Sum(p => (double)p.Intensity));
            TestContext.Progress.WriteLine($"feature {feature.Id}: aggregated total {aggregateTotal:e3}, " +
                $"summed over charge XICs {perChargeTotal:e3}");
            Assert.That(aggregateTotal, Is.EqualTo(perChargeTotal).Within(0.001).Percent);
        }

        /// <summary>
        /// The third charge mode: instead of summing the charge states, represent the feature by the single most
        /// abundant one. The resulting XIC must be exactly the per-charge XIC of the strongest charge - same
        /// intensities, one constant charge on every peak - and so can never carry more signal than the summed XIC.
        /// </summary>
        [Test]
        public static void MostAbundantChargeUsesOnlyTheStrongestChargeState()
        {
            var scans = SynthesiseScans();
            var features = TopFdFeatureXmlFile.ReadFile(FeatureXmlPath);

            var summed = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, 3, features,
                minMass: 3000, maxMass: 50000, chargeMode: FeatureChargeMode.SumAllCharges);
            var dominant = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, 3, features,
                minMass: 3000, maxMass: 50000, chargeMode: FeatureChargeMode.MostAbundantCharge);
            var perCharge = new FeatureXmlXicConstructor(new PpmTolerance(20), 2, 1, 3, features,
                minMass: 3000, maxMass: 50000, chargeMode: FeatureChargeMode.SeparateCharges);

            var summedXics = summed.GetAllXics(scans);
            var dominantXics = dominant.GetAllXics(scans);
            var perChargeXics = perCharge.GetAllXics(scans);
            TestContext.Progress.WriteLine($"summed: {summedXics.Count} XICs; most-abundant-charge: " +
                $"{dominantXics.Count} XICs; per charge: {perChargeXics.Count} XICs");
            Assert.That(dominantXics.Count, Is.GreaterThan(0));

            // one XIC per feature, as with summing - but a single charge's trace can fall below MinNumberOfPeaks
            // where the summed trace cleared it, so this mode can only ever emit fewer
            Assert.That(dominant.LastFeatureByXic.Values.Distinct().Count(), Is.EqualTo(dominantXics.Count));
            Assert.That(dominantXics.Count, Is.LessThanOrEqualTo(summedXics.Count));

            // the chosen charge is recorded for every feature that elutes within these scans, which is at least
            // every feature that produced an XIC
            Assert.That(dominant.LastMostAbundantChargeByFeature.Count, Is.GreaterThanOrEqualTo(dominantXics.Count));
            TestContext.Progress.WriteLine(
                $"charge chosen for {dominant.LastMostAbundantChargeByFeature.Count} features, " +
                $"{dominant.LastMostAbundantChargeMatchingTopFdCount} agreeing with TopFD rep_charge");

            int comparedFeatures = 0;
            foreach (var kvp in dominant.LastFeatureByXic)
            {
                var feature = kvp.Value;

                // every peak carries the one chosen charge, not a per-scan winner
                var masses = kvp.Key.Peaks.OfType<IndexedMass>().ToList();
                Assert.That(masses.Count, Is.EqualTo(kvp.Key.Peaks.Count));
                Assert.That(masses.Select(m => m.Charge).Distinct().ToList(),
                    Is.EqualTo(new List<int> { dominant.LastMostAbundantChargeByFeature[feature] }));

                var perChargeForFeature = perCharge.LastFeatureByXic
                    .Where(x => x.Value == feature).Select(x => x.Key).ToList();
                if (perChargeForFeature.Count < 2)
                    continue; // single-charge feature: nothing was selected between

                double dominantTotal = kvp.Key.Peaks.Sum(peak => (double)peak.Intensity);
                double strongestPerChargeTotal = perChargeForFeature.Max(x => x.Peaks.Sum(peak => (double)peak.Intensity));
                double summedTotal = summed.LastFeatureByXic.Where(x => x.Value == feature)
                    .Sum(x => x.Key.Peaks.Sum(peak => (double)peak.Intensity));

                // it IS the strongest charge state's own XIC ...
                Assert.That(dominantTotal, Is.EqualTo(strongestPerChargeTotal).Within(0.001).Percent);
                // ... so it is a part of the summed XIC, never the whole of it
                Assert.That(dominantTotal, Is.LessThanOrEqualTo(summedTotal * 1.000001));
                comparedFeatures++;
            }

            Assert.That(comparedFeatures, Is.GreaterThan(0),
                "no multi-charge feature produced an XIC in this fixture, so nothing was actually selected between");
            TestContext.Progress.WriteLine($"{comparedFeatures} multi-charge features compared against their " +
                "per-charge XICs");
        }
    }
}
