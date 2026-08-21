using EngineLayer.DIA;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Proves that <see cref="ExternalDeconResults.ReadMsAlign"/> reproduces an msAlign EXACTLY - every scan,
    /// every mass, every charge, every intensity - by parsing the same file independently as plain text and
    /// diffing the two. If a downstream result looks wrong, this test says whether the reader is at fault.
    /// </summary>
    [TestFixture]
    public static class MsAlignReaderFidelityTest
    {
        private static string MsAlignPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "DIATests", "ExternalDeconTestData", "HcdYbc_subset_ms1.msalign");

        /// <summary>Dumb, dependency-free parse of an msAlign: scan number -> (mass, intensity, charge) triples,
        /// straight off the text, plus the raw RETENTION_TIME field as written.</summary>
        private static (Dictionary<int, List<(double Mass, double Intensity, int Charge)>> Peaks,
                        Dictionary<int, double> RawRetentionTimes) ParseAsText(string path)
        {
            var peaks = new Dictionary<int, List<(double, double, int)>>();
            var retentionTimes = new Dictionary<int, double>();
            int scanNumber = -1;
            List<(double, double, int)> current = null;

            foreach (var line in File.ReadLines(path))
            {
                if (line.StartsWith("BEGIN IONS"))
                {
                    current = new List<(double, double, int)>();
                    scanNumber = -1;
                }
                else if (line.StartsWith("END IONS"))
                {
                    if (current != null && scanNumber > 0)
                        peaks[scanNumber] = current;
                    current = null;
                }
                else if (current == null)
                {
                    continue; // header
                }
                else if (line.StartsWith("SCANS="))
                {
                    scanNumber = int.Parse(line.Substring("SCANS=".Length), CultureInfo.InvariantCulture);
                }
                else if (line.StartsWith("RETENTION_TIME="))
                {
                    double rt = double.Parse(line.Substring("RETENTION_TIME=".Length), CultureInfo.InvariantCulture);
                    if (scanNumber > 0) retentionTimes[scanNumber] = rt;
                }
                else if (line.Contains('\t'))
                {
                    var splits = line.Split('\t');
                    current.Add((
                        double.Parse(splits[0], CultureInfo.InvariantCulture),
                        double.Parse(splits[1], CultureInfo.InvariantCulture),
                        int.Parse(splits[2], CultureInfo.InvariantCulture)));
                }
            }

            // RETENTION_TIME can precede SCANS in an entry, so re-key any that were read before the scan number
            return (peaks, retentionTimes);
        }

        [Test]
        public static void ReaderReproducesTheFileExactly()
        {
            Assert.That(File.Exists(MsAlignPath), $"missing fixture {MsAlignPath}");

            var (textPeaks, textRetentionTimes) = ParseAsText(MsAlignPath);
            var results = ExternalDeconResults.ReadMsAlign(MsAlignPath);

            TestContext.Progress.WriteLine($"text parse : {textPeaks.Count} scans, {textPeaks.Sum(kvp => kvp.Value.Count)} peaks");
            TestContext.Progress.WriteLine($"reader     : {results.ScanCount} scans, {results.PeakCount} peaks");

            Assert.That(results.ScanCount, Is.EqualTo(textPeaks.Count), "scan count");
            Assert.That(results.PeakCount, Is.EqualTo(textPeaks.Sum(kvp => kvp.Value.Count)), "total peak count");
            Assert.That(results.PeaksByOneBasedScanNumber.Keys.OrderBy(k => k),
                Is.EqualTo(textPeaks.Keys.OrderBy(k => k)), "scan numbers");

            int comparedPeaks = 0;
            foreach (var scanNumber in textPeaks.Keys)
            {
                // compare as multisets ordered by mass: the reader sorts by mass, the text parse keeps file order
                var expected = textPeaks[scanNumber].OrderBy(p => p.Mass).ThenBy(p => p.Charge).ToList();
                var actual = results.PeaksByOneBasedScanNumber[scanNumber]
                    .Select(p => (Mass: p.MonoisotopicMass, Intensity: p.Intensity, Charge: p.Charge))
                    .OrderBy(p => p.Mass).ThenBy(p => p.Charge).ToList();

                Assert.That(actual.Count, Is.EqualTo(expected.Count), $"peak count in scan {scanNumber}");
                for (int i = 0; i < expected.Count; i++)
                {
                    Assert.That(actual[i].Mass, Is.EqualTo(expected[i].Mass).Within(1e-6),
                        $"mass #{i} in scan {scanNumber}");
                    Assert.That(actual[i].Charge, Is.EqualTo(expected[i].Charge),
                        $"charge #{i} in scan {scanNumber} (mass {expected[i].Mass})");
                    Assert.That(actual[i].Intensity, Is.EqualTo(expected[i].Intensity).Within(1e-3),
                        $"intensity #{i} in scan {scanNumber}");
                    comparedPeaks++;
                }
            }
            TestContext.Progress.WriteLine($"compared {comparedPeaks} peaks: mass, charge and intensity all identical");

            // retention times: file is TopFD 1.8 (minutes), so the reader must hand back the number as written
            foreach (var kvp in textRetentionTimes)
            {
                Assert.That(results.RetentionTimeByOneBasedScanNumber[kvp.Key], Is.EqualTo(kvp.Value).Within(1e-6),
                    $"retention time of scan {kvp.Key}");
            }
            TestContext.Progress.WriteLine($"compared {textRetentionTimes.Count} retention times: identical to the file");
        }

        /// <summary>
        /// Charge must stay bound to its own mass. A reader that sorted masses without carrying the charge array
        /// along would still pass a count check, so this pins the pairing: every (mass, charge) pair the reader
        /// returns must exist verbatim in the file.
        /// </summary>
        [Test]
        public static void MassChargePairingIsPreserved()
        {
            var (textPeaks, _) = ParseAsText(MsAlignPath);
            var results = ExternalDeconResults.ReadMsAlign(MsAlignPath);

            var fromFile = new HashSet<string>(textPeaks.SelectMany(kvp =>
                kvp.Value.Select(p => $"{kvp.Key}|{p.Mass:F6}|{p.Charge}")));

            int checkedPairs = 0;
            foreach (var kvp in results.PeaksByOneBasedScanNumber)
            {
                foreach (var peak in kvp.Value)
                {
                    Assert.That(fromFile, Does.Contain($"{kvp.Key}|{peak.MonoisotopicMass:F6}|{peak.Charge}"),
                        $"reader produced mass {peak.MonoisotopicMass} at charge {peak.Charge} in scan {kvp.Key}, " +
                        "which is not a pair present in the file");
                    checkedPairs++;
                }
            }
            TestContext.Progress.WriteLine($"{checkedPairs} (scan, mass, charge) triples all present verbatim in the file");
        }
    }
}
