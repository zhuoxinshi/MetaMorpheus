using EngineLayer.DIA;
using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    /// <summary>
    /// Proves <see cref="MzrtXicWriter"/> emits a TopDIA-style mzrt.csv that round-trips through the reader
    /// <see cref="TopDiaMzrtFile"/>: a precursor and a co-eluting fragment are written as two rows on one common
    /// cycle axis (dense XIC vectors), and reading the file back recovers each XIC's cycles and intensities.
    /// </summary>
    [TestFixture]
    public static class MzrtXicWriterTests
    {
        [Test]
        public static void Write_ProducesTopDiaMzrt_ThatRoundTripsThroughReader()
        {
            // precursor (mass 10000, z10) and a co-eluting fragment (mass 2000, z2), both over cycles 5,6,7.
            var precursor = Xic(10000.0, 10, (5, 5.00, 100), (6, 5.10, 300), (7, 5.20, 150)); // apex at cycle 6
            var fragment  = Xic(2000.0, 2,  (5, 5.00, 50),  (6, 5.10, 150), (7, 5.20, 80));

            var pair = new PrecursorFragmentPair(precursor, fragment, correlation: 0.98, overlap: 1.0) { FragmentRank = 1 };
            var group = new PrecursorFragmentsGroup(precursor, new List<PrecursorFragmentPair> { pair })
            {
                PFgroupIndex = 42,
                PrecursorFeatureId = 7
            };

            string dir = FreshDir();
            string path = Path.Combine(dir, "pf.mzrt.csv");
            MzrtXicWriter.Write(new[] { group }, path);

            // ---- read it back with the existing TopDIA reader ----
            var records = TopDiaMzrtFile.ReadFile(path, out int cycleCount);
            Assert.That(cycleCount, Is.EqualTo(8), "dense vector length = largest cycle index (7) + 1");
            Assert.That(records, Has.Count.EqualTo(2), "one row for the precursor, one for the fragment");
            Assert.That(records.All(r => r.Id == 42), "both rows share the group id");

            var prec = records.Single(r => r.Charge == 10);
            Assert.That(prec.Mass, Is.EqualTo(10000.0).Within(1e-3));
            Assert.That(prec.SpectrumLow, Is.EqualTo(5));
            Assert.That(prec.SpectrumHigh, Is.EqualTo(7));
            Assert.That(prec.FirstCycle, Is.EqualTo(5), "trace starts at cycle 5");
            Assert.That(prec.Intensities.Select(v => (double)v).ToArray(),
                Is.EqualTo(new[] { 100.0, 300.0, 150.0 }).Within(1e-3), "precursor intensities at cycles 5,6,7");

            var frag = records.Single(r => r.Charge == 2);
            Assert.That(frag.Mass, Is.EqualTo(2000.0).Within(1e-3));
            Assert.That(frag.Intensities.Select(v => (double)v).ToArray(),
                Is.EqualTo(new[] { 50.0, 150.0, 80.0 }).Within(1e-3), "fragment intensities at cycles 5,6,7");

            // ---- the extra grouping columns are present in the raw file ----
            var lines = File.ReadAllLines(path);
            Assert.That(lines[0], Does.Contain("role").And.Contain("correlation").And.Contain("XIC").And.Contain("Envelope"));
            var header = lines[0].Split(',');
            int roleCol = System.Array.IndexOf(header, "role");
            int corrCol = System.Array.IndexOf(header, "correlation");
            var fragRow = lines.Skip(1).Select(l => l.Split(',')).Single(c => c[header.ToList().IndexOf("Charge")] == "2");
            Assert.That(fragRow[roleCol], Is.EqualTo("fragment"));
            Assert.That(fragRow[corrCol], Is.EqualTo("0.9800"), "fragment correlation to the precursor is written");
        }

        // Build an ExtractedIonChromatogram from synthetic mass peaks (one IndexedMass per cycle).
        private static ExtractedIonChromatogram Xic(double mass, int charge, params (int scan, double rt, double inten)[] pts)
        {
            var peaks = pts.Select(p => (IIndexedPeak)new IndexedMass(
                new IsotopicEnvelope(0, new List<(double, double)> { (mass, p.inten) }, mass, charge, p.inten, 0),
                p.rt, p.scan, 1)).ToList();
            return new ExtractedIonChromatogram(peaks);
        }

        private static string FreshDir()
        {
            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "MzrtXicWriterTests",
                TestContext.CurrentContext.Test.Name);
            if (Directory.Exists(dir)) Directory.Delete(dir, true);
            Directory.CreateDirectory(dir);
            return dir;
        }
    }
}
