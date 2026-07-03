using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// Writes ISD/DIA pseudo-MS2 scans (produced by <see cref="ISDEngine"/> / <see cref="DIAEngine"/> via
    /// consensus XIC construction + precursor-fragment grouping) to a TopPIC-compatible ms2.msalign file so
    /// they can be searched directly in TopPIC as if they were ordinary DDA MS2 spectra.
    ///
    /// The pseudo scans are already neutral-mass deconvoluted (each experimental fragment is an
    /// <see cref="IsotopicEnvelope"/> with a monoisotopic mass and charge), which is exactly what an
    /// ms2.msalign encodes. A companion ms1.msalign lists each distinct precursor as a single-peak MS1 entry
    /// so TopPIC has the precursor features it expects.
    /// </summary>
    public static class IsdMsAlignExporter
    {
        private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

        /// <summary>
        /// Write pseudo-MS2 scans to a ms2.msalign file (TopPIC format).
        /// </summary>
        /// <param name="pseudoScans">Pseudo-MS2 scans, e.g. from MetaMorpheusTask.GetMs2Scans on an ISD file.</param>
        /// <param name="filePath">Output path; conventionally ends in _ms2.msalign.</param>
        /// <param name="activation">Activation label written to each entry (ISD source-region fragmentation ≈ HCD).</param>
        public static void WriteMs2Align(IEnumerable<Ms2ScanWithSpecificMass> pseudoScans, string filePath,
            DissociationType activation = DissociationType.HCD)
        {
            if (pseudoScans == null) throw new ArgumentNullException(nameof(pseudoScans));

            using var sw = new StreamWriter(filePath, append: false);
            int id = 0;
            foreach (var scan in pseudoScans.OrderBy(s => s.OneBasedScanNumber))
            {
                id++;
                int scanNum = scan.OneBasedScanNumber > 0 ? scan.OneBasedScanNumber : id;
                int ms1Scan = scan.OneBasedPrecursorScanNumber ?? 0;

                sw.WriteLine("BEGIN IONS");
                sw.WriteLine($"ID={id}");
                sw.WriteLine("FRACTION_ID=0");
                sw.WriteLine($"SCANS={scanNum}");
                // msalign retention time is in seconds; MsDataScan.RetentionTime is in minutes
                sw.WriteLine($"RETENTION_TIME={(scan.RetentionTime * 60.0).ToString("F2", Ci)}");
                sw.WriteLine("LEVEL=2");
                sw.WriteLine($"ACTIVATION={activation.ToString().ToUpperInvariant()}");
                sw.WriteLine($"MS_ONE_ID={ms1Scan}");
                sw.WriteLine($"MS_ONE_SCAN={ms1Scan}");
                sw.WriteLine($"PRECURSOR_MZ={scan.PrecursorMonoisotopicPeakMz.ToString("F6", Ci)}");
                sw.WriteLine($"PRECURSOR_CHARGE={scan.PrecursorCharge}");
                sw.WriteLine($"PRECURSOR_MASS={scan.PrecursorMass.ToString("F6", Ci)}");
                sw.WriteLine($"PRECURSOR_INTENSITY={scan.PrecursorIntensity.ToString("F2", Ci)}");

                // fragment peaks: monoisotopic-neutral-mass  intensity  charge   (tab separated)
                foreach (var frag in (scan.ExperimentalFragments ?? Array.Empty<IsotopicEnvelope>())
                             .OrderBy(f => f.MonoisotopicMass))
                {
                    int z = frag.Charge != 0 ? Math.Abs(frag.Charge) : 1;
                    sw.WriteLine(string.Join("\t",
                        frag.MonoisotopicMass.ToString("F5", Ci),
                        frag.TotalIntensity.ToString("F2", Ci),
                        z.ToString(Ci)));
                }
                sw.WriteLine("END IONS");
                sw.WriteLine();
            }
        }

        private const double Proton = 1.007276466879;

        /// <summary>
        /// Write consensus <see cref="MassFeature"/>s to a TopFD/FlashDeconv-style `.ms1.feature` file that
        /// mzLib's <c>FromFileDeconvolutionParameters</c> can load — so MetaMorpheus can assemble precursors
        /// (and thus <c>Ms2ScanWithSpecificMass</c>) from these traced features (the "FromFile" path). Filters
        /// to intact proteoforms on the precursor side (mass ≥ minMass, ≥ minChargeCount charge states).
        /// Path MUST end in `.ms1.feature` for the reader's format auto-detection.
        /// </summary>
        public static void WriteMs1FeatureFile(IEnumerable<MassFeature> features, string filePath,
            double minMass = 3000, int minChargeCount = 3)
        {
            if (features == null) throw new ArgumentNullException(nameof(features));
            using var sw = new StreamWriter(filePath, append: false);
            sw.WriteLine("Feature_ID\tMass\tIntensity\tMin_time\tMax_time\tApex_time\tIntensity_Apex\tMin_charge\tMax_charge");
            int id = 0;
            foreach (var f in features)
            {
                if (f.ConsensusMass < minMass) continue;
                if (f.ChargeCount < minChargeCount) continue;
                id++;
                int zmin = f.Charges.Min();
                int zmax = f.Charges.Max();
                double apex = (f.RTStart + f.RTEnd) / 2.0;
                sw.WriteLine(string.Join("\t",
                    id.ToString(Ci),
                    f.ConsensusMass.ToString("F5", Ci),
                    f.SummedIntensity.ToString("F1", Ci),
                    f.RTStart.ToString("F4", Ci),
                    f.RTEnd.ToString("F4", Ci),
                    apex.ToString("F4", Ci),
                    f.SummedIntensity.ToString("F1", Ci),
                    zmin.ToString(Ci),
                    zmax.ToString(Ci)));
            }
        }

        /// <summary>
        /// Write pseudo-MS2 scans to an MGF (one BEGIN IONS block per scan) whose fragment peaks are the
        /// neutral fragment masses written as singly-charged m/z (mass + proton). This matches a top-down
        /// search configured with UseProvidedPrecursorInfo + product MaxAssumedChargeState=1 (e.g. the
        /// project's td_pseudoMS2.toml), so consensus-path pseudo scans can be searched with the exact same
        /// config as the earlier fd_profile MGFs — an apples-to-apples ID comparison.
        /// </summary>
        public static void WriteMgf(IEnumerable<Ms2ScanWithSpecificMass> pseudoScans, string filePath)
        {
            if (pseudoScans == null) throw new ArgumentNullException(nameof(pseudoScans));
            using var sw = new StreamWriter(filePath, append: false);
            int id = 0;
            foreach (var scan in pseudoScans.OrderBy(s => s.OneBasedScanNumber))
            {
                id++;
                int z = scan.PrecursorCharge != 0 ? Math.Abs(scan.PrecursorCharge) : 1;
                sw.WriteLine("BEGIN IONS");
                sw.WriteLine($"TITLE=isd.consensus.scan{id}");
                sw.WriteLine($"PEPMASS={scan.PrecursorMonoisotopicPeakMz.ToString("F5", Ci)}");
                sw.WriteLine($"CHARGE={z}+");
                sw.WriteLine($"RTINSECONDS={(scan.RetentionTime * 60.0).ToString("F1", Ci)}");
                sw.WriteLine($"SCANS={id}");
                foreach (var frag in (scan.ExperimentalFragments ?? Array.Empty<IsotopicEnvelope>())
                             .OrderBy(f => f.MonoisotopicMass))
                {
                    sw.WriteLine($"{(frag.MonoisotopicMass + Proton).ToString("F5", Ci)} {frag.TotalIntensity.ToString("F1", Ci)}");
                }
                sw.WriteLine("END IONS");
            }
        }

        /// <summary>
        /// Write a minimal companion ms1.msalign: one MS1 entry per distinct precursor (mass, intensity, charge).
        /// TopPIC pairs ms1.msalign with ms2.msalign; this provides the precursor features the pseudo scans encode.
        /// </summary>
        public static void WriteMs1Align(IEnumerable<Ms2ScanWithSpecificMass> pseudoScans, string filePath)
        {
            if (pseudoScans == null) throw new ArgumentNullException(nameof(pseudoScans));

            using var sw = new StreamWriter(filePath, append: false);
            int id = 0;
            // one MS1 entry per unique precursor scan number (fall back to per-scan if none)
            foreach (var grp in pseudoScans
                         .GroupBy(s => s.OneBasedPrecursorScanNumber ?? s.OneBasedScanNumber)
                         .OrderBy(g => g.Key))
            {
                id++;
                var rep = grp.First();
                int scanNum = grp.Key > 0 ? grp.Key : id;
                sw.WriteLine("BEGIN IONS");
                sw.WriteLine($"ID={id}");
                sw.WriteLine("FRACTION_ID=0");
                sw.WriteLine($"SCANS={scanNum}");
                sw.WriteLine($"RETENTION_TIME={(rep.RetentionTime * 60.0).ToString("F2", Ci)}");
                sw.WriteLine("LEVEL=1");
                foreach (var s in grp
                             .GroupBy(s => Math.Round(s.PrecursorMass, 2))
                             .Select(g => g.First()))
                {
                    int z = s.PrecursorCharge != 0 ? Math.Abs(s.PrecursorCharge) : 1;
                    sw.WriteLine(string.Join("\t",
                        s.PrecursorMass.ToString("F5", Ci),
                        s.PrecursorIntensity.ToString("F2", Ci),
                        z.ToString(Ci)));
                }
                sw.WriteLine("END IONS");
                sw.WriteLine();
            }
        }
    }
}
