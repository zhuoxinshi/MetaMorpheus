using Chemistry;
using MassSpectrometry;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace EngineLayer.DIA
{
    /// <summary>
    /// How an <see cref="MsDataScan"/> from the raw file is matched to its entry in an external
    /// deconvolution result file.
    /// </summary>
    public enum ExternalDeconScanMatchMode
    {
        /// <summary>Match on SCANS= (the msAlign entry's one-based scan number). Correct whenever the
        /// deconvolution software was run on the same file, or on a channel-split file that preserved the
        /// parent file's scan numbering (which is what the ISD channel splits fed to TopFD do).</summary>
        OneBasedScanNumber,

        /// <summary>Match on retention time, nearest entry within a tolerance. Use this when the file handed to
        /// the deconvolution software was renumbered from 1 (some mzML splitters do this), so scan numbers no
        /// longer line up with the raw file.</summary>
        RetentionTime
    }

    /// <summary>
    /// The unit the RETENTION_TIME field of an msAlign is written in. mzLib's reader assumes SECONDS and always
    /// divides by 60, which is right for FLASHDeconv and TopFD through 1.7.x but wrong for TopFD 1.8+, which
    /// writes MINUTES - leaving those retention times 60x too small.
    /// </summary>
    public enum MsAlignRetentionTimeUnit
    {
        /// <summary>Read the software and version out of the msAlign header and pick the unit from it: TopFD 1.8
        /// and later are minutes, everything else seconds.</summary>
        Auto,
        /// <summary>The file is in seconds (mzLib's assumption); take its conversion as-is.</summary>
        Seconds,
        /// <summary>The file is in minutes; undo mzLib's divide-by-60.</summary>
        Minutes
    }

    /// <summary>
    /// One deconvoluted peak read from an external deconvolution result file: a neutral monoisotopic mass,
    /// its intensity and its charge state. This is the msAlign peak line
    /// (mass \t intensity \t charge [\t score]) in object form.
    /// </summary>
    public class ExternalDeconPeak
    {
        public double MonoisotopicMass { get; }
        public double Intensity { get; }
        public int Charge { get; }

        public ExternalDeconPeak(double monoisotopicMass, double intensity, int charge)
        {
            MonoisotopicMass = monoisotopicMass;
            Intensity = intensity;
            Charge = charge;
        }

        /// <summary>The m/z this mass was observed at, i.e. the mass folded back onto its charge state.
        /// Used to apply an isolation-window filter, since external results carry no m/z of their own.</summary>
        public double Mz => MonoisotopicMass.ToMz(Charge == 0 ? 1 : Charge);

        /// <summary>Converts to the mzLib envelope type the indexing engine consumes. The single synthetic peak
        /// carries the whole intensity so that <see cref="IndexedMass.Intensity"/> ends up being the
        /// deconvoluted intensity reported by the external software.</summary>
        public IsotopicEnvelope ToIsotopicEnvelope()
        {
            int charge = Charge == 0 ? 1 : Charge;
            return new IsotopicEnvelope(
                new List<(double mz, double intensity)> { (MonoisotopicMass.ToMz(charge), Intensity) },
                MonoisotopicMass, charge, Intensity, 0);
        }
    }

    /// <summary>
    /// Deconvolution results produced OUTSIDE MetaMorpheus and read back in from msAlign files
    /// (TopFD *_ms1.msalign / *_ms2.msalign, and equally FLASHDeconv or IsoDec msAligns - the format is
    /// shared), indexed by the scan they came from.
    ///
    /// This is the "external decon" half of the ISD/DIA workflow: instead of running MetaMorpheus's own
    /// deconvolution per scan (<see cref="NeutralMassXicConstructor"/> -> <see cref="Deconvoluter"/>), the
    /// per-scan neutral masses come from TopFD, and everything downstream - XIC construction,
    /// precursor-fragment grouping, pseudo-MS2 generation - proceeds unchanged. See
    /// <see cref="ExternalDeconXicConstructor"/>, which consumes this class.
    ///
    /// Several files can be loaded into one instance; entries are keyed by one-based scan number, so an ISD run
    /// whose channels were deconvoluted separately (one msAlign per source voltage) merges cleanly as long as
    /// the split preserved the parent file's scan numbers. If it did not, match on retention time instead
    /// (<see cref="ExternalDeconScanMatchMode.RetentionTime"/>).
    /// </summary>
    public class ExternalDeconResults
    {
        /// <summary>Deconvoluted peaks keyed by the one-based scan number written in the msAlign SCANS= field.</summary>
        public Dictionary<int, List<ExternalDeconPeak>> PeaksByOneBasedScanNumber { get; }

        /// <summary>Retention time (MINUTES) of each msAlign entry, keyed by one-based scan number. msAlign stores
        /// seconds; the mzLib reader converts to minutes, matching <see cref="MsDataScan.RetentionTime"/>.</summary>
        public Dictionary<int, double> RetentionTimeByOneBasedScanNumber { get; }

        /// <summary>The msAlign files this instance was read from, in load order.</summary>
        public List<string> SourceFilePaths { get; }

        /// <summary>The source files as a single string. Scalar so that <see cref="DIAparameters.GetSettingsString"/>
        /// records WHICH external deconvolution went into a given search.</summary>
        public string SourceFiles => string.Join("; ", SourceFilePaths);

        /// <summary>Entries sorted by retention time, for nearest-RT matching.</summary>
        private readonly List<(double RetentionTime, int OneBasedScanNumber)> _entriesByRetentionTime;

        private ExternalDeconResults(Dictionary<int, List<ExternalDeconPeak>> peaksByScan,
            Dictionary<int, double> rtByScan, List<string> sourceFilePaths)
        {
            PeaksByOneBasedScanNumber = peaksByScan;
            RetentionTimeByOneBasedScanNumber = rtByScan;
            SourceFilePaths = sourceFilePaths;
            _entriesByRetentionTime = rtByScan
                .Select(kvp => (RetentionTime: kvp.Value, OneBasedScanNumber: kvp.Key))
                .OrderBy(e => e.RetentionTime)
                .ToList();
        }

        public int ScanCount => PeaksByOneBasedScanNumber.Count;

        public int PeakCount => PeaksByOneBasedScanNumber.Sum(kvp => kvp.Value.Count);

        /// <summary>Pattern for the ISD source voltage in a Thermo scan filter ("... sid=60.00 Full ms ..."), the
        /// same one <see cref="ISDEngine"/> uses to group ISD scans into channels.</summary>
        public static readonly Regex SourceVoltagePattern = new Regex(@"sid=(\d+)");

        /// <summary>
        /// Maps one-based scan number to ISD source voltage by reading the sid= out of each scan filter. Scans
        /// with no sid= are omitted. This is the authoritative channel assignment: it comes from the raw file
        /// itself and makes no assumption about cycle length or scan ordering.
        /// </summary>
        public static Dictionary<int, int> GetIsdVoltageByScanNumber(IEnumerable<MsDataScan> scans)
        {
            var voltageByScanNumber = new Dictionary<int, int>();
            foreach (var scan in scans ?? Enumerable.Empty<MsDataScan>())
            {
                if (scan == null || string.IsNullOrEmpty(scan.ScanFilter))
                    continue;
                var match = SourceVoltagePattern.Match(scan.ScanFilter);
                if (match.Success)
                    voltageByScanNumber[scan.OneBasedScanNumber] = int.Parse(match.Groups[1].Value);
            }
            return voltageByScanNumber;
        }

        public static ExternalDeconResults ReadMsAlign(params string[] filePaths) =>
            ReadMsAlign((IEnumerable<string>)filePaths);

        /// <summary>
        /// Reads one or more msAlign files and merges their entries into a single scan-number-keyed lookup.
        /// A scan number appearing in more than one file has its peak lists concatenated, which is what you
        /// want when the ms1 and ms2 msAligns of a run are loaded together (they never share scan numbers)
        /// but would double-count if the same file were passed twice.
        /// </summary>
        public static ExternalDeconResults ReadMsAlign(IEnumerable<string> filePaths) =>
            ReadMsAlign(filePaths, MsAlignRetentionTimeUnit.Auto);

        /// <summary>
        /// Reads one or more msAlign files, interpreting RETENTION_TIME per <paramref name="retentionTimeUnit"/>.
        /// </summary>
        public static ExternalDeconResults ReadMsAlign(IEnumerable<string> filePaths,
            MsAlignRetentionTimeUnit retentionTimeUnit)
        {
            var peaksByScan = new Dictionary<int, List<ExternalDeconPeak>>();
            var rtByScan = new Dictionary<int, double>();
            var loaded = new List<string>();

            foreach (var filePath in filePaths ?? Enumerable.Empty<string>())
            {
                if (!File.Exists(filePath))
                    throw new MetaMorpheusException("External deconvolution file not found: " + filePath);

                // Ms1Align vs Ms2Align differ only in the MS level assumed when an entry has no LEVEL= line,
                // so pick by file name and let the entries themselves carry the truth. Constructed directly
                // rather than through MsDataFileReader so that msAligns not named *_ms1/*_ms2 still load -
                // TopFD names every ISD channel split *_ms1.msalign regardless of what the channel contains.
                MsAlign msAlign = filePath.EndsWith("_ms2.msalign", StringComparison.OrdinalIgnoreCase)
                    ? new Ms2Align(filePath)
                    : new Ms1Align(filePath);
                msAlign.LoadAllStaticData();
                loaded.Add(filePath);

                // mzLib always divides RETENTION_TIME by 60 on the assumption the file is in seconds. Undo that
                // for files that were actually written in minutes, so the stored RTs are minutes either way.
                double retentionTimeScale = ResolveRetentionTimeUnit(filePath, retentionTimeUnit) ==
                    MsAlignRetentionTimeUnit.Minutes ? 60.0 : 1.0;

                foreach (var scan in msAlign.GetAllScansList())
                {
                    if (scan == null)
                        continue;

                    if (!(scan.MassSpectrum is NeutralMassSpectrum neutralMassSpectrum))
                        throw new MetaMorpheusException(
                            "Expected neutral-mass spectra in external deconvolution file " + filePath +
                            ", but scan " + scan.OneBasedScanNumber + " was not deconvoluted.");

                    if (!peaksByScan.TryGetValue(scan.OneBasedScanNumber, out var peakList))
                    {
                        peakList = new List<ExternalDeconPeak>();
                        peaksByScan[scan.OneBasedScanNumber] = peakList;
                    }

                    for (int i = 0; i < neutralMassSpectrum.XArray.Length; i++)
                    {
                        peakList.Add(new ExternalDeconPeak(
                            neutralMassSpectrum.XArray[i],
                            neutralMassSpectrum.YArray[i],
                            neutralMassSpectrum.Charges[i]));
                    }

                    rtByScan[scan.OneBasedScanNumber] = scan.RetentionTime * retentionTimeScale;
                }
            }

            if (loaded.Count == 0)
                throw new MetaMorpheusException("No external deconvolution files were provided.");

            return new ExternalDeconResults(peaksByScan, rtByScan, loaded);
        }

        /// <summary>
        /// Decides which unit an msAlign's RETENTION_TIME is in. Only TopFD 1.8 and later write minutes; TopFD
        /// 1.7.x and FLASHDeconv write seconds, which is what mzLib assumes.
        /// </summary>
        private static MsAlignRetentionTimeUnit ResolveRetentionTimeUnit(string filePath, MsAlignRetentionTimeUnit requested)
        {
            if (requested != MsAlignRetentionTimeUnit.Auto)
                return requested;

            // the software and version live on the first header line, e.g. "#TopFD 1.8.1"
            foreach (var line in File.ReadLines(filePath).Take(40))
            {
                if (!line.StartsWith("#"))
                    break;

                var match = Regex.Match(line, @"TopFD\s*v?(\d+)\.(\d+)", RegexOptions.IgnoreCase);
                if (!match.Success)
                    continue;

                int major = int.Parse(match.Groups[1].Value);
                int minor = int.Parse(match.Groups[2].Value);
                return (major > 1 || (major == 1 && minor >= 8))
                    ? MsAlignRetentionTimeUnit.Minutes
                    : MsAlignRetentionTimeUnit.Seconds;
            }

            return MsAlignRetentionTimeUnit.Seconds;
        }

        /// <summary>
        /// Splits these results into one <see cref="ExternalDeconResults"/> per ISD source voltage, using an
        /// explicit scan-number to voltage map (build one from the raw file with
        /// <see cref="GetIsdVoltageByScanNumber"/>). Scans absent from the map are dropped.
        ///
        /// This is what a SINGLE combined <c>*_ms2.msalign</c> needs: when the deconvolution was run on the whole
        /// ISD file rather than on per-channel splits, every voltage's spectra live in one file and are told apart
        /// only by their scan number.
        /// </summary>
        public Dictionary<int, ExternalDeconResults> SplitByVoltage(IReadOnlyDictionary<int, int> voltageByOneBasedScanNumber)
        {
            if (voltageByOneBasedScanNumber == null || voltageByOneBasedScanNumber.Count == 0)
                throw new MetaMorpheusException("No scan-number to ISD voltage map was supplied.");

            var peaksByVoltage = new Dictionary<int, Dictionary<int, List<ExternalDeconPeak>>>();
            var rtByVoltage = new Dictionary<int, Dictionary<int, double>>();

            foreach (var kvp in PeaksByOneBasedScanNumber)
            {
                if (!voltageByOneBasedScanNumber.TryGetValue(kvp.Key, out int voltage))
                    continue;

                if (!peaksByVoltage.TryGetValue(voltage, out var peaks))
                {
                    peaks = new Dictionary<int, List<ExternalDeconPeak>>();
                    peaksByVoltage[voltage] = peaks;
                    rtByVoltage[voltage] = new Dictionary<int, double>();
                }

                peaks[kvp.Key] = kvp.Value;
                if (RetentionTimeByOneBasedScanNumber.TryGetValue(kvp.Key, out double retentionTime))
                    rtByVoltage[voltage][kvp.Key] = retentionTime;
            }

            if (peaksByVoltage.Count == 0)
                throw new MetaMorpheusException(
                    "None of the " + ScanCount + " scans in the external deconvolution results (" +
                    string.Join(", ", SourceFilePaths) + ") appear in the supplied scan-number to voltage map.");

            return peaksByVoltage.ToDictionary(
                kvp => kvp.Key,
                kvp => new ExternalDeconResults(kvp.Value, rtByVoltage[kvp.Key], SourceFilePaths));
        }

        /// <summary>
        /// Splits these results by ISD source voltage using the repeating scan CYCLE, for when the raw file is not
        /// at hand. <paramref name="voltagesInCycleOrder"/> lists the voltage of each position in one cycle, in
        /// acquisition order - for an ISD run acquired as MS1, 60, 80, 100 that is <c>{ 0, 60, 80, 100 }</c>, with
        /// 0 marking the MS1 position. Scan n belongs to cycle position
        /// <c>(n - firstScanNumber) mod cycleLength</c>.
        ///
        /// Prefer <see cref="SplitByVoltage"/> with <see cref="GetIsdVoltageByScanNumber"/> when the raw file is
        /// available: it reads the actual scan filters and so survives a run whose cycle is interrupted.
        /// </summary>
        public Dictionary<int, ExternalDeconResults> SplitByIsdCycle(IReadOnlyList<int> voltagesInCycleOrder,
            int firstScanNumber = 1)
        {
            if (voltagesInCycleOrder == null || voltagesInCycleOrder.Count == 0)
                throw new MetaMorpheusException("No ISD cycle was supplied to split the external results by.");

            int cycleLength = voltagesInCycleOrder.Count;
            var voltageByScanNumber = new Dictionary<int, int>();
            foreach (int scanNumber in PeaksByOneBasedScanNumber.Keys)
            {
                int offset = ((scanNumber - firstScanNumber) % cycleLength + cycleLength) % cycleLength;
                voltageByScanNumber[scanNumber] = voltagesInCycleOrder[offset];
            }

            return SplitByVoltage(voltageByScanNumber);
        }

        /// <summary>
        /// Returns the external deconvolution peaks belonging to <paramref name="scan"/>, or null if that scan
        /// has no entry in the loaded files (a scan the external software skipped, or a scan belonging to a
        /// channel that was not deconvoluted).
        /// </summary>
        /// <param name="retentionTimeTolerance">Half-width, in MINUTES, of the RT window used when
        /// <paramref name="matchMode"/> is <see cref="ExternalDeconScanMatchMode.RetentionTime"/>.</param>
        public List<ExternalDeconPeak> GetPeaksForScan(MsDataScan scan,
            ExternalDeconScanMatchMode matchMode = ExternalDeconScanMatchMode.OneBasedScanNumber,
            double retentionTimeTolerance = 0.01)
        {
            if (scan == null)
                return null;

            if (matchMode == ExternalDeconScanMatchMode.OneBasedScanNumber)
            {
                return PeaksByOneBasedScanNumber.TryGetValue(scan.OneBasedScanNumber, out var peaks) ? peaks : null;
            }

            int nearest = FindNearestScanNumberByRetentionTime(scan.RetentionTime, retentionTimeTolerance);
            return nearest > 0 && PeaksByOneBasedScanNumber.TryGetValue(nearest, out var rtPeaks) ? rtPeaks : null;
        }

        /// <summary>
        /// Binary searches the RT-sorted entries for the one closest to <paramref name="retentionTime"/>.
        /// Returns -1 when the closest entry is further away than <paramref name="tolerance"/>.
        /// </summary>
        private int FindNearestScanNumberByRetentionTime(double retentionTime, double tolerance)
        {
            if (_entriesByRetentionTime.Count == 0)
                return -1;

            int low = 0;
            int high = _entriesByRetentionTime.Count - 1;
            while (low < high)
            {
                int mid = (low + high) / 2;
                if (_entriesByRetentionTime[mid].RetentionTime < retentionTime)
                    low = mid + 1;
                else
                    high = mid;
            }

            // 'low' is the first entry at or after retentionTime; the nearest is either it or its predecessor
            int best = -1;
            double bestDistance = double.MaxValue;
            for (int i = Math.Max(0, low - 1); i <= Math.Min(_entriesByRetentionTime.Count - 1, low); i++)
            {
                double distance = Math.Abs(_entriesByRetentionTime[i].RetentionTime - retentionTime);
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    best = _entriesByRetentionTime[i].OneBasedScanNumber;
                }
            }

            return bestDistance <= tolerance ? best : -1;
        }
    }
}
