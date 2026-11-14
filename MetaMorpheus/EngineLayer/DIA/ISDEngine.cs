using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace EngineLayer.DIA
{
    /// <summary>
    /// ISDEngine inherits from DIAEngine and employs a similar workflow of generating XICs, grouping XICs and make pseudo Ms2 scans for ISD data analysis. 
    /// </summary>
    public class ISDEngine : DIAEngine
    {
        private readonly MsDataFile DataFile;

        public ISDEngine(MsDataFile dataFile, CommonParameters commonParameters) : base(dataFile, commonParameters)
        {
            DataFile = dataFile;
        }
       
        public override IEnumerable<Ms2ScanWithSpecificMass> GetPseudoMs2Scans()
        {
            //read in scans and isd scan pre-process
            var allScans = DataFile.GetAllScansList().ToArray();
            var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
            ReLabelIsdScans(isdVoltageMap, allScans);

            //Get all MS1 XICs
            var allMs1Xics = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans);

            //When the precursor XICs come from TopFD features, remember each XIC's feature so the pseudo MS2 scans can
            //name it. Copied here because the constructor replaces the map on every call.
            var featureByPrecursorXic = DIAparams.Ms1XicConstructor is FeatureXmlXicConstructor ms1FeatureConstructor
                ? new Dictionary<ExtractedIonChromatogram, TopFdFeature>(ms1FeatureConstructor.LastFeatureByXic)
                : null;

            //Precursor-fragment Grouping
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            if (DIAparams.CombineFragments)
            {
                var allMs2Xics = new Dictionary<int, List<ExtractedIonChromatogram>>();
                foreach (var kvp in isdVoltageMap)
                {
                    var ms2Xics = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(kvp.Value.ToArray());
                    allMs2Xics.Add(kvp.Key, ms2Xics);
                }
                allPfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics.Values.SelectMany(p => p).ToList()).ToList();
            }
            else
            {
                foreach(var kvp in isdVoltageMap)
                {
                    var ms2Xics = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(kvp.Value.ToArray());
                    var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, ms2Xics);
                    allPfGroups.AddRange(pfGroups);
                }
            }
            PrecursorFragmentsGroup.AssignPrecursorFeatureIds(allPfGroups, featureByPrecursorXic);

            //The grouping engine runs in parallel and returns groups in no fixed order. Number the pseudo MS2 scans in the
            //order of the precursor XICs (a stable sort, so groups of one precursor keep their voltage order) so that the
            //scan numbers do not depend on thread scheduling.
            var precursorXicOrder = new Dictionary<ExtractedIonChromatogram, int>(allMs1Xics.Count);
            for (int i = 0; i < allMs1Xics.Count; i++)
                precursorXicOrder[allMs1Xics[i]] = i;
            allPfGroups = allPfGroups.OrderBy(g => precursorXicOrder[g.PrecursorXic]).ToList();

            foreach (var pfGroup in allPfGroups)
            {
                OneBasedScanNumber++;
                pfGroup.PFgroupIndex = OneBasedScanNumber;
                var pseudoScan = pfGroup.GetPseudoMs2ScanFromPfGroup(DIAparams.PseudoMs2ConstructionType, CommonParams, DataFile.FilePath,
                    DIAparams.PseudoMs2IntensityType);
                yield return pseudoScan;
            }
        }

        /// <summary>
        /// Relabel ISD scans to have correct Ms level, precursor scan number and isolation m/z/range.
        /// </summary>
        public static void ReLabelIsdScans(Dictionary<int, List<MsDataScan>> isdVoltageScanMap, MsDataScan[] ms1Scans)
        {
            for (int i = 0; i < isdVoltageScanMap.Count; i++)
            {
                foreach (var scan in isdVoltageScanMap.Values.ElementAt(i))
                {
                    scan.SetMsnOrder(2);
                    int oneBasedPrecursorScanNumber = scan.OneBasedScanNumber - i - 1;
                    var ms1Scan = ms1Scans.Where(s => s.OneBasedScanNumber == oneBasedPrecursorScanNumber).First();
                    scan.SetOneBasedPrecursorScanNumber(oneBasedPrecursorScanNumber);
                    scan.SetIsolationRange(ms1Scan.ScanWindowRange.Minimum, ms1Scan.ScanWindowRange.Maximum);
                    var scanWindowWidth = ms1Scan.ScanWindowRange.Maximum - ms1Scan.ScanWindowRange.Minimum;
                    scan.SetIsolationMz(ms1Scan.ScanWindowRange.Minimum + scanWindowWidth / 2);
                }
            }
        }

        /// <summary>
        /// Divide all ISD scans into groups based on their source energy level and extract real MS1 scans. 
        /// </summary>
        /// <returns> A dictionary where the key is ISD voltage and the values are MS2 scans. </returns>
        public static Dictionary<int, List<MsDataScan>> ConstructIsdGroups(MsDataScan[] scans, out MsDataScan[] ms1Scans)
        {
            var isdVoltageScanMap = new Dictionary<int, List<MsDataScan>>();
            //ISD voltage is written in the scan filter like "sid=15" in the raw file
            string pattern = $@"sid=(\d+)";
            foreach (var scan in scans)
            {
                int voltage = 0;
                var match = Regex.Match(scan.ScanFilter, pattern);
                if (match.Success) voltage = int.Parse(match.Groups[1].Value);
                if (!isdVoltageScanMap.ContainsKey(voltage))
                {
                    isdVoltageScanMap[voltage] = new List<MsDataScan> { scan };
                }
                else
                {
                    isdVoltageScanMap[voltage].Add(scan);
                }
            }
            var sortedMap = isdVoltageScanMap.OrderBy(kvp => kvp.Key).ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
            //ISD runs start with MS1 scans so the first kvp is removed from the dictionary
            ms1Scans = sortedMap.First().Value.ToArray();
            sortedMap.Remove(sortedMap.First().Key);
            return sortedMap;
        }

        /// <summary>
        /// Writes pseudo MS2 scans as a TopFD-style msAlign file (one BEGIN IONS ... END IONS block per scan, peaks as
        /// monoisotopic mass, intensity and charge), ordered by scan number. The scans must carry neutral fragments,
        /// i.e. come from <see cref="PseudoMs2ConstructionType.Mass"/>. The precursor feature id is left empty when
        /// unknown, as TopFD writes it.
        /// </summary>
        public static void WriteMsAlignFile(string filePath, IEnumerable<Ms2ScanWithSpecificMass> ms2ScansWithMass)
        {
            var culture = CultureInfo.InvariantCulture;
            using var writer = new StreamWriter(filePath);
            foreach (var scan in ms2ScansWithMass.OrderBy(s => s.OneBasedScanNumber))
            {
                writer.WriteLine("BEGIN IONS");
                writer.WriteLine(string.Create(culture, $"ID={scan.OneBasedScanNumber}"));
                writer.WriteLine("FRACTION_ID=0");
                writer.WriteLine($"FILE_NAME={scan.FullFilePath}");
                writer.WriteLine(string.Create(culture, $"SPECTRUM_ID={scan.OneBasedScanNumber}"));
                writer.WriteLine(string.Create(culture, $"TITLE=Scan_{scan.OneBasedScanNumber}"));
                writer.WriteLine(string.Create(culture, $"SCANS={scan.OneBasedScanNumber}"));
                writer.WriteLine(string.Create(culture, $"RETENTION_TIME={Math.Round(scan.RetentionTime, 2)}"));
                writer.WriteLine("LEVEL=2");
                writer.WriteLine(string.Create(culture, $"MS_ONE_ID={scan.OneBasedScanNumber}"));
                writer.WriteLine(string.Create(culture, $"MS_ONE_SCAN={scan.OneBasedScanNumber}"));
                writer.WriteLine(string.Create(culture, $"PRECURSOR_WINDOW_BEGIN={scan.PrecursorMonoisotopicPeakMz}"));
                writer.WriteLine(string.Create(culture, $"PRECURSOR_WINDOW_END={scan.PrecursorMonoisotopicPeakMz + 3}"));
                writer.WriteLine("ACTIVATION=HCD");
                writer.WriteLine(string.Create(culture, $"PRECURSOR_MZ={scan.PrecursorMonoisotopicPeakMz}"));
                writer.WriteLine(string.Create(culture, $"PRECURSOR_CHARGE={scan.PrecursorCharge}"));
                writer.WriteLine(string.Create(culture, $"PRECURSOR_MASS={scan.PrecursorMass}"));
                writer.WriteLine(string.Create(culture, $"PRECURSOR_INTENSITY={scan.PrecursorIntensity}"));
                writer.WriteLine(string.Create(culture, $"PRECURSOR_FEATURE_ID={scan.PrecursorFeatureId}"));

                foreach (var fragment in scan.ExperimentalFragments)
                {
                    writer.WriteLine(string.Create(culture, $"{fragment.MonoisotopicMass}\t{fragment.TotalIntensity}\t{fragment.Charge}"));
                }

                writer.WriteLine("END IONS");
                writer.WriteLine();
            }
        }
    }
}
