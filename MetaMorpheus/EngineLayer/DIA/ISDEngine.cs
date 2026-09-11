using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
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

            //Remember which deconvoluted MS1 feature each precursor XIC came from, so the pseudo MS2 scans can
            //name it. The map has to be snapshotted HERE: the constructor overwrites it on every call, and the
            //fragment channels are built below. Only a feature-file constructor knows feature ids; anything else
            //leaves the map null and the ids are written as unknown.
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

            foreach (var pfGroup in allPfGroups)
            {
                OneBasedScanNumber++;
                pfGroup.PFgroupIndex = OneBasedScanNumber;
                PrecursorFragmentsGroup.AssignPrecursorFeatureIds(new[] { pfGroup }, featureByPrecursorXic);
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
                if (scan.IsDefaultOrNull() || string.IsNullOrEmpty(scan.ScanFilter)) 
                        continue;
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

        public static void WriteMsAlignFile(string filePath, IEnumerable<Ms2ScanWithSpecificMass> ms2ScansWithMass, string envCnnMsAlignPath = null)
        {
            // When an external-decon msAlign (e.g. TopFD *_ms2.msalign) is supplied, look up each fragment's
            // EnvCNN score from it and write it as a 4th peak column (mass  intensity  charge  EnvCNN).
            var apexEnvCnn = string.IsNullOrEmpty(envCnnMsAlignPath) ? null : ReadApexEnvCnnScores(envCnnMsAlignPath);
            var sortedScans = ms2ScansWithMass.OrderBy(s => s.OneBasedScanNumber);
            using (var writer = new StreamWriter(filePath))
            {
                foreach (var scan in sortedScans)
                {
                    // Write scan entry header
                    writer.WriteLine("BEGIN IONS");
                    writer.WriteLine($"ID={scan.OneBasedScanNumber}");
                    writer.WriteLine("FRACTION_ID=0");
                    writer.WriteLine($"FILE_NAME={scan.FullFilePath}");
                    writer.WriteLine($"SPECTRUM_ID={scan.OneBasedScanNumber}");
                    writer.WriteLine($"TITLE=Scan_{scan.OneBasedScanNumber}");
                    writer.WriteLine($"SCANS={scan.OneBasedScanNumber}");
                    writer.WriteLine($"RETENTION_TIME={Math.Round(scan.RetentionTime, 2)}");
                    writer.WriteLine($"LEVEL=2");
                    writer.WriteLine($"MS_ONE_ID={scan.OneBasedScanNumber}");
                    writer.WriteLine($"MS_ONE_SCAN={scan.OneBasedScanNumber}");
                    writer.WriteLine($"PRECURSOR_WINDOW_BEGIN={scan.PrecursorMonoisotopicPeakMz}");
                    writer.WriteLine($"PRECURSOR_WINDOW_END={scan.PrecursorMonoisotopicPeakMz + 3}");
                    writer.WriteLine($"ACTIVATION=HCD");
                    writer.WriteLine($"PRECURSOR_MZ={scan.PrecursorMonoisotopicPeakMz}");
                    writer.WriteLine($"PRECURSOR_CHARGE={scan.PrecursorCharge}");
                    writer.WriteLine($"PRECURSOR_MASS={scan.PrecursorMass}");
                    writer.WriteLine($"PRECURSOR_INTENSITY={scan.PrecursorIntensity}");
                    //TopFD writes this EMPTY when no feature backs the precursor, and 0 is a real feature id, so
                    //an unknown id must not be written as 0 - that is what made every exported scan claim feature 0
                    writer.WriteLine($"PRECURSOR_FEATURE_ID={scan.PrecursorFeatureId?.ToString() ?? string.Empty}");

                    // Write peaks: monoMass, intensity, charge [, EnvCNN score at the fragment's apex]
                    for (int i = 0; i < scan.ExperimentalFragments.Length; i++)
                    {
                        var frag = scan.ExperimentalFragments[i];
                        if (apexEnvCnn == null)
                        {
                            writer.WriteLine($"{frag.MonoisotopicMass}\t{frag.TotalIntensity}\t{frag.Charge}");
                        }
                        else
                        {
                            double env = LookupApexEnvCnn(apexEnvCnn, frag.Charge, frag.MonoisotopicMass);
                            writer.WriteLine($"{frag.MonoisotopicMass}\t{frag.TotalIntensity}\t{frag.Charge}\t{(double.IsNaN(env) ? string.Empty : env.ToString())}");
                        }
                    }

                    writer.WriteLine("END IONS");
                    writer.WriteLine();
                }
            }
        }

        // For each (charge, 0.01-Da mass bin), keep the EnvCNN (4th msAlign column) of the most intense peak -
        // the peak a fragment XIC apexes on. Peak lines are the tab-delimited numeric lines (no '=') in BEGIN/END IONS.
        private static Dictionary<(int, long), double> ReadApexEnvCnnScores(string ms2MsAlignPath)
        {
            var bestIntensity = new Dictionary<(int, long), double>();
            var scoreAtBest = new Dictionary<(int, long), double>();
            foreach (var raw in File.ReadLines(ms2MsAlignPath))
            {
                if (raw.Length == 0 || raw[0] == '#' || raw.Contains("=") || raw.StartsWith("BEGIN") || raw.StartsWith("END"))
                    continue;
                var s = raw.Split('\t');
                if (s.Length < 4) continue;
                if (!double.TryParse(s[0], out double mass) || !double.TryParse(s[1], out double intensity)
                    || !int.TryParse(s[2], out int charge) || !double.TryParse(s[3], out double env))
                    continue;
                var key = (charge, (long)Math.Round(mass * 100.0));
                if (!bestIntensity.TryGetValue(key, out double bi) || intensity > bi)
                {
                    bestIntensity[key] = intensity;
                    scoreAtBest[key] = env;
                }
            }
            return scoreAtBest;
        }

        private static double LookupApexEnvCnn(Dictionary<(int, long), double> scores, int charge, double mass)
        {
            long bin = (long)Math.Round(mass * 100.0);
            for (long d = -1; d <= 1; d++) // tolerate the 0.01-Da bin boundary (XIC apex mass is stored as float)
                if (scores.TryGetValue((charge, bin + d), out double v))
                    return v;
            return double.NaN;
        }
    }
}
