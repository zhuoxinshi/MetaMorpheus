using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;
using static Proteomics.RetentionTimePrediction.SSRCalc3;

namespace EngineLayer.DIA
{
    public class ISDEngine : DIAEngine
    {
        private readonly MsDataFile DataFile;
        public ISDEngine(DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) : base(DIAparameters, dataFile, commonParameters, fileSpecificParameters, nestedIds)
        {
            DIAparams = DIAparameters;
            DataFile = dataFile;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            PseudoMs2Scans = GetPseudoMs2Scans();

            if (DIAparams.CombineFragments)
            {
                var combinedScans = new List<Ms2ScanWithSpecificMass>();
                var groupedPseudoScans = PseudoMs2Scans.GroupBy(s => new { s.PrecursorMass, s.PrecursorCharge, s.RetentionTime});
                foreach(var group in groupedPseudoScans)
                {
                    var combinedFragments = group.SelectMany(s => s.ExperimentalFragments);
                    var newScan = new Ms2ScanWithSpecificMass(group.FirstOrDefault().TheScan, group.FirstOrDefault().PrecursorMonoisotopicPeakMz, group.Key.PrecursorCharge, group.FirstOrDefault().FullFilePath, CommonParameters, combinedFragments.ToArray(), group.FirstOrDefault().PrecursorIntensity);
                    combinedScans.Add(newScan);
                }
                PseudoMs2Scans = combinedScans;
            }
            return new MetaMorpheusEngineResults(this);
        }

        public override IEnumerable<Ms2ScanWithSpecificMass> GetPseudoMs2Scans()
        {
            //read in scans and isd scan pre-process
            var allScans = DataFile.GetAllScansList().ToArray();
            var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
            ReLabelIsdScans(isdVoltageMap, allScans);

            //Get all MS1 and MS2 XICs
            var allMs1Xics = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks, out var indexingEngine);

            foreach (var ms2Group in isdVoltageMap)
            {
                var ms2Xics = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(ms2Group.Value.ToArray(), out matchedPeaks, out indexingEngine);
                var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, ms2Xics);

                foreach (var pfGroup in pfGroups)
                {
                    OneBasedScanNumber++;
                    pfGroup.PFgroupIndex = OneBasedScanNumber;
                    var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, DIAparams.PseudoMs2ConstructionType, CommonParameters, DataFile.FilePath);
                    yield return pseudoScan;
                }
            }
        }


        public static void ReLabelIsdScans(Dictionary<double, List<MsDataScan>> isdVoltageScanMap, MsDataScan[] ms1Scans)
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

        public static Dictionary<double, List<MsDataScan>> ConstructIsdGroups(MsDataScan[] scans, out MsDataScan[] ms1Scans)
        {
            var isdVoltageScanMap = new Dictionary<double, List<MsDataScan>>();
            string pattern = $@"sid=(\d+)";
            foreach (var scan in scans)
            {
                double voltage = 0;
                var match = Regex.Match(scan.ScanFilter, pattern);
                if (match.Success) voltage = double.Parse(match.Groups[1].Value);
                if (!isdVoltageScanMap.ContainsKey(voltage))
                {
                    isdVoltageScanMap[voltage] = new List<MsDataScan> { scan };
                }
                else
                {
                    isdVoltageScanMap[voltage].Add(scan);
                }
            }
            var soretedMap = isdVoltageScanMap.OrderBy(kvp => kvp.Key).ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
            ms1Scans = soretedMap.First().Value.ToArray();
            soretedMap.Remove(soretedMap.First().Key);
            return soretedMap;
        }

        public void GetMS2OnlyPseudoScans()
        {
            var allScans = DataFile.GetAllScansList().ToArray();
            var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
            ReLabelIsdScans(isdVoltageMap, allScans);
            var ms2Scans = isdVoltageMap.Values.SelectMany(p => p).ToArray();

            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            foreach(var scan in ms2Scans)
            {
                var neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(scan, CommonParameters);

                var pseudoScan = new Ms2ScanWithSpecificMass(scan, 1,
                                   1, null, CommonParameters, neutralExperimentalFragments, precursorHighestIsotopeMz: 1);
                pseudoScans.Add(pseudoScan);
            }
            PseudoMs2Scans = pseudoScans;
        }
    }
}
