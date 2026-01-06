using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
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

        //protected override MetaMorpheusEngineResults RunSpecific()
        //{
        //    //read in scans and isd scan pre-process
        //    var allScans = DataFile.GetAllScansList().ToArray();
        //    var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
        //    ReLabelIsdScans(isdVoltageMap, allScans);

        //    //Get all MS1 and MS2 XICs
        //    var allMs1Xics = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks, out var indexingEngine);
        //    var allMs2Xics = new Dictionary<double, List<ExtractedIonChromatogram>>();
        //    foreach (var ms2Group in isdVoltageMap)
        //    {
        //        allMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(ms2Group.Value.ToArray(), out matchedPeaks, out indexingEngine);
        //    }

        //    //Precursor-fragment Grouping
        //    var allPfGroups = new List<PrecursorFragmentsGroup>();
        //    if (DIAparams.CombineFragments)
        //    {
        //        allPfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics.Values.SelectMany(p => p));
        //    }
        //    else
        //    {
        //        foreach (var ms2Group in isdVoltageMap.Keys)
        //        {
        //            var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics[ms2Group]);
        //            allPfGroups.AddRange(pfGroups);
        //        }
        //    }

        //    //filtering fragments
        //    foreach (var pfGroup in allPfGroups)
        //    {
        //        pfGroup.PFpairs = pfGroup.PFpairs.OrderByDescending(pf => pf.FragmentXic.ApexPeak.Intensity).Take(200).ToList();
        //    }

        //    //Convert pfGroups to pseudo MS2 scans
        //    var pseudoScans = new List<Ms2ScanWithSpecificMass>();
        //    int pfGroupIndex = 1;
        //    foreach (var pfGroup in allPfGroups)
        //    {
        //        pfGroup.PFgroupIndex = pfGroupIndex;
        //        var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, DIAparams.PseudoMs2ConstructionType, CommonParameters, DataFile.FilePath);
        //        pseudoScans.Add(pseudoScan);
        //        pfGroupIndex++;
        //    }
        //    PseudoMs2Scans = pseudoScans;

        //    return new MetaMorpheusEngineResults(this);
        //}

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //read in scans
            var allScans = DataFile.GetAllScansList().ToArray();
            var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);

            string flashMs1FilePath = @"E:\ISD Project\FW-DIA\TestMyData\Yeast\YD_preFilter\FlashDeconvResults\09-10-25_YD_81min_ISD60-80-100_preFilter700-900-1100_rep1_ms1_ms1.tsv";
            string flash60FilePath = @"E:\ISD Project\FW-DIA\TestMyData\Yeast\YD_preFilter\FlashDeconvResults\09-10-25_YD_81min_ISD60-80-100_preFilter700-900-1100_rep1_isd60_ms1.tsv";
            string flash80FilePath = @"E:\ISD Project\FW-DIA\TestMyData\Yeast\YD_preFilter\FlashDeconvResults\09-10-25_YD_81min_ISD60-80-100_preFilter700-900-1100_rep1_isd80_ms1.tsv";
            string flash100FilePath = @"E:\ISD Project\FW-DIA\TestMyData\Yeast\YD_preFilter\FlashDeconvResults\09-10-25_YD_81min_ISD60-80-100_preFilter700-900-1100_rep1_isd100_ms1.tsv";

            var ms1XicConstructor = new ChargeEnvelopeXicConstructor(flashMs1FilePath, new PpmToleranceWithNotch(10, 2, 2), 2, 0.5, 3, DIAparams.Ms1XicConstructor.XicSplineEngine);
            var allMs1Xics = ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks, out var indexingEngine);
            var allMs2Xics = new Dictionary<double, List<ExtractedIonChromatogram>>();
            var isdFiles = new Dictionary<double, string> { { 60, flash60FilePath }, { 80, flash80FilePath }, { 100, flash100FilePath } };
            foreach (var kvp in isdFiles)
            {
                var ms2XicConstructor = new ChargeEnvelopeXicConstructor(kvp.Value, new PpmToleranceWithNotch(20, 1, 1), 2, 0.5, 3, DIAparams.Ms2XicConstructor.XicSplineEngine);
                allMs2Xics[kvp.Key] = ms2XicConstructor.GetAllXicsWithXicSpline(isdVoltageMap[kvp.Key].ToArray(), out matchedPeaks, out indexingEngine);
            }

            //Precursor-fragment Grouping
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            if (DIAparams.CombineFragments)
            {
                allPfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics.Values.SelectMany(p => p));
            }
            else
            {
                foreach (var ms2Group in isdFiles.Keys)
                {
                    var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics[ms2Group]);
                    allPfGroups.AddRange(pfGroups);
                }
            }

            //filtering fragments
            foreach (var pfGroup in allPfGroups)
            {
                pfGroup.PFpairs = pfGroup.PFpairs.OrderByDescending(pf => pf.FragmentXic.ApexPeak.Intensity).Take(100).ToList();
            }

            //Convert pfGroups to pseudo MS2 scans
            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            int pfGroupIndex = 1;
            foreach (var pfGroup in allPfGroups)
            {
                pfGroup.PFgroupIndex = pfGroupIndex;
                var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, DIAparams.PseudoMs2ConstructionType, CommonParameters, DataFile.FilePath);
                pseudoScans.Add(pseudoScan);
                pfGroupIndex++;
            }
            PseudoMs2Scans = pseudoScans;
            return new MetaMorpheusEngineResults(this);
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
    }
}
