using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Omics.Modifications;
using System.Collections.Concurrent;
using EngineLayer.Util;
using Omics;
using FlashLFQ;
using System.IO;

namespace EngineLayer.DIA
{
    /// <summary>
    /// DIAEgine defines a workflow of generating DDA-like pseudo MS2 scans for DIA data analysis. It includes the processes of extracting precursor and 
    /// fragment XICs, grouping them into PrecursorFragmentsGroup objects, and constructing pseudo MS2 scans.
    /// <summary>
    public class DIAEngine : MetaMorpheusEngine
    {
        private readonly MsDataFile DataFile;
        public DIAparameters DIAparams { get; set; } 
        public IEnumerable<Ms2ScanWithSpecificMass> PseudoMs2Scans { get; set; }
        protected override MetaMorpheusEngineResults RunSpecific()
        {
            PseudoMs2Scans = GetPseudoMs2Scans();
            return new MetaMorpheusEngineResults(this);
        }

        public virtual IEnumerable<Ms2ScanWithSpecificMass> GetPseudoMs2Scans()
        {
            //read in scans
            var ms1Scans = DataFile.GetMS1Scans().ToArray();
            var ms2Scans = DataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var DIAScanWindowMap = ConstructMs2Groups(ms2Scans);
            int oneBasedScanNumber = 0;

            foreach (var kvp in DIAScanWindowMap)
            {
                var ms1Xics = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans, out var matchedPeaks, out var indexingEngine, new MzRange(kvp.Key.min, kvp.Key.max));
                var ms2Xics = DIAparams.Ms2XicConstructor.GetAllXicsWithXicSpline(kvp.Value.ToArray(), out matchedPeaks, out indexingEngine);
                var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);

                foreach (var pfGroup in pfGroups)
                {
                    pfGroup.PFgroupIndex = oneBasedScanNumber;
                    oneBasedScanNumber++;
                    var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, DIAparams.PseudoMs2ConstructionType, CommonParameters, DataFile.FilePath);
                    yield return pseudoScan;
                }
            }
        }

        public DIAEngine(DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) :base(commonParameters, fileSpecificParameters, nestedIds)
        {
            DIAparams = DIAparameters;
            DataFile = dataFile;
        }

        public static Dictionary<(double min, double max), List<MsDataScan>> ConstructMs2Groups(MsDataScan[] ms2Scans)
        {
            var DIAScanWindowMap = new Dictionary<(double min, double max), List<MsDataScan>>();
            foreach (var ms2 in ms2Scans)
            {
                (double min, double max) range = new(Math.Round(ms2.IsolationRange.Minimum, 2), Math.Round(ms2.IsolationRange.Maximum, 2));
                if (!DIAScanWindowMap.ContainsKey(range))
                {
                    DIAScanWindowMap[range] = new List<MsDataScan> { ms2 };
                }
                else
                {
                    DIAScanWindowMap[range].Add(ms2);
                }
            }
            return DIAScanWindowMap;
        }
    }
}
