using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Collections.Concurrent;
using System.Linq;

namespace EngineLayer.DIA
{
    /// <summary>
    /// XicGroupingEngine is a specific implementation of PfGroupingEngine that groups precursor and fragment XICs into PrecursorFragmentsGroup objects 
    /// based on specified criteria including apex RT tolerance, overlap threshold, and correlation threshold.
    /// </summary>
    public class XicGroupingEngine : PfGroupingEngine
    {
        public float ApexRTTolerance { get; set; }
        public double OverlapThreshold { get; set; }
        public double CorrelationThreshold { get; set; }
        public int? PrecursorRankThreshold { get; set; } 
        public int? FragmentRankThreshold { get; set; }
        public ApexToleranceUnit ApexToleranceUnit { get; set; } //unit of ApexRTTolerance: minutes or number of scans

        /// <summary>
        /// A generic method of grouping correlated XICs into PrecursorFragmentsGroup objects based on specified criteria.
        /// </summary>
        /// <param name="apexRTTolerance"> maximum apex difference allowed between two XICs, in the unit given by apexToleranceUnit </param>
        /// <param name="overlapThreshold"> minimum overlap between two XICs </param>
        /// <param name="correlationThreshold"> minimum Pearson's correlation between two XICs </param>
        /// <param name="maxThreadsForGrouping"> the maximum threads for parallelization of grouping </param>
        /// <param name="minFragmentCountForGrouping"> minimun number of fragments needed to construct a PF group </param>
        /// <param name="precursorRankThreshold">  threshold for filtering out PF pairs ranked low for the precursor </param>
        /// <param name="fragmentRankThreshold"> threshold for filtering out PF pairs ranked low for the fragment </param>
        /// <param name="apexToleranceUnit"> whether apexRTTolerance is in minutes of retention time or in scans </param>
        public XicGroupingEngine(float apexRTTolerance, double overlapThreshold, double correlationThreshold, int maxThreadsForGrouping = 1, int minFragmentCountForGrouping = 0, int? precursorRankThreshold = null, int? fragmentRankThreshold = null,
            ApexToleranceUnit apexToleranceUnit = ApexToleranceUnit.RetentionTime)
        {
            ApexRTTolerance = apexRTTolerance;
            OverlapThreshold = overlapThreshold;
            CorrelationThreshold = correlationThreshold;
            MaxThreadsForGrouping = maxThreadsForGrouping;
            MinFragmentCountForPfGroup = minFragmentCountForGrouping;
            PrecursorRankThreshold = precursorRankThreshold;
            FragmentRankThreshold = fragmentRankThreshold;
            ApexToleranceUnit = apexToleranceUnit;
        }

        /// <summary>
        /// Given a list of precursor XICs and all eligible fragment XICs, loop over each precursor and group fragments with the precursor based on the grouping criteria.
        public override IEnumerable<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, List<ExtractedIonChromatogram> fragments)
        {
            var pfGroups = new ConcurrentBag<PrecursorFragmentsGroup>();
            var apexSortedFragmentXics = BuildApexSortedXics(fragments, ApexToleranceUnit);

            Parallel.ForEach(Partitioner.Create(0, precursors.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreadsForGrouping },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursors[i];
                        double precursorApex = ApexToleranceUnit == ApexToleranceUnit.ScanIndex ? precursor.ApexScanIndex : precursor.ApexRT;
                        var fragmentsInRange = GetXicsInRange(apexSortedFragmentXics, precursorApex, ApexRTTolerance, ApexToleranceUnit);
                        var pfGroup = GroupFragmentsForOnePrecursor(precursor, fragmentsInRange);
                        if (pfGroup != null)
                        {
                            lock (pfGroups)
                                pfGroups.Add(pfGroup);
                        }
                    }
                });

            var pfGroupsList = pfGroups.ToList();

            //Filter precursor-fragment pairs by rank if rank thresholds are set
            FilterPfPairsByRank(pfGroupsList, PrecursorRankThreshold, FragmentRankThreshold);

            //Remove groups with insufficient fragment pairs after filtering
            pfGroupsList.RemoveAll(g => g.PFpairs.Count < MinFragmentCountForPfGroup);

            return pfGroupsList;
        }

        /// <summary>
        /// Given one precursor XIC and all eligibile fragment XICs, select fragments that meet the grouping criteria and construct a precursor-fragment group for this precursor.
        /// </summary>
        public virtual PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, IEnumerable<ExtractedIonChromatogram> fragmentXics)
        {
            return GroupFragmentsForOnePrecursor(precursorXic, fragmentXics, ApexRTTolerance, OverlapThreshold, CorrelationThreshold, MinFragmentCountForPfGroup);
        }

        public static PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, IEnumerable<ExtractedIonChromatogram> fragmentXics, float apexRtTolerance, double overlapThreshold, double correlationThreshold, int minFragmentCountForGrouping)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragmentXic in fragmentXics)
            {
                double overlap = PrecursorFragmentsGroup.CalculateXicOverlapRatio(precursorXic, fragmentXic);
                if (overlap >= overlapThreshold)
                {
                    double correlation = PrecursorFragmentsGroup.CalculateXicCorrelation(precursorXic, fragmentXic);
                    if (correlation >= correlationThreshold)
                    {
                        var pfPair = new PrecursorFragmentPair(precursorXic, fragmentXic, correlation, overlap);
                        pfPairs.Add(pfPair);
                    }
                }
            }
            if (pfPairs.Count > minFragmentCountForGrouping)
            {
                var pfGroup = new PrecursorFragmentsGroup(precursorXic, pfPairs);
                return pfGroup;
            }
            return null;
        }

        /// <summary>
        /// The key a fragment XIC is sorted under: its apex retention time rounded to 0.01 min, or its apex scan index.
        /// </summary>
        protected static double GetApexKey(ExtractedIonChromatogram xic, ApexToleranceUnit apexToleranceUnit = ApexToleranceUnit.RetentionTime)
        {
            return apexToleranceUnit == ApexToleranceUnit.ScanIndex ? xic.ApexScanIndex : Math.Round(xic.ApexRT, 2);
        }

        protected static SortedDictionary<double, List<ExtractedIonChromatogram>> BuildApexSortedXics(IEnumerable<ExtractedIonChromatogram> xics, ApexToleranceUnit apexToleranceUnit = ApexToleranceUnit.RetentionTime)
        {
            var tree = new SortedDictionary<double, List<ExtractedIonChromatogram>>();
            foreach (var xic in xics)
            {
                double apexKey = GetApexKey(xic, apexToleranceUnit);
                if (!tree.ContainsKey(apexKey))
                    tree[apexKey] = new List<ExtractedIonChromatogram>();
                tree[apexKey].Add(xic);
            }
            return tree;
        }

        protected static List<ExtractedIonChromatogram> GetXicsInRange(SortedDictionary<double, List<ExtractedIonChromatogram>> apexSortedFragmentXics, double targetApex, double tolerance,
            ApexToleranceUnit apexToleranceUnit = ApexToleranceUnit.RetentionTime)
        {
            double min = targetApex - tolerance;
            double max = targetApex + tolerance;
            if (apexToleranceUnit == ApexToleranceUnit.RetentionTime)
            {
                min = Math.Round(min, 2);
                max = Math.Round(max, 2);
            }

            var subset = apexSortedFragmentXics.Where(kvp => kvp.Key >= min && kvp.Key <= max).SelectMany(kvp => kvp.Value).ToList();
            return subset;
        }
    }
}
