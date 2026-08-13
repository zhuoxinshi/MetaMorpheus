using System.Collections.Generic;

namespace EngineLayer.DIA.FeatureAlign
{
    public enum MassToleranceUnit { Ppm, Dalton }

    /// <summary>How features are grouped across runs.</summary>
    public enum GroupingMethod
    {
        /// <summary>Transitive connected components — features are linked to each other (the original behaviour).
        /// Can chain a ladder of near-mass features into one large cluster.</summary>
        UnionFind,
        /// <summary>OpenMS ConsensusFeatureGroupDetector style: greedy, each group is anchored to the most-abundant
        /// remaining feature and candidates are matched to that anchor only. The group can't be wider than the
        /// tolerance (no runaway chaining); one feature per run is kept via <see cref="DuplicateResolution"/>.</summary>
        ReferenceAnchor
    }

    /// <summary>Which feature to keep when one run has several candidates for the same group (ReferenceAnchor).
    /// Mirrors OpenMS's -when_duplicate.</summary>
    public enum DuplicateResolution { MaxAbundance, NearestMass }

    /// <summary>
    /// Options for <see cref="FeatureAligner"/>: which feature tables to match across runs, the mass and
    /// retention-time tolerances used to decide two features are the same proteoform, which quantity column
    /// to read, and where to write the output tsv. Works on FLASHQuant or FLASHDeconv feature tables.
    /// </summary>
    public class FeatureAlignParameters
    {
        /// <summary>Full paths to the input feature tables (.tsv). At least two are required.</summary>
        public List<string> InputFiles { get; set; } = new();

        /// <summary>Mass tolerance value; the unit is given by <see cref="MassUnit"/>. Applied as +/- (symmetric).</summary>
        public double MassTolerance { get; set; } = 20;

        public MassToleranceUnit MassUnit { get; set; } = MassToleranceUnit.Ppm;

        /// <summary>Retention-time tolerance in seconds, applied as +/- on the apex retention time.</summary>
        public double RtTolerance { get; set; } = 150;

        /// <summary>Full path of the tsv to write.</summary>
        public string OutputPath { get; set; }

        /// <summary>Preferred quantity column; if absent in a file, common alternatives are tried automatically.</summary>
        public string PreferredQuantColumn { get; set; } = "FeatureGroupQuantity";

        /// <summary>Number of +1 Da (and -1 Da) monoisotope "notches" allowed in the mass match, applied around the
        /// ppm tolerance (like MzLibUtil.PpmToleranceWithNotch). 0/0 = exact monoisotope only (recommended for
        /// same-tool cross-run alignment). Only used when <see cref="MassUnit"/> is Ppm.</summary>
        public int PositiveNotch { get; set; } = 0;
        public int NegativeNotch { get; set; } = 0;

        /// <summary>Grouping strategy across runs (default UnionFind = original behaviour).</summary>
        public GroupingMethod Grouping { get; set; } = GroupingMethod.UnionFind;

        /// <summary>When one run has several candidates for a group, which to keep (ReferenceAnchor only).</summary>
        public DuplicateResolution WhenDuplicate { get; set; } = DuplicateResolution.MaxAbundance;

        public FeatureAlignParameters() { }

        public FeatureAlignParameters(IEnumerable<string> inputFiles, string outputPath,
            double massTolerance = 20, MassToleranceUnit massUnit = MassToleranceUnit.Ppm,
            double rtTolerance = 150, string preferredQuantColumn = "FeatureGroupQuantity")
        {
            InputFiles = new List<string>(inputFiles);
            OutputPath = outputPath;
            MassTolerance = massTolerance;
            MassUnit = massUnit;
            RtTolerance = rtTolerance;
            PreferredQuantColumn = preferredQuantColumn;
        }
    }
}
