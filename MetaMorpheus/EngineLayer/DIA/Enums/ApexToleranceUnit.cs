namespace EngineLayer.DIA
{
    /// <summary>
    /// Unit of <see cref="XicGroupingEngine.ApexRTTolerance"/>, the largest apex difference allowed between a
    /// precursor XIC and a fragment XIC.
    /// </summary>
    public enum ApexToleranceUnit
    {
        /// <summary>
        /// Retention time in minutes; apexes are compared after rounding to 0.01 min.
        /// </summary>
        RetentionTime,

        /// <summary>
        /// Number of scans, compared on <see cref="MassSpectrometry.ExtractedIonChromatogram.ApexScanIndex"/>,
        /// the apex position in the scan array the XIC was built from. In ISD data every acquisition cycle holds
        /// one scan per source voltage, so the precursor and fragment scan arrays are aligned cycle by cycle and
        /// a tolerance of N allows the two apexes to be up to N cycles apart.
        /// </summary>
        ScanIndex,
    }
}
