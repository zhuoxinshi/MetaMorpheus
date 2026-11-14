namespace EngineLayer.DIA
{
    /// <summary>
    /// Which value of an XIC becomes the intensity of a pseudo MS2 peak, and of the pseudo MS2 precursor.
    /// The two options are not a rescaling of each other: the apex-to-total ratio depends on the peak width,
    /// so switching between them can reorder the peaks within a pseudo MS2 spectrum.
    /// </summary>
    public enum PseudoMs2IntensityType
    {
        /// <summary>
        /// Intensity of the XIC's apex point, the analogue of what one real MS2 scan would show.
        /// </summary>
        ApexPeak,

        /// <summary>
        /// Sum of the intensities of all XIC points, i.e. the total signal over the elution. For XICs read
        /// from TopFD features (<see cref="FeatureXmlXicConstructor"/>) each point is a summed isotopic
        /// envelope, so this is the feature's summed envelope intensity over the scans the XIC was built from.
        /// </summary>
        XicSum,
    }
}
