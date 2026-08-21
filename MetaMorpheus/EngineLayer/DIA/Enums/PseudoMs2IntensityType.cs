namespace EngineLayer.DIA
{
    /// <summary>
    /// Which number from an XIC becomes the intensity written into a pseudo MS2 scan - for the fragment peaks
    /// and for the precursor alike.
    ///
    /// A pseudo MS2 scan is not a scan: it summarises a whole coeluting group, so there is a real choice about
    /// what "intensity" should mean. The two options are not a rescaling of one another - on TopFD features the
    /// apex is anywhere from 2% to 17% of the total, because the ratio tracks peak width - so switching between
    /// them reorders the peaks within a spectrum, it does not just change the units.
    /// </summary>
    public enum PseudoMs2IntensityType
    {
        /// <summary>
        /// Total signal over the XIC's whole elution: every point summed. For XICs built by
        /// <see cref="FeatureXmlXicConstructor"/> each point is already a summed isotopic envelope, so this is
        /// the summed envelope intensity of the feature over the scans that were handed in.
        ///
        /// That is CLOSE TO but not identical with the feature's own <c>&lt;intensity&gt;</c> in the feature xml
        /// (the element directly below <c>&lt;mono_mass&gt;</c>): TopFD computes that as the sum of the
        /// per-charge <c>&lt;intensity&gt;</c> values, which are not themselves the raw sums of their envelopes.
        /// On well-sampled features the two agree to a fraction of a percent, but on sparse ones - a handful of
        /// envelopes spread over a few scans - this sum can run ~20% high. Use it as the same KIND of quantity
        /// as the feature xml's intensity, not as the identical number.
        ///
        /// Summing the XIC rather than reading the feature attribute is deliberate: it respects the scan array
        /// actually handed in, so on an ISD channel it counts only that channel's envelopes, and it is the same
        /// signal the grouping correlation was computed on.
        /// </summary>
        XicSum,

        /// <summary>
        /// Intensity at the apex scan only - the analogue of what one real MS2 scan would show, since each XIC
        /// point is already a summed isotopic envelope. Closer to the semantics of a genuine msAlign peak, at the
        /// cost of resting on a single scan and of under-weighting broad, well-measured fragments.
        /// </summary>
        ApexPeak,
    }
}
