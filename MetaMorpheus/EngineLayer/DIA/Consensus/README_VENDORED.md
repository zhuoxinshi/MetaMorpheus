# Vendored: mzLib consensus mass tracer

These 7 files are a **verbatim copy** of the consensus mass tracer from the mzLib source worktree
(`mzLib/MassSpectrometry/Deconvolution/Consensus/`, branch `isd-work`, commit `56a5d884ff78`), namespace
`MassSpectrometry.Deconvolution.Consensus`:

`MassTrace`, `MassTraceBuilder`, `CorrectedEnvelope`, `CorrectedTrace`, `TraceCorrector`, `MassFeature`,
`MassFeatureBuilder`.

## Why vendored
The consensus tracer (mzLib PR 1069) is NOT in the NuGet mzLib version this repo references (`mzLib 1.0.579`);
it lives only in the `isd-work` mzLib branch (based on 1.0.574). `deconscan` uses it via a direct source
reference, but MetaMorpheus consumes mzLib as a package. Rather than force a solution-wide mzLib version bump,
these files are vendored here so `ConsensusMassXicConstructor` can build **consensus feature tracing +
precursor-fragment grouping + pseudo-MS2** end to end. The tracer depends only on `IsotopicEnvelope`,
`MsDataScan`, and `Chemistry`, all present in 1.0.579, so it compiles unchanged.

## TODO (cleanup)
When MetaMorpheus is bumped to an mzLib that includes `MassSpectrometry.Deconvolution.Consensus`, **delete this
folder** and let `ConsensusMassXicConstructor` resolve the tracer from mzLib instead (its `using` is already
`MassSpectrometry.Deconvolution.Consensus`, so no code change is needed once the namespace comes from the package).
