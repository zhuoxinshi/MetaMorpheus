# ISD pseudo-MS2 spectra from TopFD features

This branch (`ISD_pub`) holds the code that builds the pseudo-MS2 spectra of the in-source dissociation (ISD)
data in the HCD_YBC benchmark: six runs (mixtures 1-1 and 1-2, three replicates each) acquired with a repeating
cycle of four full scans at source voltages 15 V (precursors) and 60, 80 and 100 V (fragments). Each run gives one
`*_pseudo_ms2.msalign` file for top-down search with TopPIC.

## Workflow

1. **Split each .raw file by source voltage** (mzLib, branch `ISD_pub`, `LabelCorrection.SeparateFilesForTopFD`):
   `{stem}_ms1.mzML` (sid=15) and `{stem}_isd60.mzML`, `{stem}_isd80.mzML`, `{stem}_isd100.mzML`. Scans keep their
   .raw scan numbers.
2. **TopFD 1.8.1** on each of the four mzML files, giving `{stem}_ms1_feature.xml`, `{stem}_isd60_feature.xml`,
   `{stem}_isd80_feature.xml` and `{stem}_isd100_feature.xml`. All 24 runs used the same settings: centroid data,
   maximum charge 50, maximum monoisotopic mass 30,000 Da, peak m/z error tolerance 0.02, MS1 signal/noise ratio 3,
   feature minimum scan number 3, intensity ratio for splitting features 2.5, feature ECScore cutoff 0.1.
3. **Pseudo-MS2 spectra** (this branch), per run, from the .raw file and the four feature files:
   - Precursor XICs are read from the 15 V features and fragment XICs from the 60, 80 and 100 V features
     (`FeatureXmlXicConstructor`). Each feature gives one XIC: per scan, the envelope intensities of all its charge
     states are summed, at the feature's monoisotopic mass. XICs are interpolated with a cubic spline on the scan
     index axis.
   - Fragment XICs of the three voltages are pooled and grouped with each precursor XIC (`XicGroupingEngine`)
     by apex distance (in acquisition cycles) and Pearson correlation of the splined profiles.
   - Each group becomes one pseudo-MS2 spectrum (`PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup`): one peak
     per fragment XIC, with the feature's monoisotopic mass, the charge of the XIC's first point (the most intense
     charge state in that scan) and the summed XIC intensity. The precursor has the precursor feature's
     monoisotopic mass, the charge at the XIC apex, the summed XIC intensity and the TopFD feature id
     (`PRECURSOR_FEATURE_ID`).
   - `ISDEngine.WriteMsAlignFile` writes the spectra in msAlign format.

## Settings

| Step | Setting | Value |
|---|---|---|
| Precursor XICs | feature mass | 4,000 to 70,000 Da |
| Fragment XICs | feature mass | below 70,000 Da |
| Both | minimum ECScore | 0 (no filter beyond TopFD's 0.1 cutoff) |
| Both | minimum XIC points (scans of the channel with an envelope) | 3 |
| Both | maximum XIC width (first to last point) | 2 min |
| Both | spline | cubic, scan-index axis, sampled every 0.05 scan, one zero point added one scan beyond each end |
| Grouping | apex tolerance | 2 scans (`ApexToleranceUnit.ScanIndex`) |
| Grouping | Pearson correlation | at least 0.5 |
| Grouping | XIC overlap threshold | 0 (no requirement) |
| Grouping | fragments kept per precursor | the 100 best correlated |
| Grouping | precursor rank filter | none |
| Grouping | fragments needed for a spectrum | more than 10 |
| Spectrum | peak and precursor intensity | summed XIC intensity (`PseudoMs2IntensityType.XicSum`) |

All inputs and settings are set in one test, `IsdPseudoMs2Reproduction.WriteHcdYbcIsdPseudoMs2`
(`Test/DIATests/IsdPseudoMs2Reproduction.cs`). It builds the `DIAparameters` of each run and passes them to
`WritePseudoMs2MsAlign`, which also writes them to `DIAparameters.txt`.

## Running it

The explicit test `IsdPseudoMs2Reproduction.WriteHcdYbcIsdPseudoMs2` runs step 3 for all six runs. Set the folders
(raw files, TopFD feature files, output) and, if needed, the settings at the top of the test, then:

```
dotnet test MetaMorpheus/Test/Test.csproj -c Release --filter "FullyQualifiedName~Test.DIATests.IsdPseudoMs2Reproduction"
```

It writes `{stem}_pseudo_ms2.msalign` for each run, an empty `{stem}_pseudo_ms1.msalign` beside it (TopPIC's
combined-file mode expects one) and `DIAparameters.txt`. `TestFeatureXmlXicConstructor` covers the feature-file
reader and XIC construction on synthetic data and runs with the normal test suite.

## Verification

This branch (MetaMorpheus master `c037f64b6` plus the changes listed below; .NET 8, mzLib 1.0.583; builds in
Visual Studio 2022) regenerated the six published files (folder `2_0.5_100_null_10_fewer_var`, written 2026-09-30 by the analysis branch
`isd_export_zs`) from the same .raw and TopFD files. The content is identical: 17,247 spectra and 876,836 peak lines,
with every header field and peak line identical as text, and each file has the same size in bytes. Only the scan
numbers (`ID`, `SCANS`, ...) differ. The original numbering followed the order in which the parallel grouping
threads finished; this branch numbers spectra in precursor-XIC order instead, and two independent runs gave
byte-identical files. The same changes on MetaMorpheus master `bafd1d819` (.NET 10, mzLib 1.0.597) gave
byte-identical files as well.

## Code

| File | Role |
|---|---|
| `EngineLayer/DIA/XicConstruction/TopFdFeatureXmlFile.cs` | Reader for TopFD `*_feature.xml` (features, charge states, per-scan envelopes). |
| `EngineLayer/DIA/XicConstruction/FeatureXmlXicConstructor.cs` | XICs from TopFD features; `BySourceVoltage` picks each fragment channel's feature file from `sid=` in the scan filter. |
| `EngineLayer/DIA/Enums/ApexToleranceUnit.cs` | Apex tolerance in minutes (default) or in scans. |
| `EngineLayer/DIA/Enums/PseudoMs2IntensityType.cs` | Apex intensity (default) or summed XIC intensity for pseudo-MS2 peaks. |
| `EngineLayer/DIA/ISDEngine.cs` | Precursor feature ids, splined fragment XICs, scan numbering in precursor order, `WriteMsAlignFile`. |
| `EngineLayer/DIA/PrecursorFragmentGroup.cs` | Intensity option, fragment envelopes with real intensity and charge, precursor intensity and feature id. |
| `EngineLayer/DIA/PrecursorFragmentGrouping/XicGroupingEngine.cs` | Apex window in scans as an option. |
| `EngineLayer/DIA/DIAparameters.cs`, `DIAEngine.cs`, `Ms2ScanWithSpecificMass.cs` | Intensity option, settings record, `PrecursorFeatureId`. |

The new options default to the existing behaviour (apex window in minutes, apex intensities); the ISD settings
above are passed explicitly.
