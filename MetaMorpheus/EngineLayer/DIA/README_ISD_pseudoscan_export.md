# ISD pseudo-scan generation + msalign export

This branch (`isd-pseudoscan-msalign-export`) contains the code that turns an all-MS1 **ISD** acquisition into
searchable **pseudo-MS2** scans, plus a static exporter that writes those scans to **TopPIC** `.msalign` files.

## What generates the pseudo scans (already in this branch)
The ISD pipeline lives in `EngineLayer/DIA/`:
- **XIC construction** (`XicConstruction/NeutralMassXicConstructor.cs`) — deconvolutes each MS1 scan and builds
  neutral-mass extracted-ion chromatograms (the "consensus feature tracing" step).
- **Precursor–fragment grouping** (`PrecursorFragmentGrouping/XicGroupingEngine.cs`) — groups fragment XICs to a
  precursor XIC by **apex-RT tolerance + RT-overlap + Pearson correlation** of the elution profiles.
- **Pseudo-MS2 generation** (`ISDEngine.cs`, `PrecursorFragmentGroup.GetPseudoMs2ScanFromPfGroup`) — turns each
  precursor–fragment group into an `Ms2ScanWithSpecificMass` (precursor + neutral-mass fragments).

Source voltage is read from the Thermo scan filter (`sid=<V>`); the lowest voltage is the intact/precursor
channel, the higher voltages are fragment channels.

## New in this branch: `IsdMsAlignExporter`
`EngineLayer/DIA/IsdMsAlignExporter.cs` — static methods:
- `WriteMs2Align(pseudoScans, path)` → a TopPIC-format `ms2.msalign` (one `BEGIN IONS…END IONS` block per
  pseudo scan: `PRECURSOR_MZ/CHARGE/MASS/INTENSITY` + `mass⇥intensity⇥charge` fragment peaks).
- `WriteMs1Align(pseudoScans, path)` → a companion `ms1.msalign` (one MS1 entry per precursor).

Each pseudo scan is already neutral-mass deconvoluted, so it maps directly onto the msalign format and can be
searched in TopPIC exactly like ordinary DDA MS2 data.

## How to run it (any computer)
```
git clone https://github.com/trishorts/MetaMorpheus
cd MetaMorpheus
git checkout isd-pseudoscan-msalign-export
```
Then, in code (see the worked example in `Test/DIATests/IsdMsAlignExportTest.cs`):
```csharp
var dataFile = MsDataFileReader.GetDataFile("your_ISD.mzML"); dataFile.LoadAllStaticData();
var dia = new DIAparameters(DIAanalysisType.ISD,
    new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1,20,4,3)),
    new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1,20,4,3)),
    new XicGroupingEngine(apexRTTolerance: 0.2f, overlapThreshold: 0.5, correlationThreshold: 0.5, maxThreadsForGrouping: 1),
    PseudoMs2ConstructionType.Mass, combineFragments: true);
var common = new CommonParameters { DIAparameters = dia };
var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, "your_ISD.mzML", common).ToArray();
IsdMsAlignExporter.WriteMs2Align(pseudoScans, "your_ISD_ms2.msalign");
IsdMsAlignExporter.WriteMs1Align(pseudoScans, "your_ISD_ms1.msalign");   // then search the pair in TopPIC
```
To verify the build end-to-end on the bundled fixture:
```
dotnet test MetaMorpheus/MetaMorpheus/Test/Test.csproj --filter FullyQualifiedName~IsdMsAlignExportTest
```

## Grouping knobs (`XicGroupingEngine`)
- `apexRTTolerance` (min) — max apex-RT difference precursor↔fragment (default used here: **0.2 min**).
- `overlapThreshold` — min RT-overlap ratio (0.5).
- `correlationThreshold` — min Pearson correlation of the two elution profiles (0.5).
- `combineFragments` — true pools 60/80/100 V fragments into one pseudo scan per precursor; false makes one per voltage.

> Note: the ISD path is currently reachable from code/tests (set `CommonParameters.DIAparameters`); it is not yet
> wired to a TOML/CLI switch. The example above is the supported entry point.
