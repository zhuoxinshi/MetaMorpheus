using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using EngineLayer;
using NUnit.Framework;
using System.IO;
using TaskLayer;
using EngineLayer.DIA;
using MzLibUtil;
using Plotly.NET;
using UsefulProteomicsDatabases;
using Nett;
using System.IO.Compression;
using Chemistry;
using MzLibUtil;
using MathNet.Numerics.Statistics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;
using Accord.Math;
using Accord.Statistics.Models.Regression.Fitting;
using Accord.Math.Optimization.Losses;
using Microsoft.VisualStudio.TestPlatform.Utilities;
using static Microsoft.FSharp.Core.ByRefKinds;
using Accord.Statistics.Models.Regression.Linear;
using Omics;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using Omics.Fragmentation;
using OxyPlot;
using Org.BouncyCastle.Asn1.Mozilla;
using NWaves.Filters;
using NWaves.Signals;

namespace Test.TestDIA
{
    public class newTestStatic
    {
        [Test]
        public static void TestScanBased()
        {
            string snip = @"E:\ISD Project\CE_241118\11-21-24_CE_6prot-ammon-acet_5AA_500nL-pHjunction_ISDrep1.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 8000, maxMass: 9000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.OverlapFirst,
            pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISD_scanBased, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc);

            var scans = ISD_scanBased.GetPseudoMs2Scans(dataFile, task.CommonParameters, task.CommonParameters.DIAparameters);

        }

        [Test]
        public static void TestXIC()
        {
            string snip = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT28.01-32.69.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, analysisType: AnalysisType.ISDEngine_static);

            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            var allMs1PeakCurves1 = ISDEngine_static.GetAllPeakCurves(ms1Scans, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.DeconHighestPeak,
                task.CommonParameters.DIAparameters.Ms1PeakFindingTolerance, task.CommonParameters.DIAparameters.MaxRTRangeMS1,
                out List<Peak>[] peaksByScan1);

            //test case ms1 peak: deconHighestPeak
            var testPeak1 = peaksByScan1[57].Where(p => Math.Abs(p.Mz - 1032.68) < 0.01).First();
            Assert.That(Math.Round(testPeak1.PeakCurve.StartRT, 2).Equals(28.62));
            Assert.That(Math.Round(testPeak1.PeakCurve.EndRT, 2).Equals(29.85));
            Assert.That(Math.Round(testPeak1.PeakCurve.ApexRT, 2).Equals(28.73));

            var testPeak2 = peaksByScan1[73].Where(p => Math.Abs(p.Mz - 1087.57) < 0.01).First();
            Assert.That(Math.Round(testPeak2.PeakCurve.StartRT, 2).Equals(29.65));
            Assert.That(Math.Round(testPeak2.PeakCurve.EndRT, 2).Equals(30.26));
            Assert.That(Math.Round(testPeak2.PeakCurve.ApexRT, 2).Equals(29.85));

            //Test ms1 peak XIC: isoEnvelopeTotal
            var allMs1PeakCurves2 = ISDEngine_static.GetAllPeakCurves(ms1Scans, task.CommonParameters, task.CommonParameters.DIAparameters, XICType.isoEnvelopeTotal,
                task.CommonParameters.DIAparameters.Ms1PeakFindingTolerance, 2, out List<Peak>[] peaksByScan2);
            var testPeak3 = peaksByScan2[57].Where(p => Math.Abs(p.Mz - 1032.68) < 0.01).First();
            Assert.That(Math.Round(testPeak3.PeakCurve.StartRT, 2).Equals(28.62));
            Assert.That(Math.Round(testPeak3.PeakCurve.EndRT, 2).Equals(29.85));
            Assert.That(Math.Round(testPeak3.PeakCurve.ApexRT, 2).Equals(28.73));
        }

        [Test]
        public static void TestMatchedFragmentIons()
        {
            //load ISD PSMs
            var psmFileISD = @"E:\ISD Project\ISD_240606\search_w-writingSpectralLib\Task1-SearchTask\AllPSMs.psmtsv";
            var psmsISD = PsmTsvReader.ReadTsv(psmFileISD, out List<string> warnings).Where(psm => psm.QValue < 0.01);
            var psmGroupedBySequence_ISD = psmsISD.GroupBy(psm => psm.FullSequence).ToList();
            var sequenceToTest = psmGroupedBySequence_ISD.Select(g => g.Key).ToArray()[5];
            var ubISD = psmGroupedBySequence_ISD.Where(g => g.Key == sequenceToTest).First();
            var massToTest = Convert.ToDouble(ubISD.First().PeptideMonoMass);
            var allMatchedIonsISD = ubISD.Select(psm => psm.MatchedIons).SelectMany(i => i).ToList();

            //load DDA PSMs
            var psmFileDDA = @"E:\ISD Project\ISD_240606\sample10-DDA_w-writingSpectralLib\Task1-SearchTask\AllPSMs.psmtsv";
            var psmsDDA = PsmTsvReader.ReadTsv(psmFileDDA, out List<string> warningsDDA).Where(psm => psm.QValue < 0.01);
            var psmGroupedBySequence_DDA = psmsDDA.GroupBy(psm => psm.FullSequence).ToList();
            var ubDDA = psmGroupedBySequence_DDA.Where(g => g.Key == sequenceToTest).First(); ;
            var sharedIons_DDA = ubDDA.Select(psm => psm.MatchedIons.Select(i => i.Annotation)).Aggregate((current, next) => current.Intersect(next).ToList());
            var allMatchedIonsDDA = ubDDA.Select(psm => psm.MatchedIons).SelectMany(i => i).ToList();

            //Target ion list
            var overlapIons = allMatchedIonsISD.Where(i => allMatchedIonsDDA.Select(i => i.Annotation).Contains(i.Annotation)).ToList();
            var targetList = overlapIons.Select(i => new { annotation = i.Annotation, mass = Math.Round(i.NeutralTheoreticalProduct.MonoisotopicMass, 0), charge = i.Charge })
                .Distinct().ToList();

            //precursor fragment grouping
            string snip = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD_RT35.01-37.96.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 8000, maxMass: 9000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.Area_correlation,
            pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISD_scanBased, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc);

            //check if target ions are included in corresponding pseudoMS2scans
            var pseudoScans = ISDEngine_static.GetPseudoMs2Scans(dataFile, task.CommonParameters, task.CommonParameters.DIAparameters);
            var scans = pseudoScans.Where(scan => Math.Round(scan.PrecursorMass, 0) == Math.Round(massToTest, 0)).ToList();

            var scores = new List<double>();
            var notFound = new Dictionary<int, List<(string, double, double, double)>>();
            foreach (var scan in scans)
            {
                double score = 0;
                var no = new List<(string, double, double, double)>();
                foreach (var ion in targetList)
                {
                    var fragments = scan.ExperimentalFragments.Where(e => Math.Round(e.MonoisotopicMass) == ion.mass && e.Charge == ion.charge);
                    if (!fragments.Any())
                    {
                        no.Add(new(ion.annotation, ion.mass, ion.charge, ion.mass.ToMz(ion.charge)));
                    }
                    else
                    {
                        score++;
                    }
                }
                scores.Add(score / (double)targetList.Count);
                notFound.Add(scan.OneBasedScanNumber, no);
            }


        }

        [Test]
        public static void TestMatchedFragmentIons2()
        {
            //load ISD PSMs
            var psmFileISD = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task4-SearchTask\Individual File Results\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100-calib-averaged_PSMs.psmtsv";
            var psmsISD = PsmTsvReader.ReadTsv(psmFileISD, out List<string> warnings).Where(psm => psm.QValue < 0.01 && psm.DecoyContamTarget == "T");
            var psmGroupedBySequence_ISD = psmsISD.GroupBy(psm => psm.FullSequence).ToList();
            var sequenceToTest = "APSAKATAAKKAVVKGTNGKKALKVRTSATFRLPKTLKLARAPKYASKAVPHYNRLDSYKVIEQPITSETAMKK[Common Biological:Methylation on K]VEDGNILVFQVSMKANKYQIKKAVKELYEVDVLKVNTLVRP" +
                "NGTKKAYVRLTADYDALDIANRIGYI";
            var ISD = psmGroupedBySequence_ISD.Where(g => g.Key == sequenceToTest).First();
            var allMatchedIonsISD = ISD.Select(psm => psm.MatchedIons).SelectMany(i => i).ToList();

            //load DDA PSMs
            var psmFileDDA = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task4-SearchTask\Individual File Results\08-12-24_PEPPI_FractionD_orbiMS1_DDA-calib-averaged_PSMs.psmtsv";
            var psmsDDA = PsmTsvReader.ReadTsv(psmFileDDA, out List<string> warningsDDA).Where(psm => psm.QValue < 0.01);
            var psmGroupedBySequence_DDA = psmsDDA.GroupBy(psm => psm.FullSequence).ToList();
            var DDA = psmGroupedBySequence_DDA.Where(g => g.Key == sequenceToTest).First();
            var allMatchedIonsDDA = DDA.Select(psm => psm.MatchedIons).SelectMany(i => i).ToList();

            //Target ion list
            var overlapIons = allMatchedIonsISD.Where(i => allMatchedIonsDDA.Select(i => i.Annotation).Contains(i.Annotation)).ToList();
            var targetList = overlapIons.Select(i => new { annotation = i.Annotation, mass = Math.Round(i.NeutralTheoreticalProduct.MonoisotopicMass, 0), charge = i.Charge })
                .Distinct().ToList();

            //search
            string snip = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT28.01-32.69.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, maxMass: 20000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.ScanCycle,
            pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc);

            //check if target ions are included in corresponding pseudoMS2scans
            var pseudoScans = ISDEngine_static.GetPseudoMs2Scans(dataFile, task.CommonParameters, task.CommonParameters.DIAparameters);
            var scans = pseudoScans.Where(scan => Math.Round(scan.PrecursorMass, 0) == Math.Round(15632.6539, 0)).ToList();
        }

        [Test]
        public static void TestUbiquitinSnip()
        {
            string snip = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD_RT35.01-37.96.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 8000, maxMass: 9000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005 , ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.DeconHighestPeak, pfGroupingType: PFGroupingType.ScanCycle,
            pseudoMs2Type: PseudoMs2ConstructionType.neutralMass, analysisType: AnalysisType.ISD_slow, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc);

            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\ISD_slow_sample1_5uL_RT35.01-37.96_DeconHighestPeakMS1&MS2_ms1Tol5ppm_maxRT2_maxMissed2_scanCycle_apexCycle3_overlap0.3_corr0.5_para";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { snip }, "test");
        }

        [Test]
        public static void TestSearch()
        {
            string snip = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT28.01-32.69.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 8, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, maxMass: 20000, type: "ISD", apexCycleTolerance: 5,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.OverlapFirst,
            pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc);

            string myDatabase = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task3-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\ISDEngine_static_FD-RT28.01-32.69_DeconHighestPeakMS1-maxRT2_mzPeakMS2-maxRT8_Tol5ppm_maxMissed2_overlapFirst_apexCycle5_overlap0.3_corr0.5";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { snip }, "test");
        }

        [Test]
        public static void AveragedISD100_deconTest()
        {
            var originalPath = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD_RT35.01-37.78_labelCorrected.mzML";
            var averagedPath = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD_RT35.01-37.78_averageISD100_noNorm_weightEvenly_labelCorrected.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myFileManagers = new MyFileManager(true);
            var originalFile = myFileManagers.LoadFile(originalPath, task.CommonParameters);
            var averagedFile = myFileManagers.LoadFile(averagedPath, task.CommonParameters);
            var isd100_orig = originalFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            var isd100_avg = averagedFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();

            var envelopes_orig = new List<IsotopicEnvelope>();
            var envelopes_avg = new List<IsotopicEnvelope>();
            foreach (var scan in isd100_orig)
            {
                envelopes_orig.AddRange(Deconvoluter.Deconvolute(scan, task.CommonParameters.PrecursorDeconvolutionParameters));
            }
            foreach (var scan in isd100_avg)
            {
                envelopes_avg.AddRange(Deconvoluter.Deconvolute(scan, task.CommonParameters.PrecursorDeconvolutionParameters));
            }

            var ms2scans = MetaMorpheusTask.GetMs2Scans(originalFile, originalPath, task.CommonParameters);
            var scan1 = ms2scans.FirstOrDefault();
            var testEnvelope = scan1.ExperimentalFragments.FirstOrDefault();
            var neutralFrag = new IsotopicEnvelope[] { testEnvelope };
            var newMs2WithMass = new Ms2ScanWithSpecificMass(scan1.TheScan, 1, 1, null, task.CommonParameters, neutralFrag);
        }

        [Test]
        public static void TestLinearRegression()
        {
            var path = @"E:\ISD Project\CE_241213\12-13-24_CE_6prot-ammon-acet_5AA_500nL-pHjunction_ISD100_micro4_60k_RT30.08-32.01_labelCorrected.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task Settings\Task4-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 8, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, maxMass: 20000, type: "ISD", apexCycleTolerance: 5,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.OverlapFirst,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, task.CommonParameters);
            var allScans = dataFile.GetAllScansList();
            var ms1scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2scans = allScans.Where(s => s.MsnOrder == 2).ToArray();
            var psmsFile = @"E:\ISD Project\CE_241213\snip_highestPeak\Task1-SearchTask\AllPSMs.psmtsv";
            var psms = PsmTsvReader.ReadTsv(psmsFile, out List<string> warnings);
            var psmsContainTheIon = psms.Where(p => p.MatchedIons.Select(i => Math.Round(i.Mz, 2)).Contains(1017.07)).ToList();
            var selectedPsms = psmsContainTheIon.OrderByDescending(psm => psm.PrecursorIntensity).GroupBy(p => p.FullSequence).Select(g => g.First()).ToList();
            var highestPeakMzs = selectedPsms.Select(psm => psm.PrecursorHighestPeakMz).ToList();

            //Find the fragment XIC and the precursor XICs
            var pc1017 = PeakCurve.PeakTracing(1017.58, 1, ms2scans, new PpmTolerance(10), 100, 2, 8);
            var precursorPCs = new List<PeakCurve>();
            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves_Peak(ms1scans, task.CommonParameters.DIAparameters, new PpmTolerance(5),
                8, out List<Peak>[] peaks);
            var allPeaks = peaks.Where(v => v!= null).SelectMany(p => p).OrderBy(p => p.Mz).ToList();
            foreach (var mz in highestPeakMzs)
            {
                var peak = allPeaks.Where(p => Math.Round(p.Mz, 5) == Math.Round(mz, 5)).First();
                if (peak.PeakCurve != null)
                {
                    precursorPCs.Add(peak.PeakCurve);
                }
            }

            //Find the RT range for all the XICs
            var allPCs = new List<PeakCurve> { pc1017 };
            allPCs.AddRange(precursorPCs);
            var startCycle = allPCs.Min(p => p.StartCycle);
            var endCycle = allPCs.Max(p => p.EndCycle);
            var rtLength = endCycle - startCycle + 1;
            var totalPreIntensity = precursorPCs.Select(p => p.Peaks.Select(p => p.Intensity).Sum()).Sum();

            //plot
            var chartList = new List<GenericChart>();
            foreach(var pc in allPCs)
            {
                var plot = pc.VisualizeRaw("line");
                chartList.Add(plot);
            }
            var combinedPlot = Chart.Combine(chartList);
            //combinedPlot.Show();

            //construct the feature matrix
            var dataMatrix = new double[precursorPCs.Count, rtLength];
            for (int i = 0; i < precursorPCs.Count; i++)
            {
                for (int j = startCycle; j <= endCycle; j++)
                {
                    if (j < precursorPCs[i].StartCycle || j > precursorPCs[i].EndCycle)
                    {
                        dataMatrix[i, j] = 0;
                        continue;
                    }
                    int index = Array.BinarySearch(precursorPCs[i].Peaks.Select(p => p.ZeroBasedScanIndex).ToArray(), j);
                    if (index >= 0)
                    {
                        dataMatrix[i, j] = precursorPCs[i].Peaks[index].Intensity;
                    } else
                    {
                        dataMatrix[i, j] = 0; 
                    }
                }
            }
            var normalizedMatrix = new double[precursorPCs.Count, rtLength];
            for (int i = 0; i < precursorPCs.Count; i++)
            {
                //var sum = precursorPCs[i].Peaks.Select(p => p.Intensity).Sum();
                for (int j = 0; j < rtLength; j++)
                {
                    normalizedMatrix[i, j] = dataMatrix[i, j] / totalPreIntensity;
                }
            }
            var A = Matrix<double>.Build.DenseOfArray(dataMatrix);
            var columnSums = A.ColumnSums();
            var precursorArray = new double[rtLength][];
            for (int j = 0; j < rtLength; j++)
            {
                var totalIntensityForTheRT = columnSums[j];
                var maxIntensityForTheRT = A.Column(j).Max();
                var minIntensityForTheRT = A.Column(j).Min();
                precursorArray[j] = new double[precursorPCs.Count];
                for (int i = 0; i < precursorPCs.Count; i++)
                {
                    precursorArray[j][i] = dataMatrix[i, j] / totalIntensityForTheRT;
                }
            }

            //construct the result vector
            var fragmentPC = new double[rtLength];
            var fragmentTotalIntensity = pc1017.Peaks.Select(p => p.Intensity).Sum();
            var fragmentMaxIntensity = pc1017.Peaks.Select(p => p.Intensity).Max();
            var fragmentMinIntensity = pc1017.Peaks.Select(p => p.Intensity).Min();
            for (int i = 0; i < fragmentPC.Length; i++)
            {
                int index = Array.BinarySearch(pc1017.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray(), i);
                if (index >= 0)
                {
                    fragmentPC[i] = pc1017.Peaks[index].Intensity / fragmentTotalIntensity;
                } else
                {
                    fragmentPC[i] = 0;
                }
            }

            //NNLS
            //TODO: try different normalization; try including all the precursor peaks instead of grouping by full sequence
            var b = Vector<double>.Build.Dense(fragmentPC);
            var nnls = new NonNegativeLeastSquares()
            {
                MaxIterations = 1000
            };
            var regression = nnls.Learn(precursorArray, fragmentPC);
            var coefficients = regression.Weights;
            double[] prediction = regression.Transform(precursorArray);
            double error = new SquareLoss(expected: fragmentPC).Loss(actual: prediction);
            //foreach (var mz in highestPeakMzs)
            //{
            //    int i = 0;
            //    var pc = PeakCurve.PeakTracing(mz, 0, ms1scans, new PpmTolerance(10), 100, 2, 8);
            //    while (pc == null)
            //    {
            //        i++;
            //        pc = PeakCurve.PeakTracing(mz, i, ms1scans, new PpmTolerance(10), 100, 2, 8);
            //    }
            //    precursorPCs.Add(pc);
            //}
            //precursorPCs = precursorPCs.OrderByDescending(p => p.AveragedIntensity).ToList();
            //var allGroups = new List<PrecursorFragmentsGroup>();
            //foreach (var pc in precursorPCs)
            //{
            //    if (pc.PFGroup != null)
            //    {
            //        continue;
            //    }
            //    var ppGroup = new PrecursorFragmentsGroup(pc);
            //    foreach (var p in precursorPCs)
            //    {
            //        var corr = PrecursorFragmentPair.CalculatePeakCurveCorr(pc, p);
            //        if (corr > 0.85)
            //        {
            //            var ppPair = new PrecursorFragmentPair(pc, p, corr);
            //            ppGroup.PFpairs.Add(ppPair);
            //            p.PFGroup = ppGroup;
            //        }
            //    }
            //    allGroups.Add(ppGroup);
            //}

            //var startCycle = precursorPCs.Min(p => p.StartCycle);
            //var endCycle = precursorPCs.Max(p => p.EndCycle);
            //var combinedPeakCurves = new List<PeakCurve>();
            //var dataMatrix = new double[allGroups.Count, endCycle - startCycle + 1];
            //foreach (var group in allGroups)
            //{
            //    var pcs = new List<PeakCurve> { group.PrecursorPeakCurve };
            //    pcs.AddRange(group.PFpairs.Select(pair => pair.FragmentPeakCurve));
            //    for (int i = 0; i < dataMatrix.Length; i++)
            //    {

            //    }
            //}

        }

        [Test]
        public static void TestNNLS()
        {
            var inputs = new double[][]
            {
                new[] { 1.0, 1.0, 1.0 },
                new[] { 2.0, 4.0, 8.0 },
                new[] { 3.0, 9.0, 27.0 },
                new[] { 4.0, 16.0, 64.0 },
            };

            var outputs = new double[] { 3, 14, 39, 84 };

            // Create a NN LS learning algorithm
            var nnls = new NonNegativeLeastSquares()
            {
                MaxIterations = 1000
            };

            // Use the algorithm to learn a multiple linear regression
            MultipleLinearRegression regression = nnls.Learn(inputs, outputs);

            // None of the regression coefficients should be negative:
            double[] coefficients = regression.Weights; // should be

            // Check the quality of the regression:
            double[] prediction = regression.Transform(inputs);

            double error = new SquareLoss(expected: outputs).Loss(actual: prediction); // should be 0
        }

        [Test]
        public static void TestBottomUpISD()
        {
            var path = @"E:\DIA\Data\DIA-bottom-up_241218\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            string tomlFile = @"E:\DIA\Data\DIA-bottom-up_241218\bu-ISD_mix1_120k\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 0, maxMass: 5000, type: "ISD", apexCycleTolerance: 5,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.ScanCycle,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc,
                cutMs1Peaks: false, cutMs2Peaks: false);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, task.CommonParameters);

            var peptide = new PeptideWithSetModifications("VLDALDSIK", new Dictionary<string, Modification>());
            var fragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, task.CommonParameters.DigestionParams.FragmentationTerminus, fragments);

            //Get ms1 XICs
            var ms1Scans = dataFile.GetMS1Scans().Where(s => s.RetentionTime >= 38 && s.RetentionTime <= 40).ToArray();
            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, task.CommonParameters, task.CommonParameters.DIAparameters, 
                task.CommonParameters.DIAparameters.Ms1XICType, task.CommonParameters.DIAparameters.Ms1PeakFindingTolerance, 
                task.CommonParameters.DIAparameters.MaxRTRangeMS1, out List<Peak>[] peaksByScan, task.CommonParameters.DIAparameters.CutMs1Peaks);

            //Get ms2 XICs
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).Where(s => s.RetentionTime >= 38 && s.RetentionTime <= 40).ToArray();
            var allMs2PeakCurves = ISDEngine_static.GetAllPeakCurves(ms2Scans, task.CommonParameters, task.CommonParameters.DIAparameters, 
                task.CommonParameters.DIAparameters.Ms2XICType, task.CommonParameters.DIAparameters.Ms2PeakFindingTolerance, 
                task.CommonParameters.DIAparameters.MaxRTRangeMS2, out List<Peak>[] peaksByScan2, task.CommonParameters.DIAparameters.CutMs2Peaks);

            //Group XICs for the peptide
            var preXIC = allMs1PeakCurves.Where(pc => Math.Abs(pc.AveragedMz - 487.28) < 0.01 && pc.Charge == 2
            ).First();

            preXIC.VisualizePeakRegions();
            preXIC.CutPeak();
            preXIC.VisualizeRaw("line").Show();
            var map = ISDEngine_static.GetRtIndexMap(ms1Scans);
            preXIC.GetSGfilterSmoothedData(map, 9, 3);
            //Find the minimum, compare the discrimination factor of the left and the right point, keep the point with higher discrimination factor side
            //SG filter first, then cut peak; wavelet or cutPeak
            //var SGSmoothed = preXIC.GetSGfilterSmoothedData(9, 2);
            //preXIC.Visualize(preXIC.Peaks.Select(p => p.RetentionTime).ToArray(), SGSmoothed).Show();

            var pfGroup = ISDEngine_static.PFgrouping(preXIC, allMs2PeakCurves, task.CommonParameters.DIAparameters);
            var pseudoMs2 = ISDEngine_static.ConstructNewMs2Scans(pfGroup, task.CommonParameters, task.CommonParameters.DIAparameters.PseudoMs2ConstructionType, dataFile.FilePath);

            //deconvolute the raw spectrum at apex
            var apexScan = ms1Scans.Where(s => s.RetentionTime == preXIC.ApexRT).First();
            var apexMs2 = ms2Scans.Where(s => s.OneBasedScanNumber == apexScan.OneBasedScanNumber + 1).First();
            var allObservedMasses = Deconvoluter.Deconvolute(apexMs2, task.CommonParameters.PrecursorDeconvolutionParameters);
            var observedTheoreticalMasses = new List<Product>();
            var roundedObservedMasses = allObservedMasses.Select(m => Math.Round(m.MonoisotopicMass, 0)).ToList();
            foreach (var fragment in fragments)
            {
                if (roundedObservedMasses.Contains(Math.Round(fragment.MonoisotopicMass, 0)))
                {
                    observedTheoreticalMasses.Add(fragment);
                }
            }
            var matchedMs2Curves = pfGroup.PFpairs.Select(p => p.FragmentPeakCurve).ToList();
            var deconvolutedMasses = pseudoMs2.ExperimentalFragments.Select(p => Math.Round(p.MonoisotopicMass, 1)).ToList();
            var matchedMasses = new List<double>();
            var unmatchedMasses = new List<double>();
            foreach (var fragment in fragments)
            {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                if (deconvolutedMasses.Contains(Math.Round(fragment.MonoisotopicMass, 1)))
                {
                    matchedMasses.Add(fragment.MonoisotopicMass);
                }
                else
                {
                    unmatchedMasses.Add(fragment.MonoisotopicMass);
                }
            }
        }

        [Test]
        public static void TestVisualization()
        {
            var path = @"E:\DIA\Data\DIA-bottom-up_241218\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            string tomlFile = @"E:\DIA\Data\DIA-bottom-up_241218\bu-ISD_mix1_120k\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 0, maxMass: 5000, type: "ISD", apexCycleTolerance: 5,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.ScanCycle,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc,
                cutMs1Peaks: false, cutMs2Peaks: false);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, task.CommonParameters);

            var ms1Scans = dataFile.GetMS1Scans().Where(s => s.RetentionTime >= 38 && s.RetentionTime <= 40).ToArray();
            var allMs1PeakCurves = ISDEngine_static.GetAllPeakCurves(ms1Scans, task.CommonParameters, task.CommonParameters.DIAparameters,
                task.CommonParameters.DIAparameters.Ms1XICType, task.CommonParameters.DIAparameters.Ms1PeakFindingTolerance,
                task.CommonParameters.DIAparameters.MaxRTRangeMS1, out List<Peak>[] peaksByScan, task.CommonParameters.DIAparameters.CutMs1Peaks);
            var preXIC = allMs1PeakCurves.Where(pc => Math.Abs(pc.AveragedMz - 487.28) < 0.01 && pc.Charge == 2
            ).First();

            //preXIC.VisualizeRaw("line").Show();
            var map = ISDEngine_static.GetRtIndexMap(ms1Scans);
            var smoothedData = preXIC.GetSGfilterSmoothedData(map, 51, 2);
            preXIC.VisualizeGeneral(smoothedData.Select(p => p.Item1).ToArray(), smoothedData.Select(p => p.Item2).ToArray()).Show();
        }

        [Test]
        public static void TestTemp()
        {
            var spectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 10, 20, 30 }, false);
            double[] intensity = { 645523.0313,  2840027.875, 3638575.93, 10811454.23, 18073620.34,
                       18725528.47, 27232004.29, 40920044.22, 42583805.95, 55198342.05, 56973095.88, 74472979.85, 81412777.7,
                       86125057.49, 72334569.94, 76413935.04, 62766233.79, 61580305.79, 63975710.61, 50340044.05, 45728131.18,
                       53059923.13, 46062358.94, 33987709.73,  28573307.18, 31451496.96, 30863762.17, 31758881.85,
                       44452680.13, 51583712.98, 65542667.65, 59104996.7, 54094069.27, 49302712.43, 47378242.51, 42609288.13,
                       39850707.92, 35733177.29, 42743239.44, 39426074.07, 34068996.89, 31745392.82, 34066514.36, 33100802.68,
                       27264968.85, 28240874.96, 28044224.33, 24225874.99, 24133904.89,  23105600.01,
                       24196788.96, 26124754.13, 23572928.42, 22178019.58, 19562688.35, 18567397.43, 17071497.79,
                       14410106.09, 15909702.51, 15499873.91, 13752143.89, 14551270.77, 13127558.88, 14314424.74, 12676548.31,
                       12269520.52, 12036858.71, 10631350.14, 11887922.79, 9704722.864,  9030259.031, 9308738.631,
                       7320634.788, 7756802.779, 7140971.624, 6266278.887, 6337805.943, 5971872.621,  5847247.232,
                       4387509.91, 4206881.253, 2682847.883,
                       2188591.858, 2596483.797, 2942508.441, 2321496.598, 2051064.086,  2538946.367, 1944507.633,
                       1850311.938, 1831861.398, 1252752.586 };
            float[] y = new float[intensity.Length];    
            for (int i = 0; i < intensity.Length; i++)
            {
                y[i] = (float)intensity[i];
            }
            var filter = new SavitzkyGolayFilter(9);
            var signal = new DiscreteSignal(1, y);
            var smoothed = filter.ApplyTo(signal);
            var x = new double[intensity.Length];
            var extra = filter.Size / 2;
            var smoothedData = smoothed.Samples.Skip(extra).Take(smoothed.Samples.Length - extra*2).ToArray();
            for (int i = 0; i < intensity.Length; i++)
            {
                x[i] = i + 1;
            }
            var plot1 = Chart2D.Chart.Line<double, double, string>(
                        x: x,
                        y: intensity).WithTraceInfo($"raw", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("blue"));
            var plot2 = Chart2D.Chart.Line<double, float, string>(
                        x: x,
                        y: smoothedData).WithTraceInfo($"smoothed", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
            var combinedPlot = Chart.Combine(new[] { plot1, plot2 });
            combinedPlot.Show();
        }

        [Test]
        public static void TestCutPeak()
        {
            var path = @"E:\DIA\Data\DIA-bottom-up_241218\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            string tomlFile = @"E:\DIA\Data\DIA-bottom-up_241218\bu-ISD_mix1_120k\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 0, maxMass: 5000, type: "ISD", apexCycleTolerance: 5,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.ScanCycle,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc,
                cutMs1Peaks: true, cutMs2Peaks: true);
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, task.CommonParameters);

            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).Where(s => s.RetentionTime >= 26.5 && s.RetentionTime <= 27.5).ToArray();
            var allMs2PeakCurves = ISDEngine_static.GetAllPeakCurves(ms2Scans, task.CommonParameters, task.CommonParameters.DIAparameters,
                task.CommonParameters.DIAparameters.Ms2XICType, task.CommonParameters.DIAparameters.Ms2PeakFindingTolerance,
                task.CommonParameters.DIAparameters.MaxRTRangeMS2, out List<Peak>[] peaksByScan2, task.CommonParameters.DIAparameters.CutMs2Peaks);
            var allMs2PeakCurves_noCut = ISDEngine_static.GetAllPeakCurves(ms2Scans, task.CommonParameters, task.CommonParameters.DIAparameters,
                task.CommonParameters.DIAparameters.Ms2XICType, task.CommonParameters.DIAparameters.Ms2PeakFindingTolerance,
                task.CommonParameters.DIAparameters.MaxRTRangeMS2, out List<Peak>[] peaksByScan3, false);
            var allMs2PeakCurves_cutPeak = ISDEngine_static.GetAllPeakCurves(ms2Scans, task.CommonParameters, task.CommonParameters.DIAparameters,
            XICType.Peak_cutPeak, task.CommonParameters.DIAparameters.Ms2PeakFindingTolerance, task.CommonParameters.DIAparameters.MaxRTRangeMS2, 
            out List<Peak>[] peaksByScan2_cutPeak);

            var originalPeakNum = ms2Scans.Sum(s => s.MassSpectrum.Size);
            var peakNumAfterXIC = allMs2PeakCurves.Sum(pc => pc.Peaks.Count);
            var peakNumAfterXIC_noCut = allMs2PeakCurves_noCut.Sum(pc => pc.Peaks.Count);
            var peakNumAfterXIC_cutPeak = allMs2PeakCurves_cutPeak.Sum(pc => pc.Peaks.Count);
        }

        [Test]
        public static void TestSearchBottomUp()
        {
            var path = @"E:\DIA\Data\DIA-bottom-up_241218\12-18-24_bu-ISD100_5pro_mix1_labelCorrected_RT13.21-16.21.mzML";
            string tomlFile = @"E:\DIA\Data\DIA-bottom-up_241218\bu-ISD_noPeakTrim\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20),
                maxNumMissedScan: 1, binSize: 100, overlapRatioCutOff: 0, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 150, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, precursorIntensityCutOff: 300000, 
                splineTimeInterval: 0.005f, minMass: 0, maxMass: 5000, apexCycleTolerance: 15, scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, 
                ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.ScanCycle, pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISDEngine_static, 
                correlationType: CorrelationType.CubicSpline_scanCycle_preCalc, cutMs1Peaks: true, cutMs2Peaks: true);

            string myDatabase = @"E:\ISD Project\ISD_240606\idmapping_2024_06_11.xml";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\ISDEngine_static_BU-RT13.21-16.21_DeconHighestPeakMS1_PeakMS2_Tol5ppm_maxRT2_cutMs1Ms2Peak_maxMissed2_scanCycle_apexCycle15_overlap0_corr0.5";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { path }, "test");
        }

        

        [Test]
        public static void TestSimpleFakeData()
        {
            var precursorPeaks = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(10, i + 1, i + 1, ZeroBasedScanNumber: i);
                precursorPeaks.Add(peak);
            }
            var precursorPC = new PeakCurve(precursorPeaks);
            precursorPC.MonoisotopicMass = 50;
            precursorPC.Charge = 10;

            //Test PeakCurve properties
            Assert.That(precursorPC.AveragedMz == 10);
            Assert.That(precursorPC.AveragedIntensity == (double)91/21);
            Assert.That(precursorPC.ApexRT == 6);

            var fragmentPeaks1 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(2, i + 1.1, i + 1, ZeroBasedScanNumber: i);
                fragmentPeaks1.Add(peak);
            }
            var fragmentPC1 = new PeakCurve(fragmentPeaks1);
            fragmentPC1.MonoisotopicMass = 10;
            fragmentPC1.Charge = 5;

            var fragmentPeaks2 = new List<Peak>();
            for (int i = 0; i < 6; i++)
            {
                var peak = new Peak(5, i + 1.2, (i + 1)*1000, ZeroBasedScanNumber: i);
                fragmentPeaks2.Add(peak);
            }
            var fragmentPC2 = new PeakCurve(fragmentPeaks2);
            fragmentPC2.MonoisotopicMass = 10;
            fragmentPC2.Charge = 2;

            //Test correlation calculations
            var corr1 = PrecursorFragmentPair.CalculatePeakCurveCorr(precursorPC, fragmentPC2);
            precursorPC.GetScanCycleSmoothedData(0.05f);
            fragmentPC1.GetScanCycleSmoothedData(0.05f);
            fragmentPC2.GetScanCycleSmoothedData(0.05f);
            var corr2 = PrecursorFragmentPair.CalculateCorr_scanCycleSpline_preCalculated(precursorPC, fragmentPC2);
            var corr3 = PrecursorFragmentPair.CalculateCorr_spline(precursorPC, fragmentPC2, "cubic", 0.005);
            Assert.That(corr1 == 1);
            Assert.That(Math.Round(corr2, 2) == 1);
            Assert.That(Math.Round(corr3, 2) == 1);

            var pfPair1 = new PrecursorFragmentPair(precursorPC, fragmentPC1, 0.9);
            var pfPair2 = new PrecursorFragmentPair(precursorPC, fragmentPC2, 0.9);
            var pfGroup = new PrecursorFragmentsGroup(precursorPC);
            pfGroup.PFpairs.AddRange(new List<PrecursorFragmentPair> { pfPair1, pfPair2 });

            //Test make new ms2withspecificmass
            var ms2WithPre_mzPeak = ISDEngine_static.ConstructNewMs2Scans(pfGroup, new CommonParameters(), PseudoMs2ConstructionType.mzPeak, "dataFile");
            Assert.That(ms2WithPre_mzPeak.PrecursorMass == precursorPC.MonoisotopicMass);
            Assert.That(ms2WithPre_mzPeak.TheScan.MassSpectrum.XArray.Length == 2);
            Assert.That(ms2WithPre_mzPeak.TheScan.MassSpectrum.SumOfAllY == fragmentPC1.AveragedIntensity + fragmentPC2.AveragedIntensity);

            var ms2WithPre_mass = ISDEngine_static.ConstructNewMs2Scans(pfGroup, new CommonParameters(), PseudoMs2ConstructionType.neutralMass, "dataFile");
            Assert.That(ms2WithPre_mass.PrecursorMass == precursorPC.MonoisotopicMass);
            Assert.That(ms2WithPre_mass.ExperimentalFragments.Length == 2);
            Assert.That(ms2WithPre_mass.ExperimentalFragments[0].MonoisotopicMass == 10);
            Assert.That(ms2WithPre_mass.ExperimentalFragments[1].MonoisotopicMass == 10);
        }

        [Test]
        public static void TestRtMap()
        {
            var path = @"E:\DIA\Data\DIA-bottom-up_241218\12-18-24_bu-ISD100_5pro_mix1_labelCorrected.mzML";
            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, new CommonParameters());
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var ms1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToArray();
            var rtMap = ISDEngine_static.GetRtMap(ms1Scans, ms2Scans);
            Assert.That(rtMap.Count == ms2Scans.Length);
        }

        [Test]
        public static void TestMzSpectrumInitialization()
        {
            var spectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 10, 20, 30 }, false);
            
        }

    }
}
