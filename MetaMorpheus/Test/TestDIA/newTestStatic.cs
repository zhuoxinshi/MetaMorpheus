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
using MathNet.Numerics.Statistics;
using static System.Net.WebRequestMethods;

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
                scanCycleSplineInterval: 0.005, cutPeaks: false, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.Area_correlation,
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
                scanCycleSplineInterval: 0.005, cutPeaks: false, analysisType:AnalysisType.ISDEngine_static);

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
                scanCycleSplineInterval: 0.005, cutPeaks: false, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.Area_correlation,
            pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.ISD_scanBased, correlationType: CorrelationType.CubicSpline_scanCycle_preCalc);

            //check if target ions are included in corresponding pseudoMS2scans
            var pseudoScans = ISDEngine_static.GetPseudoMs2Scans(dataFile, task.CommonParameters, task.CommonParameters.DIAparameters);
            var scans = pseudoScans.Where(scan => Math.Round(scan.PrecursorMass, 0) == Math.Round(massToTest, 0)).ToList();

            var scores = new List<double>();
            var notFound = new Dictionary<int, List<(string, double, double, double)>> ();
            foreach (var scan in scans)
            {
                double score = 0;
                var no = new List<(string, double, double, double)>();
                foreach (var ion in targetList)
                {
                    var fragments = scan.ExperimentalFragments.Where(e => Math.Round(e.MonoisotopicMass) == ion.mass && e.Charge == ion.charge);
                    if (!fragments.Any())
                    {
                        no.Add(new (ion.annotation, ion.mass, ion.charge, ion.mass.ToMz(ion.charge)));
                    }
                    else
                    {
                        score++;
                    }
                }
                scores.Add(score/(double)targetList.Count);
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
                scanCycleSplineInterval: 0.005, cutPeaks: false, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.ScanCycle,
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
                scanCycleSplineInterval: 0.005, cutPeaks: false, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.DeconHighestPeak, pfGroupingType: PFGroupingType.ScanCycle,
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
                scanCycleSplineInterval: 0.005, cutPeaks: false, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.OverlapFirst,
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
            foreach(var scan in isd100_orig)
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
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(path, task.CommonParameters);
            var allScans = dataFile.GetAllScansList();
            var ms1scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2scans = allScans.Where(s => s.MsnOrder == 2).ToArray();
            var psmsFile = @"E:\ISD Project\CE_241213\snip_highestPeak\Task1-SearchTask\AllPSMs.psmtsv";
            var psms = PsmTsvReader.ReadTsv(psmsFile, out List<string> warnings);
            var psmsContainTheIon = psms.Where(p => p.MatchedIons.Select(i => Math.Round(i.Mz, 2)).Contains(1017.07)).ToList();
            var highestPeakMzs = psmsContainTheIon.Select(psm => psm.PrecursorHighestPeakMz).GroupBy(p => Math.Round(p, 2)).OrderBy(p => p.Key).ToList();

            var pc1017 = PeakCurve.PeakTracing(1017.58, 1, ms2scans, new PpmTolerance(10), 100, 2, 8);
        }

    }
}
