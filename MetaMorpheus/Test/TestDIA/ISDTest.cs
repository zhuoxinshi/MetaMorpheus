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
using ThermoFisher.CommonCore.Data.Business;
using Easy.Common.Extensions;
using Readers;
using MathNet.Numerics.Statistics;

namespace Test.TestDIA
{
    public class ISDTest
    {
        [Test]
        public static void TestConstructMs2Groups()
        {
            string filePath1 = @"E:\ISD Project\TestIsdDataAnalysis\data\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string filePath2 = @"E:\ISD Project\Data for test\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100-calib-averaged.mzML";
            string snip = @"E:\ISD Project\TestIsdDataAnalysis\data\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_2mzstep_relativeToTics(1)_RT28.01-32.69.mzML";

            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            //var file = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD.mzML";

            string myDatabase = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task3-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\isdEngine_FD60only-RT29.03-33.52_corr0.5_highestPeakXIC_ms1Tol10ppm_apexRT0.3_maxMissed2_overlap0.3";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, cutPeaks: false);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { snip }, "test");
        }

        [Test]
        public static void TestEnvelopeXIC()
        {
            string snip = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT28.01-32.69.mzML";
            var myFileManagers = new MyFileManager(true);
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD");
            var myMsDataFile = myFileManagers.LoadFile(snip, task.CommonParameters);
            var isdEngine = new ISDEngine(myMsDataFile, task.CommonParameters, task.CommonParameters.DIAparameters);
            isdEngine.Ms1PeakIndexing();
            isdEngine.ConstructMs2Group();
            isdEngine.GetMs1PeakEnvelopeCurve();  
            //isdEngine.GetMs1PeakCurves();
            //isdEngine.GetPrecursorEnvelopeXIC2();
            //isdEngine.GetEnvelopeXICForMs2();
            //isdEngine.PrecursorFragmentPairing();
        }

        [Test]
        public static void TestDecon()
        {
            string file = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100.mzML";
            var myFileManagers = new MyFileManager(true);
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            var myMsDataFile = myFileManagers.LoadFile(file, task.CommonParameters);
            //var scan = myMsDataFile.GetOneBasedScan(1330);
            //var envelopes = Deconvoluter.Deconvolute(scan, task.CommonParameters.PrecursorDeconvolutionParameters).OrderBy(e => e.MostAbundantObservedIsotopicMass.ToMz(e.Charge)).ToList();

        }

        [Test]
        public static void IsolationRange()
        {
            string file = @"E:\DIA\Data\DIA_241108\11-11-24_DDA_5pro_65min_40-45-50stpHCD.raw";
            var myFileManagers = new MyFileManager(true);
            var myMsDataFile = myFileManagers.LoadFile(file, new CommonParameters());
            var ms2scans = myMsDataFile.GetAllScansList().Where(p => p.MsnOrder == 2).ToList();
            var isoRange = new List<MzRange>();
            var isoMz = new List<double>();
            var isolationWindow_5mz = new List<(double, double)>();
            var isolationWindow_10mz = new List<(double, double)>();
            var isolationWindow_30mz = new List<(double, double)>();
            var selectedScans = ms2scans.OrderByDescending(p => p.SelectedIonIntensity).Take(100).ToList();
            foreach (var scan in selectedScans)
            {
                var range = scan.IsolationRange;
                isoMz.Add(scan.IsolationMz.Value);
                isoRange.Add(range);
                isolationWindow_5mz.Add((Math.Round(scan.IsolationMz.Value - 2.5, 2), Math.Round(scan.IsolationMz.Value + 2.5, 2)));
                isolationWindow_10mz.Add((Math.Round(scan.IsolationMz.Value - 5, 2), Math.Round(scan.IsolationMz.Value + 5, 2)));
                isolationWindow_30mz.Add((Math.Round(scan.IsolationMz.Value - 15, 2), Math.Round(scan.IsolationMz.Value + 15, 2)));
            }
            var output5mz = @"B:\Users\Zhuoxin\DIA\5mzIso_1.csv";
            var output10mz = @"B:\Users\Zhuoxin\DIA\10mzIso_1.csv";
            var output30mz = @"B:\Users\Zhuoxin\DIA\30mzIso_1.csv";
            isolationWindow_5mz = isolationWindow_5mz.Distinct().ToList();
            isolationWindow_10mz = isolationWindow_10mz.Distinct().ToList();
            isolationWindow_30mz = isolationWindow_30mz.Distinct().ToList();
            using (var sw = new StreamWriter(File.Create(output5mz)))
            {
                sw.WriteLine("Calculated m/z Window");
                foreach (var range in isolationWindow_5mz.OrderBy(w => w.Item1))
                {
                    sw.WriteLine($"{Math.Round(range.Item1,1)}-{Math.Round(range.Item2,1)}");
                }
            }
            using (var sw = new StreamWriter(File.Create(output10mz)))
            {
                sw.WriteLine("Calculated m/z Window");
                foreach (var range in isolationWindow_10mz.OrderBy(w => w.Item1))
                {
                    sw.WriteLine($"{Math.Round(range.Item1, 1)}-{Math.Round(range.Item2, 1)}");
                }
            }
            using (var sw = new StreamWriter(File.Create(output30mz)))
            {
                sw.WriteLine("Calculated m/z Window");
                foreach (var range in isolationWindow_30mz.OrderBy(w => w.Item1))
                {
                    sw.WriteLine($"{Math.Round(range.Item1, 1)}-{Math.Round(range.Item2, 1)}");
                }
            }
        }

        [Test]
        public static void TestReadFeatureFile()
        {
            var file = @"E:\toppic-windows-1.6.5\ISD\test1\id_08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_MS1_ms1.feature";
            var featureFile = new Ms1FeatureFile(file);
            var feature = featureFile.First();
            featureFile.LoadResults();
        }

        [Test]
        public static void XICstandards()
        {
            var psmFileISD = @"E:\ISD Project\ISD_240606\search_w-writingSpectralLib\Task1-SearchTask\AllPSMs.psmtsv";
            var psmsISD = PsmTsvReader.ReadTsv(psmFileISD, out List<string> warnings).Where(psm => psm.QValue < 0.01);
            var sequences = psmsISD.GroupBy(psm => psm.FullSequence).ToList();
            var ubISD = sequences[3];
            var matchedIonsListISD = ubISD.Select(psm => psm.MatchedIons).ToList();
            var allMatchedIonsISD = matchedIonsListISD.SelectMany(i => i).ToList();
            var ionsGroupByLabel_ISD = allMatchedIonsISD.GroupBy(i => i.Annotation).Select(g => g.Key).ToList();
            var sd_ISD = Statistics.StandardDeviation(allMatchedIonsISD.Select(i => i.Intensity));
            var filteredIons_intensity = allMatchedIonsISD.Where(i => i.Intensity > 50000).ToList();

            var psmFileDDA = @"E:\ISD Project\ISD_240606\sample10-DDA_w-writingSpectralLib\Task1-SearchTask\AllPSMs.psmtsv";
            var psmsDDA = PsmTsvReader.ReadTsv(psmFileDDA, out List<string> warningsDDA).Where(psm => psm.QValue < 0.01);
            var sequences_DDA = psmsDDA.GroupBy(psm => psm.FullSequence).ToList();
            var ubDDA = sequences_DDA[0];
            var matchedIonsListDDA = ubDDA.Select(psm => psm.MatchedIons).ToList();
            var allMatchedIonsDDA = matchedIonsListDDA.SelectMany(i => i).ToList();
            var ionsGroupByMass_DDA = allMatchedIonsDDA.GroupBy(i => i.Annotation).Select(g => g.Key).ToList();
            var sd_DDA = Statistics.StandardDeviation(allMatchedIonsDDA.Select(i => i.Intensity));
            var filteredIons_intensity_DDA = allMatchedIonsDDA.Where(i => i.Intensity > 50000).ToList();

            var overlapIons = ionsGroupByLabel_ISD.Where(i => ionsGroupByMass_DDA.Contains(i)).ToList();
            //there are fragments with highly correlated XIC but not matched
        }

        [Test]
        public static void TestGrouping()
        {
            string snip = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT28.01-32.69.mzML";
            var myFileManagers = new MyFileManager(true);
            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD");
            var myMsDataFile = myFileManagers.LoadFile(snip, task.CommonParameters);

            var isdEngine = new ISDEngine(myMsDataFile, task.CommonParameters, task.CommonParameters.DIAparameters);
            isdEngine.Ms1PeakIndexing();
            isdEngine.ConstructMs2Group();
            isdEngine.GetMs1PeakCurves();
            isdEngine.GetMs2PeakCurves();
            var allMs2PeakCurves = isdEngine.Ms2PeakCurves.Values.SelectMany(p => p).ToList();
            var testMz = 857.47;
            var testPeak = PeakCurve.GetPeakFromScan(testMz, isdEngine.Ms1PeakTable, 0, new PpmTolerance(5), 100);
            //var testPC = PeakCurve.FindPeakCurve(testPeak, isdEngine.Ms1PeakTable, ms1scans, null, 2, new PpmTolerance(5), 100, 1.5);
            //var testGroup = ISDEngine.GroupPrecursorFragments(testPC, allMs2PeakCurves, task.CommonParameters.DIAparameters);
            //int match = 0;
            //foreach (var ion in matchedIons)
            //{
            //    if contains
            //    match += 1
            //}

            isdEngine.GetPrecursorEnvelopeXIC2();
            isdEngine.GetEnvelopeXICForMs2();
            isdEngine.PrecursorFragmentPairing();
        }

        [Test]
        public static void TestUseToppicFeatures()
        {
            string filePath1 = @"E:\ISD Project\TestIsdDataAnalysis\data\06-11-24_mix_sample1_2uL_ISD_RT45.01-48.09.mzML";
            string filePath2 = @"E:\ISD Project\Data for test\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100-calib-averaged.mzML";
            string snip = @"E:\ISD Project\TestIsdDataAnalysis\data\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_2mzstep_relativeToTics(1)_RT28.01-32.69.mzML";

            string tomlFile = @"E:\ISD Project\ISD_240606\2024-10-24-15-44-25\Task Settings\Task1-SearchTaskconfig.toml";
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            //var file = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD.mzML";

            string myDatabase = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task3-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            //string library = @"E:\ISD Project\TestIsdDataAnalysis\SpectralLibraryDDA\Task1-SearchTask\SpectralLibrary_2024-07-09-17-24-30.msp";
            DbForTask db = new DbForTask(myDatabase, false);
            string outputFolder = @"E:\ISD Project\TestSearch\TopFD_overlapRT0.3_apex0.3_overlapArea0.2";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }

            //task.CommonParameters.GetType().GetProperty("MaxThreadsToUsePerFile").SetMethod.Invoke(task.CommonParameters, new object[] { 1 });
            //var type = task.CommonParameters.GetType();
            //var property = type.GetProperty("MaxThreadsToUsePerFile");
            //property.SetMethod.Invoke(task.CommonParameters, new object[] { 1 });

            //task.CommonParameters.ProductMassTolerance = new PpmTolerance(50);
            //task.CommonParameters.PrecursorMassTolerance = new PpmTolerance(50);
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(5),
                maxNumMissedScan: 2, binSize: 100, overlapRatioCutOff: 0.3, correlationCutOff: 0.5, apexRtTolerance: 0.3,
                fragmentRankCutOff: 2000, precursorRankCutOff: 10, maxRTrangeMS1: 2, maxRTrangeMS2: 2, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f, minMass: 12000, type: "ISD", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, cutPeaks: false);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { snip }, "test");
        }
    }
}
