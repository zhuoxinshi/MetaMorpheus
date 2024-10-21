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

namespace Test.TestDIA
{
    public class Random
    {
        [Test]
        public static void Ms1DeconTest()
        {
            SearchTask task = new SearchTask();
            string outputFolder = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string DIAfile = @"E:\DIA\FragPipe\DIA\snip_data\Extract_Isolation_Window.mzML";
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            //Directory.CreateDirectory(outputFolder);
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
            
            MyFileManager myFileManager = new MyFileManager(true);
            var myMsDataFile = myFileManager.LoadFile(DIAfile, task.CommonParameters);
            var allMs1Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
            var allEnvelopes = new List<IsotopicEnvelope>[allMs1Scans.Length];
            var allPrecursors = new List<DeconvolutedMass>();
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], task.CommonParameters.PrecursorDeconvolutionParameters);
                foreach (var envelope in envelopes)
                {
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new DeconvolutedMass(envelope, charge, allMs1Scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
                    allPrecursors.Add(precursor);
                }
            }
            var pre = allPrecursors.GroupBy(p => new {mass = Math.Round(p.MonoisotopicMass, 2), p.Charge }).ToList();
            var pre2 = allPrecursors.GroupBy(p => Math.Round(p.MonoisotopicMass, 2)).ToList();
            var allPrecursorClusters = new List<PrecursorCluster>();
            foreach(var group in pre)
            {
                var precursors = group.ToList();
                double totalIntensity = group.ToList().Sum(p => p.Envelope.TotalIntensity);
                var newPC = new PrecursorCluster(group.Key.mass, group.Key.Charge, totalIntensity, precursors);
                allPrecursorClusters.Add(newPC);
            }
            var filteredPC = allPrecursorClusters.OrderByDescending(pc => pc.EnvelopeCount >= 5).ToList();
            var XICs = new List<(double rt, double intensity, double mass)>();
            var lfqPeaks = new List<(string fileName, double, int, string, string, double, string)>();
            List<string> myOutput = new();
            foreach(var cluster in allPrecursorClusters)
            {
                myOutput.Add("SmallCalibratible_Yeast" + "\t" + cluster.ApexRT + "\t" + cluster.Charge + "\t" + "" + "\t" + "" + "\t" + cluster.MonoisotopicMass + "\t" + "");
                lfqPeaks.Add(("SmallCalibratible_Yeast", cluster.ApexRT, cluster.Charge, "", "", cluster.MonoisotopicMass, ""));
            }

            System.IO.File.WriteAllLines(@"E:\DIA\output.tsv", myOutput);
            //var masses = pre.Select(p => p.Key);
            //Directory.CreateDirectory(outputFolder);
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
        }

        [Test]
        public void TestGetAllMs1PeakCurves()
        {
            var diaFile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var commonParameters = new CommonParameters();
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 2, 100, 0.5, 0.5, 1);
            var diaEngine = new DIAEngine(diaDataFile, commonParameters, diaParam);
            diaEngine.Ms1PeakIndexing();
            diaEngine.GetMs1PeakCurves();
            diaEngine.ConstructMs2Group();
            diaEngine.GetMs2PeakCurves();
            diaEngine.PrecursorFragmentPairing();
            diaEngine.ConstructNewMs2Scans();
            var pfgroup = diaEngine.PFgroups[306];
            var list = new List<(double rt, double intensity, double corr)>();
            int index = 0;
            foreach(var peak in pfgroup.PrecursorPeakCurve.Peaks)
            {
                list.Add((peak.RetentionTime, peak.Intensity, 1));
            }
            index++;
            foreach(var pair in pfgroup.PFpairs)
            {
                foreach(var fragPeak in pair.FragmentPeakCurve.Peaks)
                {
                    list.Add((fragPeak.RetentionTime, fragPeak.Intensity, pair.Correlation));
                }
                index++;
            }

            string path = @"E:\DIA\Visualization\pfXIC.csv";
            if(File.Exists(path))
            {
                File.Delete(path);
            }
            using (var sw = new StreamWriter(File.Create(path)))
            {
                sw.WriteLine("rt, intensity, corr");
                foreach (var line in list)
                {
                    sw.WriteLine($"{line.rt}, {line.intensity}, {line.corr}");
                }
            }
        }

        [Test]
        public void TestUnrelated()
        {
            var fd_dda_path = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_DDA.raw";
            MyFileManager myFileManager = new MyFileManager(true);
            var commonParameters = new CommonParameters();
            var fd_DDA = myFileManager.LoadFile(fd_dda_path, commonParameters);
            var fd_DDA_ms2 = fd_DDA.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            var isolationWindows = fd_DDA_ms2.Select(s => s.IsolationRange).ToList();

            var proteoforms_dda_path = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task4-SearchTask\Individual File Results\08-12-24_PEPPI_FractionD_orbiMS1_DDA-calib-averaged_Proteoforms.psmtsv";
            var proteoforms_DDA = PsmTsvReader.ReadTsv(proteoforms_dda_path, out List<string> errors).Where(p => p.QValue < 0.01 && p.QValueNotch < 0.01 && p.DecoyContamTarget == "T").ToList();
            var proteoforms_isd_path = @"E:\ISD Project\ISD_240812\FB-FD_lessGPTMD\Task4-SearchTask\Individual File Results\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100-calib-averaged_Proteoforms.psmtsv";
            var proteoforms_ISD = PsmTsvReader.ReadTsv(proteoforms_isd_path, out errors).Where(p => p.QValue < 0.01 && p.QValueNotch < 0.01 && p.DecoyContamTarget == "T"). ToList();
            var uniqueISD = new List<PsmFromTsv>();
            var allseqInDDA = proteoforms_DDA.Select(p => p.FullSequence).ToList();
            foreach(var proteoform in proteoforms_ISD)
            {
                if (!allseqInDDA.Contains(proteoform.FullSequence))
                {
                    uniqueISD.Add(proteoform);
                }
            }

            var unselected = new List<PsmFromTsv>();
            var selected = new List<PsmFromTsv>();
            foreach(var id in uniqueISD)
            {
                var rtRange_start = id.RetentionTime - 0.5;
                var rtRange_end = id.RetentionTime + 0.5;
                var ddaScansInRange = fd_DDA_ms2.Where(s => s.RetentionTime >= rtRange_start && s.RetentionTime <= rtRange_end).ToList();
                foreach(var scan in ddaScansInRange)
                {
                    if (id.PrecursorMz >= scan.IsolationRange.Minimum && id.PrecursorMz <= scan.IsolationRange.Maximum)
                    {
                        selected.Add(id);
                        break;
                    }
                }
            }
            var uniqueISDpreMz = uniqueISD.Select(p => p.PrecursorMz).ToList();
        }

        [Test]
        public static void TestXICplot()
        {
            int bin = 100;
            string snipFile = @"E:\DIA\FragPipe\DIA\snip_data\Extract_Isolation_Window.mzML";
            MyFileManager myFileManager = new MyFileManager(true);
            var commonParameters = new CommonParameters();
            var snipData = myFileManager.LoadFile(snipFile, commonParameters);
            var ms1scans = snipData.GetMS1Scans().ToList();
            var allEnvelopes = new List<IsotopicEnvelope>[ms1scans.Count];
            for (int i = 0; i < ms1scans.Count; i++)
            {
                allEnvelopes[i] = new List<IsotopicEnvelope>();
                var envelopes = Deconvoluter.Deconvolute(ms1scans[i], commonParameters.PrecursorDeconvolutionParameters);
                allEnvelopes[i].AddRange(envelopes);
            }
            var peaks = allEnvelopes.SelectMany(e => e).Select(e => e.Peaks).SelectMany(p => p).ToList();
            var mzList = new List<double>();
            foreach(var peak in peaks)
            {
                var roundedMz = Math.Round(peak.mz, 2);
                if (!mzList.Contains(roundedMz))
                {
                    mzList.Add(roundedMz);
                }
            }
            mzList = mzList.OrderBy(m => m).ToList();
            var tol = new PpmTolerance(10);
            var lines = new List<List<double>>();
            foreach(var ms1 in ms1scans)
            {
                var mzs = new List<double>();
                for (int i = 0; i < mzList.Count; i++)
                {
                    int index = ms1.MassSpectrum.GetClosestPeakIndex(mzList[i]);
                    if (tol.Within(mzList[i], ms1.MassSpectrum.XArray[index]))
                    {
                        mzs.Add(ms1.MassSpectrum.YArray[index]);
                    }
                    else
                    {
                        mzs.Add(0);
                    }
                }
                lines.Add(mzs);
            }
            lines = lines.Select(l => l.ToList()).ToList();
            var path = @"E:\DIA\FragPipe\DIA\snip_data\clustering.csv";
            using(var sw = new StreamWriter(File.Create(path)))
            {
                sw.WriteLine(string.Join(",", mzList));
                foreach(var line in lines)
                {
                    sw.WriteLine(string.Join(",", line));
                }
            }
        }

        [Test]
        public static void TestSpline()
        {
            var diaFile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var commonParameters = new CommonParameters();
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 2, 100, 0.5, 0.5,2);
            var diaEngine = new DIAEngine(diaDataFile, commonParameters, diaParam);
            diaEngine.Ms1PeakIndexing();
            diaEngine.GetMs1PeakCurves();
            var allPC = diaEngine.Ms1PeakCurves;
            var samplePC1 = allPC[49];
            var samplePC2 = allPC[50];
            samplePC1.GetCubicSpline();
            samplePC2.GetCubicSpline();
            var startRT = Math.Max(samplePC1.StartRT, samplePC2.StartRT);
            var endRT = Math.Min(samplePC1.EndRT, samplePC2.EndRT);

            //plot
            var rtSeq = new List<double>();
            var intensities1 = new List<double>();
            var intensities2 = new List<double>();
            for (double i = startRT; i < endRT; i += 0.005)
            {
                rtSeq.Add(i);
            }
            foreach (var rt in rtSeq)
            {
                intensities1.Add(samplePC1.CubicSpline.Interpolate(rt));
                intensities2.Add(samplePC2.CubicSpline.Interpolate(rt));
            }
            var plot1 = Chart2D.Chart.Point<double, double, string>(
                x: samplePC1.Peaks.Select(p => p.RetentionTime),
                y: samplePC1.Peaks.Select(p => p.Intensity)).WithTraceInfo("original1").WithMarkerStyle(Color: Color.fromString("red"));
            var plot2 = Chart2D.Chart.Point<double, double, string>(
                x: samplePC2.Peaks.Select(p => p.RetentionTime),
                y: samplePC2.Peaks.Select(p => p.Intensity)).WithTraceInfo("original2").WithMarkerStyle(Color: Color.fromString("blue"));
            var plot3 = Chart2D.Chart.Point<double, double, string>(
                x: rtSeq,
                y: intensities1).WithTraceInfo("interpolate1").WithMarkerStyle(Color: Color.fromString("red"));
            var plot4 = Chart2D.Chart.Point<double, double, string>(
                x: rtSeq,
                y: intensities2).WithTraceInfo("interpolate2").WithMarkerStyle(Color: Color.fromString("blue"));
            var combinedPlot1 = Chart.Combine(new[] { plot1, plot2 });
            var combinedPlot2 = Chart.Combine(new[] { plot3, plot4 });
            combinedPlot1.Show();
            combinedPlot2.Show();
            int stop = 0;
        }

        [Test]
        public static void TestNumberOfPeaks()
        {
            var diaFile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var commonParameters = new CommonParameters();
            commonParameters.TrimMsMsPeaks = true;
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 2, 100, 0.5, 0.5, 2);
            var diaEngine = new DIAEngine(diaDataFile, commonParameters, diaParam);
            var allms2scans = diaDataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            var allms1scans = diaDataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToList();
            var numMs2Peaks = allms2scans.Sum(s => s.MassSpectrum.Size);
            var numMs1Peaks = allms1scans.Sum(s => s.MassSpectrum.Size);
            diaEngine.Ms1PeakIndexing();
            diaEngine.GetMs1PeakCurves();
            diaEngine.ConstructMs2Group();
            diaEngine.GetMs2PeakCurves();
            var diaMs2Peaks = diaEngine.Ms2PeakCurves.Values.SelectMany(pc => pc).Sum(pc => pc.Peaks.Count);
        }

        [Test]
        public static void TestMs1PeakCurveDetection()
        {
            var diaFile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var commonParameters = new CommonParameters();
            commonParameters.TrimMsMsPeaks = true;
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 1, 100, 0.5, 0.5, 2);

            var ms1scans = diaDataFile.GetMS1Scans().ToArray();
            var allMs1Peaks = Peak.GetAllPeaks(ms1scans, diaParam.PeakSearchBinSize);
            var rankedMs1Peaks = allMs1Peaks.OrderByDescending(p => p.Intensity).ToList();
            var ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, diaParam.PeakSearchBinSize);
            var ms1PeakCurves = new List<PeakCurve>();
            foreach (var peak in rankedMs1Peaks)
            {
                if (peak.PeakCurve == null)
                {
                    var newPeakCurve = PeakCurve.FindPeakCurve(peak, ms1PeakTable, ms1scans, ms1scans[0].IsolationRange,
                        diaParam.MaxNumMissedScan, diaParam.Ms2PeakFindingTolerance, diaParam.PeakSearchBinSize, 2);
                    if (newPeakCurve.Peaks.Count > 0)
                    {
                        ms1PeakCurves.Add(newPeakCurve);
                    }
                }
            }
            //check
            var numOfPeaks = ms1PeakCurves.Sum(pc => pc.Peaks.Count);
        }

        [Test]
        public static void VisualizeScanXIC()
        {
            var diaFile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var commonParameters = new CommonParameters();
            commonParameters.TrimMsMsPeaks = false;
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 1, 100, 0.2, 0.5, 0.25, 100, 10, 2);
            var diaEngine2 = new DIAEngine2(diaDataFile, commonParameters, diaParam);
            diaEngine2.Ms1PeakIndexing();
            diaEngine2.ConstructMs2Group();
            diaEngine2.GetMs1PeakCurves();
            diaEngine2.GetMs2PeakCurves();
            //diaEngine2.PrecursorFragmentPairing();
            //diaEngine2.ConstructNewMs2Scans();

            string pepFile = @"E:\DIA\TestSearch\test2.0_corr0.5_noGroup_highestPeakXIC_cubicSpline_apexRT0.25_noPeakTrim_maxMissed1_overlap0.2_FragRank100_preRank10_maxRT2\AllPeptides.psmtsv";
            var psms = PsmTsvReader.ReadTsv(pepFile, out List<string> warnings);
            var seq = "";
            int scanNumber = 207;
            var precursorPeakCurve = diaEngine2.Ms1PeakCurves.Values.SelectMany(v => v).Where(pc => pc.Index == scanNumber).First();
            var ms2curves = diaEngine2.Ms2PeakCurves.Where(d => precursorPeakCurve.AveragedMz >= d.Key.min && precursorPeakCurve.AveragedMz <= d.Key.max)
                .First().Value;
            var pfgroup = DIAEngine2.GroupPrecursorFragments(precursorPeakCurve, ms2curves, diaParam);
            var psm = psms.Where(psm => psm.Ms2ScanNumber == scanNumber).First();
            var matchedIonMzs = psm.MatchedIons.Select(i => i.Mz);
            var matchedIonIntensities = psm.MatchedIons.Select(i => i.Intensity);   
            //var fragmentPCs = diaEngine2.Ms2PeakCurves.Values.SelectMany(v => v).Where(pc => matchedIonMzs.Contains(pc.AveragedMz)).ToList();

            var plots = new List<GenericChart>();
            var normalizedIntensity = precursorPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
            var precursorPlot = Chart2D.Chart.Line<double, double, string>(
                x: precursorPeakCurve.Peaks.Select(p => p.RetentionTime),
                y: normalizedIntensity).WithTraceInfo("precursor").WithMarkerStyle(Color: Color.fromString("red"));
            plots.Add(precursorPlot);
            var matchedIonPeakCurves = new List<PeakCurve>();
            foreach (double mz in matchedIonMzs)
            {
                //write a function, given the precursor peakCurve index, find the corresponding fragment peakCurve
                var fragmentPeakCurve = pfgroup.PFpairs.Where(p => Math.Abs(p.FragmentPeakCurve.AveragedMz - mz) < 0.005).First().FragmentPeakCurve;
                matchedIonPeakCurves.Add(fragmentPeakCurve);
                var norm = fragmentPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity)).ToList();
                var fragmentPlot = Chart2D.Chart.Line<double, double, string>(
                    x: fragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                    y: norm).WithTraceInfo($"fragment_{fragmentPeakCurve.AveragedMz}").WithMarkerStyle(Color: Color.fromString("blue"));
                plots.Add(fragmentPlot);
            }
            var combinedPlot = Chart.Combine(plots);
            combinedPlot.Show();

            //spline plots
            var splinePlots = new List<GenericChart>();
            var allPCs = new List<PeakCurve> { precursorPeakCurve};
            allPCs.AddRange(matchedIonPeakCurves);
            var startRT = allPCs.Max(pc => pc.StartRT);
            var endRT = allPCs.Min(pc => pc.EndRT);
            var rtSeq = new List<double>();
            for (double i = startRT; i < endRT; i += 0.005)
            {
                rtSeq.Add(i);
            }
            //precursor spline
            precursorPeakCurve.GetCubicSpline();
            var preIntensities = rtSeq.Select(rt => precursorPeakCurve.CubicSpline.Interpolate(rt));
            var pre_norm = preIntensities.Select(i => Math.Log10(i));
            var prePlot = Chart2D.Chart.Point<double, double, string>(
                                   x: rtSeq,
                                   y: pre_norm).WithTraceInfo("precursor").WithMarkerStyle(Color: Color.fromString("red"));
            splinePlots.Add(prePlot);
            foreach (var pc in matchedIonPeakCurves)
            {
                pc.GetCubicSpline();
                var intensities = rtSeq.Select(rt => pc.CubicSpline.Interpolate(rt));
                var normalized = intensities.Select(i => Math.Log10(i));
                var plot_spline = Chart2D.Chart.Point<double, double, string>(
                                       x: rtSeq,
                                       y: normalized).WithTraceInfo("fragment").WithMarkerStyle(Color: Color.fromString("blue"));
                splinePlots.Add(plot_spline);
            }
            var combinedSpline = Chart.Combine(splinePlots);
            combinedSpline.Show();

            int stop = 0;

            //all fragment PeakCurves
            var allFragmentPCs = pfgroup.PFpairs.Select(p => p.FragmentPeakCurve).ToList();
            var allplots = new List<GenericChart>();
            allplots.AddRange(plots);
            foreach(var pc in allFragmentPCs)
            {
                if (matchedIonPeakCurves.Contains(pc))
                {
                    continue;
                }
                var norm_all = pc.Peaks.Select(p => Math.Log10(p.Intensity)).ToList();
                var fragmentPlot = Chart2D.Chart.Line<double, double, string>(
                                       x: pc.Peaks.Select(p => p.RetentionTime),
                                       y: norm_all).WithTraceInfo($"fragment_{pc.AveragedMz}").WithMarkerStyle(Color: Color.fromString("green"));
                allplots.Add(fragmentPlot);
            }   
            var all = Chart.Combine(allplots);
            all.Show();

            int stop2 = 0;
        }

        [Test]
        public static void Visualize()
        {
            var diaFile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var commonParameters = new CommonParameters();
            commonParameters.TrimMsMsPeaks = false;
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20), 1, 100, 0.2, 0.5, 0.2, 100, 10, 0.5, precursorIntensityCutOff: 300000);

            string pepFile = @"E:\DIA\TestSearch\test2.0_corr0.5_highestPeakXIC_ms1Tol5ppm_cubicSpline_apexRT0.2_noPeakTrim_maxMissed0.5_overlap0.2_FragRank100_preRank10_maxRT0.5_300000\AllPeptides.psmtsv";
            var allpsms = PsmTsvReader.ReadTsv(pepFile, out List<string> warnings).Where(psm => psm.QValue < 0.01);
            var targetPsms = allpsms.Where(p => p.DecoyContamTarget == "T").ToList();  
            var decoyPsms = allpsms.Where(p => p.DecoyContamTarget == "D").ToList();
            var Tpsms = targetPsms.Where((item, index) =>
                (index >= 1 && index <= 10) ||
                (index >= 1000 && index <= 1010) ||
                (index >= 2000 && index <= 2010)).ToList();
            var Dpsms = decoyPsms.Where((item, index) =>
                (index >= 1 && index <= 10)).ToList();
            var listPsmsToLook = Tpsms.Concat(decoyPsms).ToList();
            var charts = PrecursorFragmentsGroup.VisualizePFgroups(diaDataFile, listPsmsToLook, commonParameters, diaParam);
            foreach(var chart in charts)
            {
                chart.Show();
            }
        }

        [Test]
        public static void TestDatabaseSize()
        {
            var fastaGPTMD_isd = @"E:\ISD Project\ISD_240927\allISD_oldMM_cali-avg-gptmd-fasta\Task3-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02GPTMD.xml";
            var xmlGPTMD_isd = @"E:\ISD Project\ISD_240927\allISD_oldMM_cali-avg-gptmd_xml\Task3-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            var xml = ProteinDbLoader.LoadProteinXML(xmlGPTMD_isd, true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var ok);
            var fasta = ProteinDbLoader.LoadProteinXML(fastaGPTMD_isd, true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var ok2);
        }

    }
}
