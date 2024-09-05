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
using System.IO;

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
            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            //Directory.CreateDirectory(outputFolder);
            //task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
            
            MyFileManager myFileManager = new MyFileManager(true);
            var myMsDataFile = myFileManager.LoadFile(DIAfile, task.CommonParameters);
            var allMs1Scans = myMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).ToArray();
            var allEnvelopes = new List<IsotopicEnvelope>[allMs1Scans.Length];
            var allPrecursors = new List<Precursor>();
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], task.CommonParameters.PrecursorDeconvolutionParameters);
                foreach (var envelope in envelopes)
                {
                    var charge = envelope.Charge;
                    double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
                    double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
                    var precursor = new Precursor(envelope, charge, allMs1Scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
                        allMs1Scans[i].OneBasedScanNumber, i);
                    allPrecursors.Add(precursor);
                }
            }
            var pre = allPrecursors.GroupBy(p => new {mass = Math.Round(p.MonoisotopicMass, 2), p.Charge });
            var allPrecursorClusters = new List<PrecursorCluster>();
            foreach(var group in pre)
            {
                var precursors = group.ToList();
                double totalIntensity = group.ToList().Sum(p => p.Envelope.TotalIntensity);
                var newPC = new PrecursorCluster(precursors, group.Key.mass, group.Key.Charge, group.ToList().Count, totalIntensity);
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
            var diaParam = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 2, 100, 0.5, 0.5);
            var diaEngine = new DIAEngine(diaDataFile, commonParameters, diaParam);
            diaEngine.Ms1PeakIndexing();
            diaEngine.GetMs1PeakCurves();
            diaEngine.ConstructMs2Group();
            diaEngine.GetMs2PeakCurves();
            diaEngine.PrecursorFragmentPairing();
            diaEngine.ConstructNewMs2Scans();
            var pfgroup = diaEngine.PFgroups[1242];
            var list = new List<(double rt, double intensity, int index)>();
            int index = 0;
            foreach(var peak in pfgroup.PrecursorPeakCurve.Peaks)
            {
                list.Add((peak.RetentionTime, peak.Intensity, index));
            }
            index++;
            foreach(var pair in pfgroup.PFpairs)
            {
                foreach(var fragPeak in pair.FragmentPeakCurve.Peaks)
                {
                    list.Add((fragPeak.RetentionTime, fragPeak.Intensity, index));
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
                sw.WriteLine("rt, intensity, index");
                foreach (var line in list)
                {
                    sw.WriteLine($"{line.rt}, {line.intensity}, {line.index}");
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

        }
    }
}
