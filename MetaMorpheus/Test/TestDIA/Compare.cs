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
using Org.BouncyCastle.Asn1.Cmp;

namespace Test.TestDIA
{
    public class Compare
    {
        [Test]
        public static void ComparePrecursors()
        {
            var task = new SearchTask();
            task.CommonParameters.TrimMsMsPeaks = false;
            task.CommonParameters.TrimMs1Peaks = false;
            var q1 = @"E:\DIA\FragPipe\MMsearchUmpireOutput\UmpirePseudoMS2Files\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT_Q1.mzML";
            var q2 = @"E:\DIA\FragPipe\MMsearchUmpireOutput\UmpirePseudoMS2Files\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT_Q2.mzML";
            var q3 = @"E:\DIA\FragPipe\MMsearchUmpireOutput\UmpirePseudoMS2Files\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT_Q3.mzML";
            var myFileManager = new MyFileManager(false);
            var q1Scans = myFileManager.LoadFile(q1, task.CommonParameters).GetAllScansList().ToArray();
            var q2Scans = myFileManager.LoadFile(q2, task.CommonParameters).GetAllScansList().ToArray();
            var q3Scans = myFileManager.LoadFile(q3, task.CommonParameters).GetAllScansList().ToArray();

            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var diaFile = myFileManager.LoadFile(DIAfile, task.CommonParameters);
            var ms1scan = diaFile.GetMS1Scans().First(s => s.OneBasedScanNumber == 1642);
            var envelopes = Deconvoluter.Deconvolute(ms1scan, task.CommonParameters.PrecursorDeconvolutionParameters, new MzRange(630, 650));
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20),
                maxNumMissedScan: 1, overlapRatioCutOff: 0.2, correlationCutOff: 0.5, apexRtTolerance: 0.1,
                fragmentRankCutOff: 200, precursorRankCutOff: 20, maxRTrangeMS1: 0.5, maxRTrangeMS2: 0.5, highCorrThreshold: 0.5, numHighCorrFragments: 0,
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.05f, minMass: 0, maxMass: 10000, type: "DIA", apexCycleTolerance: 3,
                scanCycleSplineInterval: 0.005, ms1XICType: XICType.DeconHighestPeak, ms2XICType: XICType.Peak, pfGroupingType: PFGroupingType.RetentionTime,
                pseudoMs2Type: PseudoMs2ConstructionType.mzPeak, analysisType: AnalysisType.DIAEngine_static, ms1SplineType: SplineType.CubicSpline, ms2SplineType: SplineType.CubicSpline,
                splineRtInterval: 0.005);

            var ms1scans = diaFile.GetMS1Scans().ToArray();
            var allMs1PCs = ISDEngine_static.GetAllPeakCurves_DeconHighestPeak(ms1scans, task.CommonParameters, task.CommonParameters.DIAparameters, new PpmTolerance(10), 0.5, 
                out List<Peak>[] allPeaksByScan, false, new MzRange(630, 650));
            var targetPC = allMs1PCs.Where(pc => Math.Abs(pc.AveragedMz - 636.36) < 0.01 && Math.Round(pc.ApexRT, 2) == 38.31 && pc.Charge == 2).First();

            var allMs1PCs630_650 = ISDEngine_static.GetAllPeakCurves_DeconHighestPeak(ms1scans, task.CommonParameters, task.CommonParameters.DIAparameters, new PpmTolerance(5), 0.5,
                                               out List<Peak>[] allPeaksByScan2, false, new MzRange(630, 650));
            var allMs1PCs613_631 = ISDEngine_static.GetAllPeakCurves_DeconHighestPeak(ms1scans, task.CommonParameters, task.CommonParameters.DIAparameters, new PpmTolerance(5), 0.5,
                                               out List<Peak>[] allPeaksByScan3, false, new MzRange(613, 631));
            var uniqueIDs = PsmTsvReader.ReadTsv(@"E:\DIA\umpire_unique.psmtsv", out List<string> warnings);
            var uniqueIDs_RetentionTime = uniqueIDs.Select(id => id.RetentionTime.Value).ToList();
            var uniqueIDs_PrecursorMz = uniqueIDs.Select(id => id.PrecursorMz).ToList();
            var uniqueIDs_PrecursorCharge = uniqueIDs.Select(id => id.PrecursorCharge).ToList();
            var decon = new bool[uniqueIDs.Count()];
            var ms1PC = new bool[uniqueIDs.Count()];
            var pcLength = new double[uniqueIDs.Count()];
            var pcApex = new double[uniqueIDs.Count()];
            var pcAveragedMz = new double[uniqueIDs.Count()];
            var ms1scanRt = ms1scans.Select(s => s.RetentionTime).ToList();
            foreach (var id in uniqueIDs)
            {
                var index = FindClosestValue(ms1scanRt, id.RetentionTime.Value);
                var envelopes2 = new List<IsotopicEnvelope>();
                var allPCs = new List<PeakCurve>();
                if (id.PrecursorMz >= 630 && id.PrecursorMz <= 650)
                {
                    envelopes2 = Deconvoluter.Deconvolute(ms1scans[index], task.CommonParameters.PrecursorDeconvolutionParameters, new MzRange(630, 650)).ToList();
                    
                } else
                {
                    envelopes2 = Deconvoluter.Deconvolute(ms1scans[index], task.CommonParameters.PrecursorDeconvolutionParameters, new MzRange(613, 630)).ToList();
                }
                foreach(var envelope in envelopes2)
                {
                    if (Math.Abs(envelope.MonoisotopicMass - id.PrecursorMass) < 1 && envelope.Charge == id.PrecursorCharge)
                    {
                        if (envelope.Peaks.Select(p => Math.Round(p.mz, 2)).Contains(Math.Round(id.PrecursorMz, 2)))
                        {
                            decon[uniqueIDs.IndexOf(id)] = true;
                            var highestPeak = envelope.Peaks.OrderByDescending(p => p.intensity).First();
                            PeakCurve pc = null;
                            if (id.PrecursorMz >= 630 && id.PrecursorMz <= 650)
                            {
                                pc = allMs1PCs630_650.Where(pc => Math.Abs(pc.AveragedMz - highestPeak.mz) < 0.01 && Math.Abs(pc.ApexRT - ms1scans[index].RetentionTime) < 0.2 && pc.Charge == envelope.Charge).FirstOrDefault();
                            }
                            else
                            {
                                pc = allMs1PCs613_631.Where(pc => Math.Abs(pc.AveragedMz - highestPeak.mz) < 0.01 && Math.Abs(pc.ApexRT - ms1scans[index].RetentionTime) < 0.2 && pc.Charge == envelope.Charge).FirstOrDefault();
                            }
                            if (pc != null)
                            {
                                ms1PC[uniqueIDs.IndexOf(id)] = true;
                                pcLength[uniqueIDs.IndexOf(id)] = pc.Peaks.Count();
                                pcApex[uniqueIDs.IndexOf(id)] = pc.ApexRT;
                                pcAveragedMz[uniqueIDs.IndexOf(id)] = pc.AveragedMz;
                            }
                        }
                    }
                }
            }

            var outPath = @"E:\DIA\PrecursorCompare.csv";
            // Ensure the file exists and is properly closed
            if (!File.Exists(outPath))
            {
                using (File.Create(outPath)) { }  // Immediately close the created file
            }
            using (StreamWriter writer = new StreamWriter(outPath, false))
            {
                writer.WriteLine("id_RetentionTime, id_PrecursorMz, id_PrecursorCharge, Deconvoluted, xicFound, xicLength, xicApex, xicAveragedMz");
                for (int i = 0; i < uniqueIDs.Count; i++)
                {
                    writer.WriteLine($"{uniqueIDs_RetentionTime[i]},{uniqueIDs_PrecursorMz[i]},{uniqueIDs_PrecursorCharge[i]}, {decon[i]}, {ms1PC[i]}, {pcLength[i]}, {pcApex[i]}, {pcAveragedMz[i]}");
                }
            }
        }

        public static int FindClosestValue(List<double> sortedList, double target)
        {
            int index = sortedList.BinarySearch(target);

            if (index >= 0) // Exact match found
                return index;

            // Convert negative index to insertion point
            int insertionPoint = ~index;

            if (insertionPoint == 0) // Target is smaller than the first element
                return 0;

            if (insertionPoint == sortedList.Count) // Target is greater than the last element
                return sortedList.Count - 1;

            // Compare the two nearest values
            double prev = sortedList[insertionPoint - 1];
            double next = sortedList[insertionPoint];

            return (Math.Abs(prev - target) <= Math.Abs(next - target)) ? insertionPoint - 1 : insertionPoint;
        }
    }
}
