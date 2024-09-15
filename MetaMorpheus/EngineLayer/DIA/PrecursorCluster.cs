using EngineLayer.ClassicSearch;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.DIA
{
    public class PrecursorCluster
    {
        public PrecursorCluster(double monoisotopicMass, int charge, double totalIntensity = 0, List<Precursor> precursors = null)
        {
            Precursors = precursors;
            MonoisotopicMass = monoisotopicMass;
            Charge = charge;
            TotalIntensity = totalIntensity;
        }

        public List<Precursor> Precursors;
        public List<PeakCurve> PeakCurves;
        public PeakCurve PeakCurveForGrouping {  get; set; }
        public double ApexRT => Precursors.OrderByDescending(p => p.Envelope.TotalIntensity).First().RetentionTime;
        public List<double> RetentionTimes => Precursors.Select(p => p.RetentionTime).ToList();
        public List<double> TotalIntensities => Precursors.Select(p => p.Envelope.TotalIntensity).ToList();
        public double MonoisotopicMass;
        public int Charge;
        public int EnvelopeCount;
        public double TotalIntensity;

        public PeakCurve GetPeakCurveForGrouping()
        {
            var fakePeakList = new List<Peak>();
            foreach (var precursor in Precursors)
            {
                var totalIntensity = precursor.Envelope.TotalIntensity;
                var fakePeak = new Peak(double.NaN, precursor.RetentionTime, totalIntensity, 1, precursor.ScanNumber, precursor.ZeroBasedScanIndex);
                fakePeakList.Add(fakePeak);
            }
            var minMzs = Precursors.Select(p => p.Envelope.Peaks.OrderBy(peak => peak.mz).First().mz).ToList();
            var startmz = minMzs.OrderBy(m => m).First();
            var maxMzs = Precursors.Select(p => p.Envelope.Peaks.OrderByDescending(peak => peak.mz).First().mz).ToList();
            var endmz = maxMzs.OrderByDescending(m => m).First();
            var newCurve = new PeakCurve(fakePeakList, 1, null, MonoisotopicMass, Charge, startMz: startmz, endMz : endmz);
            return newCurve;
        }

        //public static List<PeakCurve> GetMs1PeakCurves(MsDataScan[] allMs1Scans, DIAparameters DIAparameters, CommonParameters commonParameters)
        //{
        //    //Get all precursors
        //    var allPrecursors = new List<Precursor>();
        //    for (int i = 0; i < allMs1Scans.Length; i++)
        //    {
        //        var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], commonParameters.PrecursorDeconvolutionParameters);
        //        foreach (var envelope in envelopes)
        //        {
        //            var charge = envelope.Charge;
        //            double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().mz;
        //            double highestPeakIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).FirstOrDefault().intensity;
        //            var precursor = new Precursor(envelope, charge, allMs1Scans[i].RetentionTime, highestPeakMz, highestPeakIntensity, envelope.MonoisotopicMass,
        //                allMs1Scans[i].OneBasedScanNumber, i);
        //            allPrecursors.Add(precursor);
        //        }
        //    }
        //    var pre = allPrecursors.GroupBy(p => new { mass = Math.Round(p.MonoisotopicMass, 2), p.Charge });
        //    var allPrecursorClusters = new List<PrecursorCluster>();
        //    foreach (var group in pre)
        //    {
        //        var precursors = group.ToList();
        //        double totalIntensity = group.ToList().Sum(p => p.Envelope.TotalIntensity);
        //        var newPC = new PrecursorCluster( group.Key.mass, group.Key.Charge, totalIntensity, precursors);
        //        allPrecursorClusters.Add(newPC);
        //    }
        //    var filteredPC = allPrecursorClusters.Where(pc => pc.EnvelopeCount >= 5).ToList();
        //    var allPeakCurvesForGrouping = new List<PeakCurve>();
        //    foreach(var pc in filteredPC)
        //    {
        //        allPeakCurvesForGrouping.Add(pc.GetPeakCurveForGrouping());
        //    }
        //    return allPeakCurvesForGrouping;
        //}

        public static List<Peak> FindEnvelopeFromScan(List<Peak> targetPeaks, List<Peak>[] peakTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            var theorIntensityRatio = targetPeaks.Select(p => p.Intensity/targetPeaks.Sum(p => p.Intensity)).ToArray();
            var intensityCount = new List<int>();
            var peaksFound = new List<Peak>();
            foreach (var peak in targetPeaks)
            {
                var peakFound = PeakCurve.GetPeakFromScan(peak.Mz, peakTable, zeroBasedScanIndex, tolerance, binSize);
                if (peakFound != null)
                {
                    peaksFound.Add(peakFound);
                    intensityCount.Add(1);
                }
                else
                {
                    intensityCount.Add(0);
                }
            }
            var percentIntensityFound = theorIntensityRatio.Zip(intensityCount, (a,b) => a*b).Sum();
            if (percentIntensityFound >= 0.3)
            {
                return peaksFound;
            }
            else
            {
                return null;
            }
        }

        public static List<List<Peak>> FindIsotopePeakCurve(List<Peak> targetPeaks, List<Peak>[] peakTable, int zeroBasedScanIndex, MsDataScan[] scans, Tolerance tolerance, int binSize, int maxMissedScans)
        {
            var isotopeList = new List<List<Peak>>();

            //go right
            int missedScans = 0;
            for (int i = targetPeaks.First().ZeroBasedScanIndex + 1; i < scans.Length; i++)
            {
                var peaks = FindEnvelopeFromScan(targetPeaks, peakTable, i, tolerance, binSize);
                if (peaks != null)
                {
                    isotopeList.Add(peaks);
                }
                else
                {
                    missedScans ++;
                }
                if (missedScans > maxMissedScans)
                {
                    break;
                }
            }

            //go left
            missedScans = 0;
            for (int i = targetPeaks.First().ZeroBasedScanIndex - 1; i >=0; i--)
            {
                var peaks = FindEnvelopeFromScan(targetPeaks, peakTable, i, tolerance, binSize);
                if (peaks != null)
                {
                    isotopeList.Add(peaks);
                }
                else
                {
                    missedScans++;
                }
                if (missedScans > maxMissedScans)
                {
                    break;
                }
            }
            return isotopeList;
        }

        public static List<PeakCurve> GetMs1PeakCurves_isotope(MsDataScan[] allMs1Scans, List<Peak>[] ms1PeakTable, DIAparameters DIAparameters, CommonParameters commonParameters)
        {
            var allMs1PeakCurves = new List<PeakCurve>();
            var allPrecursors = new List<Precursor>();
            for (int i = 0; i < allMs1Scans.Length; i++)
            {
                var envelopes = Deconvoluter.Deconvolute(allMs1Scans[i], commonParameters.PrecursorDeconvolutionParameters);
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
            var preGroups = allPrecursors.GroupBy(p => new { mass = Math.Round(p.MonoisotopicMass, 2), p.Charge }).ToList();
            var referencePrecursors = preGroups.Select(g => g.OrderByDescending(p => p.Envelope.Peaks.Count).ThenBy(p => p.Envelope.TotalIntensity).First()).ToList();
            foreach (var precursor in referencePrecursors)
            {
                var targetList = precursor.Envelope.Peaks.Select(p => p.mz).ToList();
                var targetPeaks = new List<Peak>();
                foreach(var targetMz in targetList)
                {
                    var peak = PeakCurve.GetPeakFromScan(targetMz, ms1PeakTable, precursor.ZeroBasedScanIndex, new PpmTolerance(0), DIAparameters.PeakSearchBinSize);
                    targetPeaks.Add(peak);
                }
                var isotopePeakList = FindIsotopePeakCurve(targetPeaks, ms1PeakTable, precursor.ZeroBasedScanIndex, allMs1Scans, DIAparameters.Ms1PeakFindingTolerance,
                    DIAparameters.PeakSearchBinSize, DIAparameters.MaxNumMissedScan);
                var fakePeaks = new List<Peak>();
                foreach (var envelopePeaks in isotopePeakList)
                {
                    var mz = envelopePeaks.OrderByDescending(p => p.Intensity).First().Mz;
                    var rt = envelopePeaks.First().RetentionTime;
                    var intensity = envelopePeaks.Sum(p => p.Intensity);
                    var fakePeak = new Peak(mz, rt, intensity, 1, envelopePeaks.First().ScanNumber, envelopePeaks.First().ZeroBasedScanIndex);
                    fakePeaks.Add(fakePeak);
                }
                var newPC = new PeakCurve(fakePeaks, 1, null, precursor.MonoisotopicMass, precursor.Charge, startMz: targetList.Min(), endMz: targetList.Max());
                allMs1PeakCurves.Add(newPC);
            }
            return allMs1PeakCurves;
        }
    }
}
