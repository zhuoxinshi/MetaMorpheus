using EngineLayer.ClassicSearch;
using FlashLFQ;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PrecursorCluster
    {
        public PrecursorCluster(List<Precursor> precursors, double monoisotopicMass, int charge, int envelopeCount, double totalIntensity)
        {
            Precursors = precursors;
            MonoisotopicMass = monoisotopicMass;
            Charge = charge;
            EnvelopeCount = envelopeCount;
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

        public static List<PeakCurve> GetMs1PeakCurves(MsDataScan[] allMs1Scans, DIAparameters DIAparameters, CommonParameters commonParameters)
        {
            //Get all precursors
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
            var pre = allPrecursors.GroupBy(p => new { mass = Math.Round(p.MonoisotopicMass, 2), p.Charge });
            var allPrecursorClusters = new List<PrecursorCluster>();
            foreach (var group in pre)
            {
                var precursors = group.ToList();
                double totalIntensity = group.ToList().Sum(p => p.Envelope.TotalIntensity);
                var newPC = new PrecursorCluster(precursors, group.Key.mass, group.Key.Charge, group.ToList().Count, totalIntensity);
                allPrecursorClusters.Add(newPC);
            }
            var filteredPC = allPrecursorClusters.Where(pc => pc.EnvelopeCount >= 5).ToList();
            var allPeakCurvesForGrouping = new List<PeakCurve>();
            foreach(var pc in filteredPC)
            {
                allPeakCurvesForGrouping.Add(pc.GetPeakCurveForGrouping());
            }
            return allPeakCurvesForGrouping;
        }
    }
}
