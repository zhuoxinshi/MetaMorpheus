using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PeakEnvelopeCurve
    {
        public PeakEnvelopeCurve(List<PeakEnvelope> peakEnvelopes)
        {
            PeakEnvelopes = peakEnvelopes;
        }

        public List<PeakEnvelope> PeakEnvelopes { get; set; }
        public double MonoisotopicMass { get; set; }
        public int Charge { get; set; }
        public int MsLevel { get; set; }
        public MzRange IsolationRange { get; set; }
        public double ApexRT => PeakEnvelopes.OrderByDescending(p => p.TotalIntensity).First().RetentionTime;
        public int Index { get; set; }
        public PeakCurve FakePeakCurve { get; set; }

        public void GetFakePeakCurve()
        {
            List<Peak> fakePeaks = new List<Peak>();
            foreach (var peakEnvelope in PeakEnvelopes)
            {
                var fakePeak = new Peak(MonoisotopicMass.ToMz(Charge), peakEnvelope.RetentionTime, peakEnvelope.TotalIntensity, MsLevel, peakEnvelope.ScanNumber, 
                    peakEnvelope.ZeroBasedScanIndex);
                fakePeaks.Add(fakePeak);
            }
            FakePeakCurve = new PeakCurve(fakePeaks, MsLevel, IsolationRange, MonoisotopicMass, Charge, index: Index);
            var allPeaks = PeakEnvelopes.SelectMany(pe => pe.IsotopicPeaks).ToList();
            foreach(var peak in allPeaks)
            {
                peak.PeakCurve = FakePeakCurve;
            }
        }

        public static PeakEnvelopeCurve GetPeakEnvelopeCurve((double mz, double intensity)[] targetPeaks, double[] theorIntensityRatio, int highestPeakIndex, List<Peak>[] peakTable,
            MsDataScan[] scans, int zeroBasedScanIndex, DIAparameters diaParam, Tolerance peakFindingTolerance, double maxRTRange)
        {
            var peakEnvelopeList = new List<PeakEnvelope> ();
            var newPEC = new PeakEnvelopeCurve(peakEnvelopeList);

            //this scan
            var firstEnvelope = PeakEnvelope.FindEnvelopeFromScan(targetPeaks, theorIntensityRatio, highestPeakIndex, peakTable, zeroBasedScanIndex,
                peakFindingTolerance, diaParam.PeakSearchBinSize);
            if (firstEnvelope != null)
            {
                peakEnvelopeList.Add(firstEnvelope);
            }
            else
            {
                return null;
            }

            //go right
            int missedScans = 0;
            for (int i = zeroBasedScanIndex + 1; i < scans.Length; i++)
            {
                var peakEnvelope = PeakEnvelope.FindEnvelopeFromScan(targetPeaks, theorIntensityRatio, highestPeakIndex, peakTable, i, peakFindingTolerance, diaParam.PeakSearchBinSize);
                if (peakEnvelope != null && peakEnvelope.PeakEnvelopeCurve == null)
                {
                    peakEnvelopeList.Add(peakEnvelope);
                    missedScans = 0;
                }
                else
                {
                    missedScans++;
                }
                if (missedScans > diaParam.MaxNumMissedScan)
                {
                    break;
                }
                if (scans[i].RetentionTime - newPEC.ApexRT > maxRTRange)
                {
                    break;
                }
            }

            //go left
            missedScans = 0;
            for (int i = zeroBasedScanIndex - 1; i >= 0; i--)
            {
                var peakEnvelope = PeakEnvelope.FindEnvelopeFromScan(targetPeaks, theorIntensityRatio, highestPeakIndex, peakTable, i, peakFindingTolerance, diaParam.PeakSearchBinSize);
                if (peakEnvelope != null && peakEnvelope.PeakEnvelopeCurve == null)
                {
                    peakEnvelopeList.Add(peakEnvelope);
                    missedScans = 0;
                }
                else
                {
                    missedScans++;
                }
                if (missedScans > diaParam.MaxNumMissedScan)
                {
                    break;
                }
                if (newPEC.ApexRT - scans[i].RetentionTime > maxRTRange)
                {
                    break;
                }
            }
            
            if (peakEnvelopeList.Count < 5)
            {
                return null;
            }

            //cut peaks
            if (diaParam.CutMs1Peaks)
            {
                var NL = peakEnvelopeList.Select(p => p.TotalIntensity).OrderBy(p => p).First();
                //DiscriminationFactorToCutPeak is a parameter, default 0.6 in FlashLFQ
                double DiscriminationFactorToCutPeak = 0.6;

                peakEnvelopeList = peakEnvelopeList.OrderBy(p => p.RetentionTime).ToList();
                var apexPeak = peakEnvelopeList.OrderByDescending(p => p.TotalIntensity).First();
                int apexIndex = peakEnvelopeList.IndexOf(apexPeak);
                PeakEnvelope valleyPeakEnvelope = null;
                int indexOfValley = 0;

                //go left
                for (int i = apexIndex; i >= 0; i--)
                {
                    PeakEnvelope timepoint = peakEnvelopeList[i];

                    if (valleyPeakEnvelope == null || timepoint.TotalIntensity < valleyPeakEnvelope.TotalIntensity)
                    {
                        valleyPeakEnvelope = timepoint;
                        indexOfValley = peakEnvelopeList.IndexOf(valleyPeakEnvelope);
                    }

                    double discriminationFactor =
                        (timepoint.TotalIntensity - valleyPeakEnvelope.TotalIntensity) / (timepoint.TotalIntensity - NL);

                    if (discriminationFactor > DiscriminationFactorToCutPeak)
                    {
                        peakEnvelopeList.RemoveAll(p => p.RetentionTime < valleyPeakEnvelope.RetentionTime);
                        break;
                    }
                }

                //go right
                valleyPeakEnvelope = null;
                for (int i = apexIndex; i < peakEnvelopeList.Count; i++)
                {
                    PeakEnvelope timepoint = peakEnvelopeList[i];

                    if (valleyPeakEnvelope == null || timepoint.TotalIntensity < valleyPeakEnvelope.TotalIntensity)
                    {
                        valleyPeakEnvelope = timepoint;
                        indexOfValley = peakEnvelopeList.IndexOf(valleyPeakEnvelope);
                    }

                    double discriminationFactor =
                        (timepoint.TotalIntensity - valleyPeakEnvelope.TotalIntensity) / (timepoint.TotalIntensity - NL);

                    if (discriminationFactor > DiscriminationFactorToCutPeak)
                    {
                        peakEnvelopeList.RemoveAll(p => valleyPeakEnvelope.RetentionTime < p.RetentionTime);
                        break;
                    }
                }
            }

            foreach(var peakEnvelope in peakEnvelopeList)
            {
                peakEnvelope.PeakEnvelopeCurve = newPEC;
                foreach(var peak in peakEnvelope.IsotopicPeaks)
                {
                    peak.PeakEnvelope = peakEnvelope;
                }
            }
            newPEC.PeakEnvelopes = peakEnvelopeList.OrderBy(pe => pe.RetentionTime).ToList();
            return newPEC;
        }
    }
}
