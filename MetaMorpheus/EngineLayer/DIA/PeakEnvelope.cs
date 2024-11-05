using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PeakEnvelope
    {
        public PeakEnvelope(List<Peak> peaks, double monoIsotopicMass, int charge, double retentionTime, int msLevel, int scanNumber, int zeroBasedScanIndex)
        {
            IsotopicPeaks = peaks;
            MonoisotopicMass = monoIsotopicMass;
            Charge = charge;
            RetentionTime = retentionTime;
            MsLevel = msLevel;
            ScanNumber = scanNumber;
            ZeroBasedScanIndex = zeroBasedScanIndex;
        }
        
        public PeakEnvelope(List<Peak> peaks)
        {
            IsotopicPeaks = peaks;
        }

        public List<Peak> IsotopicPeaks { get; set; }
        public double MonoisotopicMass { get; set; }
        public int Charge { get; set; }
        public double RetentionTime { get; set; }
        public int MsLevel { get; set; }
        public MzRange IsolationRange { get; set; }
        public int ScanNumber { get; set; }
        public int ZeroBasedScanIndex { get; set; }
        public double TotalIntensity => IsotopicPeaks.Sum(p => p.Intensity);
        
        public PeakEnvelopeCurve PeakEnvelopeCurve { get; set; }

        public static PeakEnvelope FindEnvelopeFromScan((double mz, double intensity)[] targetPeaks, double[] theorIntensityRatio, int highestPeakIndex, List<Peak>[] peakTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            var intensityCount = new int[targetPeaks.Length];
            var peaksFound = new List<Peak>();
            var newPE = new PeakEnvelope(peaksFound);

            //go left
            int missedIsotopes = 0;
            for (int i = highestPeakIndex; i >= 0; i--)
            {
                var peak = PeakCurve.GetPeakFromScan(targetPeaks[i].mz, peakTable, zeroBasedScanIndex, tolerance, binSize);
                if (peak != null && peak.PeakEnvelope == null)
                {
                    peaksFound.Add(peak);
                    intensityCount[i] = 1;
                    missedIsotopes = 0;
                }
                else
                {
                    intensityCount[i] = 0;
                    missedIsotopes++;
                }
                if (missedIsotopes > 1)
                {
                    break;
                }
            }

            //go right
            missedIsotopes = 0;
            for (int i = highestPeakIndex + 1; i < targetPeaks.Length; i++)
            {
                var peak = PeakCurve.GetPeakFromScan(targetPeaks[i].mz, peakTable, zeroBasedScanIndex, tolerance, binSize);
                if (peak != null && peak.PeakEnvelope == null)
                {
                    peaksFound.Add(peak);
                    intensityCount[i] = 1;
                    missedIsotopes = 0;
                }
                else
                {
                    intensityCount[i] = 0;
                    missedIsotopes++;
                }
                if (missedIsotopes > 1)
                {
                    break;
                }
            }

            var percentIntensityFound = theorIntensityRatio.Zip(intensityCount, (a, b) => a * b).Sum();
            if (percentIntensityFound >= 0.5)
            {
                newPE.RetentionTime = peaksFound.First().RetentionTime;
                newPE.ScanNumber = peaksFound.First().ScanNumber;
                newPE.ZeroBasedScanIndex = zeroBasedScanIndex;
                return newPE;
            }
            else
            {
                return null;
            }
        }
    }
}
