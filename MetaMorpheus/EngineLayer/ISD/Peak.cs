using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace EngineLayer
{
    public class Peak
    {
        public Peak(double mz, double rt, double intensity, int msLevel, int scanNumber = 0, int index = 0, PeakCurve curve = null)
        {
            Mz = mz;
            Intensity = intensity;
            RT = rt;
            ScanNumber = scanNumber;
            Index = index;
            XIC = curve;
            MsLevel = msLevel;
        }

        public double Mz { get; set; }
        public double Intensity { get; set; }
        public double RT { get; set; }
        public int ScanNumber { get; set; }
        public int Index { get; set; }
        public PeakCurve XIC { get; set; }   
        public int MsLevel { get; set; }
        public List<Peak> otherPeaks {  get; set; }

        public static List<Peak> GetAllPeaks(MsDataScan[] scans)
        {
            var allPeaks = new List<Peak>();
            int index = 0;
            for(int i = 0; i < scans.Length; i++)
            {
                var spectrum = scans[i].MassSpectrum;
                for (int j = 0; j < spectrum.XArray.Length; j++)
                {
                    Peak newPeak = new Peak(spectrum.XArray[j], scans[i].RetentionTime, spectrum.YArray[j], scans[i].MsnOrder, 
                        scans[i].OneBasedScanNumber, index);
                    allPeaks.Add(newPeak);
                    index++;
                }
            }
            return allPeaks;
        }
    }
}
