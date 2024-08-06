using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer.ISD;
using MassSpectrometry;

namespace EngineLayer
{
    public class Peak
    {
        public Peak(double mz, double rt, double intensity, int msLevel, int scanNumber = 0, int index = 0, string label = null)
        {
            Mz = mz;
            Intensity = intensity;
            RT = rt;
            ScanNumber = scanNumber;
            Index = index;
            MsLevel = msLevel;
            Label = label;
        }

        public double Mz { get; set; }
        public double Intensity { get; set; }
        public double RT { get; set; }
        public int ScanNumber { get; set; }
        public int Index { get; set; }
        public PeakCurve XICforDIA { get; set; }   
        public int MsLevel { get; set; }
        public List<Peak> XICpeaks {  get; set; }
        public double ApexRT {  get; set; }
        public string Label {  get; set; }
        public XIC XIC {  get; set; }

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

        public static List<Peak>[] GetAllPeaksByScan(MsDataScan[] scans)
        {
            var allPeaks = new List<Peak>[scans.Max(s => s.OneBasedScanNumber) + 1];
            int index = 0;
            foreach (var scan in scans)
            {
                try
                {
                    allPeaks[scan.OneBasedScanNumber] = new List<Peak>();
                    var spectrum = scan.MassSpectrum;
                    for (int j = 0; j < spectrum.XArray.Length; j++)
                    {
                        Peak newPeak = new Peak(spectrum.XArray[j], scan.RetentionTime, spectrum.YArray[j], scan.MsnOrder,
                            scan.OneBasedScanNumber, index);
                        allPeaks[scan.OneBasedScanNumber].Add(newPeak);
                        index++;
                    }
                }
                catch
                {
                    var scanToLook = scan.MassSpectrum;
                    bool stop = true;
                }
                
            }
            return allPeaks;
        }
    }
}
