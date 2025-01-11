using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using EngineLayer.DIA;
using MassSpectrometry;

namespace EngineLayer
{
    public class Peak
    {
        public Peak(double mz = 0, double rt = 0, double intensity = 0, int msLevel = 0, int scanNumber = 0, int ZeroBasedScanNumber = 0, int index = 0, PeakCurve peakCurve = null)
        {
            Mz = mz;
            Intensity = intensity;
            RetentionTime = rt;
            ScanNumber = scanNumber;
            Index = index;
            MsLevel = msLevel;
            PeakCurve = peakCurve;
            ZeroBasedScanIndex = ZeroBasedScanNumber;
        }

        public double Mz { get; set; }
        public double Intensity { get; set; }
        public double RetentionTime { get; set; }
        public int ScanNumber { get; set; }
        public int ZeroBasedScanIndex {  get; set; }
        public int Index { get; set; }
        public int MsLevel { get; set; }
        public PeakCurve PeakCurve { get; set; }
        public PeakEnvelope PeakEnvelope { get; set; }

        public static List<Peak> GetAllPeaks(MsDataScan[] scans, int binsPerDalton)
        {
            var allPeaks = new List<Peak>();
            int index = 0;
            int zeroBasedScanIndex = 0;
            //scanIndexMap = new Dictionary<int, int>();
            for (int i = 0; i < scans.Length; i++)
            {
                var spectrum = scans[i].MassSpectrum;
                //scanIndexMap.Add(scans[i].OneBasedScanNumber, zeroBasedScanIndex);
                for (int j = 0; j < spectrum.XArray.Length; j++)
                {
                    Peak newPeak = new Peak(spectrum.XArray[j], scans[i].RetentionTime, spectrum.YArray[j], scans[i].MsnOrder,
                        scans[i].OneBasedScanNumber, zeroBasedScanIndex, index);
                    allPeaks.Add(newPeak);
                    index++;
                }
                zeroBasedScanIndex++;
            }
            return allPeaks;
        }

        public static List<Peak>[] GetAllPeaksByScan(MsDataScan[] scans)
        {
            var maxScanNum = scans[scans.Length - 1].OneBasedScanNumber;
            var peaksByScan = new List<Peak>[maxScanNum + 1];
            int index = 0;
            for (int i = 0; i < scans.Length; i++)
            {
                peaksByScan[scans[i].OneBasedScanNumber] = new List<Peak>();
                var spectrum = scans[i].MassSpectrum;
                for (int j = 0; j < spectrum.XArray.Length; j++)
                {
                    Peak newPeak = new Peak(spectrum.XArray[j], scans[i].RetentionTime, spectrum.YArray[j], scans[i].MsnOrder,
                        scans[i].OneBasedScanNumber, i, index, null);
                    peaksByScan[scans[i].OneBasedScanNumber].Add(newPeak);
                    index++;
                }
            }
            return peaksByScan;
        }

        public static List<Peak>[] GetAllPeaksByScan(MsDataScan[] scans, int numScansPerCycle)
        {
            var maxScanNum = scans[scans.Length - 1].OneBasedScanNumber;
            var peaksByScan = new List<Peak>[maxScanNum + 1];
            int index = 0;
            for (int i = 0; i < scans.Length; i++)
            {
                peaksByScan[scans[i].OneBasedScanNumber] = new List<Peak>();
                var spectrum = scans[i].MassSpectrum;
                var zeroBasedScanIndex = (scans[i].OneBasedScanNumber - 1)/numScansPerCycle;
                for (int j = 0; j < spectrum.XArray.Length; j++)
                {
                    Peak newPeak = new Peak(spectrum.XArray[j], scans[i].RetentionTime, spectrum.YArray[j], scans[i].MsnOrder,
                        scans[i].OneBasedScanNumber, zeroBasedScanIndex, index, null);
                    peaksByScan[scans[i].OneBasedScanNumber].Add(newPeak);
                    index++;                        
                }
            }
            return peaksByScan;
        }

        public static List<Peak>[] GetPeakTable(List<Peak>[] allPeaks, int binsPerDalton, out Dictionary<int, int> scanIndexMap)
        {
            var table = new List<Peak>[(int)Math.Ceiling(allPeaks.Where(v => v != null).SelectMany(p => p).Max(p => p.Mz) * binsPerDalton) + 1];
            int zeroBasedScanIndex = 0;
            scanIndexMap = new Dictionary<int, int>();

            for (int i = 0; i < allPeaks.Length; i++)
            {
                scanIndexMap.Add(allPeaks[i].FirstOrDefault().ScanNumber, zeroBasedScanIndex);
                for (int j = 0; j < allPeaks[i].Count; j++)
                {
                    //Label the peak with zeroBasedScanIndex
                    allPeaks[i][j].ZeroBasedScanIndex = zeroBasedScanIndex;

                    int roundedMz = (int)Math.Round(allPeaks[i][j].Mz * binsPerDalton, 0);

                    if (table[roundedMz] == null)
                    {
                        table[roundedMz] = new List<Peak>();
                    }
                    table[roundedMz].Add(allPeaks[i][j]);
                }
                zeroBasedScanIndex++;
            }
            return table;
        }

        public static List<Peak>[] GetPeakTable(List<Peak> allPeaks, int binsPerDalton)
        {
            var table = new List<Peak>[(int)Math.Ceiling(allPeaks.Max(p => p.Mz) * binsPerDalton) + 1];
            //var peaks = allPeaks.ToArray();
            for (int i = 0; i < allPeaks.Count; i++)
            {
                int roundedMz = (int)Math.Round(allPeaks[i].Mz * binsPerDalton, 0);

                if (table[roundedMz] == null)
                {
                    table[roundedMz] = new List<Peak>();
                }
                table[roundedMz].Add(allPeaks[i]);
            }
            return table;
        }
    }
}
