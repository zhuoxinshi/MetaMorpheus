using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer
{
    public class Peak
    {
        public Peak(double mz, double rt, double intensity, int msLevel, int scanNumber, int ZeroBasedScanNumber, int index = 0, PeakCurve peakCurve = null, double voltage = 15)
        {
            Mz = mz;
            Intensity = intensity;
            RetentionTime = rt;
            ScanNumber = scanNumber;
            Index = index;
            MsLevel = msLevel;
            PeakCurve = peakCurve;
            ZeroBasedScanIndex = ZeroBasedScanNumber;
            Voltage = voltage;
        }

        public double Mz { get; set; }
        public double Intensity { get; set; }
        public double RetentionTime { get; set; }
        public int ScanNumber { get; set; }
        public int ZeroBasedScanIndex {  get; set; }
        public int Index { get; set; }
        public int MsLevel { get; set; }
        public PeakCurve PeakCurve { get; set; }
        public double Voltage { get; set; }

        public static List<Peak>[] GetAllPeaks(List<Peak>[] allPeaks, MsDataScan[] scans)
        {
            int index = 0;
            string pattern = $@"sid=(\d+)";
            for (int i = 0; i < scans.Length; i++)
            {
                allPeaks[scans[i].OneBasedScanNumber] = new List<Peak>();
                var spectrum = scans[i].MassSpectrum;
                var match = Regex.Match(scans[i].ScanFilter, pattern);
                double voltage = double.Parse(match.Groups[1].Value);
                for (int j = 0; j < spectrum.XArray.Length; j++)
                {
                    Peak newPeak = new Peak(spectrum.XArray[j], scans[i].RetentionTime, spectrum.YArray[j], scans[i].MsnOrder,
                        scans[i].OneBasedScanNumber, 0, index, null, voltage);
                    allPeaks[scans[i].OneBasedScanNumber].Add(newPeak);
                    index++;
                }
            }
            return allPeaks;
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

    }
}
