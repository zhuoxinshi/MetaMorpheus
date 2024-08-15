using MassSpectrometry;
using MathNet.Numerics.RootFinding;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http.Headers;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;
using MathNet.Numerics.Statistics;
using Easy.Common.Extensions;
using Chemistry;
using System.Collections.Concurrent;
namespace EngineLayer.ISD
{
    public class XIC
    {
        public List<Peak> XICpeaks {  get; set; }
        public double ApexRT => XICpeaks.OrderByDescending(p => p.Intensity).First().RT;
        public double StartRT => XICpeaks.Min(p => p.RT);
        public double EndRT => XICpeaks.Max(p => p.RT);
        public int ApexScanNumber => XICpeaks.OrderByDescending(p => p.Intensity).First().ScanNumber;
        public double AveragedMz => AverageMz();
        public double AveragedIntensity => AverageIntensity();
        public int RoundedMz { get; set; }
        public XICgroup Group { get; set; }
        public double CorrelationToGroup {  get; set; }

        public XIC CorrelatedXIC { get; set; }

        public XIC(List<Peak> peaks, int roundedMz, XICgroup group = null, double correlation = 0 )
        {
            XICpeaks = peaks;
            RoundedMz = roundedMz;
            CorrelationToGroup = correlation;
        }
        public double AverageMz()
        {
            double sumIntensity = XICpeaks.Sum(p => p.Intensity);
            double sumMz = XICpeaks.Sum(p => p.Mz);
            double averagedMz = 0;
            foreach(var peak in XICpeaks)
            {
                double weight = peak.Intensity / sumIntensity;
                averagedMz += weight * peak.Mz;
            }
            return averagedMz;
        }

        public double AverageIntensity()
        {
            double sumIntensity = XICpeaks.Sum(p => p.Intensity);
            double averagedIntensity = 0;
            foreach (var peak in XICpeaks)
            {
                double weight = peak.Intensity / sumIntensity;
                averagedIntensity += weight * peak.Intensity;
            }
            return averagedIntensity;
        }

        //get XIC straight from the peakTable without considering tolerance
        public static XIC GetXIC(List<Peak>[] peakTable, Peak peak, int binSize)
        {
            int roundedMz = (int)Math.Round(peak.Mz * binSize, 0);
            var peaks = peakTable[roundedMz];
            return new XIC(peaks, roundedMz);
        }

        //get XIC straight from the peakTable without considering tolerance
        public static XIC GetXICwithMz(List<Peak>[] peakTable, double mz, int binSize)
        {
            int roundedMz = (int)Math.Round(mz * binSize, 0);
            var peaks = peakTable[roundedMz];
            return new XIC(peaks, roundedMz);
        }

        //get XIC straight from the peakTable without considering tolerance
        public static List<Peak> GetXICPeaksWithMz(List<Peak>[] peakTable, double mz, int binSize)
        {
            var newXIC = GetXICwithMz(peakTable, mz, binSize);
            return newXIC.XICpeaks;
        }

        //get XICs straight from the peakTable without considering tolerance
        public static List<XIC> GetAllXICs(List<Peak>[] peakTable)
        {
            var allXICs = new List<XIC>();
            for (int i = 0; i < peakTable.Length; i++)
            {
                if (peakTable[i] == null)
                {
                    continue;
                }
                var XIC = new XIC(peakTable[i], i);
                allXICs.Add(XIC);
                foreach(var peak in peakTable[i])
                {
                    peak.XIC = XIC;
                }
            }
            return allXICs;
        }

        
        public static List<XIC> GroupFragmentXIC(XIC preXIC, List<XIC> ms2XICs, double apexTolerance, double rtShift, double corrCutOff)
        {
            var filteredXICs = ms2XICs.Where(xic => Math.Abs(xic.ApexRT - preXIC.ApexRT) < apexTolerance);
            var fragmentXICs = new List<XIC>();
            foreach (var xic in filteredXICs)
            {
                double corr = XICCorr(preXIC.XICpeaks, xic.XICpeaks, rtShift);
                if (corr >= corrCutOff)
                {
                    fragmentXICs.Add(xic);
                }
            }
            return fragmentXICs;
        }

        public static Peak GetPeakFromScanWithinRTWindow(double mz, MsDataScan scan, List<Peak>[] peakTable, Tolerance tolerance, int binSize, double minRT,
            double maxRT)
        {
            Peak bestPeak = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(mz) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(mz) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < peakTable.Length && peakTable[j] != null)
                {
                    var bin = peakTable[j].Where(p => p.RT >= minRT && p.RT <= maxRT).OrderBy(p => p.ScanNumber).ToArray();
                    int index = Array.BinarySearch(bin.Select(p => p.ScanNumber).ToArray(), scan.OneBasedScanNumber);
                    if (index < 0)
                    {
                        continue;
                    }

                    for (int i = index; i < bin.Length; i++)
                    {
                        Peak peak = bin[i];

                        if (peak.ScanNumber > scan.OneBasedScanNumber)
                        {
                            break;
                        }

                        if (tolerance.Within(peak.Mz, mz) && peak.ScanNumber == scan.OneBasedScanNumber
                            && (bestPeak == null || Math.Abs(peak.Mz - mz) < Math.Abs(bestPeak.Mz - mz)))
                        {
                            bestPeak = peak;
                        }
                    }
                }
            }
            return bestPeak;
        }

        public static XIC FindMs2XICWithinRTWindow(double mz, MsDataScan[] ms2Scans, List<Peak>[] ms2Table, Tolerance tolerance, int binSize, double minRT, double maxRT)
        {
            List<Peak> XICpeaks = new List<Peak>();
            for (int i = 0; i < ms2Scans.Length; i++)
            {
                var peak = GetPeakFromScanWithinRTWindow(mz, ms2Scans[i], ms2Table, tolerance, binSize, minRT, maxRT);
                if (peak != null)
                {
                    XICpeaks.Add(peak);
                }
            }
            XIC newXIC = new XIC(XICpeaks, (int)Math.Round(mz * binSize, 0));
            return newXIC;
        }

        //no consideration on apex time, just compute Pearson's correlation
        public static double XICCorr(List<Peak> XIC1, List<Peak> XIC2, double rtShift)
        {
            var RT_1 = XIC1.OrderBy(p => p.RT).Select(p => Math.Round(p.RT, 2)).ToArray();
            var RT_2 = XIC2.OrderBy(p => p.RT).Select(p => Math.Round((p.RT + rtShift), 2)).ToArray();

            if (RT_1 == null || RT_2 == null)
            {
                return double.NaN;
            }

            var ms1Intensity = new List<double>();
            var ms2Intensity = new List<double>();
            for (int i = 0; i < RT_1.Length; i++)
            {
                int index = Array.BinarySearch(RT_2, RT_1[i]);
                if (index >= 0)
                {
                    ms1Intensity.Add(XIC1[i].Intensity);
                    ms2Intensity.Add(XIC2[index].Intensity);
                }
            }
            if (ms1Intensity.Count >= 5 && ms2Intensity.Count >= 5)
            {
                // Calculate Pearson correlation
                double correlation = Correlation.Pearson(ms1Intensity, ms2Intensity);
                return correlation;
            }
            else
            {
                return double.NaN;
            }
        }

        //used in grouping all ms1 peaks
        public static XIC GetXIC_LFQ(Peak targetPeak, MsDataScan[] scans, List<Peak>[] peakTable, Tolerance tolerance, int binSize)
        {
            if (targetPeak.XIC == null)
            {
                List<Peak> XICpeaks = new List<Peak> { targetPeak };
                XIC newXIC = new XIC(XICpeaks, 0);
                targetPeak.XIC = newXIC;
                for (int i = 0; i < scans.Length; i++)
                {
                    var peak = XICfromLFQ.GetPeakFromScan(targetPeak.Mz, scans[i], peakTable, tolerance, binSize);
                    if (peak != null && peak.XIC == null)
                    {
                        XICpeaks.Add(peak);
                        peak.XIC = newXIC;
                    }
                }
                return newXIC;
            }
            return targetPeak.XIC;
        }

        //used in grouping all ms1 peaks
        public static List<XIC> GetAllXICs_LFQ(List<Peak> allPeaks, MsDataScan[] scans, List<Peak>[] peakTable, Tolerance tolerance, int binSize)
        {
            var allXICs = new List<XIC>();
            foreach (Peak peak in allPeaks)
            {
                var xic = GetXIC_LFQ(peak, scans, peakTable, tolerance, binSize);
                if (!allXICs.Contains(xic))
                {
                    allXICs.Add(xic);
                }
            }
            return allXICs;
        }

        //used in grouping all ms1 peaks
        public static Peak GetPeakFromScan(double mz, MsDataScan scan, List<Peak>[] peakTable, Tolerance tolerance, int binSize)
        {
            Peak bestPeak = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(mz) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(mz) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < peakTable.Length && peakTable[j] != null)
                {
                    var bin = peakTable[j].Where(p => p.XIC == null).OrderBy(p => p.ScanNumber).ToList();
                    int index = Array.BinarySearch(bin.Select(p => p.ScanNumber).ToArray(), scan.OneBasedScanNumber);
                    if (index < 0)
                    {
                        continue;
                    }

                    for (int i = index; i < bin.Count; i++)
                    {
                        Peak peak = bin[i];

                        if (peak.ScanNumber > scan.OneBasedScanNumber)
                        {
                            break;
                        }

                        if (tolerance.Within(peak.Mz, mz) && peak.ScanNumber == scan.OneBasedScanNumber
                            && (bestPeak == null || Math.Abs(peak.Mz - mz) < Math.Abs(bestPeak.Mz - mz)))
                        {
                            bestPeak = peak;
                        }
                    }
                }
            }
            return bestPeak;
        }

        public static List<Ms2ScanWithSpecificMass>[] GroupFragmentIons_allXICs(List<Ms2ScanWithSpecificMass>[] scansWithPrecursors, MsDataScan[] ms1scans, MsDataScan[] ms2scans, CommonParameters commonParameters, int binSize, double rtShift, double corrCutOff)
        {
            var ms1Peaks = Peak.GetAllPeaks(ms1scans);
            var ms1Table = XICfromLFQ.GetXICTable(ms1Peaks, binSize);
            var ms2PeaksByScans = Peak.GetAllPeaksByScan(ms2scans);
            var ms2Peaks = ms2PeaksByScans.Where(s => s != null).SelectMany(s => s).ToList();
            var ms2Table = XICfromLFQ.GetXICTable(ms2Peaks, binSize);

            //var allMs2XICs = GetAllXICs_LFQ(ms2Peaks, ms2scans, ms2Table, commonParameters.ProductMassTolerance, binSize);
            var allMs2XICs = GetAllXICs(ms2Table);

            Parallel.ForEach(Partitioner.Create(0, scansWithPrecursors.Length), new ParallelOptions { MaxDegreeOfParallelism = 18 }, //max number of threads modified to use locally
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        if (scansWithPrecursors[i].Count > 0)
                        {
                            for (int j = 0; j < scansWithPrecursors[i].Count; j++)
                            {
                                var targetScan = scansWithPrecursors[i][j];
                                var preXIC = GetXICwithMz(ms1Table, targetScan.MostAbundantPrePeak, binSize);
                                if (preXIC == null)
                                {
                                    continue;
                                }
                                var peaks = new List<(double mz, double intensity)>();
                                var ms2PeakCandidates = ms2PeaksByScans[targetScan.OneBasedScanNumber].ToArray();
                                for (int k = 0; k < ms2PeakCandidates.Length; k++)
                                {
                                    var ms2XIC = ms2PeakCandidates[k].XIC;
                                    double corr = XICfromLFQ.CalculatePearsonCorr(preXIC.XICpeaks, ms2XIC.XICpeaks, rtShift);
                                    if (corr > corrCutOff)
                                    {
                                        peaks.Add((ms2PeakCandidates[k].Mz, ms2PeakCandidates[k].Intensity));
                                    }
                                }
                                MzSpectrum diaSpectrum = new MzSpectrum(peaks.Select(p => p.mz).ToArray(), peaks.Select(p => p.intensity).ToArray(), false);
                                MsDataScan newScan = new MsDataScan(diaSpectrum, targetScan.TheScan.OneBasedScanNumber, targetScan.TheScan.MsnOrder, targetScan.TheScan.IsCentroid,
                        targetScan.TheScan.Polarity, targetScan.TheScan.RetentionTime, targetScan.TheScan.ScanWindowRange, targetScan.TheScan.ScanFilter, targetScan.TheScan.MzAnalyzer,
                        targetScan.TheScan.TotalIonCurrent, targetScan.TheScan.InjectionTime, targetScan.TheScan.NoiseData, targetScan.TheScan.NativeId,
                        targetScan.TheScan.SelectedIonMZ, targetScan.TheScan.SelectedIonChargeStateGuess, targetScan.TheScan.SelectedIonIntensity,
                        targetScan.TheScan.IsolationMz, targetScan.TheScan.IsolationWidth, targetScan.TheScan.DissociationType, null,
                        targetScan.TheScan.SelectedIonMonoisotopicGuessMz, targetScan.TheScan.HcdEnergy, targetScan.TheScan.ScanDescription);
                                if (commonParameters.DeconvoluteMs2Type == "withCharge" || commonParameters.DeconvoluteMs2Type == "none")
                                {
                                    var newMs2WithPre = new Ms2ScanWithSpecificMass(newScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath,
                        commonParameters, null, targetScan.Pre_RT, mostAbundantPrePeak: targetScan.MostAbundantPrePeak);
                                    scansWithPrecursors[i][j] = newMs2WithPre;
                                }
                                if (commonParameters.DeconvoluteMs2Type == "mono")
                                {
                                    var neutralFragments = new IsotopicEnvelope[diaSpectrum.XArray.Length];
                                    for (int l = 0; l < diaSpectrum.XArray.Length; l++)
                                    {
                                        double mz = diaSpectrum.XArray[l];
                                        double intensity = diaSpectrum.YArray[l];
                                        neutralFragments[l] = new IsotopicEnvelope(new List<(double mz, double intensity)> { (mz, intensity) },
                                            mz.ToMass(1), 1, intensity, 0, 0);
                                    }
                                    var newMs2WithPre = new Ms2ScanWithSpecificMass(newScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath,
                        commonParameters, neutralFragments, targetScan.Pre_RT, mostAbundantPrePeak: targetScan.MostAbundantPrePeak);
                                    scansWithPrecursors[i][j] = newMs2WithPre;
                                }

                            }
                        }
                    }
                });
            
            return scansWithPrecursors;
        }


    }
}
