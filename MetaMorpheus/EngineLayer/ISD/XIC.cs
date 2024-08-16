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

        //construct peak table for making XIC
        public static List<Peak>[] GetXICTable(List<Peak> allPeaks, int binsPerDalton)
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
                allPeaks[i].XICpeaks = table[roundedMz];
            }
            return table;
        }

        //get XIC straight from the peakTable without considering tolerance
        public static XIC GetXICwithMz(List<Peak>[] peakTable, double mz, int binSize)
        {
            int roundedMz = (int)Math.Round(mz * binSize, 0);
            var peaks = peakTable[roundedMz];
            return new XIC(peaks, roundedMz);
        }

        //get all XICs in the run straight from the peakTable without considering tolerance
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

        //get the XIC for a peak using the method from LFQ (go to each scan, find the best peak within tolerance)
        public static XIC GetXIC_LFQ(Peak targetPeak, MsDataScan[] scans, List<Peak>[] peakTable, Tolerance tolerance, int binSize)
        {
            if (targetPeak.XIC == null)
            {
                List<Peak> XICpeaks = new List<Peak> { targetPeak };
                XIC newXIC = new XIC(XICpeaks, 0);
                targetPeak.XIC = newXIC;
                for (int i = 0; i < scans.Length; i++)
                {
                    var peak = XIC.GetPeakFromScan(targetPeak.Mz, scans[i], peakTable, tolerance, binSize);
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

        //get all XICs in the run using the method from LFQ
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

        //Find the best peak in one scan with a specific mz within tolerance
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

        //Group precursor XIC with all the fragment XICs that have high correlation and close apex RT
        public static List<XIC> GroupFragmentXIC(XIC preXIC, List<XIC> ms2XICs, double apexTolerance, double rtShift, double corrCutOff)
        {
            var filteredXICs = ms2XICs.Where(xic => Math.Abs(xic.ApexRT - preXIC.ApexRT) < apexTolerance);
            var fragmentXICs = new List<XIC>();
            foreach (var xic in filteredXICs)
            {
                double corr = CalculatePearsonCorr(preXIC.XICpeaks, xic.XICpeaks, rtShift);
                if (corr >= corrCutOff)
                {
                    fragmentXICs.Add(xic);
                }
            }
            return fragmentXICs;
        }

        //compute Pearson's correlation between two XICs
        public static double CalculatePearsonCorr(List<Peak> XIC1, List<Peak> XIC2, double rtShift)
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

        //get the correlation between two XICs (apex RT + Pearson)
        public static double XICCorr(XIC xic1, XIC xic2, double rtShift, double apexTolerance)
        {
            if (Math.Abs(xic1.ApexRT - xic2.ApexRT) > apexTolerance)
            {
                return 0;
            }
            double corr = CalculatePearsonCorr(xic1.XICpeaks, xic2.XICpeaks, rtShift);
            return corr;
        }
        
        //Group precursor ion with fragment ions
        public static List<Ms2ScanWithSpecificMass>[] GroupFragmentIons_allXICs(List<Ms2ScanWithSpecificMass>[] scansWithPrecursors, MsDataScan[] ms1scans, MsDataScan[] ms2scans, CommonParameters commonParameters, int binSize, double rtShift, double corrCutOff)
        {
            var ms1Peaks = Peak.GetAllPeaks(ms1scans);
            var ms1Table = GetXICTable(ms1Peaks, binSize);
            var ms2PeaksByScans = Peak.GetAllPeaksByScan(ms2scans);
            var ms2Peaks = ms2PeaksByScans.Where(s => s != null).SelectMany(s => s).ToList();
            var ms2Table = GetXICTable(ms2Peaks, binSize);

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
                                    double corr = XICCorr(preXIC, ms2XIC, rtShift, 0.3);
                                    if (corr > corrCutOff)
                                    {
                                        peaks.Add((ms2PeakCandidates[k].Mz, ms2PeakCandidates[k].Intensity));
                                    }
                                }
                                MzSpectrum diaSpectrum = new MzSpectrum(peaks.Select(p => p.mz).ToArray(), peaks.Select(p => p.intensity).ToArray(), false);
                                //if (commonParameters.DeconvoluteMs2Type == "withCharge" || commonParameters.DeconvoluteMs2Type == "none")
                                //{
                                //    var neutralFragments = Deconvoluter.Deconvolute(diaSpectrum, commonParameters.ProductDeconvolutionParameters, diaSpectrum.Range).ToArray();
                                //    var newMs2WithPre = new Ms2ScanWithSpecificMass(targetScan.TheScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath,
                                //        commonParameters, neutralFragments, targetScan.Pre_RT, mostAbundantPrePeak: targetScan.MostAbundantPrePeak);
                                //    scansWithPrecursors[i][j] = newMs2WithPre;
                                //}
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
                                    var newMs2WithPre = new Ms2ScanWithSpecificMass(targetScan.TheScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath,
                                        commonParameters, neutralFragments, targetScan.Pre_RT, mostAbundantPrePeak: targetScan.MostAbundantPrePeak);
                                    scansWithPrecursors[i][j] = newMs2WithPre;
                                }

                            }
                        }
                    }
                });
            
            return scansWithPrecursors;
        }

        //For deconvoluting MS2 scans
        public static List<Peak> GetDeconvolutedPeaks(MsDataScan scan, CommonParameters commonParam, string type)
        {
            var neutralExperimentalFragmentMasses =
                Deconvoluter.Deconvolute(scan, commonParam.ProductDeconvolutionParameters, scan.MassSpectrum.Range).ToList();

            if (commonParam.AssumeOrphanPeaksAreZ1Fragments)
            {
                HashSet<double> alreadyClaimedMzs = new HashSet<double>(neutralExperimentalFragmentMasses
                    .SelectMany(p => p.Peaks.Select(v => Chemistry.ClassExtensions.RoundedDouble(v.mz).Value)));

                for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                {
                    double mz = scan.MassSpectrum.XArray[i];
                    double intensity = scan.MassSpectrum.YArray[i];

                    if (!alreadyClaimedMzs.Contains(Chemistry.ClassExtensions.RoundedDouble(mz).Value))
                    {
                        neutralExperimentalFragmentMasses.Add(new IsotopicEnvelope(
                            new List<(double mz, double intensity)> { (mz, intensity) },
                            mz.ToMass(1), 1, intensity, 0, 0));
                    }
                }
            }

            var newPeaks = new Dictionary<double, double>();
            var allPeaks = new List<Peak>();

            switch (type)
            {
                case "mono":
                    foreach (var envelope in neutralExperimentalFragmentMasses)
                    {
                        double newMz = envelope.MonoisotopicMass.ToMz(1);
                        double newIntensity = envelope.Peaks.OrderByDescending(p => p.intensity).First().intensity;
                        if (!newPeaks.ContainsKey(newMz))
                        {
                            newPeaks.Add(newMz, newIntensity);
                        }
                        else
                        {
                            newPeaks[newMz] = newPeaks[newMz] + newIntensity;
                        }
                    }
                    break;

                case "cs":
                    foreach (var envelope in neutralExperimentalFragmentMasses)
                    {
                        var peaks = envelope.Peaks;
                        foreach (var peak in peaks)
                        {
                            double newMz = peak.mz.ToMass(envelope.Charge).ToMz(1);
                            double newIntensity = peak.intensity;
                            if (!newPeaks.ContainsKey(newMz))
                            {
                                newPeaks.Add(newMz, newIntensity);
                            }
                            else
                            {
                                newPeaks[newMz] = newPeaks[newMz] + newIntensity;
                            }
                        }

                    }
                    break;

                case "withCharge":
                    foreach (var envelope in neutralExperimentalFragmentMasses)
                    {
                        var peaks = envelope.Peaks;
                        foreach (var peak in peaks)
                        {
                            newPeaks.Add(peak.mz, peak.intensity);
                        }
                    }
                    break;

            }
            foreach (var peak in newPeaks)
            {
                allPeaks.Add(new Peak(peak.Key, scan.RetentionTime, peak.Value, 2, scan.OneBasedScanNumber));
            }
            return allPeaks;
        }

        //Take the deconvoluted peaks and make a new scan
        public static MsDataScan GetDeconvolutedScan(MsDataScan scan, CommonParameters commonParameters, string type)
        {
            var peaks = GetDeconvolutedPeaks(scan, commonParameters, type).OrderBy(p => p.Mz).ToList();
            var spectrum = new MzSpectrum(peaks.Select(p => p.Mz).ToArray(), peaks.Select(p => p.Intensity).ToArray(), false);
            MsDataScan newScan = new MsDataScan(spectrum, scan.OneBasedScanNumber, scan.MsnOrder, scan.IsCentroid,
                scan.Polarity, scan.RetentionTime, scan.ScanWindowRange, scan.ScanFilter, scan.MzAnalyzer,
                scan.TotalIonCurrent, scan.InjectionTime, scan.NoiseData, scan.NativeId,
                scan.SelectedIonMZ, scan.SelectedIonChargeStateGuess, scan.SelectedIonIntensity,
                scan.IsolationMz, scan.IsolationWidth, scan.DissociationType, scan.OneBasedPrecursorScanNumber,
                scan.SelectedIonMonoisotopicGuessMz, scan.HcdEnergy, scan.ScanDescription);
            return newScan;
        }

        //Calculate the average RT interval between MS1 and ISD scans to match the RTs
        public static double GetRTshift(MsDataScan[] allScans)
        {
            List<double> rtShift = new List<double>();
            for (int i = 0; i < allScans.Length - 1; i++)
            {
                double rtDiff = allScans[i].RetentionTime - allScans[i + 1].RetentionTime;
                rtShift.Add(rtDiff);
            }
            double meanRTdiff = rtShift.Average();
            return meanRTdiff;
        }
    }
}
