using Easy.Common.Extensions;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.RootFinding;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using MzLibUtil;
using Newtonsoft.Json.Linq;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using ThermoFisher.CommonCore.Data;
using ThermoFisher.CommonCore.Data.Business;
using static System.Formats.Asn1.AsnWriter;

namespace EngineLayer.ISD
{
    public class XICfromLFQ
    {
        public static List<Peak>[] GetXICTable(List<Peak> allPeaks, int binsPerDalton)
        {
            var table = new List<Peak>[(int)Math.Ceiling(allPeaks.Max(p => p.Mz) * binsPerDalton) + 1];

            var peaks = allPeaks.ToArray();
            for (int i = 0; i< peaks.Length; i++)
            {
                int roundedMz = (int)Math.Round(peaks[i].Mz * binsPerDalton, 0);
                if (table[roundedMz] == null)
                {
                    table[roundedMz] = new List<Peak>();
                }
                table[roundedMz].Add(peaks[i]);
                peaks[i].otherPeaks = table[roundedMz];
            }
            return table;
        }

        public static List<Ms2ScanWithSpecificMass>[] GroupFragmentIonsXIC(List<Ms2ScanWithSpecificMass>[] scansWithPrecursors, List<Peak>[] ms1Table, 
            List<Peak>[] ms2Table, CommonParameters commonParameters, int binSize)
        {
            for (int i = 0; i < scansWithPrecursors.Length; i++)
            {
                if (scansWithPrecursors[i].Count > 0)
                {
                    for (int j = 0; j < scansWithPrecursors[i].Count; j++)
                    {
                        var targetScan = scansWithPrecursors[i][j];
                        List<double> diaMzs = new List<double>();
                        List<double> diaIntensities = new List<double>();
                        var allMzs = targetScan.TheScan.MassSpectrum.XArray;
                        for (int k = 0; k < allMzs.Length; k++)
                        {
                            double correlation = CalculateCorrelation(targetScan.MostAbundantPrePeak, allMzs[k], ms1Table, ms2Table, binSize);
                            if (correlation > 0.5)
                            {
                                diaMzs.Add(allMzs[k]);
                                diaIntensities.Add(targetScan.TheScan.MassSpectrum.YArray[k]);
                            }
                        }
                        MzSpectrum diaSpectrum = new MzSpectrum(diaMzs.ToArray(), diaIntensities.ToArray(), false);
                        MsDataScan newScan = new MsDataScan(diaSpectrum, targetScan.OneBasedScanNumber, targetScan.TheScan.MsnOrder, targetScan.TheScan.IsCentroid,
                            targetScan.TheScan.Polarity, targetScan.TheScan.RetentionTime, targetScan.TheScan.ScanWindowRange, targetScan.TheScan.ScanFilter, targetScan.TheScan.MzAnalyzer,
                            targetScan.TheScan.TotalIonCurrent, targetScan.TheScan.InjectionTime, targetScan.TheScan.NoiseData, targetScan.TheScan.NativeId,
                            targetScan.TheScan.SelectedIonMZ, targetScan.TheScan.SelectedIonChargeStateGuess, targetScan.TheScan.SelectedIonIntensity,
                            targetScan.TheScan.IsolationMz, targetScan.TheScan.IsolationWidth, targetScan.TheScan.DissociationType, targetScan.TheScan.OneBasedPrecursorScanNumber,
                            targetScan.TheScan.SelectedIonMonoisotopicGuessMz, targetScan.TheScan.HcdEnergy, targetScan.TheScan.ScanDescription);
                        scansWithPrecursors[i][j] = new Ms2ScanWithSpecificMass(newScan, targetScan.PrecursorMonoisotopicPeakMz, targetScan.PrecursorCharge, targetScan.FullFilePath,
                commonParameters, null, targetScan.Pre_RT, mostAbundantPrePeak: targetScan.MostAbundantPrePeak);
                    }
                }
            }
            return scansWithPrecursors;
        }

        public static double CalculateCorrelation(double ms1mz, double ms2mz, List<Peak>[] ms1Table, List<Peak>[] ms2Table, int bin)
        {
            var matchedPeaks = MatchRTs(ms1mz, ms2mz, ms1Table, ms2Table, bin);
            var dataSet1 = matchedPeaks.Select(p => p.Item1).ToArray();
            var dataSet2 = matchedPeaks.Select(p => p.Item2).ToArray();
            if (matchedPeaks.Count >= 5)
            {
                // Calculate Pearson correlation
                double correlation = PeakCurve.CalculatePearsonCorrelation(dataSet1, dataSet2);
                return correlation;
            }
            else
            {
                return double.NaN;
            }
        }
        public static List<(double, double)> MatchRTs(double ms1mz, double ms2mz, List<Peak>[] ms1Table, List<Peak>[] ms2Table, int bin)
        {
            int roundedMs1mz = (int)Math.Round(ms1mz * bin, 0);
            int roundedMs2mz = (int)Math.Round(ms2mz * bin, 0);
            var RT_1 = ms1Table[roundedMs1mz].Select(p => p.RT).ToArray();
            var RT_2 = ms2Table[roundedMs2mz].Select(p => p.RT).ToArray();
            var Intensity_1 = ms1Table[roundedMs1mz].Select(p => p.Intensity).ToArray();
            var Intensity_2 = ms2Table[roundedMs2mz].Select(p => p.Intensity).ToArray();

            if (RT_1 == null || RT_2 == null)
            {
                return new List<(double, double)> { (-1.0, -1.0) };
            }

            List<(double, double)> list = new List<(double, double)>();
            List<(double, double)> list2 = new List<(double, double)>();
            List<(double, double)> list3 = new List<(double, double)>();

            for (int j = 0; j < RT_1.Length; j++)
            {
                list2.Add((RT_1[j], Intensity_1[j]));
            }

            for (int k = 0; k < RT_2.Length; k++)
            {
                list3.Add((RT_2[k], Intensity_2[k]));
            }

            list2 = list2.OrderByDescending(((double, double) i) => i.Item2).ToList();
            list3 = list3.OrderByDescending(((double, double) i) => i.Item2).ToList();
            foreach (var item in list3)
            {
                int num = 0;
                while (list2.Count > 0 && num < list2.Count)
                {
                    if (Math.Abs(list2[num].Item1 - item.Item1) / (list2[num].Item1 + item.Item1) < 0.01)
                    {
                        list.Add((list2[num].Item2, item.Item2));
                        list2.RemoveAt(num);
                        num = -1;
                        break;
                    }
                    num++;
                }
                if (list2.Count == 0)
                {
                    num++;
                }
                if (num > 0)
                {
                    list.Add((0.0, item.Item2));
                }
            }
            return list;
        }
    }
}
