using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    //public class ISDEngine_scanBased
    //{
    //    public ISDEngine_scanBased(MsDataFile myFile, CommonParameters commonParameters, DIAparameters diaParam)
    //    {
    //        MyMSDataFile = myFile;
    //        CommonParameters = commonParameters;
    //        DIAparameters = diaParam;
    //    }

    //    var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
    //        var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
    //        var ms3Scans = msNScans.Where(p => p.MsnOrder == 3).ToArray();
    //        List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];

    //        if (!ms2Scans.Any())
    //        {
    //            return scansWithPrecursors;
    //        }

    //        var averageParam = new SpectralAveragingParameters()
    //        {
    //            OutlierRejectionType = OutlierRejectionType.SigmaClipping,
    //            SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
    //            NumberOfScansToAverage = 5,
    //            ScanOverlap = 4,
    //            NormalizationType = NormalizationType.RelativeToTics,
    //            SpectralWeightingType = SpectraWeightingType.WeightEvenly,
    //        };

    //        //average ms1
    //        var allscans = myMSDataFile.GetAllScansList().OrderBy(s => s.OneBasedScanNumber).ToArray();
    //        var origMs1scans = allscans.Where(s => s.MsnOrder == 1).ToArray();
    //        var ms1scans = new MsDataScan[origMs1scans.Length];
    //        var averagedMS1 = SpectraFileAveraging.AverageSpectraFile(origMs1scans.ToList(), averageParam).ToList();
    //        for (int i = 0; i < origMs1scans.Length; i++)
    //        {
    //            //if (i < averageParam.NumberOfScansToAverage / 2 || i >= origMs1scans.Length - averageParam.NumberOfScansToAverage / 2)
    //            //{
    //            //    ms1scans[i] = origMs1scans[i];
    //            //}
    //            //else
    //            //{
    //                var newScan = new MsDataScan(averagedMS1[i].MassSpectrum, origMs1scans[i].OneBasedScanNumber, origMs1scans[i].MsnOrder, origMs1scans[i].IsCentroid,
    //                    origMs1scans[i].Polarity, origMs1scans[i].RetentionTime, origMs1scans[i].ScanWindowRange, origMs1scans[i].ScanFilter, origMs1scans[i].MzAnalyzer,
    //                    origMs1scans[i].TotalIonCurrent, origMs1scans[i].InjectionTime, origMs1scans[i].NoiseData, origMs1scans[i].NativeId,
    //                    origMs1scans[i].SelectedIonMZ, origMs1scans[i].SelectedIonChargeStateGuess, origMs1scans[i].SelectedIonIntensity,
    //                    origMs1scans[i].IsolationMz, origMs1scans[i].IsolationWidth, origMs1scans[i].DissociationType, null,
    //                    origMs1scans[i].SelectedIonMonoisotopicGuessMz, origMs1scans[i].HcdEnergy, origMs1scans[i].ScanDescription);
    //                ms1scans[i] = newScan;
    //                //ms1scans[i] = origMs1scans[i];
    //            //}
    //        }

    //        //average ms2
    //        var ms2scans = allscans.Where(s => s.MsnOrder == 2).ToArray();
    //        var ISDScanVoltageMap = PeakCurve.ConstructMs2Group(ms2scans);
    //        foreach (var ms2Group in ISDScanVoltageMap)
    //        {
    //            var origScans = ms2Group.Value.ToArray();
    //            var scans = new MsDataScan[origScans.Length];
    //            var averagedMS2 = SpectraFileAveraging.AverageSpectraFile(origScans.ToList(), averageParam).ToList();
    //            for (int i = 0; i < origScans.Length; i++)
    //            {
    //                //if (i < averageParam.NumberOfScansToAverage / 2 || i >= origScans.Length - averageParam.NumberOfScansToAverage / 2)
    //                //{
    //                //    scans[i] = origScans[i];
    //                //}
    //                //else
    //                //{
    //                    var newScan = new MsDataScan(averagedMS2[i].MassSpectrum, origScans[i].OneBasedScanNumber, origScans[i].MsnOrder, origScans[i].IsCentroid,
    //                        origScans[i].Polarity, origScans[i].RetentionTime, origScans[i].ScanWindowRange, origScans[i].ScanFilter, origScans[i].MzAnalyzer,
    //                        origScans[i].TotalIonCurrent, origScans[i].InjectionTime, origScans[i].NoiseData, origScans[i].NativeId,
    //                        origScans[i].SelectedIonMZ, origScans[i].SelectedIonChargeStateGuess, origScans[i].SelectedIonIntensity,
    //                        origScans[i].IsolationMz, origScans[i].IsolationWidth, origScans[i].DissociationType, origScans[i].OneBasedPrecursorScanNumber,
    //                        origScans[i].SelectedIonMonoisotopicGuessMz, origScans[i].HcdEnergy, origScans[i].ScanDescription);
    //                    scans[i] = newScan;
    //                    //scans[i] = origScans[i];
    //                //}
    //            }
    //            ms2Group.Value.Clear();
    //            ms2Group.Value.AddRange(scans);
    //        }

    //        var allPeaks = new List<Peak>[allscans.Length + 1];
    //        Peak.GetAllPeaks(allPeaks, ms1scans);
    //        var allMs1Peaks = allPeaks.Where(p => p != null).ToArray();
    //        var ms1PeakTable = Peak.GetPeakTable(allMs1Peaks, 100, out Dictionary<int, int> scanIndexMap);
    //        var DIAparam= commonParameters.DIAparameters;
    //        var allMs2PeakCurves = PeakCurve.GetMs2PeakCurves(ISDScanVoltageMap, allPeaks, DIAparam);
    //        //var allMs1PeakCurves = PeakCurve.GetAllMs1PeakCurves(ms1scans, allPeaks, ms1PeakTable, DIAparam);

    //        var allMs2Scans = ISDScanVoltageMap.Values.SelectMany(p => p).ToArray();
    //        Parallel.ForEach(Partitioner.Create(0, allMs2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
    //            (partitionRange, loopState) =>
    //            {
    //                var precursors = new List<(double MonoPeakMz, int Charge, double Intensity, int PeakCount, double highestPeakMz, double? FractionalIntensity)>();

    //                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
    //                {
    //                    if (GlobalVariables.StopLoops) { break; }

    //                    precursors.Clear();
    //                    MsDataScan ms2scan = allMs2Scans[i];

    //                    if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
    //                    {
    //                        MsDataScan precursorSpectrum = ms1scans.First(s => s.OneBasedScanNumber == ms2scan.OneBasedPrecursorScanNumber);

    //                        try
    //                        {
    //                            ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
    //                        }
    //                        catch (MzLibException ex)
    //                        {
    //                            Warn("Could not get precursor ion for MS2 scan #" + ms2scan.OneBasedScanNumber + "; " + ex.Message);
    //                            continue;
    //                        }

    //                        if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
    //                        {
    //                            ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
    //                        }

    //                        if (commonParameters.DoPrecursorDeconvolution)
    //                        {
    //                            foreach (IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
    //                                precursorSpectrum.MassSpectrum, commonParameters.PrecursorDeconvolutionParameters))
    //                            {
    //                                if(envelope.MonoisotopicMass < DIAparam.MinMass)
    //                                {
    //                                    continue;
    //                                }
    //                                double monoPeakMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
    //                                int peakCount = envelope.Peaks.Count();
    //                                double intensity = 1;
    //                                if (commonParameters.UseMostAbundantPrecursorIntensity)
    //                                {
    //                                    intensity = envelope.Peaks.Max(p => p.intensity);
    //                                }
    //                                else
    //                                {
    //                                    intensity = envelope.Peaks.Sum(p => p.intensity);
    //                                }

    //                                var fractionalIntensity = envelope.TotalIntensity /
    //                                      (double)precursorSpectrum.MassSpectrum.YArray
    //                                      [
    //                                          precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Minimum)..precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Maximum)
    //                                      ].Sum();

    //                                double highestPeakMz = envelope.Peaks.OrderByDescending(p => p.intensity).First().mz;

    //                                precursors.Add((monoPeakMz, envelope.Charge, intensity, peakCount, highestPeakMz,
    //                                    fractionalIntensity));
    //                            }
    //                        }
    //                    }

    //                    scansWithPrecursors[i] = new List<Ms2ScanWithSpecificMass>();
    //                    IsotopicEnvelope[] neutralExperimentalFragments = null;

    //                    foreach (var precursor in precursors)
    //                    {
    //                        // assign precursor for this MS2 scan
    //                        var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.MonoPeakMz,
    //                            precursor.Charge, fullFilePath, commonParameters, neutralExperimentalFragments,
    //                            precursor.Intensity, precursor.PeakCount, precursor.FractionalIntensity);
    //                        scan.HighestPeakMz = precursor.highestPeakMz;

    //                        scansWithPrecursors[i].Add(scan);
    //                    }
    //                }
    //            });
    //        //debug
    //        var preCount = scansWithPrecursors.SelectMany(p => p).Count();

    //        PeakCurve.GetPrecursorPeakCurve(scansWithPrecursors, ms1scans, allPeaks, ms1PeakTable, DIAparam, scanIndexMap);
    //        var newScanWithPre = PeakCurve.GetPseudoMs2Scans(scansWithPrecursors, commonParameters, DIAparam, allPeaks);
    //        //var newScanWithPre = PeakCurve.GetPseudoMs2Scans_PFgroup(scansWithPrecursors, commonParameters, DIAparam, allPeaks, allMs2PeakCurves, myMSDataFile);

    //        //debug
    //        var post = newScanWithPre.SelectMany(p => p).OrderByDescending(s => s.NumPeaks).ToList();
    //        return newScanWithPre;
//    }
}
