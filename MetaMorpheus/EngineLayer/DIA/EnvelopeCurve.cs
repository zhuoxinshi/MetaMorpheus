using FlashLFQ;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.DIA
{
    public class EnvelopeCurve
    {
        public EnvelopeCurve(List<DeconvolutedMass> envelopes, int index = 0)
        {
            Envelopes = envelopes;
            MonoisotopicMass = envelopes.First().MonoisotopicMass;
            Charge = envelopes.First().Charge;
            Index = index;
        }

        public List<DeconvolutedMass> Envelopes { get; set; }
        public int MassIndex {  get; set; }
        public List<Peak>[] Peaks { get; set; }
        public double MonoisotopicMass {  get; set; }
        public int Charge { get; set; }
        public double ApexRT { get; set; }
        public PeakCurve FakePeakCurve { get; set; }
        public int Index { get; set; }

        public PeakCurve GetEnvelopePeakCurve()
        {
            var fakePeaks = new List<Peak>();
            foreach(var envelopePeaks in Peaks)
            {
                var mz = envelopePeaks.OrderByDescending(p => p.Intensity).First().Mz;
                var rt = envelopePeaks.First().RetentionTime;
                var intensity = envelopePeaks.Sum(p => p.Intensity);
                var fakePeak = new Peak(mz, rt, intensity, 1, envelopePeaks.First().ScanNumber, envelopePeaks.First().ZeroBasedScanIndex);
                fakePeaks.Add(fakePeak);
            }
            var newPC = new PeakCurve(fakePeaks, 1, null, MonoisotopicMass, Charge = Charge);
            return newPC;
        }

        public PeakCurve GetFakePeakCurve()
        {
            var fakePeaks = new List<Peak>();
            foreach (var envelope in Envelopes)
            {
                var mz = envelope.HighestPeakMz;
                var rt = envelope.RetentionTime;
                var intensity = envelope.AdjustedTotalIntensity;
                var fakePeak = new Peak(mz, rt, intensity, envelope.MsLevel, envelope.ScanNumber, envelope.ZeroBasedScanIndex);
                fakePeaks.Add(fakePeak);
            }
            var fakePC = new PeakCurve(fakePeaks, 1, null, MonoisotopicMass, Charge = Charge);
            FakePeakCurve = fakePC;
            return fakePC;
        }

        public static EnvelopeCurve GetEnvelopeCurve(DeconvolutedMass mass, List<DeconvolutedMass>[] allMasses, List<Peak>[] peakTable, DIAparameters DIAparameters)
        {
            var list = new List<DeconvolutedMass> { mass };
            var newEC = new EnvelopeCurve(list);

            //go left
            int missedScans = 0;
            for (int i = mass.ZeroBasedScanIndex - 1; i >= 0; i--)
            {
                var envelope = FindMassFromScan(mass, allMasses[i], DIAparameters);
                if (envelope != null)
                {
                    list.Add(envelope);
                    envelope.EnvelopeCurve = newEC;
                    missedScans = 0;
                }
                else
                {
                    missedScans++;
                }
                if (missedScans > DIAparameters.MaxNumMissedScan)
                {
                    break;
                }
            }

            //go right
            missedScans = 0;
            for (int i = mass.ZeroBasedScanIndex + 1; i < allMasses.Length; i++)
            {
                var envelope = FindMassFromScan(mass, allMasses[i], DIAparameters);
                if (envelope != null && envelope.EnvelopeCurve == null)
                {
                    list.Add(envelope);
                    envelope.EnvelopeCurve = newEC;
                    missedScans = 0;
                }
                else
                {
                    missedScans++;
                }
                if (missedScans > DIAparameters.MaxNumMissedScan)
                {
                    break;
                }
            }
            newEC.Envelopes = list.OrderBy(p => p.RetentionTime).ToList();
            var targetEnvelope = newEC.Envelopes.OrderByDescending(p => p.HighestPeakIntensity).First();
            var targetPeakMzs = targetEnvelope.Envelope.Peaks.Select(p => p.mz).ToList();
            foreach(var envelope in newEC.Envelopes)
            {
                envelope.Isotopes = FindIsotopes(targetPeakMzs, peakTable, envelope.ZeroBasedScanIndex, DIAparameters);
            }
            return newEC;
        }

        public static DeconvolutedMass FindMassFromScan(DeconvolutedMass mass, List<DeconvolutedMass> massList, DIAparameters DIAparameters)
        {
            var masses = massList.Where(m => m.Charge == mass.Charge).ToArray();
            if (masses.Length == 0)
            {
                return null;
            }
            int ind = Array.BinarySearch(masses.Select(p => p.MonoisotopicMass).ToArray(), mass.MonoisotopicMass);
            if (ind >= 0)
            {
                return masses[ind]; 
            }
            else
            {
                ind = ~ind;
                if (ind == 0)
                {
                    if (DIAparameters.PrecursorMassTolerance.Within(mass.MonoisotopicMass, masses[0].MonoisotopicMass))
                    {
                        return masses[0];
                    }
                    else
                    {
                        return null;
                    }
                } 
                else if (ind == masses.Length)
                {
                    if (DIAparameters.PrecursorMassTolerance.Within(mass.MonoisotopicMass, masses[ind - 1].MonoisotopicMass))
                    {
                        return masses[ind - 1];
                    }
                    else
                    {
                        return null;
                    }
                }
                else
                {
                    var previousDiff = Math.Abs(masses[ind - 1].MonoisotopicMass - mass.MonoisotopicMass);
                    var nextDiff = Math.Abs(masses[ind].MonoisotopicMass - mass.MonoisotopicMass);
                    if (previousDiff < nextDiff)
                    {
                        if (DIAparameters.PrecursorMassTolerance.Within(mass.MonoisotopicMass, masses[ind - 1].MonoisotopicMass))
                        {
                            return masses[ind - 1];
                        }
                        else
                        {
                            return null;
                        }
                    }
                    else
                    {
                        if (DIAparameters.PrecursorMassTolerance.Within(mass.MonoisotopicMass, masses[ind].MonoisotopicMass))
                        {
                            return masses[ind];
                        }
                        else
                        {
                            return null;
                        }
                    }
                }
            }
        }

        public static List<Peak> FindIsotopes(List<double> targetPeakMzs, List<Peak>[] peakTable, int zeroBasedScanIndex, DIAparameters DIAparameters)
        {
            var isotopes = new List<Peak>();
            foreach (var mz in targetPeakMzs)
            {
                var matchedPeak = PeakCurve.GetPeakFromScan(mz, peakTable, zeroBasedScanIndex, DIAparameters.Ms1PeakFindingTolerance, DIAparameters.PeakSearchBinSize);
                if (matchedPeak != null)
                {
                    isotopes.Add(matchedPeak);
                }
            }
            return isotopes;
        }

        public void RefineEnvelopePeakCurve(List<Peak>[] peakTable, DIAparameters DIAparameters)
        {
            foreach (var envelope in Envelopes)
            {
                var isotopes = FindIsotopes(envelope.Envelope.Peaks.Select(p => p.mz).ToList(), peakTable, envelope.ZeroBasedScanIndex, DIAparameters);
                envelope.Isotopes = isotopes;
                envelope.AdjustedTotalIntensity = isotopes.Sum(p => p.Intensity);
            }
        }
    }
}
