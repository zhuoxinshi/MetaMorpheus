using Plotly.NET;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Python.Runtime.TypeSpec;
using ThermoFisher.CommonCore.Data.Business;
using MassSpectrometry;
using MzLibUtil;
using Chemistry;

namespace EngineLayer.DIA
{
    public class MassCurve : PeakCurve
    {
        //public List<DeconvolutedMass> Masses { get; set; }
        //public double ApexRtTotal => Peaks.OrderByDescending(m => (DeconvolutedMass)m.Intensity).First().RetentionTime;
        public MassCurve(List<Peak> masses)
        {
            Peaks = masses;
            Charge = masses.First().Charge;
            MonoisotopicMass = masses.First().MonoisotopicMass;
        }

        public override double AveragedMz => AverageMz();

        //public GenericChart VisualizeTotalIntensity()
        //{
        //    var plot = Chart2D.Chart.Line<double, double, string>(
        //        x: Masses.Select(m => m.RetentionTime),
        //        y: Masses.Select(m => m.TotalIntensity)).WithTraceInfo($"{Math.Round(AveragedMass, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
        //    return plot;
        //}

        //public GenericChart VisualizeHighestPeakIntensity()
        //{
        //    var plot = Chart2D.Chart.Line<double, double, string>(
        //        x: (List<DeconvolutedMass>)Peaks.Select(m => m.RetentionTime),
        //        y: (List<DeconvolutedMass>)Peaks.Select(m => (DeconvolutedMass)m.HighestPeakIntensity)).WithTraceInfo($"{Math.Round(AveragedMass, 3)}", ShowLegend: true).WithMarkerStyle(Color: Color.fromString("red"));
        //    return plot;
        //}

        public override double AverageMz()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            double averagedMz = 0;
            foreach (var peak in Peaks)
            {
                double weight = peak.Intensity / sumIntensity;
                averagedMz += weight * peak.HighestPeakMz;
            }
            return averagedMz;
        }

        public static DeconvolutedMass GetMassFromScan(DeconvolutedMass targetMass, List<DeconvolutedMass>[] massTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            DeconvolutedMass bestMass = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(targetMass.MonoisotopicMass) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(targetMass.MonoisotopicMass) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < massTable.Length && massTable[j] != null)
                {
                    List<DeconvolutedMass> bin = massTable[j];//Where(m => m.Charge == targetMass.Charge).ToList()
                    int index = BinarySearchForIndexedMass(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        DeconvolutedMass mass = bin[i];

                        if (mass.ZeroBasedScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        if (tolerance.Within(mass.MonoisotopicMass, targetMass.MonoisotopicMass) && mass.ZeroBasedScanIndex == zeroBasedScanIndex && mass.Charge == targetMass.Charge
                            && (bestMass == null || Math.Abs(mass.MonoisotopicMass - targetMass.MonoisotopicMass) < Math.Abs(bestMass.MonoisotopicMass - targetMass.MonoisotopicMass)))
                        {
                            bestMass = mass;
                        }
                    }
                }
            }
            return bestMass;
        }

        public static DeconvolutedMass GetMassFromScan(double monoMass, int charge, List<DeconvolutedMass>[] massTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            DeconvolutedMass bestMass = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(monoMass) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(monoMass) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < massTable.Length && massTable[j] != null)
                {
                    List<DeconvolutedMass> bin = massTable[j];//.Where(m => m.Charge == charge).ToList()
                    //if (bin.Count == 0) continue;
                    int index = BinarySearchForIndexedMass(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        DeconvolutedMass mass = bin[i];

                        if (mass.ZeroBasedScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        if (tolerance.Within(mass.MonoisotopicMass, monoMass) && mass.ZeroBasedScanIndex == zeroBasedScanIndex && mass.Charge == charge
                            && (bestMass == null || Math.Abs(mass.MonoisotopicMass - monoMass) < Math.Abs(bestMass.MonoisotopicMass - monoMass)))
                        {
                            bestMass = mass;
                        }
                    }
                }
            }
            return bestMass;
        }

        private static int BinarySearchForIndexedMass(List<DeconvolutedMass> massList, int zeroBasedScanIndex)
        {
            int m = 0;
            int l = 0;
            int r = massList.Count - 1;

            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                {
                    break;
                }
                if (massList[m].ZeroBasedScanIndex < zeroBasedScanIndex)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }

            for (int i = m; i >= 0; i--)
            {
                if (massList[i].ZeroBasedScanIndex < zeroBasedScanIndex)
                {
                    break;
                }
                m--;
            }

            if (m < 0)
            {
                m = 0;
            }

            return m;
        }

        public static bool ToleranceWithinNotch(double mass1, double mass2, Tolerance tolerance, int numNotches = 0)
        {
            if (tolerance.Within(mass1, mass2))
            {
                return true;
            }
            else if (numNotches > 0)
            {
                for (int notch = 1; notch <= numNotches; notch++)
                {
                    double notchStep = notch * Constants.C13MinusC12;
                    if (tolerance.Within(mass1, mass2 + notchStep) || tolerance.Within(mass1, mass2 - notchStep))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        public static MassCurve FindMassCurve(Peak targetMass, List<DeconvolutedMass>[] massTable, MsDataScan[] scans, MzRange isolationWindow, int maxMissedScans
            , Tolerance tolerance, int binSize, double maxRTrange)
        {
            var xic = new List<Peak>();
            xic.Add(targetMass);
            MassCurve newMassCurve = new MassCurve(xic); 
            targetMass.PeakCurve = newMassCurve; 

            // go right
            int missedScans = 0;
            for (int t = targetMass.ZeroBasedScanIndex + 1; t < scans.Length; t++)
            {
                var mass = GetMassFromScan((DeconvolutedMass)targetMass, massTable, t, tolerance, binSize);

                if (mass == null)
                {
                    missedScans++;
                }
                else if (mass != null)
                {
                    if (mass.PeakCurve == null) 
                    {
                        missedScans = 0;
                        xic.Add(mass);
                        mass.PeakCurve = newMassCurve;

                        if (mass.RetentionTime - targetMass.RetentionTime > maxRTrange)
                        {
                            break;
                        }
                    }
                    else
                    {
                        missedScans++;
                    }
                }

                if (missedScans > maxMissedScans)
                {
                    break;
                }
            }

            // go left
            missedScans = 0;
            for (int t = targetMass.ZeroBasedScanIndex - 1; t >= 0; t--)
            {
                var mass = GetMassFromScan((DeconvolutedMass)targetMass, massTable, t, tolerance, binSize);

                if (mass == null)
                {
                    missedScans++;
                }
                else if (mass != null)
                {
                    if (mass.PeakCurve == null) // Changed from PeakCurve to MassCurve
                    {
                        missedScans = 0;
                        xic.Add(mass);
                        mass.PeakCurve = newMassCurve; // Changed from PeakCurve to MassCurve

                        if (targetMass.RetentionTime - mass.RetentionTime > maxRTrange)
                        {
                            break;
                        }
                    }
                    else
                    {
                        missedScans++;
                    }
                }

                if (missedScans > maxMissedScans)
                {
                    break;
                }
            }
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return newMassCurve;
        }

        public static List<PeakCurve> GetAllMassCurves(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam,
            Tolerance peakFindingTolerance, double maxRTRange, double minMass, int minCharge, out List<Peak>[] allMassesByScan, bool cutPeak = false, MzRange isolationWindow = null
            )
        {
            var allMassCurves = new List<PeakCurve>();
            var deconParam = scans[0].MsnOrder == 1? commonParameters.PrecursorDeconvolutionParameters : commonParameters.ProductDeconvolutionParameters;
            allMassesByScan = DeconvolutedMass.GetAllNeutralMassesByScan(scans, deconParam, isolationWindow, minMass, minCharge);
            var allMasses = allMassesByScan.Where(v => v != null).SelectMany(p => p).ToList();
            var massTable = DeconvolutedMass.GetMassTable(allMasses, diaParam.PeakSearchBinSize);
            var sortedMasses = allMasses.OrderByDescending(p => p.Intensity).ToList();

            foreach (var mass in sortedMasses)
            {
                if (mass.PeakCurve == null)
                {
                    var newMassCurve = MassCurve.FindMassCurve((DeconvolutedMass)mass, massTable, scans, null, diaParam.MaxNumMissedScan,
                    peakFindingTolerance, diaParam.PeakSearchBinSize, maxRTRange);
                    if (cutPeak)
                    {
                        newMassCurve.CutPeak();
                    }
                    allMassCurves.Add((MassCurve)newMassCurve);
                }
            }
            return allMassCurves;
        }

        public static DeconvolutedMass GetMassFromScan(DeconvolutedMass targetMass, List<DeconvolutedMass>[,] massTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            DeconvolutedMass bestMass = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(targetMass.MonoisotopicMass) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(targetMass.MonoisotopicMass) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < massTable.Length && massTable[targetMass.Charge,j] != null)
                {
                    List<DeconvolutedMass> bin = massTable[targetMass.Charge, j];
                    int index = BinarySearchForIndexedMass(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        DeconvolutedMass mass = bin[i];

                        if (mass.ZeroBasedScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        if (tolerance.Within(mass.MonoisotopicMass, targetMass.MonoisotopicMass) && mass.ZeroBasedScanIndex == zeroBasedScanIndex && mass.Charge == targetMass.Charge
                            && (bestMass == null || Math.Abs(mass.MonoisotopicMass - targetMass.MonoisotopicMass) < Math.Abs(bestMass.MonoisotopicMass - targetMass.MonoisotopicMass)))
                        {
                            bestMass = mass;
                        }
                    }
                }
            }
            return bestMass;
        }
    }
}
