using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ChargeStateEnvelopeCurve : PeakCurve
    {
        public ChargeStateEnvelopeCurve(List<Peak> chargeEnvelopes)
        {
            Peaks = chargeEnvelopes;
            MonoisotopicMass = Peaks.First().MonoisotopicMass;
        }

        public static List<PeakCurve> GetAllEnvelopeCurves(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam,
            Tolerance peakFindingTolerance, double maxRTRange, double minMass, int minCharge, out List<Peak>[] allMassesByScan, bool cutPeak = false, MzRange isolationWindow = null
            )
        {
            var allMassCurves = new List<PeakCurve>();
            var deconParam = scans[0].MsnOrder == 1 ? commonParameters.PrecursorDeconvolutionParameters : commonParameters.ProductDeconvolutionParameters;
            allMassesByScan = ChargeStateEnvelope.GetAllChargeEnvelopeByScan(scans, deconParam, isolationWindow, minMass, minCharge);
            var allMasses = allMassesByScan.Where(v => v != null).SelectMany(p => p).ToList();
            var massTable = ChargeStateEnvelope.GetEnvelopeTable(allMasses, diaParam.PeakSearchBinSize);
            var sortedMasses = allMasses.OrderByDescending(p => p.Intensity).ToList();

            foreach (var mass in sortedMasses)
            {
                if (mass.PeakCurve == null)
                {
                    var newMassCurve = ChargeStateEnvelopeCurve.FindEnvelopeCurve((ChargeStateEnvelope)mass, massTable, scans, null, diaParam.MaxNumMissedScan,
                    peakFindingTolerance, diaParam.PeakSearchBinSize, maxRTRange);
                    if (cutPeak)
                    {
                        newMassCurve.CutPeak();
                    }
                    allMassCurves.Add((ChargeStateEnvelopeCurve)newMassCurve);
                }
            }
            return allMassCurves;
        }

        public static List<PeakCurve> GetAllEnvelopeCurvesWithNotch(MsDataScan[] scans, CommonParameters commonParameters, DIAparameters diaParam,
            Tolerance peakFindingTolerance, double maxRTRange, double minMass, int minCharge, out List<Peak>[] allMassesByScan, bool cutPeak = false, MzRange isolationWindow = null
            )
        {
            var allMassCurves = new List<PeakCurve>();
            var deconParam = scans[0].MsnOrder == 1 ? commonParameters.PrecursorDeconvolutionParameters : commonParameters.ProductDeconvolutionParameters;
            allMassesByScan = ChargeStateEnvelope.GetAllChargeEnvelopeByScan(scans, deconParam, isolationWindow, minMass, minCharge);
            var allMasses = allMassesByScan.Where(v => v != null).SelectMany(p => p).ToList();
            var sortedMasses = allMasses.OrderByDescending(p => p.Intensity).ToList();
            int maxScanNum = allMassesByScan.Length - 1;

            foreach (var mass in sortedMasses)
            {
                if (mass.PeakCurve == null)
                {
                    var newMassCurve = ChargeStateEnvelopeCurve.FindEnvelopeCurveNoTable(mass, allMassesByScan, maxScanNum, peakFindingTolerance, isolationWindow, diaParam.MaxNumMissedScan, maxRTRange);
                    if (cutPeak)
                    {
                        newMassCurve.CutPeak();
                    }
                    allMassCurves.Add((ChargeStateEnvelopeCurve)newMassCurve);
                }
            }
            return allMassCurves;
        }

        public static ChargeStateEnvelope GetEnvelopeFromScan(ChargeStateEnvelope targetMass, List<ChargeStateEnvelope>[] massTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            ChargeStateEnvelope bestMass = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(targetMass.AveragedMass) * binSize);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(targetMass.AveragedMass) * binSize);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < massTable.Length && massTable[j] != null)
                {
                    List<ChargeStateEnvelope> bin = massTable[j];//Where(m => m.Charge == targetMass.Charge).ToList()
                    int index = BinarySearchForIndexedMass(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        ChargeStateEnvelope mass = bin[i];

                        if (mass.ZeroBasedScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        if (tolerance.Within(mass.AveragedMass, targetMass.AveragedMass) && mass.ZeroBasedScanIndex == zeroBasedScanIndex 
                            && (bestMass == null || Math.Abs(mass.AveragedMass - targetMass.AveragedMass) < Math.Abs(bestMass.AveragedMass - targetMass.AveragedMass)))
                        {
                            bestMass = mass;
                        }
                    }
                }
            }
            return bestMass;
        }

        public static ChargeStateEnvelope GetEnvelopeFromScan(double targetMass, List<ChargeStateEnvelope>[] massTable, int zeroBasedScanIndex, Tolerance tolerance, int binSize)
        {
            ChargeStateEnvelope bestMass = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(targetMass * binSize));
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(targetMass * binSize));

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < massTable.Length && massTable[j] != null)
                {
                    List<ChargeStateEnvelope> bin = massTable[j];//Where(m => m.Charge == targetMass.Charge).ToList()
                    int index = BinarySearchForIndexedMass(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        ChargeStateEnvelope mass = bin[i];

                        if (mass.ZeroBasedScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        if (tolerance.Within(mass.AveragedMass, targetMass) && mass.ZeroBasedScanIndex == zeroBasedScanIndex
                            && (bestMass == null || Math.Abs(mass.AveragedMass - targetMass) < Math.Abs(bestMass.AveragedMass - targetMass)))
                        {
                            bestMass = mass;
                        }
                    }
                }
            }
            return bestMass;
        }

        public static ChargeStateEnvelope GetEnvelopeFromScanNoTable(double targetMass, Peak[] allMasses, Tolerance tolerance)
        {
            ChargeStateEnvelope bestMass = null;

            foreach(var mass in allMasses)
            {
                if (mass.PeakCurve == null)
                {
                    if (MassCurve.ToleranceWithinNotch(targetMass, mass.MonoisotopicMass, tolerance, 3)
                        && (bestMass == null || Math.Abs(bestMass.AveragedMass - targetMass) > Math.Abs(((ChargeStateEnvelope)mass).AveragedMass - targetMass)))
                    {
                        bestMass = ((ChargeStateEnvelope)mass);
                    }
                }
            }
            return bestMass;
        }

        public static ChargeStateEnvelopeCurve FindEnvelopeCurveNoTable(Peak targetMass, List<Peak>[] massesByScan, int maxScanNum, Tolerance tolerance,
            MzRange isolationWindow, int maxMissedScans, double maxRTrange)
        {
            var xic = new List<Peak>();
            xic.Add(targetMass);
            ChargeStateEnvelopeCurve newMassCurve = new ChargeStateEnvelopeCurve(xic);
            targetMass.PeakCurve = newMassCurve;

            // go right
            int missedScans = 0;
            for (int t = targetMass.ScanNumber + 4; t < maxScanNum; t+=4)
            {
                if ( massesByScan[t] == null || massesByScan[t].Count == 0)
                {
                    missedScans++;
                }
                var mass = GetEnvelopeFromScanNoTable(newMassCurve.AveragedMass, massesByScan[t].ToArray(), tolerance);

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
                if (newMassCurve.EndRT - newMassCurve.ApexRT > maxRTrange)
                {
                    break;
                }
            }

            // go left
            missedScans = 0;
            for (int t = targetMass.ScanNumber - 4; t >= 1; t-=4)
            {
                if (massesByScan[t] == null || massesByScan[t].Count == 0)
                {
                    missedScans++;
                }
                var mass = GetEnvelopeFromScanNoTable(newMassCurve.AveragedMass, massesByScan[t].ToArray(), tolerance);

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
                if (newMassCurve.ApexRT - newMassCurve.StartRT > maxRTrange)
                {
                    break;
                }
            }
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return newMassCurve;
        }

        public static ChargeStateEnvelopeCurve FindEnvelopeCurve(Peak targetMass, List<ChargeStateEnvelope>[] massTable, MsDataScan[] scans, MzRange isolationWindow, int maxMissedScans
            , Tolerance tolerance, int binSize, double maxRTrange)
        {
            var xic = new List<Peak> { targetMass };
            ChargeStateEnvelopeCurve newMassCurve = new ChargeStateEnvelopeCurve(xic);
            targetMass.PeakCurve = newMassCurve;

            // go right
            int missedScans = 0;
            for (int t = targetMass.ZeroBasedScanIndex + 1; t < scans.Length; t++)
            {
                var mass = GetEnvelopeFromScan(newMassCurve.AveragedMass, massTable, t, tolerance, binSize);

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
                if (newMassCurve.EndRT - newMassCurve.ApexRT > maxRTrange)
                {
                    break;
                }
            }

            // go left
            missedScans = 0;
            for (int t = targetMass.ZeroBasedScanIndex - 1; t >= 0; t--)
            {
                var mass = GetEnvelopeFromScan(newMassCurve.AveragedMass, massTable, t, tolerance, binSize);

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
                if (newMassCurve.ApexRT - newMassCurve.StartRT > maxRTrange)
                {
                    break;
                }
            }
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return newMassCurve;
        }

        private static int BinarySearchForIndexedMass(List<ChargeStateEnvelope> massList, int zeroBasedScanIndex)
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
    }
}
