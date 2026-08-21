using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace EngineLayer.DIA
{
    /// <summary>
    /// PrecursorFragmentsGroup represents a group of precursor-fragment pairs belonging to the same precursor XIC.
    /// </summary>
    public class PrecursorFragmentsGroup
    {
        public ExtractedIonChromatogram PrecursorXic { get; set; }

        public List<PrecursorFragmentPair> PFpairs { get; set; }
        public int PFgroupIndex { get; set; } //an index identifier for each group and also used as the scan number for the pseudo MS2 scan generated from this group

        /// <summary>
        /// Id of the deconvoluted MS1 feature this group's precursor XIC came from, or null when the XIC was not
        /// built from a feature file. When the precursor XICs come from <see cref="FeatureXmlXicConstructor"/>
        /// this is TopFD's own <c>feat_id</c>, so a pseudo MS2 scan can be joined straight back to the row in
        /// <c>*_feature.xml</c> that produced it. Note 0 is a VALID feature id, which is why "unknown" is null
        /// rather than 0 or -1.
        /// </summary>
        public int? PrecursorFeatureId { get; set; }

        public PrecursorFragmentsGroup(ExtractedIonChromatogram precursorXic, List<PrecursorFragmentPair> pfPairs)
        {
            PrecursorXic = precursorXic;
            PFpairs = pfPairs;
        }

        /// <summary>
        /// Stamps each group with the id of the TopFD feature its precursor XIC came from, using the
        /// <see cref="FeatureXmlXicConstructor.LastFeatureByXic"/> map of the constructor that built the
        /// PRECURSOR XICs. Groups whose precursor is not in the map are left unset and will export an empty
        /// feature id.
        ///
        /// Note that map is rebuilt on every <c>GetAllXics</c> call, so capture it from the precursor
        /// constructor before that same object is used for anything else.
        /// </summary>
        public static void AssignPrecursorFeatureIds(IEnumerable<PrecursorFragmentsGroup> groups,
            IReadOnlyDictionary<ExtractedIonChromatogram, TopFdFeature> featureByPrecursorXic)
        {
            if (groups == null || featureByPrecursorXic == null)
                return;

            foreach (var group in groups)
            {
                if (featureByPrecursorXic.TryGetValue(group.PrecursorXic, out var feature))
                    group.PrecursorFeatureId = feature.Id;
            }
        }

        /// <summary>
        /// The intensity an XIC contributes to a pseudo MS2 scan, under the requested convention. See
        /// <see cref="PseudoMs2IntensityType"/> for why the choice is not merely one of units.
        /// </summary>
        public static double GetXicIntensity(ExtractedIonChromatogram xic, PseudoMs2IntensityType intensityType)
        {
            switch (intensityType)
            {
                case PseudoMs2IntensityType.ApexPeak:
                    return xic.ApexPeak.Intensity;
                case PseudoMs2IntensityType.XicSum:
                    // each point is already a summed isotopic envelope, so this is the feature's total signal
                    // over the scans that were handed to the XIC constructor
                    return xic.Peaks.Sum(p => (double)p.Intensity);
                default:
                    throw new MetaMorpheusException("Invalid pseudo MS2 intensity type specified.");
            }
        }

        public static double CalculateXicCorrelation(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if (xic1.XYData == null || xic2.XYData == null)
            {
                return CalculateXicCorrelationRawPeaks(xic1, xic2);
            }
            double corr = CalculateCorrelationXYData(xic1.XYData, xic2.XYData);
            return corr;
        }

        public static double CalculateXicCorrelationRawPeaks(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            var start = Math.Max(xic1.StartScanIndex, xic2.StartScanIndex);
            var end = Math.Min(xic1.EndScanIndex, xic2.EndScanIndex);

            if (end - start + 1 < 3) // Ideally we need at least 3 points to calculate correlation
            {
                return double.NaN;
            }
            var y1 = new float[end - start + 1];
            var y2 = new float[end - start + 1];
            var peakScanIndices1 = xic1.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            var peakScanIndices2 = xic2.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            for (int i = 0; i < y1.Length; i++)
            {
                int index1 = Array.BinarySearch(peakScanIndices1, i + start);
                int index2 = Array.BinarySearch(peakScanIndices2, i + start);
                if (index1 >= 0)
                {
                    y1[i] = xic1.Peaks[index1].Intensity;
                }
                else
                {
                    y1[i] = 0;
                }
                if (index2 >= 0)
                {
                    y2[i] = xic2.Peaks[index2].Intensity;
                }
                else
                {
                    y2[i] = 0;
                }
            }
            return PearsonCorrelation(y1, y2);
        }

        public static double CalculateXicCorrXYData_Umpire(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2, int NoPointPerInterval)
        {
            double start = Math.Max(xic1.XYData[0].Item1, xic2.XYData[0].Item1);

            int num = Math.Max(xic1.XYData.Count(), xic2.XYData.Count()) / 2;
            double timeInterval = (double)2 / NoPointPerInterval;

            if (num < 6)
            {
                return double.NaN;
            }
            double[] arrayA = new double[num];
            double[] arrayB = new double[num];

            int i = 0;
            double low = start;
            double up = start + timeInterval;

            for (int j = 0; j < xic1.XYData.Length; j++)
            {
                while (xic1.XYData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (xic1.XYData[j].Item1 >= low && xic1.XYData[j].Item1 < up)
                {
                    if (xic1.XYData[j].Item2 > arrayA[i])
                    {
                        arrayA[i] = xic1.XYData[j].Item2;
                    }
                }
            }

            i = 0;
            low = start;
            up = start + timeInterval;

            for (int j = 0; j < xic2.XYData.Length; j++)
            {
                while (xic2.XYData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (xic2.XYData[j].Item1 >= low && xic2.XYData[j].Item1 < up)
                {
                    if (xic2.XYData[j].Item2 > arrayB[i])
                    {
                        arrayB[i] = xic2.XYData[j].Item2;
                    }
                }
            }

            for (int idx = 1; idx < num - 1; idx++)
            {
                if (arrayA[idx] == 0f)
                {
                    arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
                }
                if (arrayB[idx] == 0f)
                {
                    arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
                }
            }

            double corr = Correlation.Pearson(arrayA, arrayB);
            return corr;
        }

        public static double CalculateCorrelationXYData<T> ((T, T)[] xy1, (T, T)[] xy2) where T : INumber<T>
        {
            T start = T.Max(xy1[0].Item1, xy2[0].Item1);
            T end = T.Min(xy1[xy1.Length - 1].Item1, xy2[xy2.Length - 1].Item1);

            var validxy1 = xy1.Where(p => p.Item1 >= start && p.Item1 <= end).ToArray();
            var validxy2 = xy2.Where(p => p.Item1 >= start && p.Item1 <= end).ToArray();
            int numPoints = Math.Min(validxy1.Length, validxy2.Length);
            if (numPoints < 3)
            {
                return double.NaN;
            }
            var xy = validxy1.Take(numPoints).Zip(validxy2.Take(numPoints), (a, b) => (a.Item2, b.Item2)).ToArray();
            var y1 = xy.Select(p => p.Item1).ToArray();
            var y2 = xy.Select(p => p.Item2).ToArray();
            double corr = PearsonCorrelation(y1, y2);

            return corr;
        }

        public static double PearsonCorrelation<T>(T[] x, T[] y) where T : INumber<T>
        {
            if (x.Length != y.Length) throw new ArgumentException("Arrays must have the same length.");

            int n = x.Length;
            T sumX = T.Zero, sumY = T.Zero, sumXY = T.Zero, sumX2 = T.Zero, sumY2 = T.Zero;
            for (int i = 0; i < n; i++)
            {
                sumX += x[i];
                sumY += y[i];
                sumXY += x[i] * y[i];
                sumX2 += x[i] * x[i];
                sumY2 += y[i] * y[i];
            }

            double meanX = double.CreateChecked(sumX) / n;
            double meanY = double.CreateChecked(sumY) / n;

            double numerator = double.CreateChecked(sumXY) - n * meanX * meanY;
            double denominatorX = double.CreateChecked(sumX2) - n * meanX * meanX;
            double denominatorY = double.CreateChecked(sumY2) - n * meanY * meanY;
            double denominator = Math.Sqrt(denominatorX * denominatorY);

            if (denominator == 0) return 0; // Avoid division by zero

            return numerator / denominator;
        }

        public static double CalculateXicOverlapRatio(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if (xic1.EndScanIndex <= xic2.StartScanIndex || xic2.EndScanIndex <= xic1.StartScanIndex)
            {
                return 0;
            }
            var start = Math.Min(xic1.StartScanIndex, xic2.StartScanIndex);
            var end = Math.Max(xic1.EndScanIndex, xic2.EndScanIndex);
            var overlapStart = Math.Max(xic1.StartScanIndex, xic2.StartScanIndex);
            var overlapEnd = Math.Min(xic1.EndScanIndex, xic2.EndScanIndex);

            double overlap = (overlapEnd - overlapStart) / (double)(end - start);
            return overlap;
        }

        public static double CalculateXicOverlapRatio_Umpire(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            double overlap = 0;
            var ms1rtrange = xic1.EndRT - xic1.StartRT;
            var ms2rtrange = xic2.EndRT - xic2.StartRT;
            if (xic1.StartRT >= xic2.StartRT && xic1.StartRT <= xic2.EndRT && xic1.EndRT >= xic2.EndRT)
            {
                overlap = (xic2.EndRT - xic1.StartRT) / ms1rtrange;
            }
            else if (xic1.EndRT >= xic2.StartRT && xic1.EndRT <= xic2.EndRT && xic1.StartRT <= xic2.StartRT)
            {
                overlap = (xic1.EndRT - xic2.StartRT) / ms1rtrange;
            }
            else if (xic1.StartRT <= xic2.StartRT && xic1.EndRT >= xic2.EndRT)
            {
                overlap = ms2rtrange / ms1rtrange;
            }
            else if (xic1.StartRT >= xic2.StartRT && xic1.EndRT <= xic2.EndRT)
            {
                overlap = 1;
            }
            return overlap;
        }

        public Ms2ScanWithSpecificMass GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType pseudoMs2ConstructionType, CommonParameters commonParameters, string dataFilePath,
            PseudoMs2IntensityType pseudoMs2IntensityType = PseudoMs2IntensityType.XicSum)
        {
            //sort all fragment XICs by their m/z value of the apex peak
            PFpairs.Sort((a, b) => a.FragmentXic.ApexPeak.M.CompareTo(b.FragmentXic.ApexPeak.M));
            //We use the apex peak of each Xic to construct the pseudo MS2 scan (might change if other values work better)
            var mzs = PFpairs.Select(pf => (double)pf.FragmentXic.ApexPeak.M).ToArray();
            //the mass comes from the apex, but the INTENSITY follows pseudoMs2IntensityType - by default the whole
            //XIC, since a pseudo MS2 peak stands for a coeluting fragment rather than for one moment in time
            var intensities = PFpairs.Select(pf => GetXicIntensity(pf.FragmentXic, pseudoMs2IntensityType)).ToArray();
            var newMs2Scan = new MsDataScan(new MzSpectrum(mzs, intensities, false), PFgroupIndex, 2, true, Polarity.Positive, PrecursorXic.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null, MZAnalyzerType.Unknown, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: PFgroupIndex);

            // assign neutral experimental fragments for the new MS2 scan
            IsotopicEnvelope[] neutralExperimentalFragments = null;
            switch (pseudoMs2ConstructionType)
            {
                // If the fragment XIC is mzPeak-based, we deconvolute the scan as normal
                case PseudoMs2ConstructionType.MzPeak:
                    neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, commonParameters);
                    break;
                // If the fragment XIC is mass-based, we create an isotopic envelope as the neutral mass fragment for each fragment XIC
                case PseudoMs2ConstructionType.Mass:
                    //these envelopes are what WriteMsAlignFile emits, so the real intensity has to reach them -
                    //passing a constant here is what made every intensity in an exported msAlign read 1
                    neutralExperimentalFragments = PFpairs.Select((pf, i) => new IsotopicEnvelope(1,
                            new List<(double mz, double intensity)> { (pf.FragmentXic.ApexPeak.M.ToMz(GetFragmentCharge(pf)), intensities[i]) },
                            pf.FragmentXic.ApexPeak.M, GetFragmentCharge(pf), intensities[i], 0)).ToArray();
                    break;
                default:
                    throw new MetaMorpheusException("Invalid pseudo MS2 construction type specified.");
            }
            // add precursor information
            var charge = PrecursorXic.ApexPeak is IndexedMass im ? im.Charge : 1;
            var monoMz = PrecursorXic.ApexPeak.M.ToMz(charge);
            //measured on the same convention as the fragments, so precursor and fragment intensities in an
            //exported msAlign are on one scale
            var precursorIntensity = GetXicIntensity(PrecursorXic, pseudoMs2IntensityType);
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge, dataFilePath,
                commonParameters, neutralExperimentalFragments, precursorIntensity: precursorIntensity)
            {
                PrecursorFeatureId = PrecursorFeatureId
            };

            return scanWithprecursor;
        }

        /// <summary>Charge of a fragment XIC, taken from its peaks; 1 when the peaks carry no charge.</summary>
        private static int GetFragmentCharge(PrecursorFragmentPair pf)
        {
            var first = pf.FragmentXic.Peaks.FirstOrDefault() as IndexedMass;
            return first != null && first.Charge != 0 ? first.Charge : 1;
        }

        public void SetFragmentRankForPfPairs()
        {
            if (PFpairs.IsNullOrEmpty())
            {
                throw new MetaMorpheusException("The PfGroup does not contain any PFpairs.");
            }
            PFpairs.Sort((a, b) => b.Correlation.Value.CompareTo(a.Correlation.Value));
            for (int i = 0; i < PFpairs.Count; i++)
            {
                PFpairs[i].FragmentRank = i + 1;
            }
        }
    }
}
