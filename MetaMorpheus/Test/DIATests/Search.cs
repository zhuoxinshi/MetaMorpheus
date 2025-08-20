using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;
using TaskLayer;
using System.IO;
using Chemistry;
using Readers;
using MathNet.Numerics.Interpolation;
using System.Drawing.Imaging;

namespace Test.DIATests
{
    public class Search
    {
        [Test]
        public static void SearchCE()
        {
            var peakList1 = new List<IIndexedPeak>();
            double[] intensityMultipliers = { 1, 3, 6, 3, 1 };
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList1.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
            }
            var xic1 = new ExtractedIonChromatogram(peakList1);
            var peakList2 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList2.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 10 + i / 10, zeroBasedScanIndex: i + 10, mz: 501.0));
            }
            var xic2 = new ExtractedIonChromatogram(peakList2);
            var sharedXic = PrecursorFragmentPair.CalculateSharedXIC(xic1, xic2);
        }
    }
}
