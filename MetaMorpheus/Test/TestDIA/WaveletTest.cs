using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using EngineLayer;
using NUnit.Framework;
using System.IO;
using TaskLayer;
using EngineLayer.DIA;
using MzLibUtil;
using Plotly.NET;
using UsefulProteomicsDatabases;

namespace Test.TestDIA
{
    public class WaveletTest
    {
        [Test]
        public static void TestPeakSplit()
        {
            var diaFile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            var commonParameters = new CommonParameters();
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 2, 100, 0.5, 0.5, 1);
            var diaEngine2 = new DIAEngine2(diaDataFile, commonParameters, diaParam);
            diaEngine2.Ms1PeakIndexing();
            diaEngine2.GetMs1PeakCurves();

            var windows = diaEngine2.Ms1PeakCurves.Keys.ToList();
            var testPeakCurve = diaEngine2.Ms1PeakCurves[windows[0]].First(p => p.EndRT - p.StartRT > 3);
        }
    }
}
