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
            commonParameters.TrimMs1Peaks = false;
            commonParameters.TrimMsMsPeaks = false;
            MyFileManager myFileManager = new MyFileManager(true);
            var diaDataFile = myFileManager.LoadFile(diaFile, commonParameters);
            var diaParam = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20),
                maxNumMissedScan: 1, binSize: 100, overlapRatioCutOff: 0.2, correlationCutOff: 0.5, apexRtTolerance: 0.2,
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, precursorIntensityCutOff: 300000,
                splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f);
            var diaEngine2 = new DIAEngine2(diaDataFile, commonParameters, diaParam);
            diaEngine2.Ms1PeakIndexing();
            diaEngine2.ConstructMs2Group();
            diaEngine2.GetMs1PeakCurves();

            var windows = diaEngine2.Ms1PeakCurves.Keys.ToList();
            var testPeakCurve = diaEngine2.Ms1PeakCurves[windows[0]].ToList();
            //foreach (var pc in testPeakCurve)
            //{
            //    pc.VisualizePeakRegions();
            //}

            diaEngine2.GetMs2PeakCurves();
            var testMs2PeakCurve = diaEngine2.Ms2PeakCurves.Values.SelectMany(v => v).ToList();
            foreach(var pc in testMs2PeakCurve)
            {
                pc.VisualizePeakRegions();
            }
            var numPeaks = diaEngine2.Ms2PeakCurves.Values.SelectMany(v => v).Sum(v => v.Peaks.Count);
            var num = diaDataFile.GetAllScansList().Where(s => s.MsnOrder == 2).Sum(s => s.MassSpectrum.Size);
            var testPeakCurve1 = testPeakCurve[58];
            //testPeakCurve1.DetectPeakRegions();
            //testPeakCurve1.VisualizePeakRegions();
            //plot.Show();

            for (int i = 0; i < testPeakCurve.Count; i++)
            {
                var peakCurve = testPeakCurve[i];
                peakCurve.VisualizePeakRegions();
            }
        }
    }
}
