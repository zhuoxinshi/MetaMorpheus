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
    public class ISDTest
    {
        [Test]
        public static void TestConstructMs2Groups()
        {
            var file = @"E:\ISD Project\ISD_240606\06-07-24_mix_sample1_5uL_ISD.mzML";
            var commonParameters = new CommonParameters();
            commonParameters.TrimMs1Peaks = false;
            commonParameters.TrimMsMsPeaks = false;
            MyFileManager myFileManager = new MyFileManager(true);
            var myMsDataFile = myFileManager.LoadFile(file, commonParameters);
            var isdEngine = new ISDEngine(myMsDataFile, commonParameters, new DIAparameters(new PpmTolerance(5), new PpmTolerance(20), maxNumMissedScan: 1, 
                binSize: 100, overlapRatioCutOff: 0.2, correlationCutOff: 0.5, apexRtTolerance: 0.2, fragmentRankCutOff: 100, precursorRankCutOff: 10, 
                maxRTrangeMS1: 1, maxRTrangeMS2: 1, precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f));
            isdEngine.Ms1PeakIndexing();
            isdEngine.ConstructMs2Group();
            isdEngine.GetMs1PeakCurves();
            isdEngine.GetMs2PeakCurves();
            isdEngine.PrecursorFragmentPairing();
            isdEngine.PFgroupFilter();
        }
    }
}
