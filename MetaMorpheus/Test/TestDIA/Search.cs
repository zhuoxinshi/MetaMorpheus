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
using static System.Net.WebRequestMethods;
using MzLibUtil;
using Plotly.NET;

namespace Test.TestDIA
{
    public class Search
    {
        [Test]
        public static void TestDIASearch()
        {
            var task = new SearchTask();
            task.CommonParameters.TrimMsMsPeaks = false;
            task.CommonParameters.TrimMs1Peaks = false;
            string outputFolder = @"E:\DIA\TestSearch\test2.0_corr0.5_highestPeakXIC_ms1Tol5ppm_cubic_apexRT0.2_maxMissed1_overlap0.2_FragRank100_preRank10_maxRT0.5_300000_cutPeakMS1&2";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(5), new PpmTolerance(20), 
                maxNumMissedScan:1, binSize: 100, overlapRatioCutOff: 0.2, correlationCutOff: 0.5, apexRtTolerance:0.2, 
                fragmentRankCutOff: 100, precursorRankCutOff: 10, maxRTrangeMS1: 1, maxRTrangeMS2: 1, highCorrThreshold: 0.5, numHighCorrFragments: 0, 
                precursorIntensityCutOff: 300000, splitMS2Peak: false, splitMS1Peak: false, splineTimeInterval: 0.005f);
            // Use reflection to set max threads
            //task.CommonParameters.GetType().GetProperty("MaxThreadsToUsePerFile").SetMethod.Invoke(task.CommonParameters, new object[] { 1 });
            //var type = task.CommonParameters.GetType();
            //var property = type.GetProperty("MaxThreadsToUsePerFile");
            //property.SetMethod.Invoke(task.CommonParameters, new object[] { 1 });

            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string myDatabase = @"E:\ISD Project\Claire's human data\Human_9606.fasta";
            string myDatabase_xml = @"E:\ISD Project\ISD_240812\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { DIAfile }, "test");
        }
    }
}
