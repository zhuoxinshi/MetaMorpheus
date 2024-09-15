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
            task.CommonParameters.TrimMsMsPeaks = true;
            task.CommonParameters.TrimMs1Peaks = false;
            string outputFolder = @"E:\DIA\TestSearch\test_corr0.7_chargeAdded_noMassRound_envelopeXIC_cubicSpline_apexRT0.5_PeakTrim_maxMissed2";
            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 2, 100, 0.5, 0.7, 0.5);
            //string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string myDatabase = @"E:\ISD Project\Claire's human data\Human_9606.fasta";
            string myDatabase_xml = @"E:\ISD Project\ISD_240812\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { DIAfile }, "test");
        }
    }
}
