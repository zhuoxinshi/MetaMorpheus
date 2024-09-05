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
            string outputFolder = @"E:\DIA\TestSearch\test_corr0.5_chargeAdded_noMassRound_highestMZ_MS1tol10ppm";
            if (!Directory.Exists(outputFolder) )
            {
                Directory.CreateDirectory(outputFolder);
            }
            task.CommonParameters.DIAparameters = new DIAparameters(new PpmTolerance(10), new PpmTolerance(20), 2, 100, 0.5, 0.5);
            //string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string myDatabase = @"E:\ISD Project\Claire's human data\Human_9606.fasta";
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { DIAfile }, "test");
        }
    }
}
