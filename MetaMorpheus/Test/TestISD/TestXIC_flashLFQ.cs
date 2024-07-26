using NUnit.Framework;
using EngineLayer;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using MzLibUtil;
using Newtonsoft.Json.Linq;
using Proteomics.ProteolyticDigestion;
using Chemistry;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;
using System.IO;
using System.Text.RegularExpressions;
using System.Collections.Concurrent;
using EngineLayer.ISD;
using MathNet.Numerics.Distributions;
using Omics.Fragmentation;
using FlashLFQ;

namespace Test.TestISD
{
    public class TestXIC_flashLFQ
    {
        [Test]
        public static void TestConstructionXIC()
        {
            string filePath = @"E:\ISD Project\TestIsdDataAnalysis\06-07-24_mix_1pmol_5uL_ISD_RT32.16-35.59.mzML";
            string tomlFile = @"E:\ISD Project\ISD_240606\sample1-to-12_DDA&ISD_xml\Task Settings\Task1-SearchTaskconfig.toml";
            MyFileManager myFileManager = new MyFileManager(true);
            SearchTask task = Toml.ReadFile<SearchTask>(tomlFile, MetaMorpheusTask.tomlConfig);
            task.CommonParameters.DoDIA = true;
            var myMsDataFile = myFileManager.LoadFile(filePath, task.CommonParameters);

            var spectraFileInfo = new SpectraFileInfo(filePath, "", 0,0,0);
            var peaks = new PeakIndexingEngine();
        }
    }
}
