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
using System.Security.Cryptography.Xml;
using Chemistry;
using Proteomics.ProteolyticDigestion;
using Easy.Common.Extensions;
using NUnit.Framework.Constraints;

namespace Test.TestDIA
{
    public class Random
    {
        [Test]
        public static void TestXIC()
        {

        }

        [Test]
        public static void DeconTest()
        {
            var task = new SearchTask();
            var tdTask = new SearchTask()
            {
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams("top-down"))
            };
            string DIAfile = @"E:\DIA\FragPipe\DIA\CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT.mzML";
            string tdFile = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_DDA.raw";
            var fm = new MyFileManager(false);
            var dataFile = fm.LoadFile(DIAfile, task.CommonParameters);
            var ddaFile = fm.LoadFile(tdFile, tdTask.CommonParameters);
            var ms1scans = dataFile.GetMS1Scans().ToList();

            var range = ms1scans.FirstOrDefault().ScanWindowRange;
            var startMz = range.Minimum;
            var endMz = range.Maximum;
            var step = 50;
            var windowList = new List<MzRange>();
            for (double i = startMz; i < endMz; i += step)
            {
                var window = new MzRange(i - 5, i + step + 5);
                windowList.Add(window);
            }

            var fullWindowEnvelopes = new List<IsotopicEnvelope>();
            foreach(var ms1 in ms1scans)
            {
                var envelopes = Deconvoluter.Deconvolute(ms1, task.CommonParameters.PrecursorDeconvolutionParameters);
                fullWindowEnvelopes.AddRange(envelopes);
            }

            var smallWindowEnvelopes = new List<IsotopicEnvelope>();
            foreach(var ms1 in ms1scans)
            {
                foreach(var window in windowList)
                {
                    var envelopes = Deconvoluter.Deconvolute(ms1, task.CommonParameters.PrecursorDeconvolutionParameters, window);
                    smallWindowEnvelopes.AddRange(envelopes);
                }
            }

            var ms1test = ms1scans.Where(p => p.OneBasedScanNumber == 118).FirstOrDefault();
            var envelopes_full = Deconvoluter.Deconvolute(ms1test, task.CommonParameters.PrecursorDeconvolutionParameters);
            var envelopes_small = new List<IsotopicEnvelope>();
            foreach(var window in windowList)
            {
                var envelopes = Deconvoluter.Deconvolute(ms1test, task.CommonParameters.PrecursorDeconvolutionParameters, window);
                envelopes_small.AddRange(envelopes);
            }
            envelopes_small = envelopes_small.OrderBy(p => p.MonoisotopicMass.ToMz(p.Charge)).ToList();
            var diaTestScan2 = ms1scans.Where(p => p.OneBasedScanNumber == 1663).FirstOrDefault();
            var test20mz = Deconvoluter.Deconvolute(diaTestScan2, task.CommonParameters.PrecursorDeconvolutionParameters, new MzRange(613, 631));
            var test5mz = new List<IsotopicEnvelope>();
            var windows5mz = new List<MzRange> { new MzRange(613, 618), new MzRange(616, 621), new MzRange(619, 624), new MzRange(622, 627), new MzRange(625, 631) };
            var windows10mz = new List<MzRange> { new MzRange(613, 623), new MzRange(621, 631)};
            var test10mz = new List<IsotopicEnvelope>();
            foreach (var window in windows10mz)
            {
                var envelopes = Deconvoluter.Deconvolute(diaTestScan2, task.CommonParameters.PrecursorDeconvolutionParameters, window);
                test10mz.AddRange(envelopes);
            }

            var tdMs1 = ddaFile.GetMS1Scans().Where(p => p.OneBasedScanNumber == 1280).FirstOrDefault();
            var tdDeconParam = new ClassicDeconvolutionParameters(10, 40, 5, 3);
            var fullwindowTd = Deconvoluter.Deconvolute(tdMs1, tdDeconParam).OrderBy(e => e.Peaks.FirstOrDefault().mz).ToList();
            var smallwindowTd = new List<IsotopicEnvelope>();
            foreach (var window in windowList)
            {
                var envelopes = Deconvoluter.Deconvolute(tdMs1, tdDeconParam, window).OrderBy(e => e.Peaks.FirstOrDefault().mz).ToList();
                if (envelopes != null)
                {
                    smallwindowTd.AddRange(envelopes);
                }
            }

            
        }

    }
}
