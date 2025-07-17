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
using Nett;

namespace Test.TestDIA
{
    public class Other
    {
        [Test]
        public static void SimilarityCalculation()
        {
            var libraryPath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\1614_filtered_sub_sequences_predictions.msp";
            var pathList = new List<string> { libraryPath }
;           var library = new SpectralLibrary(pathList);
            var librarySpectra = library.GetAllLibrarySpectra().ToList();

            var rawPath = @"";
            var myFileManagers = new MyFileManager(true);
            var dataFile = myFileManagers.LoadFile(rawPath, new CommonParameters());
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();

            var psmTsvList = new List<PsmFromTsv>();
            foreach(var psmTsv in psmTsvList)
            {
                var librarySpectrum = librarySpectra.FirstOrDefault(s => s.Name == psmTsv.FullSequence);
            }
        }
    }
}
