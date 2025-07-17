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

namespace Test.TestDIA
{
    public class Other
    {
        [Test]
        public static void SimilarityCalculation()
        {
            var libraryPath = @"C:\Users\Zhuoxin Shi\Downloads\example_predictions.msp";
            var pathList = new List<string> { libraryPath }
;           var library = new SpectralLibrary(pathList);
            var librarySpectra = library.GetAllLibrarySpectra().ToList();
        }
    }
}
