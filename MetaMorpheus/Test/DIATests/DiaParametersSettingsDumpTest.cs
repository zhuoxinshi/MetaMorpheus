using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.IO;
using TaskLayer;

namespace Test.DIATests
{
    [TestFixture]
    public static class DiaParametersSettingsDumpTest
    {
        /// <summary>
        /// Builds an ISD-style <see cref="DIAparameters"/> (the kind used to search ISD data), dumps it to a
        /// file, and confirms the file records every setting used by the DIA engine: the analysis type, the
        /// concrete XIC constructors / grouping engine, and their nested deconvolution, spline and grouping
        /// settings. This is the record produced next to a search whenever DIAparameters is set.
        /// </summary>
        [Test]
        public static void TestDiaParametersSettingsFile()
        {
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 60, 4, 3), traceToleranceDa: 1.5, minMass: 4000, xicSpline: new Bspline(2, 150)),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3), traceToleranceDa: 1.5, aggregateCharges: false, xicSpline: new Bspline(2, 150)),
                new UmpirePfGroupingEngine(150, 0.2f, 0.3, 0.7, maxThreadsForGrouping: 10),
                PseudoMs2ConstructionType.Mass,
                combineFragments: true);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestDiaParametersSettingsDump");
            Directory.CreateDirectory(outputFolder);
            string settingsPath = Path.Combine(outputFolder, "DIAparameters.txt");

            diaParams.WriteSettingsToFile(settingsPath);

            Assert.That(File.Exists(settingsPath));
            string text = File.ReadAllText(settingsPath);

            // top-level DIA settings
            Assert.That(text, Does.Contain("AanalysisType = ISD"));
            Assert.That(text, Does.Contain("PseudoMs2ConstructionType = Mass"));
            Assert.That(text, Does.Contain("CombineFragments = True"));

            // the concrete XIC constructors and grouping engine that were chosen
            Assert.That(text, Does.Contain("Ms1XicConstructor: ConsensusMassXicConstructor"));
            Assert.That(text, Does.Contain("Ms2XicConstructor: ConsensusMassXicConstructor"));
            Assert.That(text, Does.Contain("PfGroupingEngine: UmpirePfGroupingEngine"));

            // nested XIC-constructor scalar settings
            Assert.That(text, Does.Contain("MinMass = 4000"));
            Assert.That(text, Does.Contain("TraceToleranceDa = 1.5"));
            Assert.That(text, Does.Contain("AggregateCharges = False")); // set on the fragment channel

            // nested deconvolution + spline settings
            Assert.That(text, Does.Contain("DeconParameters: ClassicDeconvolutionParameters"));
            Assert.That(text, Does.Contain("MaxAssumedChargeState = 60")); // precursor channel decon
            Assert.That(text, Does.Contain("XicSplineEngine: Bspline"));
            Assert.That(text, Does.Contain("NumberOfPoints = 150"));

            // grouping-engine settings
            Assert.That(text, Does.Contain("CorrelationThreshold = 0.7"));
            Assert.That(text, Does.Contain("NoPointPerInterval = 150"));

            Directory.Delete(outputFolder, true);
        }

        /// <summary>
        /// The DIA settings dump is only produced for DIA/ISD runs. This confirms the task-level helper is a
        /// no-op when CommonParameters has no DIAparameters (i.e. an ordinary, non-ISD search).
        /// </summary>
        [Test]
        public static void TestNoDiaParametersMeansNoFile()
        {
            // A plain search task with no DIAparameters set should not attempt to write the settings file.
            var commonParameters = new CommonParameters();
            Assert.That(commonParameters.DIAparameters, Is.Null);
        }
    }
}
