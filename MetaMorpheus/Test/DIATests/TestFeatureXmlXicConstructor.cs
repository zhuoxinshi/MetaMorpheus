using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.DIATests
{
    [TestFixture]
    public class TestFeatureXmlXicConstructor
    {
        // Two TopFD features in the layout of a TopFD 1.8.1 *_feature.xml (times in seconds):
        // feature 0 at 10,000 Da, charge 10 on scans 2, 6, 10 and charge 11 on scans 6, 10, 14;
        // feature 1 at 2,000 Da with a low ECScore, charge 2 on scans 2, 6, 10.
        private const string FeatureXml = @"<?xml version=""1.0"" encoding=""UTF-8""?>
<frac_feature_list>
<frac_feature>
  <file_name>fake_isd60.mzML</file_name>
  <frac_id>0</frac_id>
  <feat_id>0</feat_id>
  <mono_mass>10000.0</mono_mass>
  <intensity>1070.0</intensity>
  <time_begin>60.00</time_begin>
  <time_end>78.00</time_end>
  <scan_begin>2</scan_begin>
  <scan_end>14</scan_end>
  <min_charge>10</min_charge>
  <max_charge>11</max_charge>
  <apex_time>66.00</apex_time>
  <apex_scan>6</apex_scan>
  <apex_inte>700.0</apex_inte>
  <rep_charge>11</rep_charge>
  <rep_avg_mz>910.6</rep_avg_mz>
  <envelope_num>6</envelope_num>
  <ec_score>9.0e-01</ec_score>
  <single_charge_feature_list>
    <single_charge_feature>
      <charge>10</charge>
      <time_begin>60.00</time_begin>
      <time_end>72.00</time_end>
      <scan_begin>2</scan_begin>
      <scan_end>10</scan_end>
      <intensity>600.0</intensity>
      <envelope_num>3</envelope_num>
      <mono_mz>1001.0073</mono_mz>
      <average_mz>1001.6</average_mz>
      <refer_mz>1001.5</refer_mz>
      <max_intensity>300.0</max_intensity>
      <envelope_list>
        <envelope>
          <scan>2</scan>
          <rt>60.00</rt>
          <intensity_sum>100.0</intensity_sum>
          <max_intensity>20.0</max_intensity>
        </envelope>
        <envelope>
          <scan>6</scan>
          <rt>66.00</rt>
          <intensity_sum>300.0</intensity_sum>
          <max_intensity>60.0</max_intensity>
        </envelope>
        <envelope>
          <scan>10</scan>
          <rt>72.00</rt>
          <intensity_sum>200.0</intensity_sum>
          <max_intensity>40.0</max_intensity>
        </envelope>
      </envelope_list>
    </single_charge_feature>
    <single_charge_feature>
      <charge>11</charge>
      <time_begin>66.00</time_begin>
      <time_end>78.00</time_end>
      <scan_begin>6</scan_begin>
      <scan_end>14</scan_end>
      <intensity>470.0</intensity>
      <envelope_num>3</envelope_num>
      <mono_mz>910.0073</mono_mz>
      <average_mz>910.6</average_mz>
      <refer_mz>910.5</refer_mz>
      <max_intensity>400.0</max_intensity>
      <envelope_list>
        <envelope>
          <scan>6</scan>
          <rt>66.00</rt>
          <intensity_sum>400.0</intensity_sum>
          <max_intensity>80.0</max_intensity>
        </envelope>
        <envelope>
          <scan>10</scan>
          <rt>72.00</rt>
          <intensity_sum>50.0</intensity_sum>
          <max_intensity>10.0</max_intensity>
        </envelope>
        <envelope>
          <scan>14</scan>
          <rt>78.00</rt>
          <intensity_sum>20.0</intensity_sum>
          <max_intensity>4.0</max_intensity>
        </envelope>
      </envelope_list>
    </single_charge_feature>
  </single_charge_feature_list>
</frac_feature>
<frac_feature>
  <file_name>fake_isd60.mzML</file_name>
  <frac_id>0</frac_id>
  <feat_id>1</feat_id>
  <mono_mass>2000.0</mono_mass>
  <intensity>30.0</intensity>
  <time_begin>60.00</time_begin>
  <time_end>72.00</time_end>
  <scan_begin>2</scan_begin>
  <scan_end>10</scan_end>
  <min_charge>2</min_charge>
  <max_charge>2</max_charge>
  <apex_time>66.00</apex_time>
  <apex_scan>6</apex_scan>
  <apex_inte>10.0</apex_inte>
  <rep_charge>2</rep_charge>
  <rep_avg_mz>1001.0</rep_avg_mz>
  <envelope_num>3</envelope_num>
  <ec_score>5.0e-02</ec_score>
  <single_charge_feature_list>
    <single_charge_feature>
      <charge>2</charge>
      <time_begin>60.00</time_begin>
      <time_end>72.00</time_end>
      <scan_begin>2</scan_begin>
      <scan_end>10</scan_end>
      <intensity>30.0</intensity>
      <envelope_num>3</envelope_num>
      <mono_mz>1001.0073</mono_mz>
      <average_mz>1001.5</average_mz>
      <refer_mz>1001.5</refer_mz>
      <max_intensity>10.0</max_intensity>
      <envelope_list>
        <envelope>
          <scan>2</scan>
          <rt>60.00</rt>
          <intensity_sum>10.0</intensity_sum>
          <max_intensity>5.0</max_intensity>
        </envelope>
        <envelope>
          <scan>6</scan>
          <rt>66.00</rt>
          <intensity_sum>10.0</intensity_sum>
          <max_intensity>5.0</max_intensity>
        </envelope>
        <envelope>
          <scan>10</scan>
          <rt>72.00</rt>
          <intensity_sum>10.0</intensity_sum>
          <max_intensity>5.0</max_intensity>
        </envelope>
      </envelope_list>
    </single_charge_feature>
  </single_charge_feature_list>
</frac_feature>
</frac_feature_list>
";

        private string _featureXmlPath;

        [OneTimeSetUp]
        public void WriteFeatureXml()
        {
            _featureXmlPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "TestFeatureXmlXicConstructor_feature.xml");
            File.WriteAllText(_featureXmlPath, FeatureXml);
        }

        [OneTimeTearDown]
        public void DeleteFeatureXml()
        {
            if (File.Exists(_featureXmlPath))
                File.Delete(_featureXmlPath);
        }

        // one ISD channel: every 4th scan of the run, 0.1 min apart
        private static MsDataScan[] FakeChannelScans(int sourceVoltage, params int[] scanNumbers)
        {
            return scanNumbers.Select((n, i) => new MsDataScan(new MzSpectrum(new[] { 500.0 }, new[] { 1.0 }, false), n, 1, true,
                Polarity.Positive, 1.0 + 0.1 * i, new MzRange(400, 2000), $"FTMS + p NSI sid={sourceVoltage}.00 Full ms [400.0000-2000.0000]",
                MZAnalyzerType.Orbitrap, 1, 1.0, null, "scan=" + n)).ToArray();
        }

        [Test]
        public void TestReadFeatureXml()
        {
            var features = TopFdFeatureXmlFile.ReadFile(_featureXmlPath);
            Assert.That(features, Has.Count.EqualTo(2));

            var feature = features[0];
            Assert.That(feature.Id, Is.EqualTo(0));
            Assert.That(feature.Mass, Is.EqualTo(10000.0));
            Assert.That(feature.EcScore, Is.EqualTo(0.9).Within(1e-12));
            Assert.That(feature.RetentionTimeApex, Is.EqualTo(1.1).Within(1e-12)); // seconds converted to minutes
            // the same element names at three levels go to the right object
            Assert.That(feature.Intensity, Is.EqualTo(1070.0));
            Assert.That(feature.SingleChargeFeatures, Has.Count.EqualTo(2));
            Assert.That(feature.SingleChargeFeatures[0].Charge, Is.EqualTo(10));
            Assert.That(feature.SingleChargeFeatures[0].Intensity, Is.EqualTo(600.0));
            Assert.That(feature.SingleChargeFeatures[0].MaxIntensity, Is.EqualTo(300.0));
            Assert.That(feature.SingleChargeFeatures[0].Envelopes.Select(e => e.OneBasedScanNumber), Is.EqualTo(new[] { 2, 6, 10 }));
            Assert.That(feature.SingleChargeFeatures[0].Envelopes[1].IntensitySum, Is.EqualTo(300.0));
            Assert.That(feature.SingleChargeFeatures[0].Envelopes[1].MaxIntensity, Is.EqualTo(60.0));
            Assert.That(features[1].Id, Is.EqualTo(1));
            Assert.That(features[1].SingleChargeFeatures.Single().Envelopes, Has.Count.EqualTo(3));

            Assert.Throws<MetaMorpheusException>(() => TopFdFeatureXmlFile.ReadFile(_featureXmlPath + ".missing"));
        }

        [Test]
        public void TestFeatureXics()
        {
            var scans = FakeChannelScans(60, 2, 6, 10, 14);
            var constructor = new FeatureXmlXicConstructor(_featureXmlPath, minNumberOfPeaks: 3, minMass: 4000);

            // the 2,000 Da feature is below the mass limit
            var xic = constructor.GetAllXics(scans).Single();
            Assert.That(constructor.LastFeatureByXic[xic].Id, Is.EqualTo(0));

            // one point per scan: charge states summed, the more intense charge kept, all at the feature's mass
            Assert.That(xic.Peaks.Select(p => p.Intensity), Is.EqualTo(new[] { 100f, 700f, 250f, 20f }));
            Assert.That(xic.Peaks.Select(p => ((IndexedMass)p).Charge), Is.EqualTo(new[] { 10, 11, 10, 11 }));
            Assert.That(xic.Peaks.All(p => p.M == 10000f));
            Assert.That(xic.Peaks.Select(p => p.ZeroBasedScanIndex), Is.EqualTo(new[] { 0, 1, 2, 3 }));
            Assert.That(xic.ApexScanIndex, Is.EqualTo(1));
            Assert.That(xic.StartRT, Is.EqualTo(1.0).Within(1e-6)); // retention times of the scans, not of the xml
            Assert.That(xic.EndRT, Is.EqualTo(1.3).Within(1e-6));

            Assert.That(PrecursorFragmentsGroup.GetXicIntensity(xic, PseudoMs2IntensityType.XicSum), Is.EqualTo(1070.0));
            Assert.That(PrecursorFragmentsGroup.GetXicIntensity(xic, PseudoMs2IntensityType.ApexPeak), Is.EqualTo(700.0));

            // envelopes on scans outside the array are ignored, so only scans 2 and 6 count here
            Assert.That(constructor.GetAllXics(FakeChannelScans(60, 2, 6)), Is.Empty);
            Assert.That(constructor.LastFeatureByXic, Is.Empty);

            // the feature filters
            Assert.That(new FeatureXmlXicConstructor(_featureXmlPath, 3).GetAllXics(scans), Has.Count.EqualTo(2));
            Assert.That(new FeatureXmlXicConstructor(_featureXmlPath, 3, minEcScore: 0.5).GetAllXics(scans), Has.Count.EqualTo(1));
            Assert.That(new FeatureXmlXicConstructor(_featureXmlPath, 3, maxMass: 5000).GetAllXics(scans).Single().Peaks[0].M, Is.EqualTo(2000f));
            Assert.That(new FeatureXmlXicConstructor(_featureXmlPath, 5).GetAllXics(scans), Is.Empty);
            Assert.That(new FeatureXmlXicConstructor(_featureXmlPath, 3, minMass: 4000, maxXicWidth: 0.25).GetAllXics(scans), Is.Empty);
            Assert.That(new FeatureXmlXicConstructor(_featureXmlPath, 3, minMass: 4000, maxXicWidth: 0.35).GetAllXics(scans), Has.Count.EqualTo(1));
        }

        [Test]
        public void TestFeatureXicsBySourceVoltage()
        {
            var constructor = FeatureXmlXicConstructor.BySourceVoltage(new Dictionary<int, string> { [60] = _featureXmlPath }, 3, minMass: 4000);

            var xic = constructor.GetAllXics(FakeChannelScans(60, 2, 6, 10, 14)).Single();
            Assert.That(constructor.LastFeatureByXic[xic].Id, Is.EqualTo(0));

            // no feature file for this voltage, or scans from more than one voltage
            Assert.Throws<MetaMorpheusException>(() => constructor.GetAllXics(FakeChannelScans(80, 3, 7, 11, 15)));
            var mixedScans = FakeChannelScans(60, 2, 6).Concat(FakeChannelScans(80, 3, 7)).ToArray();
            Assert.Throws<MetaMorpheusException>(() => constructor.GetAllXics(mixedScans));
            Assert.Throws<MetaMorpheusException>(() => FeatureXmlXicConstructor.BySourceVoltage(new Dictionary<int, string>(), 3));
        }
    }
}
