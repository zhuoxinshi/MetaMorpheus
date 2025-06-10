using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using System.Runtime.CompilerServices;

namespace EngineLayer.DIA
{
    public class ChargeStateEnvelope : Peak
    {

        public List<IsotopicEnvelope> Envelopes { get; set; }
        public override double AveragedMass => AverageMass();
        public override double MonoisotopicMass => AverageMass();
        public override double HighestPeakMz => Envelopes.SelectMany(e => e.Peaks).OrderByDescending(p => p.intensity).First().mz;
        public override double Intensity => Envelopes.Sum(e => e.Peaks.Max(p => p.intensity));

        public ChargeStateEnvelope(List<IsotopicEnvelope> envelopes)
        {
            Envelopes = envelopes;
        }

        public ChargeStateEnvelope(IsotopicEnvelope envelope)
        {
            Envelopes = new List<IsotopicEnvelope>();
            Envelopes.Add(envelope);
            Charge = 1;
        }

        public double AverageMass()
        {
            if (Envelopes.Count == 0)
                return 0;
            double averageMass = 0;
            double sumIntensity = Envelopes.Sum(e => e.TotalIntensity);
            foreach (var envelope in Envelopes)
            {
                averageMass += envelope.MonoisotopicMass * envelope.TotalIntensity/ sumIntensity;
            }
            return averageMass;
        }

        public static List<Peak>[] GetAllChargeEnvelopeByScan(MsDataScan[] scans, DeconvolutionParameters deconParameters, MzRange mzRange = null, double minMass = 0, int minCharge = 1)
        {
            var maxScanNum = scans[scans.Length - 1].OneBasedScanNumber;
            var chargeEnvelopesByScan = new List<Peak>[maxScanNum + 1];
            int index = 0;
            for (int i = 0; i < scans.Length; i++)
            {
                chargeEnvelopesByScan[scans[i].OneBasedScanNumber] = new List<Peak>();
                var envelopes = Deconvoluter.Deconvolute(scans[i].MassSpectrum, deconParameters, mzRange).Where(e => e.MonoisotopicMass >= minMass && e.Charge >= minCharge).ToList();
                var chargeEnvelopes = Cluster(envelopes, new PpmTolerance(50), 3, 4);
                foreach(var chargeEnvelope in chargeEnvelopes)
                {
                    chargeEnvelope.RetentionTime = scans[i].RetentionTime;
                    chargeEnvelope.ScanNumber = scans[i].OneBasedScanNumber;
                    chargeEnvelope.ZeroBasedScanIndex = i;
                    chargeEnvelope.MsLevel = scans[i].MsnOrder;
                    chargeEnvelope.Index = index++;
                    chargeEnvelopesByScan[scans[i].OneBasedScanNumber].Add(chargeEnvelope);
                }
            }
            return chargeEnvelopesByScan;
        }

        public static List<ChargeStateEnvelope>[] GetEnvelopeTable(List<Peak> allMasses, int binsPerDalton)
        {
            var maxMass = allMasses.Max(p => p.MonoisotopicMass);
            var table = new List<ChargeStateEnvelope>[(int)Math.Ceiling(allMasses.Max(p => p.MonoisotopicMass) * binsPerDalton) + 1];
            for (int i = 0; i < allMasses.Count; i++)
            {
                int roundedMass = (int)Math.Round(allMasses[i].MonoisotopicMass * binsPerDalton, 0);

                if (table[roundedMass] == null)
                {
                    table[roundedMass] = new List<ChargeStateEnvelope>();
                }
                table[roundedMass].Add((ChargeStateEnvelope)allMasses[i]);
            }
            return table;
        }

        public static List<ChargeStateEnvelope> Cluster(List<IsotopicEnvelope> envelopes, PpmTolerance massTolerance, int numNotches, int minNumberInCluster)
        {
            var orderedEnvelopes = envelopes.OrderByDescending(e => e.TotalIntensity).ToArray();
            var clusters = new List<ChargeStateEnvelope>();
            var visited = new HashSet<int>();
            for (int i = 0; i < orderedEnvelopes.Length; i++)
            {
                if (visited.Count == envelopes.Count)
                    break;
                if (visited.Contains(i))
                    continue;
                var cluster = new ChargeStateEnvelope(orderedEnvelopes[i]);
                for (int j = i + 1; j < orderedEnvelopes.Length; j++)
                {
                    if (visited.Contains(j))
                        continue;
                    if (cluster.AddEnvelope(orderedEnvelopes[j], massTolerance, numNotches))
                    {
                        visited.Add(j);
                    }
                }
                if (cluster.Envelopes.Count > minNumberInCluster)
                {
                    clusters.Add(cluster);
                    visited.Add(i);
                }
            }
            return clusters;
        }

        private static bool isNeighbor(double mass1, double mass2, PpmTolerance massTolerance, int numNotches)
        {
            if (massTolerance.Within(mass1, mass2))
            {
                return true;
            }
            else
            {
                for (int notch = 1; notch <= numNotches; notch++)
                {
                    double notchStep = notch * Constants.C13MinusC12;
                    if (massTolerance.Within(mass1, mass2 + notchStep) || massTolerance.Within(mass1, mass2 - notchStep))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        public bool AddEnvelope(IsotopicEnvelope envelope, PpmTolerance massTolerance, int numNotches)
        {
            if (Envelopes.Select(e => e.Charge).Contains(envelope.Charge))
            {
                return false;
            }
            if (isNeighbor(AveragedMass, envelope.MonoisotopicMass, massTolerance, numNotches))
            {
                Envelopes.Add(envelope);
                return true;
            }
            else
            {
                return false;
            }
        }
    }

    
}
