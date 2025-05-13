using EngineLayer.DIA.Other;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;


namespace EngineLayer.DIA
{
    public class DeconvolutedMass : Peak
    {
        //public double RetentionTime { get; set; }
        public int MsLevel { get; set; }
        public override double HighestPeakMz => Envelope.Peaks.OrderByDescending(p => p.intensity).First().mz;   
        public override double HighestPeakIntensity => Envelope.Peaks.Max(p => p.intensity); 

        public IsotopicEnvelope Envelope { get; set; }
        public EnvelopeCurve EnvelopeCurve { get; set; }
        public override double TotalIntensity => Envelope.Peaks.Sum(p => p.intensity);
        public List<Peak> Isotopes { get; set; }
        public double AdjustedTotalIntensity { get; set; }
        public MassCurve MassCurve { get; set; }

        public DeconvolutedMass(IsotopicEnvelope envelope, double rt, int msLevel, int scanNumber, int zeroBasedScanNum, double basePeakIntensity = 1)
        {
            Envelope = envelope;
            RetentionTime = rt;
            Charge = envelope.Charge;
            MsLevel = msLevel;
            ScanNumber = scanNumber;
            ZeroBasedScanIndex = zeroBasedScanNum;
            MonoisotopicMass = envelope.MonoisotopicMass;
            Intensity = envelope.Peaks.Max(p => p.intensity);
            SN = Intensity / basePeakIntensity;
        }

        public static List<Peak>[] GetAllNeutralMassesByScan(MsDataScan[] scans, DeconvolutionParameters deconParameters, MzRange mzRange = null, double minMass = 0, int minCharge = 1)
        {
            var maxScanNum = scans[scans.Length - 1].OneBasedScanNumber;
            var massesByScan = new List<Peak>[maxScanNum + 1];
            int index = 0;
            for (int i = 0; i < scans.Length; i++)
            {
                massesByScan[scans[i].OneBasedScanNumber] = new List<Peak>();
                var envelopes = Deconvoluter.Deconvolute(scans[i].MassSpectrum, deconParameters, mzRange);
                foreach(var envelope in envelopes)
                {
                    if (envelope.MonoisotopicMass < minMass || envelope.Charge < minCharge)
                        continue;
                    DeconvolutedMass newMass = new DeconvolutedMass(envelope, scans[i].RetentionTime, scans[i].MsnOrder, scans[i].OneBasedScanNumber, i, scans[i].MassSpectrum.YofPeakWithHighestY.Value);
                    massesByScan[scans[i].OneBasedScanNumber].Add(newMass);
                    index++;
                }
            }
            return massesByScan;
        }

        public static List<DeconvolutedMass>[] GetMassTable(List<Peak> allMasses, int binsPerDalton)
        {
            var maxMass = allMasses.Max(p => p.MonoisotopicMass);
            var table = new List<DeconvolutedMass>[(int)Math.Ceiling(allMasses.Max(p => p.MonoisotopicMass) * binsPerDalton) + 1];
            for (int i = 0; i < allMasses.Count; i++)
            {
                int roundedMass = (int)Math.Round(allMasses[i].MonoisotopicMass * binsPerDalton, 0);

                if (table[roundedMass] == null)
                {
                    table[roundedMass] = new List<DeconvolutedMass>();
                }
                table[roundedMass].Add((DeconvolutedMass)allMasses[i]);
            }
            return table;
        }

        //public static List<Precursor> FindPrecursor(Precursor precursor, List<Precursor> allPrecursors, int maxMissedScans, double massTolerance)
        //{
        //    var precursorList = new List<Precursor>(); 
        //    precursorList.Add(precursor);
        //    int numScans = allPrecursors.Max(p => p.ZeroBasedScanIndex) + 1;


    }
}
