using EngineLayer.DIA.Other;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;


namespace EngineLayer.DIA
{
    public class DeconvolutedMass : Peak
    {
        //public double RetentionTime { get; set; }
        public override double HighestPeakMz => Envelope.Peaks.OrderByDescending(p => p.intensity).First().mz;   
        public override double HighestPeakIntensity => Envelope.Peaks.Max(p => p.intensity); 

        public IsotopicEnvelope Envelope { get; set; }
        public EnvelopeCurve EnvelopeCurve { get; set; }
        public override double TotalIntensity => Envelope.TotalIntensity;
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
                var envelopes = Deconvoluter.Deconvolute(scans[i].MassSpectrum, deconParameters, mzRange).OrderBy(e => e.MonoisotopicMass);
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

        public static List<DeconvolutedMass>[,] Get2DMassTable(List<Peak> allMasses, int binsPerDalton)
        {
            var maxMass = allMasses.Max(p => p.MonoisotopicMass);
            var maxCharge = allMasses.Max(p => p.Charge);
            var table = new List<DeconvolutedMass>[maxCharge + 1, (int)Math.Ceiling(allMasses.Max(p => p.MonoisotopicMass) * binsPerDalton) + 1];
            for (int i = 0; i < allMasses.Count; i++)
            {
                int roundedMass = (int)Math.Round(allMasses[i].MonoisotopicMass * binsPerDalton, 0);

                if (table[allMasses[i].Charge, roundedMass] == null)
                {
                    table[allMasses[i].Charge, roundedMass] = new List<DeconvolutedMass>();
                }
                table[allMasses[i].Charge, roundedMass].Add((DeconvolutedMass)allMasses[i]);
            }
            return table;
        }

        public static Peak MassLookUp(List<Peak>[] allMassesByScan, int scanNumber, double highestPeakMz)
        {
            var mass = allMassesByScan[scanNumber].FirstOrDefault(p => Math.Round(p.HighestPeakMz, 2) == highestPeakMz);
            return mass;
        }
    }
}

public class DBSCAN
{
    public static List<List<IsotopicEnvelope>> ClusterEnvelopes(
        List<IsotopicEnvelope> envelopes,
        double maximumNeutralMassDistancePerCluster,
        int minimumChargeStatesPerCluster,
        int maxNotches // How many +C13 notches to consider when clustering envelopes
        )
    {
        var clusters = new List<List<IsotopicEnvelope>>();
        var visited = new HashSet<IsotopicEnvelope>();
        var noise = new HashSet<IsotopicEnvelope>();

        foreach (var envelope in envelopes)
        {
            if (visited.Contains(envelope)) continue;

            visited.Add(envelope);
            var neighbors = GetNeighbors(envelope, envelopes, maximumNeutralMassDistancePerCluster, maxNotches);

            if (neighbors.Count < minimumChargeStatesPerCluster)
            {
                noise.Add(envelope);
            }
            else
            {
                var cluster = new List<IsotopicEnvelope>();
                ExpandCluster(envelope, neighbors, cluster, envelopes, visited, maximumNeutralMassDistancePerCluster, minimumChargeStatesPerCluster, maxNotches);
                clusters.Add(cluster);
            }
        }

        return clusters;
    }

    private static void ExpandCluster(
        IsotopicEnvelope envelope,
        List<IsotopicEnvelope> neighbors,
        List<IsotopicEnvelope> cluster,
        List<IsotopicEnvelope> envelopes,
        HashSet<IsotopicEnvelope> visited,
        double eps,
        int minPts, int maxNotches)
    {
        cluster.Add(envelope);

        for (int i = 0; i < neighbors.Count; i++)
        {
            var neighbor = neighbors[i];
            if (!visited.Contains(neighbor))
            {
                visited.Add(neighbor);
                var newNeighbors = GetNeighbors(neighbor, envelopes, eps, maxNotches);
                if (newNeighbors.Count >= minPts)
                {
                    neighbors.AddRange(newNeighbors);
                }
            }

            if (!cluster.Contains(neighbor))
            {
                cluster.Add(neighbor);
            }
        }
    }

    private static List<IsotopicEnvelope> GetNeighbors(
        IsotopicEnvelope envelope,
        List<IsotopicEnvelope> envelopes,
        double eps, int maxNotches)
    {
        double c13MinusC12 = Constants.C13MinusC12;
        List<IsotopicEnvelope> neighbors = new List<IsotopicEnvelope>();

        foreach (var e in envelopes)
        {
            double massDifference = Math.Abs(e.MonoisotopicMass - envelope.MonoisotopicMass);
            bool isMassMatch = massDifference <= eps;

            if (isMassMatch)
            {
                neighbors.Add(e);
            }
            else
            {
                for (int notch = 1; notch <= maxNotches; notch++)
                {
                    double c13Step = notch * c13MinusC12;
                    if (massDifference <= eps + c13Step || massDifference <= eps - c13Step)
                    {
                        neighbors.Add(e);
                        break;
                    }
                }
            }
        }
        return neighbors;
    }

//    public static (List<List<TChromosome>> Clusters, double TotalIntraClusterDistance)
//KMeans<TChromosome>(List<TChromosome> population, int k, int iterations = 10)
//where TChromosome : Chromosome
//    {
//        if (population == null || population.Count == 0)
//            throw new ArgumentException("Population cannot be null or empty.", nameof(population));
//        if (k <= 0)
//            throw new ArgumentException("Number of clusters (k) must be greater than 0.", nameof(k));
//        if (k > population.Count)
//            throw new ArgumentException("Number of clusters (k) cannot exceed the population size.", nameof(k));
//        var rnd = new Random();
//        var centroids = population.OrderBy(_ => rnd.Next()).Take(k).ToList();
//        List<List<TChromosome>> species = Enumerable.Range(0, k).Select(_ => new List<TChromosome>()).ToList();
//        for (int iter = 0; iter < iterations; iter++)
//        {
//            // Reset species
//            species = Enumerable.Range(0, k).Select(_ => new List<TChromosome>()).ToList();
//            // Assign each chromosome to the nearest centroid
//            foreach (var chr in population)
//            {
//                int bestIndex = 0;
//                double bestDist = double.MaxValue;
//                for (int i = 0; i < k; i++)
//                {
//                    double dist = EuclideanDistance(chr, centroids[i]);
//                    if (dist < bestDist)
//                    {
//                        bestDist = dist;
//                        bestIndex = i;
//                    }
//                }
//                species[bestIndex].Add(chr);
//            }
//            // Update centroids for non-empty clusters
//            for (int i = 0; i < k; i++)
//            {
//                if (species[i].Count > 0) // Only update centroid if the cluster is not empty
//                {
//                    centroids[i] = GetAverageChromosome(species[i]);
//                }
//            }
//        }
//        // Compute total intra-cluster distance
//        double totalDistance = 0.0;
//        for (int i = 0; i < k; i++)
//        {
//            if (species[i].Count > 0) // Skip empty clusters
//            {
//                var centroid = GetAverageChromosome(species[i]);
//                foreach (var chr in species[i])
//                {
//                    totalDistance += EuclideanDistance(chr, centroid);
//                }
//            }
//        }
//        return (species.Where(p => p.Count > 0).ToList(), totalDistance);
//    }

    //public static List<List<TChromosome>> KMeansSpeciationElbow<TChromosome>(List<TChromosome> population, int kMin = 2, int kMax = 10)
    //    where TChromosome : Chromosome
    //{
    //    if (population == null || population.Count == 0)
    //        throw new ArgumentException("Population cannot be null or empty.", nameof(population));
    //    if (kMin <= 0 || kMax <= 0 || kMin > kMax)
    //        throw new ArgumentException("Invalid range for kMin and kMax.");
    //    kMax = Math.Min(kMax, population.Count); // Ensure kMax does not exceed population size
    //    var distances = new List<(int K, double TotalDistance)>();
    //    for (int k = kMin; k <= kMax; k++)
    //    {
    //        var (_, totalDist) = KMeans(population, k);
    //        distances.Add((k, totalDist));
    //    }
    //    // Use a simple heuristic: find the "knee" where relative improvement drops
    //    int bestK = distances[0].K;
    //    double bestImprovement = 0;
    //    for (int i = 1; i < distances.Count - 1; i++)
    //    {
    //        double prev = distances[i - 1].TotalDistance;
    //        double curr = distances[i].TotalDistance;
    //        double next = distances[i + 1].TotalDistance;
    //        double improvement = prev - curr;
    //        double nextImprovement = curr - next;
    //        if (nextImprovement < 0.5 * improvement) // Elbow heuristic
    //        {
    //            bestK = distances[i].K;
    //            break;
    //        }
    //    }
    //    var (finalSpecies, _) = KMeans(population, bestK);
    //    for (int i = 0; i < finalSpecies.Count; i++)
    //    {
    //        var species = finalSpecies[i];
    //        foreach (var chrom in species)
    //        {
    //            chrom.SpeciesIndex = i;
    //        }
    //    }
    //    return finalSpecies;
    //}
}
