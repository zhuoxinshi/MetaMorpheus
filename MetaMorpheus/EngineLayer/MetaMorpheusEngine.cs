using Chemistry;
using MassSpectrometry;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using Transcriptomics.Digestion;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {
        protected static readonly Dictionary<DissociationType, List<double>> complementaryIonConversionDictionary = new Dictionary<DissociationType, List<double>>
        {
            { DissociationType.LowCID, new List<double>(){ Constants.ProtonMass } },
            { DissociationType.HCD, new List<double>(){ Constants.ProtonMass } },
            { DissociationType.ETD,new List<double>() {2 * Constants.ProtonMass } }, //presence of zplusone (zdot) makes this two instead of one
            { DissociationType.CID,new List<double>() {Constants.ProtonMass } },
            { DissociationType.EThcD,new List<double>() {Constants.ProtonMass, 2 * Constants.ProtonMass } },
            //TODO: refactor such that complementary ions are generated specifically for their complementary pair.
            //TODO: create a method to auto-determine the conversion
        };

        public readonly CommonParameters CommonParameters;
        protected readonly List<(string FileName, CommonParameters Parameters)> FileSpecificParameters;
        protected readonly List<string> NestedIds;

        // ISD weighted-PSM-scorer experiment (doc 18): when set, the classic Morpheus score credits each
        // matched ion by its pre-computed quality/association weight (encoded as the peak intensity in the
        // pseudo-MS2 MGF) instead of a flat +1. Off by default; flip per process via the env var so the same
        // binary can A/B weighted-vs-count without a rebuild. The weight is dominated by chromatographic
        // profile coherence (fragment XIC vs precursor 15 V XIC); a true fragment ~0.85, a chimeric one ~0.
        public static readonly bool UseWeightedFragmentScore =
            Environment.GetEnvironmentVariable("MM_WEIGHTED_FRAGMENT_SCORE") == "1";

        // ISD ladder-bonus scorer (Route 2, docs 20-21): add a SEQUENCE-COUPLED term to the count score —
        // the longest run of CONSECUTIVE matched backbone ions (a residue ladder a decoy cannot fake; at
        // q<=0.01 a run >=4 is 0% of decoys vs 39% of targets). Applied inside the per-scan competition so a
        // correct long-ladder sequence can WIN scans it loses on raw count. Off by default; MM_LADDER_BONUS=1
        // turns it on, MM_LADDER_LAMBDA sets the weight (default 3, the safe re-rank optimum).
        public static readonly bool UseLadderBonus =
            Environment.GetEnvironmentVariable("MM_LADDER_BONUS") == "1";
        private static readonly double LadderLambda =
            double.TryParse(Environment.GetEnvironmentVariable("MM_LADDER_LAMBDA"), out var _lam) ? _lam : 3.0;

        // backbone series whose consecutive FragmentNumbers form a residue ladder (matches the doc-20/21 b/y/c/z analysis)
        private static readonly HashSet<ProductType> LadderSeries = new HashSet<ProductType>
        {
            ProductType.a, ProductType.b, ProductType.c, ProductType.x, ProductType.y, ProductType.zDot, ProductType.zPlusOne
        };

        /// <summary>Longest run of consecutive matched ions within a single backbone series (the residue tag).</summary>
        public static int LongestResidueLadder(List<MatchedFragmentIon> matchedFragmentIons)
        {
            // group matched fragment numbers by series; only backbone types contribute to a ladder
            var bySeries = new Dictionary<ProductType, SortedSet<int>>();
            foreach (var f in matchedFragmentIons)
            {
                ProductType pt = f.NeutralTheoreticalProduct.ProductType;
                if (!LadderSeries.Contains(pt)) continue;
                if (!bySeries.TryGetValue(pt, out var set)) { set = new SortedSet<int>(); bySeries[pt] = set; }
                set.Add(f.NeutralTheoreticalProduct.FragmentNumber);
            }
            int best = 0;
            foreach (var set in bySeries.Values)
            {
                int run = 0, prev = int.MinValue;
                foreach (int n in set) // SortedSet enumerates ascending
                {
                    run = (n == prev + 1) ? run + 1 : 1;
                    prev = n;
                    if (run > best) best = run;
                }
            }
            return best;
        }

        protected MetaMorpheusEngine(CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds)
        {
            CommonParameters = commonParameters;
            FileSpecificParameters = fileSpecificParameters;
            NestedIds = nestedIds;
        }

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        public static double CalculatePeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons, bool fragmentsCanHaveDifferentCharges = false)
        {
            if(fragmentsCanHaveDifferentCharges)
            {
                return CalculateAllChargesPeptideScore(thisScan, matchedFragmentIons);
            }

            double score = 0;

            if (thisScan.MassSpectrum.XcorrProcessed)
            {
                // XCorr
                foreach (var fragment in matchedFragmentIons)
                {
                    switch (fragment.NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.aDegree:
                        case ProductType.aStar:
                        case ProductType.bWaterLoss:
                        case ProductType.bAmmoniaLoss:
                        case ProductType.yWaterLoss:
                        case ProductType.yAmmoniaLoss:
                            score += 0.01 * fragment.Intensity;
                            break;
                        case ProductType.D: //count nothing for diagnostic ions.
                            break;
                        default:
                            score += 1 * fragment.Intensity;
                            break;
                    }
                }
            }
            else
            {
                // Morpheus score
                for (int i = 0; i < matchedFragmentIons.Count; i++)
                {
                    switch (matchedFragmentIons[i].NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.D:
                            break;
                        default:
                            // weighted: the peak intensity already encodes the per-fragment quality/association
                            // weight, so credit it directly instead of a flat +1 (suppresses chimera-inflated
                            // false matches that a count rewards). Otherwise the standard count + intensity tiebreak.
                            score += UseWeightedFragmentScore
                                ? matchedFragmentIons[i].Intensity
                                : 1 + matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent;
                            break;
                    }

                }
            }

            // Route-2 sequence-coupled term: reward a correct candidate's consecutive residue ladder, which a
            // decoy structurally cannot reproduce (docs 20-21). Additive and non-negative, so it never demotes
            // an existing match below ScoreCutoff; it lifts long-ladder candidates within the per-scan competition.
            if (UseLadderBonus)
            {
                score += LadderLambda * LongestResidueLadder(matchedFragmentIons);
            }

            return score;
        }

        //Used only when user wants to generate spectral library.
        //Normal search only looks for one match ion for one fragment, and if it accepts it then it doesn't try to look for different charge states of that same fragment. 
        //The score will be the number of matched ions and plus some fraction calculated by intensity(matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent).
        //Like b1, b2, b3 will have score 3.xxx;But when generating library, we need look for match ions with all charges.So we will have b1,b2,b3, b1^2, b2^3. If using 
        //the normal scoring function, the score will be 5.xxxx, which is not proper. The score for b1 and b1^2 should also be 1 plus some some fraction calculated by intensity, 
        //because they are matching the same fragment ion just with different charges. So b1, b2, b3, b1^2, b2^3 should be also 3.xxx(but a little higher than b1, b2, b3 as 
        //the fraction part) rather than 5.xxx.
        private static double CalculateAllChargesPeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons)
        {
            double score = 0;

            // Morpheus score
            List<String> ions = new List<String>();
            for (int i = 0; i < matchedFragmentIons.Count; i++)
            {
                String ion = $"{ matchedFragmentIons[i].NeutralTheoreticalProduct.ProductType.ToString()}{  matchedFragmentIons[i].NeutralTheoreticalProduct.FragmentNumber}";
                if (ions.Contains(ion))
                {
                    score += matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent;
                }
                else
                {
                    score += 1 + matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent;
                    ions.Add(ion);
                }
            }

            return score;

        }

        public static List<MatchedFragmentIon> MatchFragmentIons(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters, bool matchAllCharges = false, bool includeExperimentalEnvelope = false, bool isLowRes = false)
        {
            // if this is a child scan and it's an ion trap 2D scan, we want to use the wider tolerance for matching
            var productMassTolerance = isLowRes? commonParameters.ProductMassTolerance_LowRes : commonParameters.ProductMassTolerance;
            if (matchAllCharges)
            {
                return MatchFragmentIonsOfAllCharges(scan, theoreticalProducts, commonParameters, isLowRes);
            }

            var matchedFragmentIons = new List<MatchedFragmentIon>();

            if (scan.TheScan.MassSpectrum.XcorrProcessed && scan.TheScan.MassSpectrum.XArray.Length != 0)
            {

                for (int i = 0; i < theoreticalProducts.Count; i++)
                {
                    var product = theoreticalProducts[i];
                    // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                    if (double.IsNaN(product.NeutralMass))
                    {
                        continue;
                    }

                    // Magic number represents mzbinning space. 
                    double theoreticalFragmentMz = Math.Round(product.NeutralMass.ToMz(1) / 1.0005079, 0) * 1.0005079;
                    var closestMzIndex = scan.TheScan.MassSpectrum.GetClosestPeakIndex(theoreticalFragmentMz);


                    if (productMassTolerance.Within(scan.TheScan.MassSpectrum.XArray[closestMzIndex], theoreticalFragmentMz))
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, theoreticalFragmentMz, scan.TheScan.MassSpectrum.YArray[closestMzIndex], 1));
                    }
                }

                return matchedFragmentIons;
            }

            // if the spectrum has no peaks
            if (scan.ExperimentalFragments != null && !scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum
            for (int i = 0; i < theoreticalProducts.Count; i++)
            {
                var product = theoreticalProducts[i];
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // get the closest peak in the spectrum to the theoretical peak
                var closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(product.NeutralMass);

                // is the mass error acceptable?
                if (closestExperimentalMass != null
                    && productMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, product.NeutralMass)
                    && Math.Abs(closestExperimentalMass.Charge) <= Math.Abs(scan.PrecursorCharge))
                {
                    if (includeExperimentalEnvelope)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIonWithEnvelope(product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                            closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge)
                        {
                            Envelope = closestExperimentalMass
                        });
                    }
                    else
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                            closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge));
                    }
                }
            }
            if (commonParameters.AddCompIons)
            {
                foreach (double massShift in complementaryIonConversionDictionary[commonParameters.DissociationType])
                {
                    double protonMassShift = massShift.ToMass(1);

                    for (int i = 0; i < theoreticalProducts.Count; i++)
                    {
                        var product = theoreticalProducts[i];
                        // unknown fragment mass or diagnostic ion or precursor; skip those
                        if (double.IsNaN(product.NeutralMass) || product.ProductType == ProductType.D || product.ProductType == ProductType.M)
                        {
                            continue;
                        }

                        double compIonMass = scan.PrecursorMass + protonMassShift - product.NeutralMass;

                        // get the closest peak in the spectrum to the theoretical peak
                        IsotopicEnvelope closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(compIonMass);

                        // is the mass error acceptable?
                        if (closestExperimentalMass != null && productMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, compIonMass) && closestExperimentalMass.Charge <= scan.PrecursorCharge)
                        {
                            //found the peak, but we don't want to save that m/z because it's the complementary of the observed ion that we "added". Need to create a fake ion instead.
                            double mz = (scan.PrecursorMass + protonMassShift - closestExperimentalMass.MonoisotopicMass).ToMz(closestExperimentalMass.Charge);

                            if (includeExperimentalEnvelope)
                            {
                                matchedFragmentIons.Add(new MatchedFragmentIonWithEnvelope(product, mz, closestExperimentalMass.TotalIntensity, closestExperimentalMass.Charge)
                                {
                                    Envelope = closestExperimentalMass
                                });
                            }
                            else
                            {
                                matchedFragmentIons.Add(new MatchedFragmentIon(product, mz, closestExperimentalMass.TotalIntensity, closestExperimentalMass.Charge));
                            }
                        }
                    }
                }
            }

            return matchedFragmentIons;
        }
        
        //Used only when user wants to generate spectral library.
        //Normal search only looks for one match ion for one fragment, and if it accepts it then it doesn't try to look for different charge states of that same fragment. 
        //But for library generation, we need find all the matched peaks with all the different charges.
        private static List<MatchedFragmentIon> MatchFragmentIonsOfAllCharges(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters, bool isLowRes = false)
        {
            var productMassTolerance = isLowRes ? commonParameters.ProductMassTolerance_LowRes : commonParameters.ProductMassTolerance;
            var matchedFragmentIons = new List<MatchedFragmentIon>();
            var ions = new List<string>();

            // if the spectrum has no peaks
            if (scan.ExperimentalFragments != null && !scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum
            foreach (Product product in theoreticalProducts)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                //get the range we can accept 
                var minMass = productMassTolerance.GetMinimumValue(product.NeutralMass);
                var maxMass = productMassTolerance.GetMaximumValue(product.NeutralMass);
                var closestExperimentalMassList = scan.GetClosestExperimentalIsotopicEnvelopeList(minMass, maxMass);
                if (closestExperimentalMassList != null)
                {
                    foreach (var x in closestExperimentalMassList)
                    {
                        String ion = $"{product.ProductType.ToString()}{ product.FragmentNumber}^{x.Charge}-{product.NeutralLoss}";
                        if (x != null 
                            && !ions.Contains(ion) 
                            && productMassTolerance.Within(x.MonoisotopicMass, product.NeutralMass) 
                            && Math.Abs(x.Charge) <= Math.Abs(scan.PrecursorCharge))//TODO apply this filter before picking the envelope
                        {
                            Product temProduct = product;
                            matchedFragmentIons.Add(new MatchedFragmentIon(temProduct, x.MonoisotopicMass.ToMz(x.Charge),
                                x.Peaks.First().intensity, x.Charge));

                            ions.Add(ion);
                        }
                    }
                }
            }

            return matchedFragmentIons;
        }
        protected abstract MetaMorpheusEngineResults RunSpecific();

        public MetaMorpheusEngineResults Run()
        {
            DetermineAnalyteType(CommonParameters);
            StartingSingleEngine();
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            this.CommonParameters.SetCustomProductTypes();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            FinishedSingleEngine(myResults);
            return myResults;
        }

        public Task<MetaMorpheusEngineResults> RunAsync() => Task.Run(Run);

        /// <summary>
        /// Determines and sets the analyte type based on CommonParameters digestion settings.
        /// This method is called automatically by MetaMorpheusEngine.Run() to handle:
        /// - RNA mode (RnaDigestionParams → Oligo)
        /// - Top-down mode (protease == "top-down" → Proteoform)  
        /// - Bottom-up/default mode (→ Peptide)
        /// 
        /// IMPORTANT: This recalculates the analyte type at runtime and may differ from the GUI mode
        /// set by GuiGlobalParamsViewModel.IsRnaMode. This is intentional to support:
        /// - File-specific parameters with different modes
        /// - Mixed mode workflows
        /// 
        /// For GUI initialization, rely on GuiGlobalParamsViewModel.IsRnaMode which sets 
        /// GlobalVariables.AnalyteType. Do NOT call this method during GUI task window initialization.
        /// </summary>
        /// <param name="commonParameters"></param>
        public static void DetermineAnalyteType(CommonParameters commonParameters)
        {
            // Comment made while DetermineAnalyteType happened at the task layer
            // TODO: note that this will not function well if the user is using file-specific settings, but it's assumed
            // that bottom-up and top-down data is not being searched in the same task. 

            // Update: Now that it is in the engine layer, analyte type specific operations will be okay at the engine layer, meaning searching top-down and bottom-up with file specific params will execute the proper control flow. However, a problem still exists in PostSearchAnalysis where that analyte type will be set to whatever the main parameters are. 

            if (commonParameters == null || commonParameters.DigestionParams == null)
                return;

            GlobalVariables.AnalyteType = commonParameters.DetermineAnalyteType();
        }

        #region Event Helpers

        public string GetId()
        {
            return string.Join(",", NestedIds);
        }

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, NestedIds));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, NestedIds));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        private void StartingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void FinishedSingleEngine(MetaMorpheusEngineResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }

        #endregion
    }
}
