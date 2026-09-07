using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DIA.FeatureAlign
{
    /// <summary>
    /// Verified C# port of the numeric core of the ISD-vs-DDA manuscript analyses
    /// (Python scripts in the ISD_quant repo, paper/scripts):
    ///
    ///   Section 1 precision  -> pool_3file.py / pool_paired.py
    ///       size factor SF_j = mean(loading)/loading_j      (MS1 TIC-area normalization)
    ///       replicate CV     = std(ddof=1)/mean  on SF-normalized triplicates
    ///
    ///   Section 2 diff. expr -> s2_yeast_de.py
    ///       anchor size factors = median(perRunMedian)/perRunMedian
    ///       log2FC = mean(log2 group2) - mean(log2 group1)
    ///       Welch t (unequal variance, 3 vs 3) on the log2 values
    ///       Benjamini-Hochberg q-values
    ///
    /// Every method reproduces the exact numpy/scipy computation, so the two
    /// implementations agree to ~1e-9 on identical inputs (see FeatureQuantStatsTest).
    /// The Student-t p-value uses a self-contained regularized incomplete beta
    /// (Numerical Recipes betai), matching scipy.stats.ttest_ind(equal_var=False);
    /// there is no external numeric dependency.
    /// </summary>
    public static class FeatureQuantStats
    {
        // ---------------------------------------------------------------- moments

        /// <summary>Arithmetic mean.</summary>
        public static double Mean(IReadOnlyList<double> x)
        {
            if (x.Count == 0) return double.NaN;
            double s = 0.0;
            for (int i = 0; i < x.Count; i++) s += x[i];
            return s / x.Count;
        }

        /// <summary>Sample variance (ddof = 1), i.e. numpy var(ddof=1).</summary>
        public static double SampleVariance(IReadOnlyList<double> x)
        {
            int n = x.Count;
            if (n < 2) return double.NaN;
            double m = Mean(x), s = 0.0;
            for (int i = 0; i < n; i++) { double d = x[i] - m; s += d * d; }
            return s / (n - 1);
        }

        /// <summary>Sample standard deviation (ddof = 1), i.e. numpy std(ddof=1).</summary>
        public static double SampleStdDev(IReadOnlyList<double> x) => Math.Sqrt(SampleVariance(x));

        /// <summary>numpy-compatible median (average of the two middle values for even n).</summary>
        public static double Median(IReadOnlyList<double> x)
        {
            int n = x.Count;
            if (n == 0) return double.NaN;
            double[] s = x.ToArray();
            Array.Sort(s);
            return (n % 2 == 1) ? s[n / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]);
        }

        /// <summary>numpy-compatible percentile with the default 'linear' interpolation. q in [0, 100].</summary>
        public static double Percentile(IReadOnlyList<double> x, double q)
        {
            int n = x.Count;
            if (n == 0) return double.NaN;
            double[] s = x.ToArray();
            Array.Sort(s);
            if (n == 1) return s[0];
            double rank = (q / 100.0) * (n - 1);
            int lo = (int)Math.Floor(rank);
            if (lo >= n - 1) return s[n - 1];
            double frac = rank - lo;
            return s[lo] + frac * (s[lo + 1] - s[lo]);
        }

        /// <summary>Interquartile range (p75 - p25), numpy 'linear'.</summary>
        public static double InterquartileRange(IReadOnlyList<double> x) => Percentile(x, 75) - Percentile(x, 25);

        // ---------------------------------------------------- Section 1: precision

        /// <summary>
        /// Per-run size factors SF_j = mean(loading)/loading_j. With MS1 TIC areas as the loading,
        /// this is the TIC-area normalization used for the replicate-CV analysis; the common
        /// reference (mean) cancels out of a CV, so any constant multiple of the loading is equivalent.
        /// </summary>
        public static double[] SizeFactors(IReadOnlyList<double> loading)
        {
            double mean = Mean(loading);
            var sf = new double[loading.Count];
            for (int j = 0; j < loading.Count; j++) sf[j] = mean / loading[j];
            return sf;
        }

        /// <summary>Coefficient of variation = std(ddof=1)/mean. Zero for a constant vector.</summary>
        public static double CoefficientOfVariation(IReadOnlyList<double> values)
        {
            double m = Mean(values);
            if (m == 0.0) return double.NaN;
            double sd = SampleStdDev(values);
            return sd / m;
        }

        /// <summary>
        /// Replicate CV after applying per-run size factors: CV( values_j * sizeFactors_j ).
        /// This is exactly the per-feature precision metric in pool_3file.py / pool_paired.py.
        /// </summary>
        public static double ReplicateCv(IReadOnlyList<double> values, IReadOnlyList<double> sizeFactors)
        {
            if (values.Count != sizeFactors.Count)
                throw new ArgumentException("values and sizeFactors must have the same length.");
            var norm = new double[values.Count];
            for (int j = 0; j < values.Count; j++) norm[j] = values[j] * sizeFactors[j];
            return CoefficientOfVariation(norm);
        }

        // --------------------------------------------- Section 2: differential expr.

        /// <summary>
        /// Anchor "median-of-medians" size factors: for a features x runs matrix, take the median
        /// intensity within each run (over features), then SF_j = median(runMedians)/runMedian_j.
        /// Mirrors s2_yeast_de.size_factors (numpy median over axis 0, then median-of-those / m).
        /// </summary>
        public static double[] MedianRatioSizeFactors(IReadOnlyList<IReadOnlyList<double>> featuresByRun)
        {
            if (featuresByRun.Count == 0) throw new ArgumentException("empty matrix.");
            int nRun = featuresByRun[0].Count;
            var runMedians = new double[nRun];
            for (int j = 0; j < nRun; j++)
            {
                var col = new double[featuresByRun.Count];
                for (int i = 0; i < featuresByRun.Count; i++) col[i] = featuresByRun[i][j];
                runMedians[j] = Median(col);
            }
            double gm = Median(runMedians);
            var sf = new double[nRun];
            for (int j = 0; j < nRun; j++) sf[j] = gm / runMedians[j];
            return sf;
        }

        /// <summary>
        /// Welch t (unequal-variance two-sample). Returns (t, df) with
        /// t = (mean(a) - mean(b)) / sqrt(var(a)/n_a + var(b)/n_b) and the Welch-Satterthwaite df.
        /// Sign follows the argument order; the two-sided p-value is sign-independent.
        /// </summary>
        public static (double t, double df) WelchT(IReadOnlyList<double> a, IReadOnlyList<double> b)
        {
            int na = a.Count, nb = b.Count;
            double va = SampleVariance(a), vb = SampleVariance(b);
            double sa = va / na, sb = vb / nb;
            double t = (Mean(a) - Mean(b)) / Math.Sqrt(sa + sb);
            double df = (sa + sb) * (sa + sb) / (sa * sa / (na - 1) + sb * sb / (nb - 1));
            return (t, df);
        }

        /// <summary>
        /// Two-sided Welch t-test, matching scipy.stats.ttest_ind(a, b, equal_var=False).
        /// Returns (t, df, p).
        /// </summary>
        public static (double t, double df, double p) WelchTTest(IReadOnlyList<double> a, IReadOnlyList<double> b)
        {
            var (t, df) = WelchT(a, b);
            return (t, df, StudentTwoSidedP(t, df));
        }

        /// <summary>
        /// Per-feature differential expression exactly as s2_yeast_de.de does it:
        /// log2FC = mean(log2 group2) - mean(log2 group1), and the two-sided Welch p on the log2 values.
        /// Both groups are the SF-normalized (raw-intensity) replicate values.
        /// </summary>
        public static (double log2Fc, double p) LogFcWelch(IReadOnlyList<double> group1, IReadOnlyList<double> group2)
        {
            var l1 = group1.Select(v => Math.Log2(v)).ToArray();
            var l2 = group2.Select(v => Math.Log2(v)).ToArray();
            double fc = Mean(l2) - Mean(l1);
            double p = WelchTTest(l2, l1).p;
            return (fc, p);
        }

        /// <summary>
        /// Benjamini-Hochberg q-values, reproducing the running-min-from-largest-p loop in the scripts:
        /// process p ascending-rank i from largest to smallest, q_i = min(prev, p_i * m / i), clamped to 1.
        /// </summary>
        public static double[] BenjaminiHochberg(IReadOnlyList<double> pValues)
        {
            int m = pValues.Count;
            var q = new double[m];
            if (m == 0) return q;
            // ascending order of p (stable)
            var order = Enumerable.Range(0, m).OrderBy(i => pValues[i]).ToArray();
            double prev = 1.0;
            for (int rank = 0; rank < m; rank++)
            {
                int idx = order[m - 1 - rank];             // largest p first
                double factor = (double)m / (m - rank);    // m / ascending-rank
                prev = Math.Min(prev, pValues[idx] * factor);
                q[idx] = prev;
            }
            return q;
        }

        // --------------------------------------------------- special functions

        /// <summary>
        /// Two-sided Student-t p-value P(|T| >= |t|) for the given degrees of freedom, computed as
        /// I_x(df/2, 1/2) with x = df/(df + t^2). Equals scipy's 2 * t.sf(|t|, df).
        /// </summary>
        public static double StudentTwoSidedP(double t, double df)
        {
            if (double.IsNaN(t) || df <= 0.0) return double.NaN;
            double x = df / (df + t * t);
            return RegularizedIncompleteBeta(df / 2.0, 0.5, x);
        }

        /// <summary>Regularized incomplete beta I_x(a, b). Numerical Recipes "betai".</summary>
        internal static double RegularizedIncompleteBeta(double a, double b, double x)
        {
            if (x <= 0.0) return 0.0;
            if (x >= 1.0) return 1.0;
            double lbeta = LogGamma(a + b) - LogGamma(a) - LogGamma(b);
            double bt = Math.Exp(lbeta + a * Math.Log(x) + b * Math.Log(1.0 - x));
            if (x < (a + 1.0) / (a + b + 2.0))
                return bt * BetaContinuedFraction(a, b, x) / a;
            return 1.0 - bt * BetaContinuedFraction(b, a, 1.0 - x) / b;
        }

        /// <summary>Lentz continued fraction for the incomplete beta (Numerical Recipes "betacf").</summary>
        private static double BetaContinuedFraction(double a, double b, double x)
        {
            const int MAXIT = 300;
            const double EPS = 3.0e-15;
            const double FPMIN = 1.0e-300;
            double qab = a + b, qap = a + 1.0, qam = a - 1.0;
            double c = 1.0;
            double d = 1.0 - qab * x / qap;
            if (Math.Abs(d) < FPMIN) d = FPMIN;
            d = 1.0 / d;
            double h = d;
            for (int m = 1; m <= MAXIT; m++)
            {
                int m2 = 2 * m;
                double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1.0 + aa * d; if (Math.Abs(d) < FPMIN) d = FPMIN;
                c = 1.0 + aa / c; if (Math.Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d; h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = 1.0 + aa * d; if (Math.Abs(d) < FPMIN) d = FPMIN;
                c = 1.0 + aa / c; if (Math.Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                double del = d * c; h *= del;
                if (Math.Abs(del - 1.0) < EPS) break;
            }
            return h;
        }

        /// <summary>Log gamma (Lanczos approximation, Numerical Recipes "gammln").</summary>
        private static double LogGamma(double x)
        {
            double[] cof =
            {
                76.18009172947146, -86.50532032941677, 24.01409824083091,
                -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
            };
            double y = x;
            double tmp = x + 5.5;
            tmp -= (x + 0.5) * Math.Log(tmp);
            double ser = 1.000000000190015;
            for (int j = 0; j < 6; j++) { y += 1.0; ser += cof[j] / y; }
            return -tmp + Math.Log(2.5066282746310005 * ser / x);
        }
    }
}
