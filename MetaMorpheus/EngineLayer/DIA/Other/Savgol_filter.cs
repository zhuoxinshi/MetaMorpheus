using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace EngineLayer.DIA
{
    public class Savgol_filter
    {
        public static double[] Apply(double[] x, double[] y, int windowLength, int polyOrder, int deriv = 0)
        {
            if (windowLength % 2 == 0 || windowLength <= 0)
                throw new ArgumentException("window_length must be a positive odd integer.");
            if (polyOrder >= windowLength)
                throw new ArgumentException("polyorder must be less than window_length.");
            if (x.Length != y.Length)
                throw new ArgumentException("Input arrays x and y must have the same length.");
            int halfWindow = (windowLength - 1) / 2;
            double delta = (x.Max() - x.Min()) / (x.Length - 1);
            var coeffs = ComputeCoefficients(windowLength, polyOrder, deriv, delta);
            double[] extendedY = PadSignal(y, halfWindow);
            double[] result = new double[y.Length];
            for (int i = 0; i < y.Length; i++)
            {
                result[i] = coeffs.Zip(extendedY.Skip(i).Take(windowLength), (c, s) => c * s).Sum();
            }
            return result;
        }
        private static double[] ComputeCoefficients(int windowLength, int polyOrder, int deriv, double delta)
        {
            int halfWindow = (windowLength - 1) / 2;
            var A = Matrix<double>.Build.Dense(windowLength, polyOrder + 1, (i, j) => Math.Pow(i - halfWindow, j));
            var ATA = A.Transpose() * A;
            var ATAInv = ATA.Inverse();
            var ATAAT = ATAInv * A.Transpose();
            double[] coeffs = new double[windowLength];
            for (int i = 0; i < windowLength; i++)
            {
                coeffs[i] = ATAAT.Row(deriv)[i] * Factorial(deriv) / Math.Pow(delta, deriv);
            }
            return coeffs;
        }
        private static double[] PadSignal(double[] y, int halfWindow)
        {
            double[] padded = new double[y.Length + 2 * halfWindow];
            Array.Copy(y, 0, padded, halfWindow, y.Length);
            for (int i = 0; i < halfWindow; i++)
            {
                padded[i] = y[0];
                padded[padded.Length - i - 1] = y[y.Length - 1];
            }
            return padded;
        }
        private static int Factorial(int n)
        {
            if (n == 0) return 1;
            return n * Factorial(n - 1);
        }
    }
}

