using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using ThermoFisher.CommonCore.BackgroundSubtraction;
using static Plotly.NET.StyleParam;

namespace EngineLayer.DIA
{
    public class WaveletMassDetector
    {
        private double NPOINTS; //number of points in the data
        private int WAVELET_ESL = -5;
        private int WAVELET_ESR = 5;
        public float[] DataPoint;
        double WaveletWindow = 0.3;
        private double[] MEXHAT;
        public int d;
        public double NPOINTS_half;
        public double MaxCurveRTRange = 2;
        public int NoPeakPerMin = 150;
        public double SymThreshold = 0.3;
        public List<(float rt, float intensity, int index)>[] PeakRidge;

        public WaveletMassDetector(float[] DataPoint, double NoPoints)
        {
            this.DataPoint = DataPoint;
            this.NPOINTS = NoPoints;
            double wstep = ((WAVELET_ESR - WAVELET_ESL) / (NPOINTS - 1)); //range of the wavelet[-5, 5]/number of points; I changed NPOINTS to NPOINTS - 1
            MEXHAT = new double[(int)NPOINTS];

            double waveletIndex = WAVELET_ESL;
            for (int j = 0; j < NPOINTS; j++)
            {
                // Pre calculate the values of the wavelet
                MEXHAT[j] = cwtMEXHATreal(waveletIndex, WaveletWindow, 0.0); //construct the theoretical mother MEXHAT wavelet distribution: calculate the relative intensity at each point
                waveletIndex += wstep;
            }

            NPOINTS_half = NPOINTS / 2;
            d = (int)NPOINTS / (WAVELET_ESR - WAVELET_ESL); //density of the points in the wavelet
        }

        public void Run()
        {

            //"Intensities less than this value are interpreted as noise",                
            //"Scale level",
            //"Number of wavelet'scale (coeficients) to use in m/z peak detection"
            //"Wavelet window size (%)",
            //"Size in % of wavelet window to apply in m/z peak detection");        
            //        int maxscale = (int) (Math.max(Math.min((DataPoint.get(DataPoint.size() - 1).getX() - DataPoint.get(0).getX()), parameter.MaxCurveRTRange), 0.5f) * parameter.NoPeakPerMin / (WAVELET_ESR + WAVELET_ESR));
            int maxscale = (int)(Math.Max(Math.Min((DataPoint[2 * (DataPoint.Length / 2 - 1)] - DataPoint[0]), MaxCurveRTRange), 0.5f) * NoPeakPerMin / (WAVELET_ESR + WAVELET_ESR));

            PeakRidge = new List<(float rt, float intensity, int index)>[maxscale];
            for (int scaleLevel = 0; scaleLevel < maxscale; scaleLevel++)
            {
                float[] wavelet = performCWT(scaleLevel * 2 + 5); //the cwt coefficient calculated at each point
                PeakRidge[scaleLevel] = new List<(float rt, float intensity, int index)>();
                int lastptidx = 0;
                int localmaxidx = -1;
                int startptidx = 0;

                var increasing = false;
                var decreasing = false;
                int localmaxintidx = -1;

                for (int cwtidx = 1; cwtidx < wavelet.Length / 2; cwtidx++)
                {
                    //debug
                    if (Math.Abs(wavelet[2 * cwtidx] - 68.1) < 0.1)
                    {
                        bool stop = true;
                        var t = wavelet[2 * cwtidx];
                    }
                    float CurrentPointY = wavelet[2 * cwtidx + 1],
                            lastptY = wavelet[2 * lastptidx + 1],
                            startptY = wavelet[2 * startptidx + 1],
                            localmaxY = localmaxidx == -1 ? float.NaN : wavelet[2 * localmaxidx + 1];
                    if (CurrentPointY > lastptY)
                    {
                        //the peak is increasing
                        if (decreasing)
                        {//first increasing point, last point was a possible local minimum
                         //check if the peak was symetric
                            if (localmaxidx != -1 && (lastptY <= startptY || Math.Abs(lastptY - startptY) / localmaxY < SymThreshold))
                            {
                                PeakRidge[scaleLevel].Add(new (wavelet[2 * localmaxidx], wavelet[2 * localmaxidx + 1], localmaxidx));
                                localmaxidx = cwtidx;
                                startptidx = lastptidx;
                            }
                        }
                        increasing = true;
                        decreasing = false;
                    }
                    else if (CurrentPointY < lastptY)
                    {
                        //peak decreasing
                        if (increasing)
                        {
                            //first point decreasing, last point was a possible local maximum
                            if (localmaxidx == -1 || localmaxY < lastptY)
                            {
                                localmaxidx = lastptidx;
                            }
                        }
                        decreasing = true;
                        increasing = false;
                    }
                    lastptidx = cwtidx;
                    float localmaxintY = localmaxintidx == -1 ? float.NaN : wavelet[2 * localmaxintidx + 1];
                    if (localmaxintidx == -1 || CurrentPointY > localmaxintY)
                    {
                        localmaxintidx = cwtidx;
                    }
                    if (cwtidx == wavelet.Length / 2 - 1 && decreasing)
                    {
                        if (localmaxidx != -1 && (CurrentPointY <= startptY || Math.Abs(CurrentPointY - startptY) / localmaxY < SymThreshold))
                        {
                            var localmax = (wavelet[2 * localmaxidx], wavelet[2 * localmaxidx + 1], localmaxidx);
                            PeakRidge[scaleLevel].Add(localmax);
                        }
                    }
                }
            }
        }

        /**
        * Perform the CWT over raw data points in the selected scale level
        */
        private float[] performCWT(int scaleLevel)
        {
            int length = DataPoint.Length / 2;
            float[] cwtDataPoints = new float[length * 2];

            int a_esl = scaleLevel * WAVELET_ESL; //stretch or shrink the wavelet range
            int a_esr = scaleLevel * WAVELET_ESR;
            double sqrtScaleLevel = Math.Sqrt(scaleLevel);

            //This for loop is the process of sliding the wavelet along the raw data points, each loop represents the snapshot of the wavelet at a certain position
            //dx is the index of the raw data points, the position where the theoretical wavelet is currently centered on 
            //The second loop inside this loop sums the intensity at each overlapped point to calculate the coefficient C(a,b) at this snapshot
            for (int dx = 0; dx < length; dx++) 
            {
                /*
                 * Compute wavelet boundaries
                 */
                int t1 = a_esl + dx;
                if (t1 < 0)
                {
                    t1 = 0;
                }
                int t2 = a_esr + dx;
                if (t2 >= length)
                {
                    t2 = (length - 1);
                } // this gives the range of raw data points that overlap with the theoretical wavelet wavelet

                /*
                 * Perform convolution
                 */
                float intensity = 0f;
                for (int i = t1; i <= t2; i++)
                {
                    int ind = (int)(NPOINTS_half) + (d * (i - dx) / scaleLevel);
                    //(d * (i - dx) / scaleLevel) is the relative position of the raw data point to the center of the wavelet in index space
                    //Because the wavelet is symmetric, MEXHAT[] is in index space and cannot take negative values
                    //(t-b)/a is transformed to index space by adding half of number of points
                    if (ind < 0)
                    {
                        ind = 0;
                    }
                    if (ind >= NPOINTS)
                    {
                        ind = (int)NPOINTS - 1;
                    }
                    intensity += DataPoint[2 * i + 1] * (float)MEXHAT[ind];
                }
                intensity /= (float)sqrtScaleLevel;
                // Eliminate the negative part of the wavelet map
                if (intensity < 0)
                {
                    intensity = 0;
                }
                cwtDataPoints[2 * dx] = DataPoint[2 * dx];
                cwtDataPoints[2 * dx + 1] = intensity;
            }
            return cwtDataPoints;
        }

        //This function calculates the value phi(t) at each point t; x is t (ESL + step, x axis of the wavelet, not in real time unit), window is sigma,
        //b is the center of the wavelet which is set to default 0
        private double cwtMEXHATreal(double x, double window, double b)
        {
            /*
             * c = 2 / ( sqrt(3) * pi^(1/4) )
             */
            double c = 0.8673250705840776;
            double TINY = 1E-200;
            double x2;

            if (window == 0.0)
            {
                window = TINY;
            }
            //x-b=t
            //window=delta
            x = (x - b) / window;
            x2 = x * x;
            return c * (1.0 - x2) * Math.Exp(-x2 / 2);
        }
    }
}
