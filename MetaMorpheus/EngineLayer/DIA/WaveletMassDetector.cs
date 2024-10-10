using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using static Plotly.NET.StyleParam;

namespace EngineLayer.DIA
{
    public class WaveletMassDetector
    {
        private int NPOINTS;
        private int WAVELET_ESL = -5;
        private int WAVELET_ESR = 5;
        public List<double>[] PeakRidges { get; set; }
        public float[] DataPoint;
        double waveletWindow = 0.3;
        private double[] MEXHAT;
        public int d;
        public double NPOINTS_half;
        public double MaxCurveRTRange = 2;
        public int NoPeakPerMin = 150;
        public double SymThreshold = 0.3;
        public List<(float rt, float intensity)>[] PeakRidge;

        public WaveletMassDetector(float[] DataPoint, int NoPoints)
        {
            this.DataPoint = DataPoint;
            this.NPOINTS = NoPoints;
            double wstep = ((WAVELET_ESR - WAVELET_ESL) / NPOINTS);
            MEXHAT = new double[(int)NPOINTS];

            double waveletIndex = WAVELET_ESL;
            for (int j = 0; j < NPOINTS; j++)
            {
                // Pre calculate the values of the wavelet
                MEXHAT[j] = cwtMEXHATreal(waveletIndex, waveletWindow, 0.0);
                waveletIndex += wstep;
            }

            NPOINTS_half = NPOINTS / 2;
            d = (int)NPOINTS / (WAVELET_ESR - WAVELET_ESL);
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

            //waveletCWT = new ArrayList[15];
            PeakRidge = new List<(float rt, float intensity)>[maxscale];
            //XYData maxint = new XYData(0f, 0f);
            for (int scaleLevel = 0; scaleLevel < maxscale; scaleLevel++)
            {
                //            ArrayList<XYData> wavelet = performCWT(scaleLevel * 2 + 5);
                float[] wavelet = performCWT(scaleLevel * 2 + 5);
                PeakRidge[scaleLevel] = new List<(float rt, float intensity)>();
                //waveletCWT[scaleLevel] = wavelet;
                //            XYData lastpt = wavelet.get(0);
                int lastptidx = 0;
                //            XYData localmax = null;
                int localmaxidx = -1;
                //            XYData startpt = wavelet.get(0);
                int startptidx = 0;

                var increasing = false;
                var decreasing = false;
                //            XYData localmaxint = null;
                int localmaxintidx = -1;

                //            for (int cwtidx = 1; cwtidx < wavelet.size(); cwtidx++) {
                for (int cwtidx = 1; cwtidx < wavelet.Length / 2; cwtidx++)
                {
                    //                XYData CurrentPoint = wavelet.get(cwtidx);
                    float CurrentPointY = wavelet[2 * cwtidx + 1],
                            lastptY = wavelet[2 * lastptidx + 1],
                            startptY = wavelet[2 * startptidx + 1],
                            localmaxY = localmaxidx == -1 ? float.NaN : wavelet[2 * localmaxidx + 1];
                    //                if (CurrentPoint.getY() > lastpt.getY()) {//the peak is increasing
                    if (CurrentPointY > lastptY)
                    {//the peak is increasing
                        if (decreasing)
                        {//first increasing point, last point was a possible local minimum
                         //check if the peak was symetric
                         //                        if (localmax != null && (lastpt.getY() <= startpt.getY() || Math.abs(lastpt.getY() - startpt.getY()) / localmax.getY() < parameter.SymThreshold)) {
                            if (localmaxidx != -1 && (lastptY <= startptY || Math.Abs(lastptY - startptY) / localmaxY < SymThreshold))
                            {
                                //                            PeakRidge[scaleLevel].add(localmax);
                                PeakRidge[scaleLevel].Add(new (wavelet[2 * localmaxidx], wavelet[2 * localmaxidx + 1]));
                                //                            localmax = CurrentPoint;
                                localmaxidx = cwtidx;
                                //                            startpt = lastpt;
                                startptidx = lastptidx;
                            }
                        }
                        increasing = true;
                        decreasing = false;
                        //                } else if (CurrentPoint.getY() < lastpt.getY()) {//peak decreasing
                    }
                    else if (CurrentPointY < lastptY)
                    {//peak decreasing
                        if (increasing)
                        {//first point decreasing, last point was a possible local maximum
                         //                        if (localmax == null || localmax.getY() < lastpt.getY()) {
                            if (localmaxidx == -1 || localmaxY < lastptY)
                            {
                                //                            localmax = lastpt;
                                localmaxidx = lastptidx;
                            }
                        }
                        decreasing = true;
                        increasing = false;
                    }
                    //                lastpt = CurrentPoint;
                    lastptidx = cwtidx;
                    float localmaxintY = localmaxintidx == -1 ? float.NaN : wavelet[2 * localmaxintidx + 1];
                    //                if (localmaxint == null || CurrentPoint.getY() > localmaxint.getY()) {
                    if (localmaxintidx == -1 || CurrentPointY > localmaxintY)
                    {
                        //                    localmaxint = CurrentPoint;
                        localmaxintidx = cwtidx;
                    }
                    //                if (cwtidx == wavelet.size() - 1 && decreasing) {
                    if (cwtidx == wavelet.Length / 2 - 1 && decreasing)
                    {
                        //                    if (localmax != null && (CurrentPoint.getY() <= startpt.getY() || Math.abs(CurrentPoint.getY() - startpt.getY()) / localmax.getY() < parameter.SymThreshold)) {
                        //                    final float startptY=wavelet[2*startptidx+1];
                        if (localmaxidx != -1 && (CurrentPointY <= startptY || Math.Abs(CurrentPointY - startptY) / localmaxY < SymThreshold))
                        {
                            var localmax = (wavelet[2 * localmaxidx], wavelet[2 * localmaxidx + 1]);
                            PeakRidge[scaleLevel].Add(localmax);
                        }
                    }
                }
            }
        }

        private float[] performCWT(int scaleLevel)
        {
            //        int length = DataPoint.size();
            int length = DataPoint.Length / 2;
            //        ArrayList<XYData> cwtDataPoints = new ArrayList<XYData>();
            float[] cwtDataPoints = new float[length * 2];

            int a_esl = scaleLevel * WAVELET_ESL;
            int a_esr = scaleLevel * WAVELET_ESR;
            double sqrtScaleLevel = Math.Sqrt(scaleLevel);
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
                }

                /*
                 * Perform convolution
                 */
                float intensity = 0f;
                for (int i = t1; i <= t2; i++)
                {
                    int ind = (int)(NPOINTS_half) + (d * (i - dx) / scaleLevel);
                    if (ind < 0)
                    {
                        ind = 0;
                    }
                    if (ind >= NPOINTS)
                    {
                        ind = (int)NPOINTS - 1;
                    }
                    //                if(i<0 || ind<0){
                    //                    System.out.print("");
                    //                }
                    //                intensity += DataPoint.get(i).getY() * MEXHAT[ind];
                    intensity += DataPoint[2 * i + 1] * (float)MEXHAT[ind];
                }
                intensity /= (float)sqrtScaleLevel;
                // Eliminate the negative part of the wavelet map
                if (intensity < 0)
                {
                    intensity = 0;
                }
                //            cwtDataPoints.add(new XYData(DataPoint.get(dx).getX(), intensity));
                cwtDataPoints[2 * dx] = DataPoint[2 * dx];
                cwtDataPoints[2 * dx + 1] = intensity;
            }
            return cwtDataPoints;
        }

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
