using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using static Plotly.NET.StyleParam;

namespace EngineLayer.DIA
{
    public class Bspline
    {
        private float[] bspline_T_ = null;

        public List<(float, float)> Run(List<(float, float)> data, int PtNum, int smoothDegree)
        {
            List<(float, float)> bsplineCollection = new List<(float, float)>();
            int p = smoothDegree;
            int n = data.Count() - 1;
            int m = data.Count() + p;
            bspline_T_ = new float[m + p];

            if (data.Count() <= p)
            {
                return data;
            }

            for (int i = 0; i <= n; i++)
            {
                bspline_T_[i] = 0;
                bspline_T_[m - i] = 1;
            }
            float intv = 1.0f / (m - 2 * p);
            for (int i = 1; i <= m - 1; i++)
            {
                bspline_T_[p + i] = bspline_T_[p + i - 1] + intv;
            }


            for (int i = 0; i <= PtNum; i++)
            {
                float t = (float)i / PtNum;
                var pt = getbspline(data, t, n, p);
                bsplineCollection.Add(pt);
            }
            if (bsplineCollection[bsplineCollection.Count() - 1].Item1 < data[data.Count() - 1].Item1)
            {
                bsplineCollection.Add(data[data.Count() - 1]);
            }
            if (bsplineCollection[0].Item1 > data[0].Item1)
            {
                bsplineCollection.Add(data[0]);
            }
            return bsplineCollection;
        }

        public (float, float) getbspline(List<(float, float)> data, float t, int n, int p)
        {
            float x = 0, y = 0;
            for (int i = 0; i <= n; i++)
            {
                float a = bspline_base(i, p, t);
                var pt = data[i];
                x += pt.Item1 * a;
                y += pt.Item2 * a;
            }
            return new(x, y);
        }

        float bspline_base(int i, int p, float t)
        {
            float n, c1, c2;
            float tn1 = 0;
            float tn2 = 0;
            if (p == 0)
            {
                if (bspline_T_[i] <= t && t < bspline_T_[i + 1] && bspline_T_[i] < bspline_T_[i + 1])
                {
                    n = 1;
                }
                else
                {
                    n = 0;
                }
            }
            else
            {
                if (bspline_T_[i + p] - bspline_T_[i] == 0)
                {
                    c1 = 0;
                }
                else
                {
                    tn1 = bspline_base(i, p - 1, t);
                    c1 = (t - bspline_T_[i]) / (bspline_T_[i + p] - bspline_T_[i]);
                }
                if (bspline_T_[i + p + 1] - bspline_T_[i + 1] == 0)
                {
                    c2 = 0;
                }
                else
                {
                    tn2 = bspline_base(i + 1, p - 1, t);
                    c2 = (bspline_T_[i + p + 1] - t) / (bspline_T_[i + p + 1] - bspline_T_[i + 1]);
                }
                n = c1 * tn1 + c2 * tn2;
            }
            return n;
        }
    }
}
