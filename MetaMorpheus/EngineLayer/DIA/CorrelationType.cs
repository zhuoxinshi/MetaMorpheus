using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public enum CorrelationType
    {
        NoSpline,
        LinearSpline_RT,
        LinearSpline_scanCycle,
        CubicSpline_RT,
        CubicSpline_scanCycle,
        CubicSpline_scanCycle_preCalc,
        Normalized
    }
}
