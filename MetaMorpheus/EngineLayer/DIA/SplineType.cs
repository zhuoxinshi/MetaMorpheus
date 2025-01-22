using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public enum SplineType
    {
        NoSpline,
        LinearSpline,
        CubicSpline,
        ScanCycleCubicSpline,
        Normalized,
        Ms1SpaceCubicSpline,
        SavgolSmoothedCubicSpline,
        CubicSplineSavgolSmoothed,
        ScanCycleCubicSplineSavgolSmoothed,
        Ms1SpaceCubicSplineSavgolSmoothed
    }
}
