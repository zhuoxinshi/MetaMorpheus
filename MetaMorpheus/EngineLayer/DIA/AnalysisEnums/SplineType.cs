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
        BSpline,
        UmpireBSpline,
        Ms1SpaceBSpline,
        ScanCycleCubicSpline,
        Normalized,
        Ms1SpaceCubicSpline,
        SavgolSmoothed,
        SavgolSmoothedCubicSpline,
        Ms1SpaceSavgolSmoothedCubicSpline,
        CubicSplineSavgolSmoothed,
        ScanCycleCubicSplineSavgolSmoothed,
        Ms1SpaceCubicSplineSavgolSmoothed,
        GaussianFit,
        SimpleGaussian,
        SimpleGaussianSpline,
        NormalizedLinearSpline,
        ExtendedCubicSpline,
    }
}
