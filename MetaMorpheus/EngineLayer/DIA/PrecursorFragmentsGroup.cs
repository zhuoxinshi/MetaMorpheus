using MassSpectrometry;
using Plotly.NET;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PrecursorFragmentsGroup
    {
        public PrecursorFragmentsGroup(PeakCurve precursorPeakCurve)
        {
            PrecursorPeakCurve = precursorPeakCurve;
            PFpairs = new List<PrecursorFragmentPair>();
        }

        public PeakCurve PrecursorPeakCurve { get; set; }
        public List<PrecursorFragmentPair> PFpairs { get; set; }
        public int Index => PrecursorPeakCurve.Index;
        public int NumHighCorrFragments { get; set; }

        //TODO: finish this
        public static PrecursorFragmentsGroup GroupPF(PeakCurve prePeakCurve, List<PeakCurve> fragPeakCurves)
        {
            return new PrecursorFragmentsGroup(prePeakCurve);
        }
        
        public void GetNumberOfHighCorrFragments(DIAparameters diaParam)
        {
            NumHighCorrFragments = PFpairs.Where(p => p.Correlation >= diaParam.HighCorrThreshold).Count();
        }
        
        public void Visualize()
        {
            var plots = new List<GenericChart>();
            var normalizedIntensity = PrecursorPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
            var precursorPlot = Chart2D.Chart.Line<double, double, string>(
                    x: PrecursorPeakCurve.Peaks.Select(p => p.RetentionTime),
                    y: normalizedIntensity).WithTraceInfo($"precursor_{Math.Round(PrecursorPeakCurve.AveragedMz, 3)}").WithMarkerStyle(Color: Color.fromString("red"));
            foreach(var pf in PFpairs)
            {
                var norm2 = pf.FragmentPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                        x: pf.FragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                        y: norm2).WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{pf.Correlation}_{pf.FragmentRank}")
                        .WithMarkerStyle(Color: Color.fromString("blue"));
                plots.Add(fragmentPlot2);
            }
            var combinedPlot = Chart.Combine(plots);
            combinedPlot.Show();
        }

        public static List<GenericChart> VisualizePFgroups(MsDataFile dataFile, List<PsmFromTsv> psms, CommonParameters commonParameters, DIAparameters diaParam)
        {
            var diaEngine2 = new DIAEngine2(dataFile, commonParameters, diaParam);
            diaEngine2.Ms1PeakIndexing();
            diaEngine2.ConstructMs2Group();
            diaEngine2.GetMs1PeakCurves();
            diaEngine2.GetMs2PeakCurves();

            var chartList = new List<GenericChart>();
            foreach (var psm in psms)
            {
                var precursorPeakCurve = diaEngine2.Ms1PeakCurves.Values.SelectMany(v => v).Where(pc => pc.Index == psm.Ms2ScanNumber).First();
                var ms2curves = diaEngine2.Ms2PeakCurves.Where(d => precursorPeakCurve.AveragedMz + 1 > d.Key.min && precursorPeakCurve.AveragedMz - 1 < d.Key.max)
                .First().Value;
                var pfgroup = DIAEngine2.GroupPrecursorFragments(precursorPeakCurve, ms2curves, diaParam);
                var plots = new List<GenericChart>();
                var normalizedIntensity = precursorPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                var precursorPlot = Chart2D.Chart.Line<double, double, string>(
                    x: precursorPeakCurve.Peaks.Select(p => p.RetentionTime),
                    y: normalizedIntensity).WithTraceInfo($"precursor_{Math.Round(precursorPeakCurve.AveragedMz, 3)}_{psm.DecoyContamTarget}").WithMarkerStyle(Color: Color.fromString("red"));
                plots.Add(precursorPlot);

                var matchedIonPeakCurves = new List<PeakCurve>();
                var matchedIonMzs = psm.MatchedIons.Select(i => i.Mz);
                foreach (double mz in matchedIonMzs)
                {
                    try
                    {
                        var fragmentPeakCurve = pfgroup.PFpairs.Where(p => Math.Abs(p.FragmentPeakCurve.AveragedMz - mz) < 0.005).First().FragmentPeakCurve;
                        matchedIonPeakCurves.Add(fragmentPeakCurve);
                    }
                    catch
                    {

                    }
                }
                foreach (var pf in pfgroup.PFpairs)
                {
                    if (!matchedIonPeakCurves.Contains(pf.FragmentPeakCurve))
                    {
                        var norm2 = pf.FragmentPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                        var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                                                       x: pf.FragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                                                       y: norm2).WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{pf.Correlation}_{pf.FragmentRank}")
                                                       .WithMarkerStyle(Color: Color.fromString("green"));
                        plots.Add(fragmentPlot2);
                    }
                    else
                    {
                        var norm2 = pf.FragmentPeakCurve.Peaks.Select(p => Math.Log10(p.Intensity));
                        var fragmentPlot2 = Chart2D.Chart.Line<double, double, string>(
                                                       x: pf.FragmentPeakCurve.Peaks.Select(p => p.RetentionTime),
                                                       y: norm2).WithTraceInfo($"fragment_{Math.Round(pf.FragmentPeakCurve.AveragedMz, 3)}_{pf.Correlation}_{pf.FragmentRank}")
                                                       .WithMarkerStyle(Color: Color.fromString("blue"));
                        plots.Add(fragmentPlot2);
                    }
                }
                var combinedPlot = Chart.Combine(plots).WithLayout(Layout.init<string>(Title: Title.init($"{psm.Ms2ScanNumber}_{psm.Score}_{psm.DecoyContamTarget}")));
                chartList.Add(combinedPlot);
            }
            return chartList;
        }
    }
}
