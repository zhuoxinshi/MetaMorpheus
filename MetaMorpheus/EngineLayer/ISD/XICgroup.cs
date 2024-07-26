using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Threading.Tasks.Dataflow;
using Easy.Common.Extensions;
using System.IO;

namespace EngineLayer.ISD
{
    public class XICgroup
    {
        public List<XIC> XIClist { get; set; }
        public XIC PrecursorXIC { get; set; }
        public XIC ReferenceXIC => XIClist.OrderByDescending(x => x.XICpeaks.Sum(p => p.Intensity)).First();

        public XICgroup(List<XIC> group)
        {
            XIClist = group;
        }

        public XICgroup()
        {
            XIClist = new List<XIC>();
        }

        public static XICgroup FindGroup(XIC target, List<XIC> XICs, double apexTolerance)
        {
            var group = new XICgroup();
            var xicMatchedByApex = XICs.Where(xic => Math.Abs(xic.ApexRT - target.ApexRT) < apexTolerance);
            foreach (XIC xic in xicMatchedByApex)
            {
                if(xic == target)
                {
                    group.XIClist.Add(xic);
                    xic.Group = group;
                    continue;
                }
                if (xic.Group.XIClist.Count == 0)
                {
                    double corr = GetCorr(target.XICpeaks, xic.XICpeaks, 0);
                    if (corr > 0.8)
                    {
                        group.XIClist.Add(xic);
                        xic.Group = group;
                    }
                }
            }
            target.Group = group;
            return group;
        }

        public static List<XICgroup> FindAllGroups(List<XIC> allXICs, double apexTolerance)
        {
            var allXICgroups = new List<XICgroup>();    
            foreach(XIC xic in allXICs)
            {
                if (xic.Group.XIClist.Count != 0)
                {
                    continue;
                }
                var group = FindGroup(xic, allXICs, apexTolerance);
                allXICgroups.Add(group);
            }
            return allXICgroups;
        }
        
        //used in grouping ms1
        public static List<XICgroup> FindAllGroups2(List<XIC> allXICs, int apexScanTolerance, double corrCutOff)
        {
            var allGroups = new List<XICgroup>();
            var sortedXICs = allXICs.OrderBy(x => x.ApexScanNumber).ThenBy(x => x.AveragedIntensity).ToArray();
            for(int i = 0; i < sortedXICs.Length; i++)
            {
                if (sortedXICs[i].Group == null)
                {
                    var newGroup = new XICgroup(new List<XIC> { allXICs[i] });
                    allXICs[i].Group = newGroup;
                    int j = i;
                    //go left
                    while (j >= 0 && allXICs[j].ApexScanNumber >= allXICs[i].ApexScanNumber - apexScanTolerance)
                    {
                        double corr = GetCorr(allXICs[i].Group.ReferenceXIC.XICpeaks, allXICs[j].XICpeaks, 0);
                        if (corr > corrCutOff && corr > allXICs[j].CorrelationToGroup)
                        {
                            newGroup.XIClist.Add(allXICs[j]);
                            allXICs[j].Group = newGroup;
                            allXICs[j].CorrelationToGroup = corr;
                            allXICs[j].CorrelatedXIC = allXICs[i].Group.ReferenceXIC;
                        }
                        j--;
                    }
                    //go right
                    int k = i;
                    while (k < sortedXICs.Length && allXICs[k].ApexScanNumber <= allXICs[i].ApexScanNumber + apexScanTolerance)
                    {
                        double corr = GetCorr(allXICs[i].Group.ReferenceXIC.XICpeaks, allXICs[k].XICpeaks, 0);
                        if (corr > corrCutOff && corr > allXICs[k].CorrelationToGroup)
                        {
                            newGroup.XIClist.Add(allXICs[k]);
                            allXICs[k].Group = newGroup;
                            allXICs[k].CorrelationToGroup = corr;
                            allXICs[k].CorrelatedXIC = allXICs[i].Group.ReferenceXIC;
                        }
                        k++;
                    }
                    allGroups.Add(newGroup);
                }
            }
            return allGroups;
        }

        public static double GetCorr(List<Peak> XIC1, List<Peak> XIC2, double rtShift)
        {
            var RT_1 = XIC1.OrderBy(p => p.RT).Select(p => Math.Round(p.RT, 2)).ToArray();
            var RT_2 = XIC2.OrderBy(p => p.RT).Select(p => Math.Round((p.RT + rtShift), 2)).ToArray();

            if (RT_1 == null || RT_2 == null)
            {
                return double.NaN;
            }

            var ms1Intensity = new List<double>();
            var ms2Intensity = new List<double>();
            for (int i = 0; i < RT_1.Length; i++)
            {
                int index = Array.BinarySearch(RT_2, RT_1[i]);
                if (index >= 0)
                {
                    ms1Intensity.Add(XIC1[i].Intensity);
                    ms2Intensity.Add(XIC2[index].Intensity);
                }
            }
            if (ms1Intensity.Count >= 5 && ms2Intensity.Count >= 5)
            {
                // Calculate Pearson correlation
                double correlation = Correlation.Pearson(ms1Intensity, ms2Intensity);
                return correlation;
            }
            else
            {
                return double.NaN;
            }
        }

        public static List<XICgroup> GroupWithApex(List<XIC> allXICs, double corrCutOff)
        {
            var filteredXICs = allXICs.Where(x => x.XICpeaks.Count > 0);
            var groups = filteredXICs.GroupBy(x => x.ApexRT);
            var allXICgroups = new List<XICgroup>();
            foreach(var group in groups)
            {
                var XIClist = filteredXICs.Where(x => x.ApexRT == group.Key).ToList();
                var highestXIC = XIClist.OrderByDescending(x => x.XICpeaks.Sum(p => p.Index)).First();
                var xicGroup = new XICgroup();
                xicGroup.XIClist.Add(highestXIC);
                highestXIC.Group = xicGroup;
                foreach (var xic in XIClist)
                {
                    double corr = GetCorr(highestXIC.XICpeaks, xic.XICpeaks, 0);
                    if (corr > corrCutOff)
                    {
                        xicGroup.XIClist.Add(xic);
                        xic.Group = xicGroup;
                    }
                }
                allXICgroups.Add(xicGroup);
            }
            return allXICgroups;
        }

        public static List<MzSpectrum> GenerateNewMs1 (List<XICgroup> allGroups)
        {
            var allSpectrum = new List<MzSpectrum>();
            foreach (var group in allGroups)
            {
                var sortedXICs = group.XIClist.OrderBy(x => x.AveragedMz);
                var mz = sortedXICs.Select(x => x.AveragedMz).ToArray();
                var intensities = sortedXICs.Select(x => x.AveragedIntensity).ToArray();
                MzSpectrum ms1spectrum = new MzSpectrum(mz, intensities, false);
                allSpectrum.Add(ms1spectrum);
            }
            return allSpectrum;
        }

        public static List<IsotopicEnvelope>[] DeconvoluteNewMs1Scans(XICgroup[] allGroups, CommonParameters commonParameters)
        {
            var allNewMs1 = GenerateNewMs1(allGroups.ToList()).ToArray();
            var deconResults = new List<IsotopicEnvelope>[allGroups.Length];
            for(int i = 0; i < allGroups.Length; i++)
            {
                deconResults[i] = new List<IsotopicEnvelope>();
                var ms1 = allNewMs1[i];
                deconResults[i].AddRange(Deconvoluter.Deconvolute(ms1, commonParameters.PrecursorDeconvolutionParameters));
            }
            return deconResults;
        }

        public static void VisualizeXICgroups(List<(double rt, double intensity, double mz, double corr)> XICs, string outputPath)
        {
            using (var sw = new StreamWriter(File.Create(outputPath)))
            {
                sw.WriteLine("Retention Time,Intensity,rounded_mz,corr");
                foreach (var xic in XICs)
                {
                    sw.WriteLine($"{xic.rt},{xic.intensity},{xic.mz},{xic.corr}");
                }
            }
        }
    }
}
