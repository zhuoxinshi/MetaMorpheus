using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PfGroupResult
    {
        public int PfGroupIndex { get; set; }
        public double MonoisotopicMass { get; set; }
        public double PrecursorIntensity { get; set; }
        public double RetentionTime { get; set; }
        public int Charge { get; set; }
        public List<double> NeutralExperimentalFragments { get; set; }
        
        public PfGroupResult(PrecursorFragmentsGroup pfGroup)
        {
            PfGroupIndex = pfGroup.PFgroupIndex;
            MonoisotopicMass = pfGroup.PrecursorXic.ApexPeak.M;
            PrecursorIntensity = pfGroup.PrecursorXic.ApexPeak.Intensity;
            RetentionTime = pfGroup.PrecursorXic.ApexPeak.RetentionTime;
            Charge = pfGroup.PrecursorXic.Peaks.First() is IndexedMass im ? im.Charge : 1;
        }
    }
}
