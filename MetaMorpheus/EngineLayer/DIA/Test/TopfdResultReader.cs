using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class TopfdResultReader
    {
        public static IEnumerable<IndexedMass> ReadMs1AlignFile(string filePath)
        {
            using (var sr = new StreamReader(filePath))
            {
                int zeroBasedScanIndex = 0;
                double rt = 0;

                bool start = false;
                while (sr.ReadLine() is { } line)
                {
                    if (line.Length == 0) continue;
                    if (line.Contains("BEGIN IONS"))
                    {
                        start = true;
                    }
                    if (!start)
                        continue;

                    if (line.Contains("SPECTRUM_ID"))
                    {
                        var splits = line.Split('=');
                        zeroBasedScanIndex = int.Parse(splits[1]); 
                    }
                    if (line.Contains("RETENTION_TIME"))
                    {
                        var splits = line.Split('=');
                        rt = double.Parse(splits[1]);
                    }
                    
                    if (line.Contains("\t"))
                    {
                        var splits = line.Split('\t');
                        var monoMass = double.Parse(splits[0]);
                        var intensity = double.Parse(splits[1]);
                        var charge = int.Parse(splits[2]);

                        var envelope = new IsotopicEnvelope(monoMass, intensity, charge);
                        var indexedMass = new IndexedMass(envelope, rt, zeroBasedScanIndex, 1);
                        yield return indexedMass;
                    }
                }
            }
        }
    }
}
