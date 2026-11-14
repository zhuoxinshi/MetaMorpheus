using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DIAparameters
    {
        public DIAanalysisType AanalysisType { get; set; }
        public XicConstructor Ms1XicConstructor { get; set; }
        public XicConstructor Ms2XicConstructor { get; set; } 
        public PfGroupingEngine PfGroupingEngine { get; set; }
        public PseudoMs2ConstructionType PseudoMs2ConstructionType { get; set; }
        public bool CombineFragments { get; set; } //If true, combine fragment XICs from different groups into one prior to PF grouping.
        public PseudoMs2IntensityType PseudoMs2IntensityType { get; set; } //Which XIC value becomes the intensity of a pseudo MS2 peak.

        public DIAparameters(DIAanalysisType analysisType, XicConstructor ms1XicConstructor, XicConstructor ms2XicConstructor, PfGroupingEngine pfGroupingEngine, PseudoMs2ConstructionType pseudoMs2ConstructionType, bool combineFragments = false,
            PseudoMs2IntensityType pseudoMs2IntensityType = PseudoMs2IntensityType.ApexPeak)
        {
            AanalysisType = analysisType;
            Ms1XicConstructor = ms1XicConstructor;
            Ms2XicConstructor = ms2XicConstructor;
            PfGroupingEngine = pfGroupingEngine;
            PseudoMs2ConstructionType = pseudoMs2ConstructionType;
            CombineFragments = combineFragments;
            PseudoMs2IntensityType = pseudoMs2IntensityType;
        }

        /// <summary>
        /// A readable record of every DIA/ISD setting: the scalar settings of these parameters and of the nested
        /// XIC constructors, XIC spline and grouping engine, with the concrete type of each nested object.
        /// </summary>
        public string GetSettingsString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("DIA/ISD Engine Settings");
            sb.AppendLine("=======================");
            sb.AppendLine("Written: " + DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss", CultureInfo.InvariantCulture));
            sb.AppendLine();
            AppendSettings(this, sb, indent: 0, depth: 0, maxDepth: 2);
            return sb.ToString();
        }

        /// <summary>
        /// Writes <see cref="GetSettingsString"/> to <paramref name="filePath"/>, creating the folder if needed.
        /// </summary>
        public void WriteSettingsToFile(string filePath)
        {
            var directory = Path.GetDirectoryName(filePath);
            if (!string.IsNullOrEmpty(directory))
                Directory.CreateDirectory(directory);
            File.WriteAllText(filePath, GetSettingsString());
        }

        // Scalar properties first, then nested settings objects as indented sections. Collections are skipped.
        private static void AppendSettings(object obj, StringBuilder sb, int indent, int depth, int maxDepth)
        {
            string pad = new string(' ', indent * 4);
            var properties = obj.GetType()
                .GetProperties(BindingFlags.Public | BindingFlags.Instance)
                .Where(p => p.GetIndexParameters().Length == 0)
                .ToList();

            foreach (var property in properties.Where(p => IsScalarType(p.PropertyType)))
            {
                sb.AppendLine($"{pad}{property.Name} = {FormatScalar(property.GetValue(obj))}");
            }

            foreach (var property in properties.Where(p => !IsScalarType(p.PropertyType) && IsSettingsType(p.PropertyType)))
            {
                object value = property.GetValue(obj);
                if (value == null)
                {
                    sb.AppendLine($"{pad}{property.Name}: (null)");
                    continue;
                }
                // the concrete type records which implementation (e.g. which XicConstructor) was used
                sb.AppendLine($"{pad}{property.Name}: {value.GetType().Name}");
                if (depth < maxDepth)
                    AppendSettings(value, sb, indent + 1, depth + 1, maxDepth);
            }
        }

        private static bool IsScalarType(Type type)
        {
            type = Nullable.GetUnderlyingType(type) ?? type;
            return type.IsPrimitive || type.IsEnum || type == typeof(string) || type == typeof(decimal)
                || typeof(Tolerance).IsAssignableFrom(type);
        }

        private static bool IsSettingsType(Type type)
        {
            if (!type.IsClass || type == typeof(string) || typeof(IEnumerable).IsAssignableFrom(type))
                return false;
            string ns = type.Namespace ?? string.Empty;
            return ns.StartsWith("EngineLayer") || ns.StartsWith("MassSpectrometry") || ns.StartsWith("Chemistry");
        }

        private static string FormatScalar(object value)
        {
            return value switch
            {
                null => "(null)",
                IFormattable formattable => formattable.ToString(null, CultureInfo.InvariantCulture),
                _ => value.ToString(),
            };
        }
    }
}
