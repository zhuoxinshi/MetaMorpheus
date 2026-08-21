using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections;
using System.Collections.Generic;
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
        public PseudoMs2IntensityType PseudoMs2IntensityType { get; set; } //Which XIC number becomes the intensity of a pseudo MS2 peak.

        public DIAparameters(DIAanalysisType analysisType, XicConstructor ms1XicConstructor, XicConstructor ms2XicConstructor, PfGroupingEngine pfGroupingEngine, PseudoMs2ConstructionType pseudoMs2ConstructionType, bool combineFragments = false, PseudoMs2IntensityType pseudoMs2IntensityType = PseudoMs2IntensityType.XicSum)
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
        /// Builds a human-readable dump of every setting used by the DIA/ISD engine, including the
        /// concrete types and all scalar settings of the nested XIC constructors, grouping engine,
        /// deconvolution parameters and XIC spline. This is intended as a record of exactly which DIA
        /// settings produced a given search, written next to the search output.
        /// </summary>
        public string GetSettingsString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("DIA/ISD Engine Settings");
            sb.AppendLine("=======================");
            sb.AppendLine("Written: " + DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss"));
            sb.AppendLine();
            AppendObject(this, sb, indent: 0, depth: 0, maxDepth: 2);
            return sb.ToString();
        }

        /// <summary>
        /// Writes <see cref="GetSettingsString"/> to <paramref name="filePath"/>, creating the directory if needed.
        /// </summary>
        public void WriteSettingsToFile(string filePath)
        {
            var directory = Path.GetDirectoryName(filePath);
            if (!string.IsNullOrEmpty(directory) && !Directory.Exists(directory))
                Directory.CreateDirectory(directory);

            File.WriteAllText(filePath, GetSettingsString());
        }

        /// <summary>
        /// Recursively appends the public scalar properties of <paramref name="obj"/> (and, up to
        /// <paramref name="maxDepth"/>, its nested settings objects) to <paramref name="sb"/>. Scalars are
        /// printed first, then nested settings objects as labelled sub-sections. Data collections are skipped.
        /// </summary>
        private static void AppendObject(object obj, StringBuilder sb, int indent, int depth, int maxDepth)
        {
            if (obj == null)
                return;

            string pad = new string(' ', indent * 4);
            var properties = obj.GetType()
                .GetProperties(BindingFlags.Public | BindingFlags.Instance)
                .Where(p => p.GetIndexParameters().Length == 0)
                .ToList();

            // pass 1: scalar settings
            foreach (var property in properties)
            {
                Type underlyingType = Nullable.GetUnderlyingType(property.PropertyType) ?? property.PropertyType;
                if (!IsScalarType(underlyingType))
                    continue;

                object value = TryGetValue(obj, property);
                sb.AppendLine($"{pad}{property.Name} = {FormatScalar(value)}");
            }

            // pass 2: nested settings objects
            foreach (var property in properties)
            {
                Type declaredType = property.PropertyType;
                if (IsScalarType(Nullable.GetUnderlyingType(declaredType) ?? declaredType) || !IsSettingsType(declaredType))
                    continue;

                object value = TryGetValue(obj, property);
                if (value == null)
                {
                    sb.AppendLine($"{pad}{property.Name}: (null)");
                    continue;
                }

                // print the CONCRETE type so polymorphic choices (e.g. which XicConstructor) are recorded
                sb.AppendLine($"{pad}{property.Name}: {value.GetType().Name}");
                if (depth < maxDepth)
                    AppendObject(value, sb, indent + 1, depth + 1, maxDepth);
            }
        }

        private static bool IsScalarType(Type type)
        {
            return type.IsPrimitive
                || type.IsEnum
                || type == typeof(string)
                || type == typeof(decimal)
                || typeof(Tolerance).IsAssignableFrom(type);
        }

        private static bool IsSettingsType(Type type)
        {
            if (type == typeof(string) || typeof(IEnumerable).IsAssignableFrom(type))
                return false;

            string ns = type.Namespace ?? string.Empty;
            return type.IsClass
                && (ns.StartsWith("EngineLayer") || ns.StartsWith("MassSpectrometry") || ns.StartsWith("Chemistry"));
        }

        private static object TryGetValue(object obj, PropertyInfo property)
        {
            try
            {
                return property.GetValue(obj);
            }
            catch
            {
                return null;
            }
        }

        private static string FormatScalar(object value)
        {
            if (value == null)
                return "(null)";
            if (value is IFormattable formattable)
                return formattable.ToString(null, System.Globalization.CultureInfo.InvariantCulture);
            return value.ToString();
        }
    }
}
