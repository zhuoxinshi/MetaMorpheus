using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EngineLayer.DIA.FeatureAlign
{
    /// <summary>
    /// Minimal reader for a FeatureAlign consensus TSV, reproducing the column selection and row
    /// filtering that the manuscript Python scripts perform (pool_3file.py, pool_paired.py, s2_yeast_de.py).
    /// Kept independent of <see cref="FeatureAligner"/>'s own writer so the real-data cross-check tests
    /// two independent implementations against each other.
    ///
    /// Quant columns are named "quant_&lt;runKey&gt;_ms1"; the &lt;runKey&gt; is exactly the key used in the
    /// MS1 TIC-area table, so <see cref="RunKey"/> strips the prefix/suffix to join the two.
    /// </summary>
    public class FeatureAlignQuantTable
    {
        private const string QuantPrefix = "quant_";
        private const string QuantSuffix = "_ms1";

        public string[] Header { get; }
        private readonly List<string[]> _rows;
        private readonly int _massCol;

        public FeatureAlignQuantTable(string path)
        {
            var lines = File.ReadAllLines(path);
            if (lines.Length == 0) throw new InvalidDataException($"empty file: {path}");
            Header = lines[0].Split('\t');
            _massCol = Array.IndexOf(Header, "mass_mean");
            if (_massCol < 0) throw new InvalidDataException($"no mass_mean column in {path}");
            _rows = lines.Skip(1).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
        }

        public int RowCount => _rows.Count;

        /// <summary>Indices of quant_ columns whose header contains <paramref name="mustContain"/>, in file order.</summary>
        public List<int> QuantColumns(string mustContain) =>
            Enumerable.Range(0, Header.Length)
                      .Where(i => Header[i].StartsWith(QuantPrefix, StringComparison.Ordinal)
                                  && Header[i].Contains(mustContain))
                      .ToList();

        /// <summary>The TIC-table key for a quant column: "quant_KEY_ms1" -> "KEY".</summary>
        public string RunKey(int col)
        {
            string h = Header[col];
            if (!h.StartsWith(QuantPrefix, StringComparison.Ordinal) || !h.EndsWith(QuantSuffix, StringComparison.Ordinal))
                throw new InvalidOperationException($"column '{h}' is not a quant_..._ms1 column");
            return h.Substring(QuantPrefix.Length, h.Length - QuantPrefix.Length - QuantSuffix.Length);
        }

        /// <summary>
        /// Feature matrix over <paramref name="cols"/>: every row whose mass_mean parses and is
        /// &gt;= <paramref name="minMass"/> and whose selected values all parse and are &gt; 0.
        /// Mirrors the Python fval + "all(v &gt; 0)" filtering; row order is file order.
        /// </summary>
        public List<double[]> FeatureMatrix(IReadOnlyList<int> cols, double minMass)
        {
            var outRows = new List<double[]>();
            foreach (var r in _rows)
            {
                if (!TryCell(r, _massCol, out double mass) || mass < minMass) continue;
                var v = new double[cols.Count];
                bool ok = true;
                for (int k = 0; k < cols.Count; k++)
                {
                    if (!TryCell(r, cols[k], out double x) || x <= 0.0) { ok = false; break; }
                    v[k] = x;
                }
                if (ok) outRows.Add(v);
            }
            return outRows;
        }

        private static bool TryCell(string[] row, int col, out double val)
        {
            val = 0.0;
            if (col < 0 || col >= row.Length) return false;
            return double.TryParse(row[col], NumberStyles.Float, CultureInfo.InvariantCulture, out val);
        }
    }
}
