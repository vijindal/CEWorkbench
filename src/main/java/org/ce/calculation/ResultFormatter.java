package org.ce.calculation;

import org.ce.model.result.ThermodynamicResult;

import java.util.List;

/**
 * Shared text-formatting utilities for {@link ThermodynamicResult} objects.
 *
 * <p>Used by both the CLI ({@code Main}) and the GUI ({@code LineScanPanel},
 * {@code OutputPanel}) to produce consistent, human-readable output without
 * duplicating format strings across layers.</p>
 *
 * <p>All methods iterate {@link QuantityDescriptor#values()} so that adding a
 * new quantity automatically appears in all formatted output.</p>
 */
public final class ResultFormatter {

    private ResultFormatter() {}   // utility class

    // =========================================================================
    // Single-result formats
    // =========================================================================

    /**
     * Returns a compact one-liner suitable for streaming to the GUI log or CLI
     * verbose output during a scan.
     * Example: {@code "T=1000.0 K  G=-12345.6789  H=-10000.0000  S= 2.3456"}
     */
    public static String oneLine(ThermodynamicResult r) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("T=%7.1f K", r.temperature));
        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            if (q.available(r)) {
                sb.append(String.format("  %s=%12.4f", q.symbol, q.extract(r)));
            }
        }
        return sb.toString();
    }

    /**
     * Returns a verbose multi-line block for a single equilibrium result,
     * suitable for CLI {@code calc_min} output and the GUI result log.
     * Example:
     * <pre>
     * EQUILIBRIUM PROPERTIES
     * ─────────────────────────────────────────────────────
     *   T         :     1000.0 K
     *   x         :  [0.5000, 0.5000]
     *   G (J/mol) :   -12345.6789012345
     *   H (J/mol) :   -10000.0000000000
     *   S (J/mol·K):        2.3456789012
     * </pre>
     */
    public static String fullBlock(ThermodynamicResult r) {
        StringBuilder sb = new StringBuilder();
        sb.append("EQUILIBRIUM PROPERTIES\n");
        sb.append("─".repeat(56)).append('\n');
        sb.append(String.format("  %-12s: %14.1f K%n", "T", r.temperature));

        // Composition
        StringBuilder xStr = new StringBuilder("[");
        for (int i = 0; i < r.composition.length; i++) {
            if (i > 0) xStr.append(", ");
            xStr.append(String.format("%.4f", r.composition[i]));
        }
        xStr.append(']');
        sb.append(String.format("  %-12s: %s%n", "x", xStr));

        // Dynamic quantity rows
        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            if (q.available(r)) {
                String label = q.symbol + " (" + q.unit + ")";
                sb.append(String.format("  %-16s: %20.10f%n", label, q.extract(r)));
            }
        }
        return sb.toString();
    }

    // =========================================================================
    // Table formats (for line scans)
    // =========================================================================

    /**
     * Returns the table header row matching {@link #tableRow(ThermodynamicResult)}.
     * Column widths are fixed so headers and rows align.
     */
    public static String tableHeader() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("  %-10s  %-22s", "T (K)", "x"));
        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            // include all potential columns so header is stable regardless of result
            sb.append(String.format("  %-16s", q.symbol + " (" + q.unit + ")"));
        }
        return sb.toString();
    }

    /**
     * Returns a separator line for use between the table header and data rows.
     */
    public static String tableSeparator() {
        return "  " + "-".repeat(10 + 2 + 22 + (18 * QuantityDescriptor.values().length));
    }

    /**
     * Returns one fixed-width data row for the given result.
     * Columns whose values are {@code NaN} (not produced by this engine) show {@code "—"}.
     */
    public static String tableRow(ThermodynamicResult r) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("  %-10.2f", r.temperature));

        // Composition: x₂ for binary, full array for ternary+
        if (r.composition.length <= 2) {
            double x2 = r.composition.length > 1 ? r.composition[1] : r.composition[0];
            sb.append(String.format("  %-22.6f", x2));
        } else {
            StringBuilder xStr = new StringBuilder();
            for (int i = 0; i < r.composition.length; i++) {
                if (i > 0) xStr.append(',');
                xStr.append(String.format("%.4f", r.composition[i]));
            }
            sb.append(String.format("  %-22s", xStr));
        }

        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            if (q.available(r)) {
                sb.append(String.format("  %-16.6f", q.extract(r)));
            } else {
                sb.append(String.format("  %-16s", "\u2014"));
            }
        }
        return sb.toString();
    }

    /**
     * Returns a full formatted table (header + separator + rows) for a list of
     * results, ready to print to stdout or append to the GUI log.
     */
    public static String table(List<ThermodynamicResult> results) {
        StringBuilder sb = new StringBuilder();
        sb.append(tableHeader()).append('\n');
        sb.append(tableSeparator()).append('\n');
        for (ThermodynamicResult r : results) {
            sb.append(tableRow(r)).append('\n');
        }
        return sb.toString();
    }
}
