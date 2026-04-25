package org.ce.debug;

import java.util.Arrays;
import java.util.function.Consumer;

/**
 * Centralised debug toggle and formatted printing for MCS pipeline diagnostics.
 *
 * <p>Set {@link #ENABLED} to {@code false} to silence all MCS debug output
 * without removing the instrumentation calls.</p>
 *
 * <p>Output is routed through an optional {@link Consumer}{@code <String>} sink
 * (the same mechanism the GUI output panel uses). When no sink is set, output
 * falls back to {@code System.out}.</p>
 */
public final class MCSDebug {

    /** Master switch — set to {@code false} to suppress all MCS debug output. */
    public static boolean ENABLED = false;

    /**
     * Thread-local sink so concurrent calculations don't clash.
     * Set via {@link #setSink(Consumer)} before a calculation starts.
     */
    private static final ThreadLocal<Consumer<String>> SINK = new ThreadLocal<>();

    private MCSDebug() {}

    // ── Sink management ───────────────────────────────────────────────────────

    /** Route all debug output through this sink (typically the GUI progressSink). */
    public static void setSink(Consumer<String> sink) {
        SINK.set(sink);
    }

    /** Clear the sink (restores System.out fallback). */
    public static void clearSink() {
        SINK.remove();
    }

    // ── Core output ───────────────────────────────────────────────────────────

    private static void emit(String line) {
        Consumer<String> sink = SINK.get();
        if (sink != null) {
            sink.accept(line);
        } else {
            System.out.println(line);
        }
    }

    // ── Formatted output helpers ──────────────────────────────────────────────

    public static void log(String tag, String msg) {
        if (ENABLED) emit(String.format("[MCS-DBG][%s] %s", tag, msg));
    }

    public static void log(String tag, String fmt, Object... args) {
        if (ENABLED) emit(String.format("[MCS-DBG][" + tag + "] " + fmt, args));
    }

    public static void separator(String title) {
        if (!ENABLED) return;
        emit("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        emit("  [MCS-DBG] " + title);
        emit("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    }

    public static void vector(String tag, String label, double[] v) {
        if (!ENABLED || v == null) return;
        StringBuilder sb = new StringBuilder();
        sb.append(label).append(" [len=").append(v.length).append("]: ");
        for (int i = 0; i < v.length; i++) {
            if (i > 0) sb.append(", ");
            sb.append(String.format("%+.8f", v[i]));
        }
        log(tag, sb.toString());
    }

    public static void vector(String tag, String label, int[] v) {
        if (!ENABLED || v == null) return;
        log(tag, "%s %s", label, Arrays.toString(v));
    }

    public static void matrix(String tag, String label, double[][] m) {
        if (!ENABLED || m == null) return;
        log(tag, "%s [%dx%d]:", label, m.length, m[0].length);
        for (int r = 0; r < m.length; r++) {
            StringBuilder sb = new StringBuilder("  row ").append(r).append(": ");
            for (int c = 0; c < m[r].length; c++) {
                if (c > 0) sb.append(", ");
                sb.append(String.format("%+.6f", m[r][c]));
            }
            log(tag, sb.toString());
        }
    }
}
