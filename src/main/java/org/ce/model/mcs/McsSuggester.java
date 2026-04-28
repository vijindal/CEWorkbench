package org.ce.model.mcs;

/**
 * Suggests MCS run parameters (L, nEquil, nAvg) for statistically meaningful results.
 *
 * <p>Rules are based on standard alloy Monte Carlo practice:
 * <ul>
 *   <li><b>L</b> — supercell edge length. BCC gives N = 2·L³ sites. Minimum L=6 (N=432)
 *       is recommended for production; smaller L suffers from finite-size effects and
 *       periodic-image correlation. L=8 (N=1024) gives better statistics at ~4× cost.</li>
 *   <li><b>nEquil</b> — equilibration sweeps (each sweep = N trial moves). Must exceed
 *       the slowest correlation time. Heuristic: at least 5 sweeps per site, minimum 500.</li>
 *   <li><b>nAvg</b> — averaging sweeps. Statistical error scales as 1/√nAvg (after
 *       accounting for autocorrelation). Heuristic: at least 10 sweeps per site, minimum 1000.</li>
 * </ul>
 *
 * <p>Users may override any value — these are starting points, not hard constraints.
 */
public final class McsSuggester {

    private McsSuggester() {}

    /** Immutable suggestion result. */
    public record Suggestion(int L, int nEquil, int nAvg, String rationale) {}

    /**
     * Returns suggested parameters for the given L.
     * Call this whenever L changes so nEquil/nAvg scale with system size.
     */
    public static Suggestion suggest(int L) {
        int validL  = Math.max(4, L);
        int N       = 2 * validL * validL * validL;   // BCC sites

        // Equilibration: at least 5 sweeps per site, minimum 500
        int nEquil  = Math.max(500, 5 * N);
        // Averaging: at least 10 sweeps per site, minimum 1000
        int nAvg    = Math.max(1000, 10 * N);

        String rationale = String.format(
                "L=%d → N=%d sites (BCC 2·L³).  "
              + "nEquil=%d (≥5N, min 500) ensures thermal equilibration.  "
              + "nAvg=%d (≥10N, min 1000) gives ~1%% statistical error on ⟨H⟩.",
                validL, N, nEquil, nAvg);

        return new Suggestion(validL, nEquil, nAvg, rationale);
    }

    /**
     * Returns the default suggested parameters for a new MCS session (L=6).
     * L=6 (N=432) is the recommended minimum for production runs —
     * large enough to avoid finite-size bias, small enough to run in seconds.
     */
    public static Suggestion defaultSuggestion() {
        return suggest(6);
    }

    /**
     * Returns a human-readable one-line summary suitable for a tooltip or hint label.
     */
    public static String hint(int L) {
        int N      = 2 * L * L * L;
        int nEquil = Math.max(500,  5 * N);
        int nAvg   = Math.max(1000, 10 * N);
        return String.format("N=%d sites · equil %d · avg %d  (adjust L to rescale)", N, nEquil, nAvg);
    }
}
