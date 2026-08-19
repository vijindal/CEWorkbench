package org.ce.model.equilibrium;

import java.util.List;

/**
 * Pure-element reference energy ("lattice stability"), the {@code G0m} term
 * of {@code G = G0m + Gm}.
 *
 * <p>{@code Gm} — the CVM mixing energy computed by
 * {@link org.ce.model.cvm.CVMGibbsModel} — has no absolute reference zero of
 * its own; two phases' {@code Gm} values are only comparable once each is
 * anchored to the same pure-element baseline. Single-phase validation never
 * needed this (only one phase's relative mixing energy at one composition
 * was ever checked), but multi-phase equilibrium does: chemical potentials
 * are compared/equalized <em>across</em> phases, so every phase's {@code G}
 * must share the same zero.</p>
 *
 * <h2>Now a thin façade over {@link SgteDatabase}</h2>
 *
 * <p>This class formerly hardcoded ~60 SGTE/Dinsdale unary polynomials, one
 * piecewise function per (element, structure), hand-transcribed from a
 * Mathematica reference. Those are gone: every method here now delegates to
 * {@link SgteDatabase}, which parses the same polynomials directly from the
 * SGTE Unary v4.4 database at {@code inputs/unary.dat}. The API is unchanged,
 * so callers need not know the difference.</p>
 *
 * <p>Two reasons for the switch:</p>
 *
 * <ul>
 *   <li><b>Coverage.</b> The hardcoded set was eight elements
 *       (Mo, Nb, Re, Ta, Ti, V, W, Zr) across four phases, and not every
 *       combination existed. The database covers the whole SGTE unary set,
 *       so a new element costs nothing.</li>
 *   <li><b>Transcription risk.</b> A hand-copied coefficient can silently
 *       lose a digit, and no internal check can detect it — every consumer
 *       reads the same wrong number. Reading the source removes that class of
 *       error entirely.</li>
 * </ul>
 *
 * <p>The two implementations were cross-checked before the hardcoded copy was
 * removed: {@code org.ce.scratch.SgteCrossCheck} compared {@code g0} and
 * {@code dG0Dt} across 8 elements × 4 phases × 6 temperatures and found 176
 * exact agreements with zero real disagreements, plus a finite-difference
 * check of each side's derivative against its own energy (180 points per side,
 * all passing). {@code org.ce.scratch.LatticeStabilityVerification} additionally
 * gates both against externally supplied reference values.</p>
 *
 * <p><b>One deliberate behaviour change.</b> SGTE temperature ranges are
 * half-open {@code (lower, upper]}, so {@code T} exactly equal to a stated
 * upper bound belongs to the range <em>ending</em> there. The hardcoded version
 * tested {@code t < upper} and fell through to the next range instead. At an
 * exact breakpoint the two therefore differ by the polynomials' continuity gap
 * — a few mJ/mol, e.g. 3.5e-3 J/mol for Ta at 2500 K. The database convention
 * is the correct one; results at exactly-on-boundary temperatures shift
 * accordingly.</p>
 *
 * <p>Coverage is still not universal, and a missing (element, phase) pair
 * remains a programming error rather than a silently-zero fallback:
 * {@link #g0}/{@link #g0m} throw rather than guess.</p>
 *
 * <p><b>Every derivative here is analytical</b> — term-by-term differentiation
 * of the same closed-form polynomial {@code g0} evaluates, never
 * finite-differenced. That is required throughout the Hillert solver, and it
 * holds under delegation: {@link SgteDatabase} differentiates the SGTE term set
 * in closed form and selects the branch with the same range logic {@code g0}
 * uses at that {@code T}.</p>
 */
public final class LatticeStability {

    private LatticeStability() {}

    /**
     * Pure-element reference energy {@code G0(element, phase, T)} in J/mol.
     *
     * @param element element symbol, e.g. {@code "Ti"}
     * @param phase   phase name, e.g. {@code "BCC_A2"}; an ordered phase maps
     *                to its disordered parent (a pure element has no ordering)
     * @param t       temperature in K
     * @throws IllegalArgumentException if the pair is absent from the database
     *                                  or {@code t} is outside its range
     */
    public static double g0(String element, String phase, double t) {
        return SgteDatabase.g0(element, phase, t);
    }

    /**
     * Analytic {@code dG0/dT} in J/mol/K — the pure-element contribution to
     * {@code GT} in the Hillert solver.
     */
    public static double dG0Dt(String element, String phase, double t) {
        return SgteDatabase.dG0Dt(element, phase, t);
    }

    /** Analytic {@code d2G0/dT2} in J/mol/K^2. */
    public static double d2G0Dt2(String element, String phase, double t) {
        return SgteDatabase.d2G0Dt2(element, phase, t);
    }

    /**
     * Composition-weighted reference energy
     * {@code G0m = Σ x_i * g0(element_i, phase, t)} — the mechanical mixture of
     * pure elements in {@code G = G0m + Gm}. Carries no entropy of mixing.
     */
    public static double g0m(List<String> elements, String phase, double[] composition, double t) {
        return SgteDatabase.g0m(elements, phase, composition, t);
    }

    /**
     * Composition-weighted {@code d(G0m)/dT}, analytic (term-by-term
     * differentiation via {@link #dG0Dt}, same composition weights {@link #g0m}
     * uses) — feeds {@code GT}'s pure-element contribution in the Hillert
     * solver.
     */
    public static double dG0mDt(List<String> elements, String phase, double[] composition, double t) {
        return SgteDatabase.dG0mDt(elements, phase, composition, t);
    }
}
