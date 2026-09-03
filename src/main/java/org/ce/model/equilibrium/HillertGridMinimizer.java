package org.ce.model.equilibrium;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cvm.CVMGibbsModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Consumer;

/**
 * V2 STEP 7 -- the Sundman grid minimizer as a <b>standalone initial-estimate
 * component</b> for a single fixed-{@code T}/{@code P} closed-system equilibrium.
 *
 * <p>This is the finite-grid global search of
 * <b>Sundman, Lu, Ohtani (2015) &sect;5.1</b> ("Obtaining start values by grid
 * minimizer") and <b>Sundman et&nbsp;al. (2021) &sect;2.3.3</b> ("Global
 * equilibrium calculation") -- <em>not</em> a generic nonlinear optimizer, and
 * <em>not</em> part of Algorithm&nbsp;A. Its job is to hand
 * {@link HillertSolver#solve} a good global starting point:
 *
 * <pre>
 *   phases + T + P + overall component amounts
 *        &darr;
 *   a finite grid of fixed-composition pseudo-phases per phase
 *        &darr;
 *   the &le;K grid points whose non-negative mixture reaches the target
 *   with the lowest total Gibbs energy         (discrete convex hull)
 *        &darr;
 *   selected stable phase set + amounts + constitution Y for EVERY phase
 *   (stable and metastable) + discrete chemical potentials
 * </pre>
 *
 * <h2>What a grid point is (STEP 7 PART 4)</h2>
 * <p>For a V1 disordered CVM phase, a grid point is a <b>fixed bulk composition
 * {@code x}</b> with its <b>{@code ncf} internal non-point CFs relaxed</b> by
 * {@link CvmNewtonSolver} (the hardened fixed-composition minimiser). It is
 * <em>not</em> a composition-only point (that would ignore the internal CV
 * degrees of freedom and misrepresent the Gibbs surface), and <em>not</em> a
 * full-{@code Y}-at-arbitrary-fixed-CFs point (that over-expands the grid
 * dimension). The K&minus;1 composition coordinates parametrise the grid; the
 * {@code ncf} internal CFs are <em>solved</em>, not gridded. The stored
 * {@link GridPoint#y} is the width-{@code ncf+K} vector {@code [relaxedU ; x]} --
 * exactly the {@link HillertSolver.Phase#uFull} layout, so a selected point
 * initialises a {@code Phase} with no transformation.</p>
 *
 * <h2>The discrete problem (STEP 7 PART 1 / PART 8)</h2>
 * <p>Every grid point {@code q} is a stoichiometric pseudo-phase with amount
 * {@code n_q}. With the target overall composition normalised to one mole
 * ({@code sum_A b_A = 1}) and {@code sum_A M_A^q = 1} for every registered CVM
 * phase (verified in {@link org.ce.model.cvm.CvmGeometry} -- site multiplicities
 * are normalised so {@code R = 1}), the finite problem is
 *
 * <pre>
 *   minimise  sum_q n_q G_q
 *   s.t.      sum_q n_q M_q,A = b_A     for A = 0..K-1
 *             sum_q n_q       = 1       (implied by the two above; kept explicit
 *                                        as the augmenting row of each solve)
 *             n_q &gt;= 0
 * </pre>
 *
 * <p>Solved by <b>exact vertex enumeration</b> (no LP library): every 1-point,
 * 2-point (K=2) or up-to-3-point (K=3) subset whose {@code M} vectors can
 * bracket / hull the target is solved as a small square linear system via
 * {@link LinearAlgebra#solveChecked}; infeasible ({@code n_q < -AMOUNT_TOL}),
 * non-conservative or near-singular subsets are rejected; the feasible subset
 * with the lowest {@code sum n_q G_q} wins, ties broken by the
 * lexicographically-smallest sorted grid-index tuple (STEP 7 PART 20). A
 * non-degenerate fixed-{@code T}/{@code P} optimum has at most {@code K}
 * independent grid points (Sundman &sect;5.6, Gibbs phase rule).</p>
 *
 * <h2>Miscibility gaps (STEP 7 PART 11)</h2>
 * <p>Two selected grid points with the <em>same</em> parent phase but distinct
 * compositions are kept as <b>two distinct phase instances</b> ({@code A#1},
 * {@code A#2}) -- never merged on parent identity alone. They are merged only
 * when their compositions agree within {@link #GAP_MERGE_TOL}.</p>
 *
 * <h2>Metastable initialisation (STEP 7 PART 12-14)</h2>
 * <p>Discrete chemical potentials {@code mu} are solved from
 * {@code sum_A M_q,A mu_A = G_q} over the selected points (same absolute / SER
 * convention as {@link HillertSolver}). For every <em>non-selected</em> phase the
 * starting constitution is its grid point of minimum tangent-plane distance
 * {@code |G_q - sum_A mu_A M_q,A|}. The result therefore carries an initial
 * {@code Y} for <b>every</b> candidate phase, not only the stable ones.</p>
 *
 * <h2>Boundaries of this step</h2>
 * <ul>
 *   <li><b>Not wired into {@link HillertSolver}.</b> A future step decides how a
 *       {@link GridMinimizationResult} seeds {@code Phase} objects and
 *       {@code solve(...)}. This class has no dependency on the Algorithm-A
 *       iteration, {@link HillertSolver.PhaseStep} or
 *       {@link HillertSolver.EquilibriumMatrix}, and they have none on it.</li>
 *   <li><b>K = 2 and K = 3 only.</b> {@code K >= 4} throws.</li>
 *   <li><b>Disordered V1 phases only.</b> {@link GridPoint} is shaped to hold a
 *       full ordered {@code Y}, but no ordered grid is built -- an ordered phase
 *       ({@code geometry.tcf > ncf + K}) is reported as skipped-unsupported. The
 *       V1 {@code UNSUPPORTED_PHASE_MODEL} guard in {@link HillertSolver} is
 *       untouched.</li>
 *   <li><b>No external LP / optimiser library.</b> Only {@link LinearAlgebra}.</li>
 * </ul>
 *
 * <p>Deterministic: every loop is over an explicit index order, every tie is
 * broken explicitly, no {@code HashMap} iteration or object identity is
 * observed.</p>
 */
public final class HillertGridMinimizer {

    /**
     * Grid density {@code N}: the number of equal subdivisions of a composition
     * edge. A separate constant from anything in {@link HillertSolver}'s
     * in-loop candidate search ({@code candidateCompositionGrid}) -- this is a
     * different algorithm (a global pre-processing grid, not an in-iteration
     * probe). Point counts: {@code K=2 -> N+1}; {@code K=3 -> (N+1)(N+2)/2}.
     * {@code N = 20} gives 21 (binary) / 231 (ternary) points per phase.
     */
    public static final int GRID_DIVISIONS = 20;

    /** A grid-point amount below {@code -AMOUNT_TOL} makes a subset infeasible. */
    public static final double AMOUNT_TOL = 1.0e-9;

    /**
     * Mass-balance check on an accepted subset: {@code max_A |sum_q n_q M_q,A -
     * b_A|} must not exceed this (guards against an ill-conditioned solve that
     * {@link LinearAlgebra#solveChecked} nonetheless returned).
     */
    public static final double CONSERVATION_TOL = 1.0e-7;

    /**
     * Relative-residual ceiling for a subset's linear solve. A subset whose
     * augmented {@code M} system comes back from
     * {@link LinearAlgebra#solveChecked} with a relative residual worse than
     * this is rejected as a degenerate simplex even though elimination returned
     * a number (STEP 7 PART 9). {@link LinearAlgebra#solveChecked} throws on an
     * exactly-singular matrix, which is caught and also treated as rejection.
     */
    public static final double MAX_SOLVE_RESIDUAL = 1.0e-8;

    /**
     * Two selected grid points of the same parent phase are one instance if
     * their compositions agree to within this {@code inf}-norm; otherwise they
     * are a miscibility-gap pair kept as {@code A#1}, {@code A#2}.
     */
    public static final double GAP_MERGE_TOL = 1.0e-6;

    /** Objectives within this (absolute, J/mol) are a tie -> lexicographic tie-break. */
    public static final double OBJECTIVE_TIE_TOL = 1.0e-6;

    /**
     * Convergence tolerance handed to the per-grid-point {@link CvmNewtonSolver}
     * (the {@code sum |dGm/du|} threshold). Matches the order of the production
     * single-point path ({@code ThermodynamicWorkflow} uses {@code 1e-5}); a
     * much tighter value ({@code 1e-10}) is near machine precision for this
     * objective and makes {@link CvmNewtonSolver} report non-convergence at
     * points whose Gibbs energy is in fact perfectly usable.
     */
    public static final double GRID_POINT_SOLVE_TOL = 1.0e-6;

    /**
     * Same-parent merge energy tolerance (V2 STEP 10). When a same-parent
     * selected group's single-phase state at the combined composition has
     * {@code G_single <= G_grid_group + MERGE_ENERGY_TOL}, the group is
     * discretisation of one continuous phase and is merged; otherwise it is a
     * genuine miscibility gap and stays split (see {@link #mergeSameParentGroup}).
     *
     * <p><b>Derivation.</b> The comparison is {@code G_single - sum_q n_q G_q},
     * a difference of absolute Gibbs energies of order {@code 1e4-1e5 J/mol}. The
     * error budget it must absorb:</p>
     * <ul>
     *   <li><b>CVM Newton convergence.</b> Both {@code G_single} and each
     *       {@code G_q} come from a {@link CvmNewtonSolver} solve to
     *       {@link #GRID_POINT_SOLVE_TOL} {@code = 1e-6} on {@code sum |dGm/du|}.
     *       Near a stationary point {@code G} error {@code ~ (grad)^2 / (2 curv)};
     *       with {@code |grad| <= 1e-6} and CVM curvature {@code O(RT) ~ 1e4}
     *       this is {@code ~1e-16 J/mol} -- negligible. The realistic residual
     *       is the last accepted step, {@code ~1e-6} relative in {@code G}, i.e.
     *       {@code ~0.1 J/mol} on a {@code 1e5 J/mol} energy.</li>
     *   <li><b>Finite grid evaluation.</b> {@code G_single} is evaluated at the
     *       exact combined composition; the {@code G_q} at grid nodes. Both use
     *       the same evaluator, so there is no grid-spacing bias in the
     *       <em>difference</em> -- only each term's own {@code ~0.1 J/mol}
     *       Newton residual.</li>
     *   <li><b>Floating-point subtraction.</b> {@code |G| * 2^-52 ~ 1e5 * 2e-16
     *       ~ 2e-11 J/mol} -- negligible.</li>
     * </ul>
     * <p>So {@code 1e-2 J/mol} would already cover the budget; {@code 1e-1 J/mol}
     * ({@code MERGE_ENERGY_ABS}) is used for margin. It is also floored
     * relatively at {@code MERGE_ENERGY_REL * |G_grid_group|} {@code = 1e-9 *
     * |G|} for the (degenerate) large-{@code |G|} case -- the same relative
     * backward-error scale {@link HillertSolver}'s {@code GAMMA_REL} uses for
     * its {@code mu.M - G} comparison. Any real miscibility-gap energy in the
     * synthetic tests is {@code O(10-1000 J/mol)} -- three-to-five orders above
     * this tolerance -- so a genuine gap is never merged.</p>
     */
    public static final double MERGE_ENERGY_ABS = 1.0e-1;   // J/mol

    /** Relative floor for {@link #MERGE_ENERGY_ABS} at large {@code |G|}. */
    public static final double MERGE_ENERGY_REL = 1.0e-9;

    private HillertGridMinimizer() {
    }

    /**
     * Evaluator for the V2 STEP 10 same-parent merge test: relax ONE parent
     * phase at a given combined composition and report its minimum Gibbs energy
     * and constitution. Phrased on the composition coordinate the merged
     * conserved amount {@code M_group} maps to -- for the currently supported V1
     * disordered phase family {@code M == x} (an established model invariant),
     * so {@code xGroup == M_group}; there is no generalised
     * "evaluate parent at arbitrary M" API and this interface does not fake one.
     *
     * <p>{@link #minimize} supplies a {@link CvmNewtonSolver}-backed
     * implementation over the caller's candidate phases;
     * {@link #minimizeDiscreteForTest(List, double[], int, MergeEvaluator, Consumer)}
     * lets a synthetic test supply an analytic {@code G(x)} so the merge
     * decision can be exercised against closed-form miscibility-gap / convex
     * cases.</p>
     */
    @FunctionalInterface
    public interface MergeEvaluator {

        /**
         * Relax parent phase {@code phaseIndex} at combined composition
         * {@code xGroup} (length K, sums to 1).
         *
         * @return the relaxed parent state, or {@code null} if it could not be
         *         evaluated (non-convergence, invalid state) -- the caller then
         *         does <b>not</b> merge (keeps the group split, the safe choice)
         */
        RelaxedParentState evaluateParentAt(int phaseIndex, double[] xGroup);
    }

    /**
     * Result of {@link MergeEvaluator#evaluateParentAt}: the independently
     * relaxed single-parent state at the combined composition.
     *
     * @param y constitution vector, width {@code ncf+K}
     *          ({@link HillertSolver.Phase#uFull} layout); owned by the caller
     * @param m {@code M_A} at {@code y}, length K
     * @param g absolute {@code G} per formula unit at {@code y}
     */
    public record RelaxedParentState(double[] y, double[] m, double g) {
        public RelaxedParentState {
            y = y.clone();
            m = m.clone();
        }

        @Override
        public double[] y() { return y.clone(); }

        @Override
        public double[] m() { return m.clone(); }
    }

    // =========================================================================
    // Public data types
    // =========================================================================

    /**
     * One fixed-composition pseudo-phase on the grid: parent phase, grid index,
     * composition, component amounts {@code M}, full constitution {@code Y}
     * (width {@code ncf+K}, the {@link HillertSolver.Phase#uFull} layout), and
     * relaxed absolute Gibbs energy {@code G = G0m + Gm} per formula unit.
     *
     * <p>Immutable; every array accessor returns a fresh copy.</p>
     */
    public static final class GridPoint {
        private final int phaseIndex;
        private final String phaseLabel;
        private final int gridIndex;
        private final double[] x;
        private final double[] m;
        private final double[] y;
        private final double g;
        private final boolean converged;

        GridPoint(int phaseIndex, String phaseLabel, int gridIndex,
                double[] x, double[] m, double[] y, double g, boolean converged) {
            this.phaseIndex = phaseIndex;
            this.phaseLabel = phaseLabel;
            this.gridIndex = gridIndex;
            this.x = x.clone();
            this.m = m.clone();
            this.y = y.clone();
            this.g = g;
            this.converged = converged;
        }

        public int phaseIndex()      { return phaseIndex; }
        public String phaseLabel()   { return phaseLabel; }
        public int gridIndex()       { return gridIndex; }
        public double[] composition(){ return x.clone(); }
        public double[] m()          { return m.clone(); }
        public double[] y()          { return y.clone(); }
        public double g()            { return g; }
        public boolean converged()   { return converged; }

        /**
         * Test-only factory for a synthetic grid point with an explicit
         * {@code (M, G)} and no CVM model behind it -- lets the discrete
         * mixture / tangent solver be exercised against analytic
         * {@code G(M)} surfaces (STEP 7 PART 16-17) whose answers are known in
         * closed form. {@code y} is set equal to {@code M} (irrelevant to the
         * mixture math).
         */
        static GridPoint synthetic(int phaseIndex, String label, int gridIndex,
                double[] m, double g) {
            return new GridPoint(phaseIndex, label, gridIndex, m, m, m, g, true);
        }

        @Override
        public String toString() {
            return String.format("GridPoint[%s#g%d x=%s G=%.6f]",
                    phaseLabel, gridIndex, Arrays.toString(round(x)), g);
        }
    }

    /**
     * The starting constitution the grid minimizer assigns to one phase
     * instance -- selected (stable, {@code amount > 0}) or metastable
     * ({@code amount == 0}). A future step turns these into
     * {@link HillertSolver.Phase} objects.
     *
     * @param phaseIndex    index into the input {@code phases} list
     * @param phaseLabel    that phase's label, with a {@code #n} suffix when the
     *                      same parent phase supplies more than one instance
     *                      (miscibility gap)
     * @param instance      1-based instance number within the parent phase
     *                      (always 1 unless a miscibility gap split it)
     * @param amount        formula-unit amount {@code N} (moles); {@code 0} for a
     *                      metastable instance
     * @param y             constitution vector, width {@code ncf+K}
     *                      ({@link HillertSolver.Phase#uFull} layout); owned by
     *                      the caller after {@link #y()}
     * @param m             {@code M_A} at {@code y}, length {@code K}
     * @param g             absolute {@code G} per formula unit at {@code y}
     * @param sourceGridIndex the grid index this constitution came from
     * @param tangentDistance {@code |G - sum_A mu_A M_A|} at {@code y} -- {@code ~0}
     *                        for a selected instance, the minimised distance for
     *                        a metastable one
     */
    public record PhaseInitialState(
            int phaseIndex,
            String phaseLabel,
            int instance,
            double amount,
            double[] y,
            double[] m,
            double g,
            int sourceGridIndex,
            double tangentDistance) {

        public PhaseInitialState {
            y = y.clone();
            m = m.clone();
        }

        @Override
        public double[] y() { return y.clone(); }

        @Override
        public double[] m() { return m.clone(); }

        public boolean selected() { return amount > 0.0; }
    }

    /**
     * The grid minimizer's output. All collections are unmodifiable and all
     * arrays are defensively copied.
     *
     * @param converged           whether a feasible discrete minimum was found
     *                            for the target
     * @param failureReason       {@code null} on success; else a short reason
     *                            (e.g. {@code "infeasible target"},
     *                            {@code "no valid grid points"})
     * @param selectedGridPoints  the &le;K grid points of the discrete minimum,
     *                            in ascending {@code (phaseIndex, gridIndex)}
     *                            order
     * @param selectedPhaseAmounts amounts {@code n_q} aligned with
     *                            {@code selectedGridPoints} (moles;
     *                            {@code sum = 1} for a normalised target)
     * @param stableInitialStates one {@link PhaseInitialState} per selected phase
     *                            instance (miscibility-gap parents contribute
     *                            two), {@code amount > 0}
     * @param metastableInitialStates one {@link PhaseInitialState} per
     *                            non-selected candidate phase, {@code amount == 0},
     *                            constitution = the tangent-closest grid point
     * @param discreteChemicalPotentials {@code mu_A} solved from the selected
     *                            points' common-tangent equations, length
     *                            {@code K}; SER / absolute convention. May be
     *                            {@code null} if underdetermined and no fallback
     *                            applied (see {@code muMethod})
     * @param muMethod            how {@code discreteChemicalPotentials} was
     *                            obtained ({@code "common tangent (K points)"},
     *                            {@code "least-squares (<K points)"},
     *                            {@code "single-phase dG/dx tangent"} ...)
     * @param totalGibbsEnergy    {@code sum_q n_q G_q} at the discrete minimum
     * @param representedAmounts  {@code sum_q n_q M_q,A}, length {@code K} --
     *                            should match the (normalised) target to
     *                            {@link #CONSERVATION_TOL}
     * @param objectiveGap        {@code totalGibbsEnergy} minus the lowest single
     *                            grid-point {@code G} that alone reaches the
     *                            target (or {@code NaN} if none does); negative
     *                            means the mixture genuinely beats every
     *                            one-point solution (a real tie line)
     * @param gridPointCount      total valid grid points across all phases
     * @param skippedGridPointCount points rejected (invalid state / non-finite /
     *                            negative M / non-converged / ordered-unsupported)
     * @param diagnostics         human-readable trace lines
     */
    public record GridMinimizationResult(
            boolean converged,
            String failureReason,
            List<GridPoint> selectedGridPoints,
            double[] selectedPhaseAmounts,
            List<PhaseInitialState> stableInitialStates,
            List<PhaseInitialState> metastableInitialStates,
            double[] discreteChemicalPotentials,
            String muMethod,
            double totalGibbsEnergy,
            double[] representedAmounts,
            double objectiveGap,
            int gridPointCount,
            int skippedGridPointCount,
            List<String> diagnostics) {

        public GridMinimizationResult {
            selectedGridPoints = List.copyOf(selectedGridPoints);
            stableInitialStates = List.copyOf(stableInitialStates);
            metastableInitialStates = List.copyOf(metastableInitialStates);
            diagnostics = List.copyOf(diagnostics);
            selectedPhaseAmounts = selectedPhaseAmounts == null ? null : selectedPhaseAmounts.clone();
            discreteChemicalPotentials =
                    discreteChemicalPotentials == null ? null : discreteChemicalPotentials.clone();
            representedAmounts = representedAmounts == null ? null : representedAmounts.clone();
        }

        @Override
        public double[] selectedPhaseAmounts() {
            return selectedPhaseAmounts == null ? null : selectedPhaseAmounts.clone();
        }

        @Override
        public double[] discreteChemicalPotentials() {
            return discreteChemicalPotentials == null ? null : discreteChemicalPotentials.clone();
        }

        @Override
        public double[] representedAmounts() {
            return representedAmounts == null ? null : representedAmounts.clone();
        }

        /** All initial states (stable then metastable), one per phase instance. */
        public List<PhaseInitialState> allInitialStates() {
            List<PhaseInitialState> all = new ArrayList<>(stableInitialStates);
            all.addAll(metastableInitialStates);
            return Collections.unmodifiableList(all);
        }
    }

    // =========================================================================
    // Entry point
    // =========================================================================

    /**
     * Run the grid minimizer.
     *
     * @param phases          the candidate phases (order defines
     *                        {@code phaseIndex}); not mutated
     * @param temperature     T in kelvin
     * @param overallAmounts  overall component amounts {@code b_A}, length
     *                        {@code K}; internally normalised to sum 1 (the
     *                        result's amounts and represented composition are on
     *                        that normalised scale)
     * @param progress        optional line sink for diagnostics; may be
     *                        {@code null}
     * @return the discrete minimum plus a starting constitution for every phase
     * @throws IllegalArgumentException if {@code K >= 4}, K is inconsistent
     *         across phases, or {@code overallAmounts} is the wrong length / not
     *         strictly positive
     */
    public static GridMinimizationResult minimize(
            List<HillertSolver.Phase> phases,
            double temperature,
            double[] overallAmounts,
            Consumer<String> progress) {

        if (phases == null || phases.isEmpty()) {
            throw new IllegalArgumentException("grid minimizer: no phases");
        }
        int k = phases.get(0).numComponents;
        for (HillertSolver.Phase p : phases) {
            if (p.numComponents != k) {
                throw new IllegalArgumentException(
                        "grid minimizer: mixed component counts (" + p.label + " has "
                                + p.numComponents + ", expected " + k + ")");
            }
        }
        if (k < 2 || k > 3) {
            throw new IllegalArgumentException(
                    "grid minimizer: only K=2 and K=3 are implemented (got K=" + k + ")");
        }
        if (overallAmounts == null || overallAmounts.length != k) {
            throw new IllegalArgumentException(
                    "grid minimizer: overallAmounts length " + safeLen(overallAmounts)
                            + " != K=" + k);
        }
        double bSum = 0.0;
        for (double b : overallAmounts) {
            if (!(b > 0.0) || !Double.isFinite(b)) {
                throw new IllegalArgumentException(
                        "grid minimizer: overallAmounts must be strictly positive and finite, got "
                                + Arrays.toString(overallAmounts));
            }
            bSum += b;
        }
        final double[] target = new double[k];
        for (int a = 0; a < k; a++) {
            target[a] = overallAmounts[a] / bSum;   // normalised: sum_A target = 1
        }

        List<String> diag = new ArrayList<>();
        emit(progress, diag, String.format(
                "grid minimizer: K=%d  T=%.2f K  target(normalised)=%s  N=%d",
                k, temperature, Arrays.toString(round(target)), GRID_DIVISIONS));

        // ---- 1. build + evaluate the grid for every phase ----
        List<GridPoint> grid = new ArrayList<>();
        int skipped = 0;
        for (int pi = 0; pi < phases.size(); pi++) {
            HillertSolver.Phase phase = phases.get(pi);
            int tcf = phase.model.geometry().tcf;
            int expectedDisorderedWidth = phase.ncf + k;
            if (tcf != expectedDisorderedWidth) {
                skipped += countGridNodes(k);
                emit(progress, diag, String.format(
                        "  phase '%s' SKIPPED (ordered: tcf=%d > ncf+K=%d; not supported by the "
                                + "grid minimizer in STEP 7)", phase.label, tcf, expectedDisorderedWidth));
                continue;
            }
            int before = grid.size();
            int localSkipped = evaluatePhaseGrid(pi, phase, temperature, k, grid, diag, progress);
            skipped += localSkipped;
            emit(progress, diag, String.format(
                    "  phase '%s': %d valid grid points, %d skipped",
                    phase.label, grid.size() - before, localSkipped));
        }

        if (grid.isEmpty()) {
            emit(progress, diag, "grid minimizer: no valid grid points -> FAIL");
            return failure("no valid grid points", 0, skipped, diag);
        }

        // Real merge evaluator: relax the parent CVM phase at a combined
        // composition via CvmNewtonSolver (V1: x_group == M_group).
        MergeEvaluator mergeEval = (phaseIndex, xGroup) ->
                relaxParentAt(phases.get(phaseIndex), temperature, xGroup, k);

        return runDiscrete(grid, target, k, phases, temperature, mergeEval, skipped, diag, progress);
    }

    /**
     * Relax parent phase at a fixed combined composition via the existing
     * fixed-composition {@link CvmNewtonSolver} -- the STEP 10 merge test's
     * "evaluate ONE parent phase at the combined composition" operation. Returns
     * {@code null} on non-convergence / invalid state (caller then does not
     * merge).
     */
    private static RelaxedParentState relaxParentAt(
            HillertSolver.Phase phase, double temperature, double[] xGroup, int k) {
        try {
            CVMGibbsModel.State state;
            if (isPureElement(xGroup)) {
                double[] u = phase.model.randomStateU(xGroup);
                state = phase.model.at(temperature, xGroup, u);
            } else {
                CvmNewtonSolver.Result r = new CvmNewtonSolver(phase.model)
                        .solve(temperature, xGroup, GRID_POINT_SOLVE_TOL, null, null);
                if (!r.converged()) {
                    return null;
                }
                state = r.state();
            }
            if (state == null || !state.isValidIncludingPoints()) {
                return null;
            }
            double g = state.g();
            double[] m = state.componentAmountsPerFormulaUnit();
            double[] y = state.cfs();
            if (!Double.isFinite(g) || !allFinite(m) || !allFinite(y)
                    || y.length != phase.ncf + k) {
                return null;
            }
            return new RelaxedParentState(y, m, g);
        } catch (RuntimeException ex) {
            return null;
        }
    }

    /**
     * The discrete stage: convex-hull mixture solve, discrete {@code mu},
     * stable + metastable initial states. Split out from {@link #minimize} so a
     * test can drive it directly against a synthetic grid (STEP 7 PART 16-17)
     * with a known analytic {@code G(M)}.
     *
     * @param phases used only for the single-phase {@code dG/dx} {@code mu}
     *               fallback; may be empty for a synthetic grid whose minimum
     *               selects {@code >= K} distinct points
     * @param mergeEval the V2 STEP 10 same-parent merge-test evaluator; may be
     *               {@code null} to skip the merge test entirely (a synthetic
     *               grid with no analytic {@code G} supplied -- then same-parent
     *               groups keep the pre-STEP-10 {@link #GAP_MERGE_TOL}-only
     *               behaviour)
     */
    static GridMinimizationResult runDiscrete(
            List<GridPoint> grid, double[] target, int k,
            List<HillertSolver.Phase> phases, double temperature,
            MergeEvaluator mergeEval,
            int skipped, List<String> diag, Consumer<String> progress) {

        // ---- 2. discrete convex-hull / linear-mixture solve ----
        MixtureSolution mix = findMinimumMixture(grid, target, k, diag, progress);
        if (mix == null) {
            emit(progress, diag, "grid minimizer: target not in the convex hull of any "
                    + "feasible subset -> infeasible");
            return failure("infeasible target", grid.size(), skipped, diag);
        }

        // ---- 3. group selected points by parent phase (miscibility gap) ----
        List<GridPoint> selPts = mix.points;
        double[] selAmt = mix.amounts;
        emit(progress, diag, String.format(
                "grid minimizer: discrete minimum uses %d grid point(s), total G = %.8f",
                selPts.size(), mix.totalG));
        for (int i = 0; i < selPts.size(); i++) {
            emit(progress, diag, String.format(
                    "    selected %s  n=%.8f", selPts.get(i), selAmt[i]));
        }

        // ---- 4. discrete chemical potentials from the selected common tangent ----
        MuSolution muSol = computeDiscreteMu(selPts, k, phases, temperature, diag, progress);

        // ---- 5. build stable + metastable initial states ----
        // (STEP 10) same-parent selected groups get a direct thermodynamic
        // merge test inside buildStableStates: relax ONE parent at the combined
        // composition and compare G_single vs the grid mixture's G.
        List<PhaseInitialState> stable = buildStableStates(
                selPts, selAmt, muSol.mu, k, mergeEval, diag, progress);
        // A parent phase merged into a single stable instance must NOT also be
        // emitted as metastable (STEP 10 PART 17): assignMetastableInitialStates
        // already filters on "parent not in the selected set", and a merged
        // parent IS in the selected set, so this is automatic -- but pass the
        // built stable list so the filter uses the post-merge parent set.
        List<PhaseInitialState> meta = assignMetastableInitialStates(
                phases, grid, selPts, muSol.mu, k, diag, progress);

        // ---- 6. represented composition + objective gap ----
        double[] represented = new double[k];
        for (int i = 0; i < selPts.size(); i++) {
            double[] m = selPts.get(i).m;
            for (int a = 0; a < k; a++) {
                represented[a] += selAmt[i] * m[a];
            }
        }
        double singleBest = bestSinglePointObjective(grid, target, k);
        double objectiveGap = Double.isNaN(singleBest) ? Double.NaN : mix.totalG - singleBest;
        emit(progress, diag, String.format(
                "grid minimizer: represented composition = %s  (target %s)",
                Arrays.toString(round(represented)), Arrays.toString(round(target))));
        emit(progress, diag, String.format(
                "grid minimizer: objective gap vs best single grid point = %s",
                Double.isNaN(objectiveGap) ? "n/a (no single point reaches target)"
                        : String.format("%.8f", objectiveGap)));

        boolean gap = hasMiscibilityGap(stable);
        if (gap) {
            emit(progress, diag, "grid minimizer: MISCIBILITY GAP detected "
                    + "(one parent phase supplies >1 distinct-composition instance)");
        }

        return new GridMinimizationResult(
                true, null,
                selPts, selAmt,
                stable, meta,
                muSol.mu, muSol.method,
                mix.totalG, represented, objectiveGap,
                grid.size(), skipped,
                diag);
    }

    /**
     * Test-only: run the discrete stage against a caller-supplied grid and
     * target (already normalised so {@code sum target == 1}). No CVM model, no
     * grid construction -- for exercising the mixture / tangent math on
     * analytic {@code G(M)} surfaces. No merge evaluator -> same-parent groups
     * keep the pre-STEP-10 {@link #GAP_MERGE_TOL}-only behaviour.
     */
    static GridMinimizationResult minimizeDiscreteForTest(
            List<GridPoint> grid, double[] target, int k, Consumer<String> progress) {
        List<String> diag = new ArrayList<>();
        return runDiscrete(grid, target, k, List.of(), Double.NaN, null, 0, diag, progress);
    }

    /**
     * Test-only: like {@link #minimizeDiscreteForTest(List, double[], int, Consumer)}
     * but with a caller-supplied {@link MergeEvaluator} so the STEP-10
     * same-parent merge decision can be exercised against a closed-form
     * {@code G(x)} (convex single phase vs. double-well miscibility gap).
     */
    static GridMinimizationResult minimizeDiscreteForTest(
            List<GridPoint> grid, double[] target, int k,
            MergeEvaluator mergeEval, Consumer<String> progress) {
        List<String> diag = new ArrayList<>();
        return runDiscrete(grid, target, k, List.of(), Double.NaN, mergeEval, 0, diag, progress);
    }

    // =========================================================================
    // 1. Grid construction + evaluation
    // =========================================================================

    /**
     * Composition grid nodes for {@code K} components at density
     * {@link #GRID_DIVISIONS}. K=2: {@code {(i/N, 1-i/N)}}. K=3: the barycentric
     * simplex {@code {(i/N, j/N, 1-(i+j)/N) : i,j>=0, i+j<=N}}. Deterministic
     * ascending order.
     */
    static List<double[]> compositionGrid(int k) {
        int n = GRID_DIVISIONS;
        List<double[]> nodes = new ArrayList<>();
        if (k == 2) {
            for (int i = 0; i <= n; i++) {
                double xa = (double) i / n;
                nodes.add(new double[] { xa, 1.0 - xa });
            }
        } else { // k == 3
            for (int i = 0; i <= n; i++) {
                for (int j = 0; j <= n - i; j++) {
                    double xa = (double) i / n;
                    double xb = (double) j / n;
                    nodes.add(new double[] { xa, xb, 1.0 - xa - xb });
                }
            }
        }
        return nodes;
    }

    static int countGridNodes(int k) {
        int n = GRID_DIVISIONS;
        return k == 2 ? (n + 1) : (n + 1) * (n + 2) / 2;
    }

    /**
     * Evaluate one phase's grid. Each interior node's internal CFs are relaxed
     * by {@link CvmNewtonSolver} at fixed composition; pure-element nodes
     * ({@code x_A = 1}) are evaluated directly (the CVM solve is trivial /
     * singular there). Invalid points are skipped and counted (STEP 7 PART 6);
     * no state is repaired.
     *
     * @return number of skipped grid points for this phase
     */
    private static int evaluatePhaseGrid(
            int phaseIndex, HillertSolver.Phase phase, double temperature, int k,
            List<GridPoint> out, List<String> diag, Consumer<String> progress) {

        CvmNewtonSolver solver = new CvmNewtonSolver(phase.model);
        int ncf = phase.ncf;
        List<double[]> nodes = compositionGrid(k);
        int skipped = 0;

        for (int gi = 0; gi < nodes.size(); gi++) {
            double[] x = nodes.get(gi);

            // clamp tiny negative round-off on the dependent coordinate
            for (int a = 0; a < k; a++) {
                if (x[a] < 0 && x[a] > -1e-12) {
                    x[a] = 0.0;
                }
            }

            CVMGibbsModel.State state;
            boolean converged;
            try {
                if (isPureElement(x)) {
                    // internal CFs of a pure element are the disordered ones; no
                    // minimisation to do (and the fixed-x Hessian is singular).
                    double[] u = phase.model.randomStateU(x);
                    state = phase.model.at(temperature, x, u);
                    converged = true;
                } else {
                    CvmNewtonSolver.Result r = solver.solve(
                            temperature, x, GRID_POINT_SOLVE_TOL, null, null);
                    state = r.state();
                    converged = r.converged();
                }
            } catch (RuntimeException ex) {
                skipped++;
                continue;
            }

            if (state == null || !converged || !state.isValidIncludingPoints()) {
                skipped++;
                continue;
            }
            double g = state.g();
            double[] m = state.componentAmountsPerFormulaUnit();
            if (!Double.isFinite(g) || !allFinite(m) || anyNegative(m, 1e-9)) {
                skipped++;
                continue;
            }
            // Y in the Phase.uFull layout [u ; x] (disordered: cfs() width == ncf+K)
            double[] y = state.cfs();
            if (y.length != ncf + k || !allFinite(y)) {
                skipped++;
                continue;
            }
            out.add(new GridPoint(phaseIndex, phase.label, gi, x, m, y, g, true));
        }
        return skipped;
    }

    // =========================================================================
    // 2. Discrete linear-mixture solve (vertex enumeration)
    // =========================================================================

    private static final class MixtureSolution {
        final List<GridPoint> points;
        final double[] amounts;
        final double totalG;
        final int[] sortedGridKey;   // for tie-break: sorted (phaseIndex,gridIndex) flattened

        MixtureSolution(List<GridPoint> points, double[] amounts, double totalG, int[] key) {
            this.points = points;
            this.amounts = amounts;
            this.totalG = totalG;
            this.sortedGridKey = key;
        }
    }

    /**
     * Exact vertex enumeration of the discrete mixture problem. Considers every
     * 1-point subset, every 2-point subset (K=2 and K=3), and every 3-point
     * subset (K=3). Each subset is solved as a square system
     * {@code [M_q ... ; 1 ... 1] n = [target ; 1]}; rejected if near-singular
     * ({@link #RCOND_FLOOR}), if any {@code n_q < -AMOUNT_TOL}, or if
     * conservation fails ({@link #CONSERVATION_TOL}). Lowest {@code sum n_q G_q}
     * wins; ties (within {@link #OBJECTIVE_TIE_TOL}) go to the
     * lexicographically-smallest sorted grid-key.
     */
    static MixtureSolution findMinimumMixture(
            List<GridPoint> grid, double[] target, int k,
            List<String> diag, Consumer<String> progress) {

        int q = grid.size();
        MixtureSolution best = null;
        long subsetsTried = 0;
        long subsetsFeasible = 0;

        // ---- 1-point ----
        for (int i = 0; i < q; i++) {
            subsetsTried++;
            MixtureSolution s = solveSubset(grid, new int[] { i }, target, k);
            if (s != null) {
                subsetsFeasible++;
                best = pick(best, s);
            }
        }
        // ---- 2-point ----
        for (int i = 0; i < q; i++) {
            for (int j = i + 1; j < q; j++) {
                subsetsTried++;
                MixtureSolution s = solveSubset(grid, new int[] { i, j }, target, k);
                if (s != null) {
                    subsetsFeasible++;
                    best = pick(best, s);
                }
            }
        }
        // ---- 3-point (K=3 only) ----
        if (k == 3) {
            for (int i = 0; i < q; i++) {
                for (int j = i + 1; j < q; j++) {
                    // cheap prune: need the target inside the (i,j,*) hull; skip
                    // only when i and j share the *same* composition (degenerate edge)
                    for (int l = j + 1; l < q; l++) {
                        subsetsTried++;
                        MixtureSolution s = solveSubset(grid, new int[] { i, j, l }, target, k);
                        if (s != null) {
                            subsetsFeasible++;
                            best = pick(best, s);
                        }
                    }
                }
            }
        }

        emit(progress, diag, String.format(
                "grid minimizer: enumerated %d subsets, %d feasible", subsetsTried, subsetsFeasible));
        return best;
    }

    /**
     * Solve one subset {@code idx} (size 1..K) for its grid-point amounts.
     * Square system of size {@code m+1} where {@code m = idx.length}:
     * {@code K} conservation rows (one per component) plus the
     * {@code sum n = 1} row -- overdetermined when {@code m < K}, so a
     * least-squares (normal-equation) solve is used there; exactly determined
     * when {@code m == K}. Returns {@code null} if infeasible / singular /
     * non-conservative.
     */
    private static MixtureSolution solveSubset(
            List<GridPoint> grid, int[] idx, double[] target, int k) {

        int m = idx.length;
        // rows: K conservation + 1 normalisation ; cols: m amounts
        double[][] a = new double[k + 1][m];
        double[] rhs = new double[k + 1];
        for (int c = 0; c < m; c++) {
            double[] mq = grid.get(idx[c]).m;
            for (int r = 0; r < k; r++) {
                a[r][c] = mq[r];
            }
            a[k][c] = 1.0;
        }
        for (int r = 0; r < k; r++) {
            rhs[r] = target[r];
        }
        rhs[k] = 1.0;

        double[] n;
        try {
            if (m == k + 1) {
                // square (K+1)x(K+1)
                LinearAlgebra.Solution sol = LinearAlgebra.solveChecked(a, rhs);
                if (sol.relativeResidual() > MAX_SOLVE_RESIDUAL) {
                    return null;
                }
                n = sol.x();
            } else {
                // overdetermined (m < K+1): normal equations AtA n = At rhs
                double[][] ata = new double[m][m];
                double[] atb = new double[m];
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < m; c++) {
                        double s = 0.0;
                        for (int t = 0; t < k + 1; t++) {
                            s += a[t][r] * a[t][c];
                        }
                        ata[r][c] = s;
                    }
                    double s = 0.0;
                    for (int t = 0; t < k + 1; t++) {
                        s += a[t][r] * rhs[t];
                    }
                    atb[r] = s;
                }
                LinearAlgebra.Solution sol = LinearAlgebra.solveChecked(ata, atb);
                if (sol.relativeResidual() > MAX_SOLVE_RESIDUAL) {
                    return null;
                }
                n = sol.x();
            }
        } catch (RuntimeException singular) {
            return null;   // exactly singular subset -- degenerate, reject
        }

        // feasibility: non-negative amounts
        for (double v : n) {
            if (!Double.isFinite(v) || v < -AMOUNT_TOL) {
                return null;
            }
        }
        // clamp tiny negatives to 0, drop exact-zero points from the reported set
        List<GridPoint> pts = new ArrayList<>();
        List<Double> amts = new ArrayList<>();
        for (int c = 0; c < m; c++) {
            double v = Math.max(0.0, n[c]);
            if (v <= AMOUNT_TOL) {
                continue;   // this vertex is not actually used
            }
            pts.add(grid.get(idx[c]));
            amts.add(v);
        }
        if (pts.isEmpty()) {
            return null;
        }
        // conservation check on the *reduced* set
        double[] rep = new double[k];
        double nSum = 0.0;
        for (int c = 0; c < pts.size(); c++) {
            double v = amts.get(c);
            nSum += v;
            double[] mq = pts.get(c).m;
            for (int r = 0; r < k; r++) {
                rep[r] += v * mq[r];
            }
        }
        for (int r = 0; r < k; r++) {
            if (Math.abs(rep[r] - target[r]) > CONSERVATION_TOL) {
                return null;
            }
        }
        if (Math.abs(nSum - 1.0) > CONSERVATION_TOL) {
            return null;
        }

        // sort the reduced set by (phaseIndex, gridIndex) for a stable report + tie-key
        Integer[] order = new Integer[pts.size()];
        for (int c = 0; c < order.length; c++) {
            order[c] = c;
        }
        final List<GridPoint> fp = pts;
        Arrays.sort(order, (o1, o2) -> {
            int c = Integer.compare(fp.get(o1).phaseIndex, fp.get(o2).phaseIndex);
            if (c != 0) return c;
            return Integer.compare(fp.get(o1).gridIndex, fp.get(o2).gridIndex);
        });
        List<GridPoint> sortedPts = new ArrayList<>();
        double[] sortedAmt = new double[order.length];
        int[] key = new int[order.length * 2];
        double totalG = 0.0;
        for (int c = 0; c < order.length; c++) {
            GridPoint gp = fp.get(order[c]);
            double v = amts.get(order[c]);
            sortedPts.add(gp);
            sortedAmt[c] = v;
            key[2 * c] = gp.phaseIndex;
            key[2 * c + 1] = gp.gridIndex;
            totalG += v * gp.g;
        }
        return new MixtureSolution(sortedPts, sortedAmt, totalG, key);
    }

    /** Lower objective wins; tie within {@link #OBJECTIVE_TIE_TOL} -> smaller sorted grid key. */
    private static MixtureSolution pick(MixtureSolution best, MixtureSolution cand) {
        if (best == null) {
            return cand;
        }
        double d = cand.totalG - best.totalG;
        if (d < -OBJECTIVE_TIE_TOL) {
            return cand;
        }
        if (d > OBJECTIVE_TIE_TOL) {
            return best;
        }
        return lexLess(cand.sortedGridKey, best.sortedGridKey) ? cand : best;
    }

    private static boolean lexLess(int[] a, int[] b) {
        int n = Math.min(a.length, b.length);
        for (int i = 0; i < n; i++) {
            if (a[i] != b[i]) {
                return a[i] < b[i];
            }
        }
        return a.length < b.length;
    }

    /** Lowest G among single grid points that alone satisfy the target; NaN if none. */
    private static double bestSinglePointObjective(List<GridPoint> grid, double[] target, int k) {
        double best = Double.NaN;
        for (GridPoint gp : grid) {
            boolean reaches = true;
            for (int a = 0; a < k; a++) {
                if (Math.abs(gp.m[a] - target[a]) > CONSERVATION_TOL) {
                    reaches = false;
                    break;
                }
            }
            if (reaches && (Double.isNaN(best) || gp.g < best)) {
                best = gp.g;
            }
        }
        return best;
    }

    // =========================================================================
    // 3-4. Discrete chemical potentials
    // =========================================================================

    private static final class MuSolution {
        final double[] mu;      // length K, or null
        final String method;

        MuSolution(double[] mu, String method) {
            this.mu = mu;
            this.method = method;
        }
    }

    /**
     * Solve {@code sum_A M_q,A mu_A = G_q} over the selected grid points.
     * <ul>
     *   <li>{@code #distinct-M points == K}: square solve -> the common tangent
     *       hyperplane through those K points.</li>
     *   <li>{@code < K} distinct points (single-phase target, or a degenerate
     *       tie): least-squares through the available equation(s); for the
     *       common single-phase case ({@code 1} point) this reduces to a
     *       one-equation, K-unknown system and we fall back to the phase's own
     *       {@code dG/dx} tangent (documented) so {@code mu} is not arbitrary.</li>
     * </ul>
     */
    static MuSolution computeDiscreteMu(
            List<GridPoint> selPts, int k, List<HillertSolver.Phase> phases,
            double temperature, List<String> diag, Consumer<String> progress) {

        // distinct-composition selected points
        List<GridPoint> distinct = new ArrayList<>();
        for (GridPoint gp : selPts) {
            boolean dup = false;
            for (GridPoint d : distinct) {
                if (infNorm(diff(d.m, gp.m)) <= GAP_MERGE_TOL) {
                    dup = true;
                    break;
                }
            }
            if (!dup) {
                distinct.add(gp);
            }
        }

        if (distinct.size() >= k) {
            // take the first K distinct points, square solve
            double[][] a = new double[k][k];
            double[] rhs = new double[k];
            for (int r = 0; r < k; r++) {
                double[] m = distinct.get(r).m;
                for (int c = 0; c < k; c++) {
                    a[r][c] = m[c];
                }
                rhs[r] = distinct.get(r).g;
            }
            try {
                LinearAlgebra.Solution sol = LinearAlgebra.solveChecked(a, rhs);
                if (sol.relativeResidual() <= MAX_SOLVE_RESIDUAL) {
                    emit(progress, diag, String.format(
                            "grid minimizer: discrete mu = %s  (common tangent, %d points)",
                            Arrays.toString(round(sol.x())), k));
                    return new MuSolution(sol.x(), "common tangent (K points)");
                }
                emit(progress, diag, "grid minimizer: K-point tangent solve ill-conditioned "
                        + "(res=" + sol.relativeResidual() + "); falling back to single-phase dG/dx");
            } catch (RuntimeException singular) {
                emit(progress, diag, "grid minimizer: K-point tangent solve was singular; "
                        + "falling back to single-phase dG/dx");
            }
        }

        // ---- fewer than K distinct points: single-phase dG/dx tangent ----
        // mu_A = G + dG/dx_A - sum_B x_B dG/dx_B     (standard CALPHAD tangent
        // intercepts for a disordered phase; M == x for V1). Uses the selected
        // point's own parent phase and its relaxed state's gmu. Needs the real
        // Phase (and its CVM model) -- unavailable on the synthetic-grid test
        // path, where phases is empty; report mu as unavailable there.
        GridPoint anchor = selPts.get(0);
        if (anchor.phaseIndex >= phases.size()) {
            emit(progress, diag, "grid minimizer: <K distinct selected points and no CVM model "
                    + "available for the dG/dx fallback -> discrete mu unavailable (underdetermined)");
            return new MuSolution(null, "unavailable (<K points, no model)");
        }
        HillertSolver.Phase phase = phases.get(anchor.phaseIndex);
        try {
            double[] u = Arrays.copyOfRange(anchor.y, 0, phase.ncf);
            double[] x = anchor.composition();
            CVMGibbsModel.State st = phase.model.at(temperature, x, u);
            double g = st.g();
            // For a disordered V1 phase the composition derivatives we need are
            // the point-block gradient of the FULL absolute G wrt x. gmuFull is
            // dG/dY over [u ; x]; its trailing-K slice is dG/dx (the internal
            // block is stationary at the relaxed grid-point state, so this is
            // the total composition derivative).
            double[] gmuFull = st.gmuFull();
            int ncf = phase.ncf;
            double[] dGdx = Arrays.copyOfRange(gmuFull, ncf, ncf + k);
            double xdot = 0.0;
            for (int a = 0; a < k; a++) {
                xdot += x[a] * dGdx[a];
            }
            double[] mu = new double[k];
            for (int a = 0; a < k; a++) {
                mu[a] = g + dGdx[a] - xdot;
            }
            emit(progress, diag, String.format(
                    "grid minimizer: discrete mu = %s  (single-phase dG/dx tangent at %s)",
                    Arrays.toString(round(mu)), anchor.phaseLabel));
            return new MuSolution(mu, "single-phase dG/dx tangent");
        } catch (RuntimeException ex) {
            emit(progress, diag, "grid minimizer: single-phase mu fallback failed (" + ex + "); "
                    + "mu left null");
            return new MuSolution(null, "unavailable (underdetermined)");
        }
    }

    // =========================================================================
    // 5. Stable + metastable initial states
    // =========================================================================

    private static List<PhaseInitialState> buildStableStates(
            List<GridPoint> selPts, double[] selAmt, double[] mu, int k,
            MergeEvaluator mergeEval, List<String> diag, Consumer<String> progress) {

        // group by parent phase, split into distinct-composition instances
        List<PhaseInitialState> out = new ArrayList<>();
        // stable sort by (phaseIndex, gridIndex) already holds from findMinimumMixture
        int i = 0;
        while (i < selPts.size()) {
            int pi = selPts.get(i).phaseIndex;
            int j = i;
            List<Integer> group = new ArrayList<>();
            while (j < selPts.size() && selPts.get(j).phaseIndex == pi) {
                group.add(j);
                j++;
            }

            // ---- (a) collapse essentially-identical points within the group ----
            List<List<Integer>> clusters = new ArrayList<>();
            for (int gIdx : group) {
                boolean placed = false;
                for (List<Integer> cl : clusters) {
                    if (infNorm(diff(selPts.get(cl.get(0)).m, selPts.get(gIdx).m)) <= GAP_MERGE_TOL) {
                        cl.add(gIdx);
                        placed = true;
                        break;
                    }
                }
                if (!placed) {
                    List<Integer> cl = new ArrayList<>();
                    cl.add(gIdx);
                    clusters.add(cl);
                }
            }

            // ---- (b) STEP 10: same-parent thermodynamic merge test ----
            // If more than one distinct-composition cluster remains, ask
            // Sundman's 2015 sec.5.1 question directly: can this same-parent grid
            // mixture be represented by ONE equilibrium state of that parent at
            // the combined composition without increasing Gibbs energy?
            if (clusters.size() > 1 && mergeEval != null) {
                MergeOutcome mo = mergeSameParentGroup(
                        selPts, selAmt, group, pi, k, mergeEval, diag, progress);
                if (mo != null && mo.merged) {
                    double td = (mu != null) ? Math.abs(mo.g - dot(mu, mo.m, k)) : Double.NaN;
                    out.add(new PhaseInitialState(
                            pi, selPts.get(group.get(0)).phaseLabel, 1, mo.amount,
                            mo.y, mo.m, mo.g,
                            selPts.get(group.get(0)).gridIndex, td));
                    i = j;
                    continue;
                }
            }

            // ---- (c) not merged: emit one instance per distinct cluster ----
            boolean split = clusters.size() > 1;
            for (int c = 0; c < clusters.size(); c++) {
                List<Integer> cl = clusters.get(c);
                // representative = lowest-G member; amount = sum
                int rep = cl.get(0);
                double amt = 0.0;
                for (int gIdx : cl) {
                    amt += selAmt[gIdx];
                    if (selPts.get(gIdx).g < selPts.get(rep).g) {
                        rep = gIdx;
                    }
                }
                GridPoint gp = selPts.get(rep);
                String label = split ? gp.phaseLabel + "#" + (c + 1) : gp.phaseLabel;
                double td = tangentDistance(gp, mu, k);
                out.add(new PhaseInitialState(
                        gp.phaseIndex, label, c + 1, amt,
                        gp.y, gp.m, gp.g, gp.gridIndex, td));
            }
            i = j;
        }
        return out;
    }

    /** Outcome of {@link #mergeSameParentGroup}. */
    private static final class MergeOutcome {
        final boolean merged;
        final double amount;   // sum_q n_q
        final double[] y;      // Y_single (the independently relaxed parent state)
        final double[] m;      // M_single
        final double g;        // G_single

        MergeOutcome(boolean merged, double amount, double[] y, double[] m, double g) {
            this.merged = merged;
            this.amount = amount;
            this.y = y;
            this.m = m;
            this.g = g;
        }

        static MergeOutcome keepSplit() {
            return new MergeOutcome(false, 0.0, null, null, Double.NaN);
        }
    }

    /**
     * V2 STEP 10 same-parent merge test (Sundman 2015 sec.5.1). For the whole
     * same-parent selected {@code group}:
     *
     * <pre>
     *   M_group      = sum_q n_q * M_q          (combined conserved composition)
     *   G_grid_group = sum_q n_q * G_q          (the grid mixture's Gibbs energy)
     *   G_single     = min Gibbs of ONE parent phase at M_group   (relaxed, real)
     *   deltaG_merge = G_single - G_grid_group
     * </pre>
     *
     * <ul>
     *   <li>{@code deltaG_merge <= tol}: the single parent at the combined
     *       composition is no more expensive than the grid mixture -- the
     *       selected points are discretisation of ONE continuous phase.
     *       <b>MERGE</b>: return one state with {@code amount = sum n_q},
     *       {@code Y = Y_single} (the independently relaxed constitution, NOT an
     *       averaged / lowest-G-node Y), {@code G = G_single}.</li>
     *   <li>{@code deltaG_merge > tol}: the mixture of distinct compositions is
     *       genuinely cheaper -- a real miscibility gap. <b>KEEP SPLIT</b>.</li>
     * </ul>
     *
     * <p>{@code tol = max(MERGE_ENERGY_ABS, MERGE_ENERGY_REL * |G_grid_group|)}
     * -- see {@link #MERGE_ENERGY_ABS}. If {@code mergeEval} cannot evaluate the
     * parent at {@code M_group} (returns {@code null}), the group is NOT merged
     * (the conservative choice).</p>
     *
     * <p><b>Per-mole normalisation.</b> {@code G_single} from the evaluator is
     * per formula unit (one mole of the phase). {@code sum_q n_q G_q} is the
     * group's <em>total</em> Gibbs energy, where {@code sum_q n_q = amountSum}
     * is the group's total amount -- <b>not necessarily 1</b> when the group is
     * one parent among several selected phases. Both sides of {@code deltaG_merge}
     * are therefore put on the same per-mole basis:
     * {@code G_grid_group_perMole = (sum_q n_q G_q) / amountSum}. For the current
     * single-candidate-phase grid problems {@code amountSum == 1} and this is a
     * no-op; it matters only if a future multi-phase grid problem gives one
     * parent a sub-unit share.</p>
     */
    private static MergeOutcome mergeSameParentGroup(
            List<GridPoint> selPts, double[] selAmt, List<Integer> group,
            int phaseIndex, int k, MergeEvaluator mergeEval,
            List<String> diag, Consumer<String> progress) {

        double amountSum = 0.0;
        double[] mGroup = new double[k];
        double gGridGroupTotal = 0.0;
        for (int gIdx : group) {
            double n = selAmt[gIdx];
            amountSum += n;
            double[] mq = selPts.get(gIdx).m;
            for (int a = 0; a < k; a++) {
                mGroup[a] += n * mq[a];
            }
            gGridGroupTotal += n * selPts.get(gIdx).g;
        }
        // per-mole grid-mixture Gibbs energy, comparable to the per-formula-unit
        // G_single the evaluator returns (amountSum == 1 for a single-candidate
        // grid problem, so this is a no-op there).
        double gGridGroup = (amountSum > 0.0) ? gGridGroupTotal / amountSum : gGridGroupTotal;
        // normalise the combined composition (n_q sum to the group's amount, not
        // necessarily 1 -- a same-parent group is usually a subset of the whole
        // selected set); M_group as a composition vector must sum to 1.
        double mSum = 0.0;
        for (double v : mGroup) {
            mSum += v;
        }
        double[] xGroup = new double[k];
        for (int a = 0; a < k; a++) {
            xGroup[a] = mGroup[a] / mSum;
        }
        // independent check: normalised M_group sums to 1
        double check = 0.0;
        for (double v : xGroup) {
            check += v;
        }
        if (Math.abs(check - 1.0) > CONSERVATION_TOL) {
            emit(progress, diag, String.format(
                    "grid minimizer: merge test for parent %d -- combined composition does not "
                            + "normalise (sum=%.9f); keeping split", phaseIndex, check));
            return MergeOutcome.keepSplit();
        }

        RelaxedParentState single = mergeEval.evaluateParentAt(phaseIndex, xGroup);
        if (single == null) {
            emit(progress, diag, String.format(
                    "grid minimizer: merge test for parent %d -- single-phase relax at combined "
                            + "composition %s did not converge; keeping split",
                    phaseIndex, Arrays.toString(round(xGroup))));
            return MergeOutcome.keepSplit();
        }

        double deltaGMerge = single.g() - gGridGroup;
        double tol = Math.max(MERGE_ENERGY_ABS, MERGE_ENERGY_REL * Math.abs(gGridGroup));

        emit(progress, diag, String.format(
                "grid minimizer: same-parent merge test (parent %d, %d selected points): "
                        + "M_group=%s  G_grid_group=%.8f  G_single=%.8f  deltaG_merge=%.6f  tol=%.6f",
                phaseIndex, group.size(), Arrays.toString(round(xGroup)),
                gGridGroup, single.g(), deltaGMerge, tol));

        if (deltaGMerge <= tol) {
            emit(progress, diag, String.format(
                    "grid minimizer:   -> MERGE parent %d into ONE instance "
                            + "(single-phase state at the combined composition is no more expensive)",
                    phaseIndex));
            return new MergeOutcome(true, amountSum, single.y(), single.m(), single.g());
        }
        emit(progress, diag, String.format(
                "grid minimizer:   -> KEEP SPLIT for parent %d "
                        + "(mixture is genuinely cheaper -> miscibility gap: deltaG_merge %.4f > tol %.4f)",
                phaseIndex, deltaGMerge, tol));
        return MergeOutcome.keepSplit();
    }

    /**
     * For every phase NOT in the selected set, pick the grid point whose Gibbs
     * value is closest to the discrete tangent plane
     * {@code |G_q - sum_A mu_A M_q,A|} (STEP 7 PART 13-14). If {@code mu} is
     * null (underdetermined), fall back to the phase's lowest-G grid point.
     */
    static List<PhaseInitialState> assignMetastableInitialStates(
            List<HillertSolver.Phase> phases, List<GridPoint> grid,
            List<GridPoint> selPts, double[] mu, int k,
            List<String> diag, Consumer<String> progress) {

        // Parent-phase count is the max phaseIndex present (grid or selection),
        // not phases.size() -- the synthetic test path passes an empty phases
        // list but real grid points still carry phaseIndex 0..P-1.
        int maxPi = phases.size() - 1;
        for (GridPoint gp : grid) {
            maxPi = Math.max(maxPi, gp.phaseIndex);
        }
        for (GridPoint gp : selPts) {
            maxPi = Math.max(maxPi, gp.phaseIndex);
        }
        int pCount = maxPi + 1;
        boolean[] selectedPhase = new boolean[pCount];
        for (GridPoint gp : selPts) {
            selectedPhase[gp.phaseIndex] = true;
        }

        List<PhaseInitialState> out = new ArrayList<>();
        for (int pi = 0; pi < pCount; pi++) {
            if (selectedPhase[pi]) {
                continue;
            }
            String piLabel = pi < phases.size() ? phases.get(pi).label : ("phase#" + pi);
            GridPoint bestPt = null;
            double bestScore = Double.POSITIVE_INFINITY;
            for (GridPoint gp : grid) {
                if (gp.phaseIndex != pi) {
                    continue;
                }
                double score = (mu != null) ? tangentDistance(gp, mu, k) : gp.g;
                if (score < bestScore
                        || (score == bestScore && bestPt != null && gp.gridIndex < bestPt.gridIndex)) {
                    bestScore = score;
                    bestPt = gp;
                }
            }
            if (bestPt == null) {
                emit(progress, diag, String.format(
                        "grid minimizer: metastable phase '%s' has no valid grid point -> no initial state",
                        piLabel));
                continue;
            }
            double td = (mu != null) ? tangentDistance(bestPt, mu, k) : Double.NaN;
            out.add(new PhaseInitialState(
                    pi, bestPt.phaseLabel, 1, 0.0,
                    bestPt.y, bestPt.m, bestPt.g, bestPt.gridIndex, td));
            emit(progress, diag, String.format(
                    "grid minimizer: metastable '%s' initialised at grid point %d "
                            + "(tangent distance %.6f, gamma=%.6f)",
                    bestPt.phaseLabel, bestPt.gridIndex, td,
                    mu != null ? gamma(bestPt, mu, k) : Double.NaN));
        }
        return out;
    }

    // =========================================================================
    // small helpers
    // =========================================================================

    private static double tangentDistance(GridPoint gp, double[] mu, int k) {
        if (mu == null) {
            return Double.NaN;
        }
        return Math.abs(gp.g - dot(mu, gp.m, k));
    }

    /** gamma_q = sum_A mu_A M_q,A - G_q  (>0 => grid-unstable, would be added). */
    private static double gamma(GridPoint gp, double[] mu, int k) {
        return dot(mu, gp.m, k) - gp.g;
    }

    private static boolean hasMiscibilityGap(List<PhaseInitialState> stable) {
        for (int i = 0; i < stable.size(); i++) {
            for (int j = i + 1; j < stable.size(); j++) {
                if (stable.get(i).phaseIndex() == stable.get(j).phaseIndex()) {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean isPureElement(double[] x) {
        for (double v : x) {
            if (v >= 1.0 - 1e-12) {
                return true;
            }
        }
        return false;
    }

    private static boolean allFinite(double[] v) {
        for (double d : v) {
            if (!Double.isFinite(d)) {
                return false;
            }
        }
        return true;
    }

    private static boolean anyNegative(double[] v, double tol) {
        for (double d : v) {
            if (d < -tol) {
                return true;
            }
        }
        return false;
    }

    private static double dot(double[] a, double[] b, int n) {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += a[i] * b[i];
        }
        return s;
    }

    private static double[] diff(double[] a, double[] b) {
        double[] d = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            d[i] = a[i] - b[i];
        }
        return d;
    }

    private static double infNorm(double[] v) {
        double m = 0.0;
        for (double d : v) {
            m = Math.max(m, Math.abs(d));
        }
        return m;
    }

    private static int safeLen(double[] v) {
        return v == null ? -1 : v.length;
    }

    private static double[] round(double[] v) {
        double[] r = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            r[i] = Math.round(v[i] * 1e6) / 1e6;
        }
        return r;
    }

    private static void emit(Consumer<String> sink, List<String> diag, String msg) {
        diag.add(msg);
        if (sink != null) {
            sink.accept(msg);
        }
    }

    private static GridMinimizationResult failure(
            String reason, int gridCount, int skipped, List<String> diag) {
        return new GridMinimizationResult(
                false, reason,
                List.of(), null,
                List.of(), List.of(),
                null, "n/a",
                Double.NaN, null, Double.NaN,
                gridCount, skipped,
                diag);
    }
}
