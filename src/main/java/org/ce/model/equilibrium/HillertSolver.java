package org.ce.model.equilibrium;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import org.ce.model.ModelSession;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertPhaseStepSolver.PerPhaseStepResult;

/**
 * N-phase equilibrium solver (Hillert's method), ported from the reference
 * {@code phaseq} Mathematica implementation -- see HILLERT_SOLVER_PLAN.md.
 *
 * <p>Outer loop, once per iteration:
 * <ol>
 *   <li>{@link HillertPhaseStepSolver} for every phase -- each
 *       phase's Newton step {@code deltaY(mu)} as an affine function of the
 *       shared trial chemical potential {@code mu}.</li>
 *   <li>Assemble {@link EquilibriumMatrix.PhaseContribution} for every
 *       currently-<em>stable</em> phase ({@code amount > 0}) using
 *       {@code G = G0m + Gm} (absolute reference energy -- see
 *       {@link LatticeStability}) and solve for numeric {@code mu}/{@code
 *       deltaN}.</li>
 *   <li>Backtracking line search: starting at {@code lambda=1}, evaluate the
 *       trial step for every phase and halve {@code lambda} until all
 *       phases' trial state is valid ({@code State.isValidIncludingPoints}),
 *       up to {@code jtrMax} tries.</li>
 *   <li>Accept the trial state, then check convergence via the per-phase
 *       step-norm.</li>
 * </ol>
 *
 * <p>Convergence follows the reference exactly: {@code Min} of the per-phase
 * {@code deltaY(mu)} L2 norms, not {@code Max} -- ported faithfully even
 * though this is an unusual choice, per session precedent of preserving
 * reference behavior rather than "fixing" it.</p>
 */
public final class HillertSolver {

    private HillertSolver() {}

    // =========================================================================
    // Inputs and outputs
    //
    // Nested here for the same reason CvmNewtonSolver.Result is nested in its
    // solver: these are this solver's working state and its output, used
    // nowhere else. Keeping them alongside the loop that drives them means the
    // whole multi-phase contract reads in one place.
    // =========================================================================

    /**
     * One candidate phase's mutable state during the solve.
     *
     * <p>The {@link CVMGibbsModel} it carries is a pure evaluator and holds no
     * per-point state, so it may safely be shared between phases of the same
     * system; what is mutable here is {@link #amount} and {@link #uFull}, which
     * the outer loop updates. (An earlier version required a separate model
     * instance per phase, because the model then carried a current
     * {@code (T, x, u)} internally.)</p>
     *
     * <p>Distinct from a single-phase {@link ModelSession}-driven calculation:
     * a Hillert phase's composition is itself an unknown, solved for jointly
     * with its internal CVM parameters by {@link HillertPhaseStepSolver} -- not
     * a fixed input the way {@link org.ce.calculation.Conditions} treats it for
     * {@code CalculationService.calculate}.</p>
     */
    public static final class Phase {

        /**
         * Amount below which a phase is treated as unstable and excluded from
         * the outer equilibrium assembly. Matches the reference's
         * amount-sign-only check ({@code amount > 0}), not a rigorous Gibbs
         * phase rule.
         */
        private static final double STABILITY_THRESHOLD = 0.0;

        public final String label;
        public final ModelSession session;
        public final CVMGibbsModel model;
        public final int ncf;
        public final int numComponents;

        /** Current phase amount (moles of formula units, "N" in the reference). */
        public double amount;

        /** Current joint internal-parameter vector {@code uFull = [u ; x]}, length {@code ncf+K}. */
        public double[] uFull;

        public Phase(String label, ModelSession session, CVMGibbsModel model,
                double initialAmount, double[] initialUFull) {
            this.label = label;
            this.session = session;
            this.model = model;
            this.ncf = model.ncf();
            this.numComponents = initialUFull.length - ncf;
            this.amount = initialAmount;
            this.uFull = initialUFull.clone();
        }

        /**
         * Current composition -- the trailing {@code K} entries of
         * {@link #uFull}. Port of the reference's {@code updateComp}:
         * composition is always exactly this slice, never a separate inversion.
         */
        public double[] composition() {
            double[] x = new double[numComponents];
            System.arraycopy(uFull, ncf, x, 0, numComponents);
            return x;
        }

        /**
         * True if this phase is currently treated as stable (amount strictly
         * positive) -- the reference's amount-sign-only check, not a rigorous
         * Gibbs phase rule.
         */
        public boolean isStable() {
            return amount > STABILITY_THRESHOLD;
        }
    }

    /**
     * Immutable output of {@link #solve} -- the multi-phase counterpart to
     * {@code ThermodynamicResult} for the single-phase path.
     *
     * <p><b>Check {@link #overallConverged} before using any value.</b> A
     * non-converged run still returns plausible-looking numbers.</p>
     */
    public record Result(
            List<PhaseResult> phases,
            double[] mu,
            boolean overallConverged,
            int outerIterations,
            double finalResidualNorm) {
    }

    /**
     * One phase's outcome: amount, composition, and energetics at the final
     * iterate.
     *
     * <p>{@link #state} is that phase's model evaluated at its converged joint
     * point {@code uFull = [u ; x]} -- the same object {@code model.atFull(T,
     * uFull)} would produce, retained rather than discarded so every other
     * property is reachable without re-solving or re-evaluating:</p>
     *
     * <pre>
     *   for (PhaseResult p : eq.phases()) {
     *       double s    = p.state().sm();              // entropy of this phase
     *       double[] dg = p.state().gmuFull();         // its widened gradient
     *       var sro     = p.state().pairSroByShell();  // its short-range order
     *   }
     * </pre>
     *
     * <p>{@link #g} is kept as its own field because it is the quantity the
     * outer equilibrium assembly actually solved with -- the absolute
     * {@code G = G0m + Gm}, which must share one zero across phases for
     * chemical potentials to be comparable. It equals {@code state().g()}; the
     * field records what the solve used, the state offers everything else.</p>
     */
    public record PhaseResult(
            String label,
            double amount,
            double[] composition,
            double g,
            CVMGibbsModel.State state,
            boolean phaseConverged) {

        /** The model this phase was evaluated against. */
        public CVMGibbsModel model() {
            return state.model();
        }
    }


    /**
     * Runs the outer/inner Hillert iteration to equilibrium.
     *
     * @param phases per-phase mutable state (updated in place on convergence path)
     * @param temperature fixed temperature (K) -- GxT/GxP not yet supported, v1 is fixed-T,P
     * @param maxOuterIterations outer Newton iteration cap
     * @param innerBacktrackTries max lambda-halving tries per outer iteration
     * @param tol convergence tolerance on the min per-phase step norm
     */
    public static Result solve(
            List<Phase> phases,
            double temperature,
            int maxOuterIterations,
            int innerBacktrackTries,
            double tol,
            Consumer<String> progressSink) {

        int numComponents = phases.get(0).numComponents;
        double[] mu = new double[numComponents];
        double finalResidualNorm = Double.POSITIVE_INFINITY;
        int outerIter = 0;
        boolean converged = false;

        for (outerIter = 1; outerIter <= maxOuterIterations; outerIter++) {
            List<PerPhaseStepResult> steps = new ArrayList<>(phases.size());
            for (Phase phase : phases) {
                steps.add(new HillertPhaseStepSolver(phase.model).step(phase.uFull, temperature));
            }

            List<EquilibriumMatrix.PhaseContribution> contributions = new ArrayList<>();
            List<Integer> stableIndices = new ArrayList<>();
            for (int p = 0; p < phases.size(); p++) {
                Phase phase = phases.get(p);
                if (!phase.isStable()) {
                    continue;
                }
                PerPhaseStepResult step = steps.get(p);
                double g = currentG(phase, temperature);
                contributions.add(new EquilibriumMatrix.PhaseContribution(
                        phase.amount, phase.composition(), g,
                        step.deltaComposition0(), step.deltaCompositionSensitivity()));
                stableIndices.add(p);
            }

            EquilibriumMatrix.EquilibriumStepResult eqStep =
                    EquilibriumMatrix.solve(contributions, numComponents);
            mu = eqStep.mu();
            double[] deltaN = eqStep.deltaN();

            double lambda = 1.0;
            double[][] trialUFull = new double[phases.size()][];
            double[] trialAmount = new double[phases.size()];
            boolean accepted = false;

            for (int tries = 0; tries < innerBacktrackTries; tries++) {
                boolean allValid = true;
                for (int p = 0; p < phases.size(); p++) {
                    Phase phase = phases.get(p);
                    double[] deltaY = steps.get(p).deltaYAt(mu);
                    double[] u = new double[phase.uFull.length];
                    for (int i = 0; i < u.length; i++) {
                        u[i] = phase.uFull[i] + lambda * deltaY[i];
                    }
                    trialUFull[p] = u;

                    int stableSlot = stableIndices.indexOf(p);
                    double amount = phase.amount;
                    if (stableSlot >= 0) {
                        amount = phase.amount + lambda * deltaN[stableSlot];
                    }
                    trialAmount[p] = amount;

                    double[] uOnly = java.util.Arrays.copyOfRange(u, 0, phase.ncf);
                    double[] xOnly = java.util.Arrays.copyOfRange(u, phase.ncf, u.length);
                    if (!phase.model.at(temperature, xOnly, uOnly).isValidIncludingPoints()) {
                        allValid = false;
                    }
                }
                if (allValid) {
                    accepted = true;
                    break;
                }
                lambda *= 0.5;
            }

            // Reference (phaseq) has no special handling for exhausting all
            // itr=1..10 backtracking tries without success: the state is
            // simply left unchanged (u2ListN/NListN untouched) and the outer
            // loop proceeds to jtr+1, which recomputes a fresh Newton step
            // from the same point. Port that literally -- do not abort the
            // outer loop here.
            if (accepted) {
                for (int p = 0; p < phases.size(); p++) {
                    phases.get(p).uFull = trialUFull[p];
                    phases.get(p).amount = trialAmount[p];
                }
            } else if (progressSink != null) {
                progressSink.accept("Hillert outer iter " + outerIter
                        + ": backtracking exhausted, state unchanged, continuing.");
            }

            double minNorm = Double.POSITIVE_INFINITY;
            for (int p = 0; p < phases.size(); p++) {
                double[] deltaY = steps.get(p).deltaYAt(mu);
                double norm = l2Norm(deltaY);
                minNorm = Math.min(minNorm, norm);
            }

            finalResidualNorm = minNorm;
            if (progressSink != null) {
                progressSink.accept("Hillert outer iter " + outerIter
                        + ": lambda=" + lambda + " accepted=" + accepted + " minNorm=" + minNorm);
            }

            if (minNorm <= tol) {
                converged = true;
                outerIter++;
                break;
            }
        }

        // Evaluate each phase once at its final point and keep the state, so a
        // caller can read any further property (entropy, gradients, SRO) on
        // demand rather than re-evaluating. G comes from that same state, so
        // the reported energy and anything derived from it cannot disagree.
        List<PhaseResult> entries = new ArrayList<>();
        for (Phase phase : phases) {
            CVMGibbsModel.State state = phase.model.atFull(temperature, phase.uFull);
            entries.add(new PhaseResult(
                    phase.label, phase.amount, phase.composition(), state.g(), state, converged));
        }

        return new Result(entries, mu, converged, outerIter - 1, finalResidualNorm);
    }

    /**
     * {@code G = G0m + Gm} at this phase's current joint state -- the absolute,
     * pure-element-anchored energy the outer equilibrium assembly needs, since
     * chemical potentials are equalised across phases and every phase's G must
     * therefore share one zero.
     *
     * <p>Evaluated through {@link org.ce.model.cvm.CVMGibbsModel}, which composes
     * the reference and mixing terms itself. This previously called
     * {@code model.evaluate(...).G} for Gm and added
     * {@code LatticeStability.g0m} separately -- a second place where
     * {@code G = G0m + Gm} was spelled out, and one that mutated the phase's
     * model as a side effect of a read.</p>
     */
    private static double currentG(Phase phase, double temperature) {
        return phase.model.atFull(temperature, phase.uFull).g();
    }

    private static double l2Norm(double[] v) {
        double sum = 0.0;
        for (double value : v) {
            sum += value * value;
        }
        return Math.sqrt(sum);
    }
}
