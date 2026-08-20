package org.ce.model.equilibrium;

import java.util.List;

/**
 * Immutable output of {@link HillertSolver#solve} — mirrors {@code
 * ThermodynamicResult}'s role for the single-phase path, but for a
 * multi-phase equilibrium (HILLERT_SOLVER_PLAN.md).
 */
public record PhaseEquilibriumResult(
        List<PhaseResultEntry> phases,
        double[] mu,
        boolean overallConverged,
        int outerIterations,
        double finalResidualNorm) {

    /**
     * One phase's outcome: amount, composition, and energetics at the final
     * iterate.
     *
     * <p>{@link #state} is that phase's model evaluated at its converged joint
     * point {@code uFull = [u ; x]} — the same object {@code model.atFull(T,
     * uFull)} would produce, retained rather than discarded so every other
     * property is reachable without re-solving or re-evaluating:</p>
     *
     * <pre>
     *   for (PhaseResultEntry p : eq.phases()) {
     *       double s   = p.state().sm();               // entropy of this phase
     *       double[] dg = p.state().gmuFull();         // its widened gradient
     *       var sro    = p.state().pairSroByShell();   // its short-range order
     *   }
     * </pre>
     *
     * <p>{@link #g} is kept as its own field because it is the quantity the
     * outer equilibrium assembly actually solved with — the absolute
     * {@code G = G0m + Gm}, which must share one zero across phases for
     * chemical potentials to be comparable. It equals {@code state().g()};
     * the field records what the solve used, the state offers everything
     * else.</p>
     */
    public record PhaseResultEntry(
            String label,
            double amount,
            double[] composition,
            double g,
            org.ce.model.cvm.CVMGibbsModel.State state,
            boolean phaseConverged) {

        /** The model this phase was evaluated against. */
        public org.ce.model.cvm.CVMGibbsModel model() {
            return state.model();
        }
    }
}
