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

    /** One phase's outcome: amount, composition, and energetics at the final iterate. */
    public record PhaseResultEntry(
            String label,
            double amount,
            double[] composition,
            double g,
            boolean phaseConverged) {}
}
