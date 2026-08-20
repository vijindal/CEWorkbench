package org.ce.model.equilibrium;

import org.ce.model.ModelSession;
import org.ce.model.cvm.CVMGibbsModel;

/**
 * Per-phase mutable state for the Hillert multi-phase equilibrium solver
 * (HILLERT_SOLVER_PLAN.md) — one instance per candidate phase.
 *
 * <p>The {@link CVMGibbsModel} it carries is a pure evaluator and holds no
 * per-point state, so it may safely be shared between phases of the same
 * system; the mutable state of a phase is this object's {@code amount} and
 * {@code uFull}, not the model's. (An earlier version required a separate
 * model instance per phase because the model then carried a current
 * {@code (T, x, u)} internally.)</p>
 *
 * <p>Distinct from a single-phase {@link ModelSession}-driven calculation:
 * a Hillert phase's composition is itself an unknown solved for jointly
 * with its internal CVM parameters (see {@code HillertPhaseStepSolver}), not a fixed input the way {@link
 * org.ce.calculation.Conditions} treats it for {@code
 * CalculationService.calculate}.</p>
 */
public final class PhaseState {

    /** Amount below which a phase is treated as unstable and excluded from the outer equilibrium assembly (matches the reference's amount-sign-only stability check: {@code amount > 0}, not a rigorous Gibbs phase rule). */
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

    public PhaseState(String label, ModelSession session, CVMGibbsModel model,
            double initialAmount, double[] initialUFull) {
        this.label = label;
        this.session = session;
        this.model = model;
        this.ncf = model.ncf();
        this.numComponents = initialUFull.length - ncf;
        this.amount = initialAmount;
        this.uFull = initialUFull.clone();
    }

    /** Current composition — the trailing {@code K} entries of {@link #uFull} (port of the reference's {@code updateComp}: composition is always exactly this slice, never a separate inversion). */
    public double[] composition() {
        double[] x = new double[numComponents];
        System.arraycopy(uFull, ncf, x, 0, numComponents);
        return x;
    }

    /** True if this phase is currently treated as stable (amount strictly positive) — matches the reference's amount-sign-only check, not a rigorous Gibbs phase rule. */
    public boolean isStable() {
        return amount > STABILITY_THRESHOLD;
    }
}
