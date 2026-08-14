package org.ce.model.equilibrium;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CVMGibbsModel.PerPhaseStepResult;

/**
 * N-phase equilibrium solver (Hillert's method), ported from the reference
 * {@code phaseq} Mathematica implementation -- see HILLERT_SOLVER_PLAN.md.
 *
 * <p>Outer loop, once per iteration:
 * <ol>
 *   <li>{@link CVMGibbsModel#solvePerPhaseStep} for every phase -- each
 *       phase's Newton step {@code deltaY(mu)} as an affine function of the
 *       shared trial chemical potential {@code mu}.</li>
 *   <li>Assemble {@link EquilibriumMatrix.PhaseContribution} for every
 *       currently-<em>stable</em> phase ({@code amount > 0}) using
 *       {@code G = G0m + Gm} (absolute reference energy -- see
 *       {@link LatticeStability}) and solve for numeric {@code mu}/{@code
 *       deltaN}.</li>
 *   <li>Backtracking line search: starting at {@code lambda=1}, evaluate the
 *       trial step for every phase and halve {@code lambda} until all
 *       phases' trial state is valid ({@link CVMGibbsModel#isValidParams}),
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

    /**
     * Runs the outer/inner Hillert iteration to equilibrium.
     *
     * @param phases per-phase mutable state (updated in place on convergence path)
     * @param temperature fixed temperature (K) -- GxT/GxP not yet supported, v1 is fixed-T,P
     * @param maxOuterIterations outer Newton iteration cap
     * @param innerBacktrackTries max lambda-halving tries per outer iteration
     * @param tol convergence tolerance on the min per-phase step norm
     */
    public static PhaseEquilibriumResult solve(
            List<PhaseState> phases,
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
            for (PhaseState phase : phases) {
                steps.add(phase.model.solvePerPhaseStep(phase.uFull, temperature));
            }

            List<EquilibriumMatrix.PhaseContribution> contributions = new ArrayList<>();
            List<Integer> stableIndices = new ArrayList<>();
            for (int p = 0; p < phases.size(); p++) {
                PhaseState phase = phases.get(p);
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
                    PhaseState phase = phases.get(p);
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
                    if (!phase.model.isValidParams(uOnly, xOnly)) {
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

        List<PhaseEquilibriumResult.PhaseResultEntry> entries = new ArrayList<>();
        for (PhaseState phase : phases) {
            double g = currentG(phase, temperature);
            entries.add(new PhaseEquilibriumResult.PhaseResultEntry(
                    phase.label, phase.amount, phase.composition(), g, converged));
        }

        return new PhaseEquilibriumResult(entries, mu, converged, outerIter - 1, finalResidualNorm);
    }

    /** G = G0m (absolute lattice-stability reference) + Gm (mixing energy from the CVM evaluator). */
    private static double currentG(PhaseState phase, double temperature) {
        double[] u = java.util.Arrays.copyOfRange(phase.uFull, 0, phase.ncf);
        double[] x = java.util.Arrays.copyOfRange(phase.uFull, phase.ncf, phase.uFull.length);
        double gm = phase.model.evaluate(u, x, temperature).G;
        List<String> elements = List.of(phase.session.systemId.elements().split("-"));
        double g0m = LatticeStability.g0m(elements, phase.session.systemId.structure(), x, temperature);
        return g0m + gm;
    }

    private static double l2Norm(double[] v) {
        double sum = 0.0;
        for (double value : v) {
            sum += value * value;
        }
        return Math.sqrt(sum);
    }
}
