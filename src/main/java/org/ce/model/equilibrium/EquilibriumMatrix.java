package org.ce.model.equilibrium;

import org.ce.model.cluster.LinearAlgebra;

import java.util.List;

/**
 * Outer-loop mass-balance + chemical-potential-equality assembly for
 * Hillert's method — port of the reference Mathematica implementation's
 * {@code genEqMat}, as actually called from {@code phaseq}'s outer
 * iteration (not a reconstruction from first principles; the full
 * {@code genEqMat} source was obtained and ported directly).
 *
 * <p><b>Where this sits relative to {@code HillertSolver.PhaseStep}
 * (HILLERT_SOLVER_PLAN.md):</b> {@code HillertSolver.PhaseStep.step} solves one
 * phase's local Newton step as an <em>affine function of the trial
 * chemical potential {@code mu}</em>, not at one fixed numeric {@code mu}
 * — tracing {@code phaseq}'s actual outer loop showed {@code delxGCVM} is
 * called with {@code mu} still symbolic, and {@code genEqMat}'s mass-
 * balance equations substitute that symbolic {@code delnN} in directly
 * before solving, so {@code mu} and {@code deltaN} are solved
 * <em>simultaneously</em> in one combined system. This class is the
 * numeric equivalent of that symbolic substitution: it takes each phase's
 * affine {@code deltaComposition(mu) = deltaComposition0 +
 * Σ_k mu[k]*deltaCompositionSensitivity[k]} and folds the {@code mu}
 * coefficients directly into the mass-balance rows' left-hand side,
 * producing one {@code (K+np)×(K+np)} system solved for {@code mu} and
 * {@code deltaN} together — exactly mirroring {@code genEqMat}'s
 * substitution, with no symbolic algebra.</p>
 *
 * <p><b>The linear system</b>, from {@code genEqMat}'s {@code gExpr}/
 * {@code nExpr} (fixed T, P, so the {@code GT*ΔT}/{@code GP*ΔP} terms
 * vanish): unknowns are {@code mu[0..K-1]} (chemical potentials) and
 * {@code deltaN[0..np-1]} (phase amount changes), {@code K+np} equations:</p>
 * <ul>
 *   <li><b>Gibbs-Duhem rows</b> (one per phase, {@code np} rows):
 *       {@code sum_i mu[i] * composition[phase][i] = G[phase]} — a phase's
 *       Gibbs energy must equal the composition-weighted sum of chemical
 *       potentials, the standard partial-molar-quantity relation. {@code G}
 *       doesn't depend on {@code mu}, so these rows are unaffected by the
 *       affine substitution below.</li>
 *   <li><b>Mass-balance rows</b> (one per component, {@code K} rows):
 *       {@code sum_phase deltaN[phase]*composition[phase][i] +
 *       sum_phase amount[phase]*deltaComposition[phase](mu)[i] = 0}.
 *       Substituting the affine form and collecting {@code mu}'s
 *       coefficients onto the left-hand side:
 *       {@code sum_phase deltaN[phase]*x[phase][i]
 *       + Σ_k mu[k]*(Σ_phase amount[phase]*deltaCompositionSensitivity[phase][k][i])
 *       = -Σ_phase amount[phase]*deltaComposition0[phase][i]}.</li>
 * </ul>
 *
 * <p>Deliberately a separate class from {@code CVMGibbsModel} (unlike
 * {@code HillertSolver.PhaseStep.step}): its inputs are aggregated quantities from
 * <em>multiple</em> phases at once, not one {@code CVMGibbsModel}
 * instance's own state — there is no single natural owner among the
 * per-phase models for this cross-phase assembly.</p>
 */
public final class EquilibriumMatrix {

    private EquilibriumMatrix() {}

    /**
     * One phase's contribution to the outer system: current amount,
     * composition, {@code G}, and its affine {@code deltaComposition(mu)}
     * (from {@link HillertSolver.PhaseStep#step}) --
     * {@code deltaComposition0} is the {@code mu=0} value,
     * {@code deltaCompositionSensitivity[k]} is {@code d(deltaComposition)/d(mu_k)}.
     */
    public record PhaseContribution(
            double amount, double[] composition, double g,
            double[] deltaComposition0, double[][] deltaCompositionSensitivity) {}

    /** Solved outer step: updated chemical potentials and each phase's amount change. */
    public record EquilibriumStepResult(double[] mu, double[] deltaN) {}

    /**
     * Solves the combined outer Gibbs-Duhem + mass-balance system for one
     * outer iteration, with each phase's {@code mu}-dependence folded in
     * algebraically (see class doc).
     *
     * @param phases K-component contributions from every phase currently
     *               treated as stable (the reference's {@code
     *               unStablePhaseRules} excludes non-positive-amount phases
     *               from this assembly entirely — callers must filter
     *               before calling, not pass all phases and expect
     *               filtering here)
     * @param numComponents K, the number of system components
     */
    public static EquilibriumStepResult solve(List<PhaseContribution> phases, int numComponents) {
        int np = phases.size();
        int n = numComponents + np; // unknowns: mu[0..K-1], deltaN[0..np-1]

        double[][] A = new double[n][n];
        double[] b = new double[n];

        // Gibbs-Duhem rows (0..np-1): sum_i mu[i]*x[phase][i] = G[phase]
        // (G has no mu-dependence, so this block is unaffected by the affine substitution below)
        for (int p = 0; p < np; p++) {
            PhaseContribution phase = phases.get(p);
            for (int i = 0; i < numComponents; i++) {
                A[p][i] = phase.composition()[i];
            }
            b[p] = phase.g();
        }

        // Mass-balance rows (np..np+K-1), with each phase's affine
        // deltaComposition(mu) substituted in: mu's coefficients move onto
        // the left-hand side (columns 0..K-1), deltaN's coefficients stay
        // as before (columns K..K+np-1), and only the mu=0 remainder stays
        // on the right-hand side.
        for (int i = 0; i < numComponents; i++) {
            int row = np + i;
            double rhs = 0.0;
            for (int p = 0; p < np; p++) {
                PhaseContribution phase = phases.get(p);
                A[row][numComponents + p] = phase.composition()[i];
                for (int k = 0; k < numComponents; k++) {
                    A[row][k] += phase.amount() * phase.deltaCompositionSensitivity()[k][i];
                }
                rhs -= phase.amount() * phase.deltaComposition0()[i];
            }
            b[row] = rhs;
        }

        double[] sol = LinearAlgebra.solve(A, b);

        double[] mu = new double[numComponents];
        System.arraycopy(sol, 0, mu, 0, numComponents);
        double[] deltaN = new double[np];
        System.arraycopy(sol, numComponents, deltaN, 0, np);

        return new EquilibriumStepResult(mu, deltaN);
    }
}
