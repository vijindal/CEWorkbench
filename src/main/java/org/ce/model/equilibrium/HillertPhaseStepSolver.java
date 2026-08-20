package org.ce.model.equilibrium;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cvm.CVMGibbsModel;

/**
 * One phase's inner Newton step for the Hillert multi-phase equilibrium
 * solver.
 *
 * <p>Holds a {@link CVMGibbsModel} and drives it from the outside, exactly as
 * {@code CvmNewtonSolver} does for the fixed-composition minimisation. The two
 * differ only in which unknowns they solve for -- and therefore in which block
 * of the same evaluated state they read:</p>
 *
 * <table border="1">
 *   <caption>The two solvers over one model</caption>
 *   <tr><th></th><th>CvmNewtonSolver</th><th>this</th></tr>
 *   <tr><td>composition</td><td>fixed constraint</td><td>an unknown</td></tr>
 *   <tr><td>reads</td><td>{@code gmu} / {@code gmuu} ({@code ncf})</td>
 *       <td>{@code gmuFull} / {@code gmuuFull} ({@code ncf+K})</td></tr>
 *   <tr><td>solves for</td><td>stationary G at fixed x</td>
 *       <td>stationary G relative to a trial mu</td></tr>
 * </table>
 *
 * <p>This previously lived on {@code CVMGibbsModel} itself, on the reasoning
 * that its only inputs were that class's own widened derivatives. That was
 * true, but it made the model both an evaluator and two different solvers.
 * Now the model evaluates and the solvers solve.</p>
 */
public final class HillertPhaseStepSolver {

    private final CVMGibbsModel model;

    public HillertPhaseStepSolver(CVMGibbsModel model) {
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        this.model = model;
    }

    /** The model this solver drives. */
    public CVMGibbsModel model() {
        return model;
    }

    /**
     * Result of {@link #step}: the joint Newton step
     * {@code deltaY(mu)}, expressed as an <b>affine function of the trial
     * chemical-potential vector {@code mu}</b> rather than a value at one
     * fixed {@code mu} -- see the note on {@link #step} for why
     * this shape, not a single numeric result, is what the outer Hillert
     * solver actually needs.
     *
     * <p>{@code deltaY(mu) = deltaY0 + Σ_k mu[k]*deltaYSensitivity[k]}, and
     * likewise for {@code deltaComposition}/{@code lambda}.</p>
     */
    public record PerPhaseStepResult(
            double[] deltaY0, double[][] deltaYSensitivity,
            double[] deltaComposition0, double[][] deltaCompositionSensitivity,
            double lambda0, double[] lambdaSensitivity) {

        /** Evaluates this affine result at a specific numeric {@code mu}. */
        public double[] deltaCompositionAt(double[] mu) {
            double[] result = deltaComposition0.clone();
            for (int k = 0; k < mu.length; k++) {
                for (int i = 0; i < result.length; i++) {
                    result[i] += mu[k] * deltaCompositionSensitivity[k][i];
                }
            }
            return result;
        }

        /** Evaluates the full joint deltaY (length ncf+K) at a specific numeric {@code mu}. */
        public double[] deltaYAt(double[] mu) {
            double[] result = deltaY0.clone();
            for (int k = 0; k < mu.length; k++) {
                for (int i = 0; i < result.length; i++) {
                    result[i] += mu[k] * deltaYSensitivity[k][i];
                }
            }
            return result;
        }
    }

    /**
     * One Hillert multi-phase equilibrium Newton step for this phase,
     * expressed as an <b>affine function of the trial chemical-potential
     * vector {@code mu}</b> -- port of the reference Mathematica
     * implementation's {@code delxGCVM}.
     *
     * <p><b>Why affine-in-mu, not a fixed-mu numeric result:</b> tracing
     * {@code phaseq}'s outer loop shows {@code delxGCVM} is called with
     * {@code mu} still <em>symbolic</em>, and {@code genEqMat}'s outer
     * mass-balance equations substitute that symbolic result in directly, so
     * {@code mu} and {@code deltaN} are solved <em>simultaneously</em> in one
     * combined system -- the per-phase step is never evaluated at a numeric
     * {@code mu} on its own. Since {@code deltaY} is provably affine in
     * {@code mu} (the matrix {@code A} does not depend on it; only the
     * right-hand side does, and only in the x-block rows), the numeric
     * equivalent is to solve the same system {@code K+1} times against basis
     * right-hand sides -- once for {@code mu=0}, once per unit vector -- and
     * let {@link EquilibriumMatrix} fold the affine
     * form into its own equations, mirroring {@code genEqMat}'s substitution
     * entirely numerically.</p>
     *
     * <p>Unrelated to {@code CvmNewtonSolver}'s single-phase Newton-Raphson loop
     * (fixed composition, stationary {@code G}) and must not be confused with
     * it: this solves for a stationary point of {@code G} <em>relative to a
     * trial {@code mu}</em>, with composition itself among the unknowns.</p>
     *
     * <p><b>The linear system</b> (at fixed T/P, so the {@code GxT*ΔT} and
     * {@code GxP*ΔP} terms vanish): unknowns are {@code deltaY[0..ncf+K-1]}
     * and {@code lambda}, over {@code ncf+K+1} equations:</p>
     * <ul>
     *   <li>Rows {@code 0..ncf-1} (u-block): {@code Guu[i,:] . deltaY = -Gu[i]}
     *       -- ordinary stationarity on the internal CFs, unconstrained by
     *       {@code mu}.</li>
     *   <li>Rows {@code ncf..ncf+K-1} (x-block): {@code Guu[i,:] . deltaY -
     *       lambda = mu[i-ncf] - Gu[i]} -- the only rows where {@code mu}
     *       appears, always with coefficient exactly {@code +1} on its own
     *       row, which is why one basis solve per component suffices.</li>
     *   <li>Row {@code ncf+K}: {@code sum(deltaY[ncf..]) = 0} -- the
     *       composition change stays on the simplex.</li>
     * </ul>
     *
     * <p>Built entirely from analytic widened derivatives -- no
     * finite-differencing anywhere.</p>
     *
     * @param uFull current joint state {@code [u ; x]}, length {@code ncf+K}
     * @param temperature current temperature, K
     */
    public PerPhaseStepResult step(double[] uFull, double temperature) {
        int ncf = model.ncf();
        int numComponents = model.numComponents();
        int width = ncf + numComponents;
        if (uFull.length != width) {
            throw new IllegalArgumentException(
                    "uFull.length=" + uFull.length + " != ncf+K=" + width);
        }

        CVMGibbsModel.State state = model.atFull(temperature, uFull);
        double[] Gu = state.gmuFull();
        double[][] Guu = state.gmuuFull();

        int n = width + 1; // + lambda

        // Matrix A is the same for every right-hand side (mu does not appear
        // in it) -- build once, reuse for all K+1 solves.
        double[][] A = new double[n][n];
        for (int i = 0; i < width; i++) {
            System.arraycopy(Guu[i], 0, A[i], 0, width);
        }
        for (int i = ncf; i < width; i++) {
            A[i][width] = -1.0; // -lambda
        }
        for (int i = ncf; i < width; i++) {
            A[width][i] = 1.0; // sum(deltaX) = 0
        }

        // b0: the mu=0 right-hand side.
        double[] b0 = new double[n];
        for (int i = 0; i < width; i++) b0[i] = -Gu[i];
        double[] sol0 = LinearAlgebra.solve(A, b0);

        // Solving A*z = e_{ncf+k} directly gives d(deltaY)/d(mu_k), since the
        // system is linear and A is shared across right-hand sides.
        double[][] deltaYSens = new double[numComponents][];
        double[] lambdaSens = new double[numComponents];
        double[][] deltaCompSens = new double[numComponents][];
        for (int k = 0; k < numComponents; k++) {
            double[] ek = new double[n];
            ek[ncf + k] = 1.0;
            double[] solK = LinearAlgebra.solve(A, ek);
            double[] deltaYk = new double[width];
            System.arraycopy(solK, 0, deltaYk, 0, width);
            deltaYSens[k] = deltaYk;
            lambdaSens[k] = solK[width];
            double[] deltaCompK = new double[numComponents];
            System.arraycopy(deltaYk, ncf, deltaCompK, 0, numComponents);
            deltaCompSens[k] = deltaCompK;
        }

        double[] deltaY0 = new double[width];
        System.arraycopy(sol0, 0, deltaY0, 0, width);
        double[] deltaComposition0 = new double[numComponents];
        System.arraycopy(deltaY0, ncf, deltaComposition0, 0, numComponents);
        double lambda0 = sol0[width];

        return new PerPhaseStepResult(deltaY0, deltaYSens, deltaComposition0, deltaCompSens, lambda0, lambdaSens);
    }

}
