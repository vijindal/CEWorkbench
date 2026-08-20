package org.ce.model.cvm;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.ProgressEvent;

import java.util.Arrays;
import java.util.concurrent.CancellationException;
import java.util.function.Consumer;

/**
 * Newton-Raphson minimisation of the CVM Gibbs free energy at <b>fixed
 * composition</b>: given {@code (T, x)}, finds the {@code u} that makes
 * {@code dGm/du = 0}.
 *
 * <p>Extracted from {@code CVMGibbsModel}, which was both the evaluator and
 * this optimiser. It now drives a {@link CVMGibbsModel} from the outside,
 * holding no physics of its own -- every energy, gradient and Hessian comes
 * from a {@link CVMGibbsModel.State} evaluated at the current iterate. What lives here is
 * only the algorithm: convergence criteria, step limiting, and the iteration
 * bookkeeping.</p>
 *
 * <p><b>Not to be confused with the Hillert per-phase step.</b> This solves at
 * fixed composition, with {@code x} a constraint and {@code u} the unknown.
 * {@code HillertPhaseStepSolver} treats composition as an unknown too
 * and solves for a stationary point relative to a trial chemical potential.
 * Both read the same {@link CVMGibbsModel.State}; they differ only in which block of it
 * they use -- the {@code ncf}-wide gradient here, the {@code (ncf+K)}-wide one
 * there.</p>
 *
 * <p>Stateless and safe to share: each {@link #solve} call carries its own
 * iterate, so one instance may serve many points. The {@code (T, x)} result
 * cache that {@code CVMGibbsModel.getEquilibriumState} keeps is memoisation of
 * this solver's output, not part of the algorithm, and stays there.</p>
 */
public final class CvmNewtonSolver {

    /** Iteration cap. A run that hits this is reported as non-converged. */
    public static final int MAX_ITER = 100;

    /**
     * Convergence threshold on the raw (unclamped) Newton step. A tiny step is
     * a genuine sign of convergence, unlike a small step produced purely by the
     * boundary clamp.
     */
    public static final double TOLX = 1.0e-12;

    private final CVMGibbsModel model;

    public CvmNewtonSolver(CVMGibbsModel model) {
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        this.model = model;
    }

    /** Outcome of a minimisation: the converged point and how it was reached. */
    public record Result(
            CVMGibbsModel.State state,
            double[] u,
            boolean converged,
            int iterations,
            double finalGradientNorm) {
    }

    /**
     * Minimises {@code Gm} over {@code u} at fixed {@code (T, x)}.
     *
     * <p>Loop structure follows the original reference solver,
     * {@code CVMBINCE.minimize()} (which solves in the <em>orthogonal</em> CF
     * basis): an early exit the moment any cluster variable is non-positive at
     * the <em>current</em> point (not just a step-limit check on the trial
     * point), a gradient-norm convergence check, then a Newton step whose size
     * is limited by {@link #stepLimit} (port of the reference's {@code stpmx})
     * so no cluster variable leaves {@code [0, 1]}, followed by a
     * raw-Newton-step-size check ({@code errx}).</p>
     *
     * <p>The reference's additional pre-clamp on the trial correlation
     * functions themselves ({@code Utils.normalU}, restricting them to a fixed
     * {@code [-1, 1]}) was tried and deliberately dropped: that bound is
     * meaningful for the orthogonal basis's Chebyshev-like CFs, but this solver
     * works in the CVCF basis, whose {@code u} components have a different,
     * non-uniform natural range (observed roughly 1e-4 to 0.2 in practice) -- a
     * blind &plusmn;1 clamp there is a no-op, not a safeguard, and did not
     * change behaviour when tested. The near-edge convergence stall this
     * porting effort was investigating (a dilute-composition Newton direction
     * that repeatedly re-approaches, but never crosses, a cluster-variable
     * boundary) remains open; see CLAUDE.md's note on near-edge ternary solver
     * fragility.</p>
     *
     * @param temperature  temperature in K
     * @param moleFractions fixed composition, length K
     * @param tolerance    convergence threshold on {@code sum |dGm/du|}
     */
    public Result solve(
            double temperature,
            double[] moleFractions,
            double tolerance,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink) {

        if (eventSink != null) {
            eventSink.accept(new ProgressEvent.EngineStart("CVM", 0));
        }

        int ncf = model.ncf();
        double[] u = initialU(moleFractions);
        double errf = 0;

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted()) {
                throw new CancellationException();
            }

            CVMGibbsModel.State state = model.at(temperature, moleFractions, u);

            // Early exit if the *current* point already has a non-positive
            // cluster variable -- mirrors the reference solver's cvMin<=0
            // check. This is NOT a failure: a cluster variable at or below zero
            // (e.g. a configuration genuinely disallowed at this
            // composition/order) means there is nowhere further for the Newton
            // step to usefully go without dividing by a near-zero probability,
            // so the reference accepts the current point as converged rather
            // than attempting a step.
            if (minClusterVariable(state) <= 0) {
                double[] gu = state.gmu();
                double errf0 = 0;
                for (double g : gu) errf0 += Math.abs(g);
                return new Result(state, u.clone(), true, its, errf0);
            }

            double[] Gu = state.gmu();
            errf = 0;
            for (double g : Gu) errf += Math.abs(g);

            if (eventSink != null) {
                eventSink.accept(new ProgressEvent.CvmIteration(
                        its, state.gm(), errf, state.hm(), state.sm(), u));
            }

            if (errf <= tolerance) {
                return new Result(state, u.clone(), true, its, errf);
            }

            try {
                double[] negGu = new double[ncf];
                for (int i = 0; i < ncf; i++) {
                    negGu[i] = -Gu[i];
                }

                double[] p = LinearAlgebra.solve(state.gmuu(), negGu);

                double errx = 0;
                for (double v : p) errx += Math.abs(v);

                // Only the cluster-variable-space clamp (stpmx) is applied,
                // directly on the unclamped u+p trial point -- see the class
                // note on why the reference's normalU pre-clamp is omitted.
                double[] uTrial = new double[ncf];
                for (int i = 0; i < ncf; i++) {
                    uTrial[i] = u[i] + p[i];
                }
                double alpha = stepLimit(u, uTrial, moleFractions);

                for (int i = 0; i < ncf; i++) {
                    u[i] += alpha * p[i];
                }

                if (errx <= TOLX) {
                    CVMGibbsModel.State finalState = model.at(temperature, moleFractions, u);
                    return new Result(finalState, u.clone(), true, its, errf);
                }

            } catch (CancellationException e) {
                // Cancellation must propagate. The original catch was on
                // Exception, which would have swallowed this and reported a
                // cancelled run as a non-converged result instead; in practice
                // the interrupt check at the top of the loop caught it first,
                // so the behaviour is unchanged, but the intent is now explicit.
                throw e;
            } catch (RuntimeException e) {
                // A singular Hessian (or any other failure of the linear
                // solve) ends the run at the last good point rather than
                // propagating: the caller checks converged and reports it.
                return new Result(state, u.clone(), false, its, errf);
            }
        }

        CVMGibbsModel.State finalState = model.at(temperature, moleFractions, u);
        return new Result(finalState, u.clone(), false, MAX_ITER, errf);
    }

    /**
     * Starting iterate: the correlation functions of the fully disordered
     * (random) state at this composition.
     */
    public double[] initialU(double[] moleFractions) {
        double[] full = model.geometry().basis
                .computeRandomCvcfCFs(moleFractions, model.geometry().pipelineResult);
        return Arrays.copyOf(full, model.ncf());
    }

    /**
     * Minimum cluster variable across all <em>non-point</em> cluster types --
     * port of the reference solver's {@code findMin}, which explicitly excludes
     * the last (point) type. The point type holds the mole fractions
     * themselves, not derived cluster probabilities, so including it would make
     * this check spuriously sensitive to composition rather than to how far the
     * solve is from a degenerate cluster configuration.
     *
     * <p>Deliberately narrower than {@code CVMGibbsModel.isValidParams}, which
     * checks the full range including the point block. Both test the same
     * physical condition; they differ in which cluster types are in scope, and
     * that difference is intentional on both sides.</p>
     */
    private double minClusterVariable(CVMGibbsModel.State state) {
        double[][][] cv = state.clusterVariables();
        CvmGeometry geo = model.geometry();
        double minCv = Double.POSITIVE_INFINITY;
        for (int t = 0; t < geo.tcdis - 1; t++) {
            double[][] tt = cv[t];
            if (tt == null) continue;
            for (double[] jj : tt) {
                if (jj == null) continue;
                for (double v : jj) minCv = Math.min(minCv, v);
            }
        }
        return minCv;
    }

    /**
     * Largest fraction of the step from {@code uOld} to {@code uTrial} that
     * keeps every cluster variable within {@code [0, 1]} -- port of the
     * reference solver's {@code stpmx(uold, unew)}.
     *
     * <p>Returns a full step when no cluster variable would leave the range;
     * otherwise {@code 0.1 * fmin}, backing well off the boundary rather than
     * landing exactly on it.</p>
     */
    public double stepLimit(double[] uOld, double[] uTrial, double[] moleFractions) {
        CvmGeometry geo = model.geometry();
        double fmin = 1.0;

        double[][][] cvOld = geo.evaluateCVs(uOld, moleFractions);
        double[][][] cvNew = geo.evaluateCVs(uTrial, moleFractions);

        for (int i = 0; i < geo.tcdis - 1; i++) {
            for (int j = 0; j < geo.lc[i]; j++) {
                for (int v = 0; v < geo.lcv[i][j]; v++) {
                    double vO = cvOld[i][j][v];
                    double vN = cvNew[i][j][v];
                    if (vN <= 0) {
                        fmin = Math.min(fmin, Math.abs(vO / (vN - vO)));
                    }
                    if (vN >= 1) {
                        fmin = Math.min(fmin, Math.abs((1.0 - vO) / (vN - vO)));
                    }
                }
            }
        }
        return (fmin >= 1.0) ? 1.0 : (0.1 * fmin);
    }

    /** The model this solver drives. */
    public CVMGibbsModel model() {
        return model;
    }
}
