package org.ce.model.equilibrium;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.ProgressEvent;

import java.util.Arrays;
import java.util.concurrent.CancellationException;
import java.util.function.Consumer;

/**
 * Newton-Raphson minimisation of the CVM Gibbs free energy at <b>fixed
 * composition</b>: given {@code (T, x)}, finds the {@code u} that makes
 * {@code dGm/du = 0}.
 *
 * <p>Holds a {@link CVMGibbsModel} and drives it from the outside, holding no
 * physics of its own -- every energy, gradient and Hessian comes from a
 * {@link CVMGibbsModel.State} evaluated at the current iterate. What lives here
 * is only the algorithm: convergence criteria, step limiting, and the iteration
 * bookkeeping.</p>
 *
 * <p><b>The companion solver.</b> {@link HillertPhaseStepSolver} sits beside
 * this one over the same model, differing only in which unknowns it solves for
 * and therefore which block of the same evaluated state it reads:</p>
 *
 * <table border="1">
 *   <caption>The two solvers over one model</caption>
 *   <tr><th></th><th>this</th><th>HillertPhaseStepSolver</th></tr>
 *   <tr><td>composition</td><td>fixed constraint</td><td>an unknown</td></tr>
 *   <tr><td>reads</td><td>{@code gmu} / {@code gmuu} ({@code ncf})</td>
 *       <td>{@code gmuFull} / {@code gmuuFull} ({@code ncf+K})</td></tr>
 *   <tr><td>solves for</td><td>stationary G at fixed x</td>
 *       <td>stationary G relative to a trial mu</td></tr>
 * </table>
 *
 * <h2>The nine steps</h2>
 *
 * <p>{@link #solve} implements this loop, ported from the reference solver
 * {@code CVMBINCE.minimize()} (which works in the <em>orthogonal</em> CF basis;
 * this one works in CVCF). Each step is marked in the code:</p>
 *
 * <ol>
 *   <li><b>Initialise</b> -- start from the fully disordered (random) state at
 *       this composition, {@link #initialU}. It is always physical, so the
 *       first iteration cannot begin outside the valid region.</li>
 *   <li><b>Evaluate</b> -- ask the model for the state at the current iterate.
 *       Everything the rest of the iteration needs comes from this one object.</li>
 *   <li><b>Degeneracy check</b> -- if any cluster variable is already
 *       non-positive, stop and accept the current point. Not a failure; see
 *       {@link #minClusterVariable}.</li>
 *   <li><b>Gradient</b> -- read {@code dGm/du} and form the convergence measure
 *       {@code errf = sum |dGm/du|} (an L1 norm, matching the reference; not
 *       the L2 norm).</li>
 *   <li><b>Convergence test</b> -- if {@code errf <= tolerance}, the point is
 *       stationary and the solve is done.</li>
 *   <li><b>Newton direction</b> -- solve {@code Guu p = -Gu} for the step
 *       {@code p}. This is the only linear algebra in the loop, and the
 *       dominant cost at {@code O(ncf^3)}.</li>
 *   <li><b>Step limit</b> -- find the largest fraction {@code alpha} of
 *       {@code p} that keeps every cluster variable inside {@code [0, 1]}
 *       ({@link #stepLimit}, the reference's {@code stpmx}). An unclamped
 *       Newton step can easily leave the physical region.</li>
 *   <li><b>Update</b> -- move the iterate, {@code u += alpha * p}.</li>
 *   <li><b>Step-size test</b> -- if the <em>raw</em> step {@code p} was already
 *       negligible ({@code sum |p| <= TOLX}), converge. Deliberately measured
 *       on {@code p} and not on {@code alpha * p}: a tiny <em>clamped</em> step
 *       means the boundary is in the way, which is not convergence.</li>
 * </ol>
 *
 * <p>The reference's additional pre-clamp on the trial correlation functions
 * themselves ({@code Utils.normalU}, restricting them to a fixed
 * {@code [-1, 1]}) was tried and deliberately dropped: that bound is meaningful
 * for the orthogonal basis's Chebyshev-like CFs, but this solver works in the
 * CVCF basis, whose {@code u} components have a different, non-uniform natural
 * range (observed roughly 1e-4 to 0.2 in practice) -- a blind &plusmn;1 clamp
 * there is a no-op, not a safeguard, and did not change behaviour when tested.
 * The near-edge convergence stall this porting effort was investigating (a
 * dilute-composition Newton direction that repeatedly re-approaches, but never
 * crosses, a cluster-variable boundary) remains open; see CLAUDE.md's note on
 * near-edge ternary solver fragility.</p>
 *
 * <p>Stateless and safe to share: each {@link #solve} call carries its own
 * iterate, so one instance may serve many points. Memoising a result for a
 * repeated {@code (T, x)} is a caller's concern -- {@code ThermodynamicWorkflow}
 * does it -- not part of the algorithm.</p>
 */
public final class CvmNewtonSolver {

    /** Iteration cap. A run that hits this is reported as non-converged. */
    public static final int MAX_ITER = 100;

    /**
     * Convergence threshold on the raw (unclamped) Newton step -- step 9. A
     * tiny step is a genuine sign of convergence, unlike a small step produced
     * purely by the boundary clamp of step 7.
     */
    public static final double TOLX = 1.0e-12;

    private final CVMGibbsModel model;

    public CvmNewtonSolver(CVMGibbsModel model) {
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        this.model = model;
    }

    /**
     * Outcome of a minimisation: the converged point and how it was reached.
     *
     * <p>{@link #state} is the model already evaluated at the converged
     * {@link #u} -- the solver built it on the final iteration and hands it
     * back rather than making the caller repeat the work. Reading a property
     * off it needs no further evaluation:</p>
     *
     * <pre>
     *   Result eq = new CvmNewtonSolver(model).solve(T, x, tol, null, null);
     *   double g = eq.state().gm();          // no re-evaluation
     *   // equivalently, but wastefully: model.at(T, x, eq.u()).gm()
     * </pre>
     *
     * <p><b>Always check {@link #converged} first.</b> A non-converged CVM run
     * still returns plausible-looking numbers; the flag is the only thing that
     * distinguishes them.</p>
     */
    public record Result(
            CVMGibbsModel.State state,
            double[] u,
            boolean converged,
            int iterations,
            double finalGradientNorm) {

        /**
         * The model this result was produced by, taken from the state rather
         * than stored separately so the two cannot disagree.
         */
        public CVMGibbsModel model() {
            return state.model();
        }
    }

    /**
     * Minimises {@code Gm} over {@code u} at fixed {@code (T, x)}. See the
     * class documentation for the nine steps this implements; they are marked
     * inline below.
     *
     * @param temperature   temperature in K
     * @param moleFractions fixed composition, length K
     * @param tolerance     convergence threshold on {@code sum |dGm/du|}
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

        // ── Step 1 ── Initialise the iterate at the fully disordered state.
        int ncf = model.ncf();
        double[] u = initialU(moleFractions);
        double errf = 0;

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted()) {
                throw new CancellationException();
            }

            // ── Step 2 ── Evaluate the model at the current iterate. Every
            // quantity the rest of this iteration needs comes from this state.
            CVMGibbsModel.State state = model.at(temperature, moleFractions, u);

            // ── Step 3 ── Degeneracy check on the *current* point (not on a
            // trial step). A cluster variable at or below zero means a
            // configuration that is genuinely disallowed at this
            // composition/order, so there is nowhere further for a Newton step
            // to go without dividing by a near-zero probability. The reference
            // solver accepts the current point as converged rather than
            // attempting a step, and so does this.
            if (minClusterVariable(state) <= 0) {
                double[] gu = state.gmu();
                double errf0 = 0;
                for (double g : gu) errf0 += Math.abs(g);
                return new Result(state, u.clone(), true, its, errf0);
            }

            // ── Step 4 ── Gradient and its L1 norm. The reference uses
            // sum |dGm/du|, not the Euclidean norm; kept as-is so the
            // convergence threshold means the same thing here as there.
            double[] Gu = state.gmu();
            errf = 0;
            for (double g : Gu) errf += Math.abs(g);

            if (eventSink != null) {
                eventSink.accept(new ProgressEvent.CvmIteration(
                        its, state.gm(), errf, state.hm(), state.sm(), u));
            }

            // ── Step 5 ── Convergence: the point is stationary.
            if (errf <= tolerance) {
                return new Result(state, u.clone(), true, its, errf);
            }

            try {
                // ── Step 6 ── Newton direction: solve Guu p = -Gu. The only
                // linear algebra in the loop, and its dominant cost.
                double[] negGu = new double[ncf];
                for (int i = 0; i < ncf; i++) {
                    negGu[i] = -Gu[i];
                }
                double[] p = LinearAlgebra.solve(state.gmuu(), negGu);

                double errx = 0;
                for (double v : p) errx += Math.abs(v);

                // ── Step 7 ── Step limit. Only the cluster-variable-space
                // clamp (stpmx) is applied, directly on the unclamped u+p trial
                // point -- see the class note on why the reference's normalU
                // pre-clamp is omitted for the CVCF basis.
                double[] uTrial = new double[ncf];
                for (int i = 0; i < ncf; i++) {
                    uTrial[i] = u[i] + p[i];
                }
                double alpha = stepLimit(u, uTrial, moleFractions);

                // ── Step 8 ── Update the iterate.
                for (int i = 0; i < ncf; i++) {
                    u[i] += alpha * p[i];
                }

                // ── Step 9 ── Step-size convergence, measured on the *raw*
                // step p rather than the clamped alpha*p: a small clamped step
                // means the boundary is in the way, which is not convergence.
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
                // A singular Hessian (or any other failure of the linear solve)
                // ends the run at the last good point rather than propagating:
                // the caller checks converged and reports it.
                return new Result(state, u.clone(), false, its, errf);
            }
        }

        // Iteration cap reached without satisfying step 5 or step 9.
        CVMGibbsModel.State finalState = model.at(temperature, moleFractions, u);
        return new Result(finalState, u.clone(), false, MAX_ITER, errf);
    }

    /**
     * <b>Step 1.</b> Starting iterate: the correlation functions of the fully
     * disordered (random) state at this composition.
     *
     * <p>Always physical -- every cluster probability is a product of mole
     * fractions and so lies strictly inside {@code (0, 1)} for any real
     * mixture -- so the first iteration cannot begin outside the valid
     * region.</p>
     */
    public double[] initialU(double[] moleFractions) {
        double[] full = model.geometry().basis
                .computeRandomCvcfCFs(moleFractions, model.geometry().pipelineResult);
        return Arrays.copyOf(full, model.ncf());
    }

    /**
     * <b>Step 3.</b> Minimum cluster variable across all <em>non-point</em>
     * cluster types -- port of the reference solver's {@code findMin}, which
     * explicitly excludes the last (point) type.
     *
     * <p>The exclusion matters: the point type holds the mole fractions
     * themselves, not derived cluster probabilities, so including it would make
     * this check sensitive to how dilute the composition is rather than to how
     * close the <em>solve</em> is to a degenerate cluster configuration. At
     * {@code x = 0.05} the point block is legitimately small while the solve is
     * perfectly healthy.</p>
     *
     * <p>Deliberately narrower than {@link CVMGibbsModel.State#isValidIncludingPoints},
     * which spans the full range and which the Hillert backtracking uses
     * precisely because there composition <em>is</em> an unknown that can drift
     * to a pure-element boundary. Both test the same physical condition; the
     * difference in scope is intentional on both sides.</p>
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
     * <b>Step 7.</b> Largest fraction of the step from {@code uOld} to
     * {@code uTrial} that keeps every cluster variable within {@code [0, 1]} --
     * port of the reference solver's {@code stpmx(uold, unew)}.
     *
     * <p>Cluster variables are linear in {@code u}, so for each one the
     * crossing fraction is exact: if the trial value {@code vN} leaves the
     * range, the fraction that lands exactly on the boundary is
     * {@code vO / (vO - vN)} for the lower bound and
     * {@code (1 - vO) / (vN - vO)} for the upper. The smallest such fraction
     * over all cluster variables is the first boundary the step would hit.</p>
     *
     * <p>Returns a full step when nothing would leave the range; otherwise
     * {@code 0.1 * fmin}, backing well off the boundary rather than landing on
     * it -- a point exactly at zero would make the next iteration's entropy
     * gradient singular.</p>
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
