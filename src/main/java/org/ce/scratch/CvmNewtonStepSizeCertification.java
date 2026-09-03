package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

/**
 * Regression gate for the step-9 (raw-Newton-step) convergence exit of
 * {@link CvmNewtonSolver}.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.CvmNewtonStepSizeCertification
 * </pre>
 *
 * <h2>The flaw this guards against</h2>
 *
 * <p>The step-size branch used to be:</p>
 * <pre>
 *   if (errx &lt;= TOLX) {
 *       State finalState = model.at(T, x, u);   // u already updated in step 8
 *       return new Result(finalState, u, true, its, errf);   // errf is STALE
 *   }
 * </pre>
 *
 * <p>{@code errf} was the gradient L1 norm at the iterate <em>before</em> the
 * step-8 update, while {@code finalState} is evaluated <em>after</em> it. The
 * step may also have been clamped by {@code alpha}. So the returned
 * {@code Result} could report {@code converged = true} with a
 * {@code finalGradientNorm} that (a) belonged to a different point than
 * {@code state()}, and (b) was never checked against {@code tolerance} for the
 * point actually returned. A tiny raw Newton direction that got clamped hard by
 * a cluster-variable boundary could therefore be certified as a converged
 * minimum when the returned state is not stationary at all.</p>
 *
 * <h2>Two checks</h2>
 *
 * <ol>
 *   <li><b>Unit-level logic demonstration</b> ({@link #demonstrateLogic}) --
 *       a deterministic table of {@code (errx, finalErrf, tolerance)} triples
 *       driven through both the OLD and NEW decision rules, proving the NEW
 *       rule refuses to certify convergence when the post-update gradient is
 *       not actually small, while the OLD rule would have.</li>
 *   <li><b>Result-consistency invariant on real solves</b>
 *       ({@link #checkRealSolves}) -- for every {@code (T, x)} in a spread of
 *       physical points, at a range of tolerances, assert that any
 *       {@code converged == true} result satisfies
 *       {@code finalGradientNorm == sum |state().gmu()|} to floating-point
 *       tolerance, and {@code finalGradientNorm <= tolerance}. Under the OLD
 *       code this invariant is violated exactly when the step-9 branch fires
 *       with a clamped step; under the NEW code it holds by construction for
 *       both the step-5 and step-9 exits.</li>
 * </ol>
 */
public final class CvmNewtonStepSizeCertification {

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  CvmNewtonSolver step-size (step 9) convergence certification");
        System.out.println("=".repeat(78));

        demonstrateLogic();
        checkRealSolves();
        checkIterationCapResultConsistency();

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " step-size certification checks failed");
        }
    }

    // ------------------------------------------------------------------------
    // Check 1 -- the decision rule, in isolation
    // ------------------------------------------------------------------------

    /** OLD step-9 rule: a tiny raw step alone declares convergence. */
    private static boolean oldRuleConverged(double errx, double tolx) {
        return errx <= tolx;
    }

    /**
     * NEW step-9 rule: a tiny raw step declares convergence only if the
     * gradient L1 norm at the post-update point is itself within tolerance.
     */
    private static boolean newRuleConverged(double errx, double tolx,
                                            double finalErrf, double tolerance) {
        return errx <= tolx && finalErrf <= tolerance;
    }

    private static void demonstrateLogic() {
        System.out.println("\n--- Check 1: OLD vs NEW step-9 decision rule " + "-".repeat(32));

        final double TOLX = CvmNewtonSolver.TOLX;
        final double tolerance = 1.0e-6;

        // {errx, finalErrf, expectedOld, expectedNew}
        double[][] rows = {
                // tiny raw step, stationary after update  -> both converge
                { 1.0e-14, 1.0e-9,  1, 1 },
                // tiny raw step, NOT stationary after update (clamped step)
                //   -> OLD wrongly converges, NEW does not
                { 1.0e-14, 1.0e-2,  1, 0 },
                { 1.0e-13, 3.5e-4,  1, 0 },
                // raw step not tiny -> neither rule fires here
                { 1.0e-3,  1.0e-9,  0, 0 },
                { 1.0e-3,  1.0e-2,  0, 0 },
                // exactly on the TOLX boundary, stationary
                { TOLX,    1.0e-8,  1, 1 },
                // exactly on the TOLX boundary, not stationary
                { TOLX,    1.0e-1,  1, 0 },
        };

        boolean sawDivergence = false;
        for (double[] r : rows) {
            double errx = r[0];
            double finalErrf = r[1];
            boolean expOld = r[2] != 0;
            boolean expNew = r[3] != 0;

            boolean gotOld = oldRuleConverged(errx, TOLX);
            boolean gotNew = newRuleConverged(errx, TOLX, finalErrf, tolerance);

            check(String.format("errx=%.1e finalErrf=%.1e  OLD", errx, finalErrf),
                    gotOld == expOld, "expected " + expOld + " got " + gotOld);
            check(String.format("errx=%.1e finalErrf=%.1e  NEW", errx, finalErrf),
                    gotNew == expNew, "expected " + expNew + " got " + gotNew);

            if (gotOld && !gotNew) sawDivergence = true;
        }

        check("at least one row where OLD converges but NEW does not",
                sawDivergence, "table must exercise the behaviour change");
    }

    // ------------------------------------------------------------------------
    // Check 2 -- result consistency on real solves
    // ------------------------------------------------------------------------

    private static final String[][] CASES = {
            { "Nb-Ti", "BCC_A2", "T" },
            { "Nb-Ti-V", "BCC_A2", "T" },
            { "Nb-Ti-V-Zr", "BCC_A2", "T" },
    };

    private static final double[] TEMPERATURES = { 800.0, 1000.0, 1273.0 };

    /** A spread of tolerances, including some tight enough to reach step 9. */
    private static final double[] TOLERANCES = { 1.0e-5, 1.0e-7, 1.0e-9, 1.0e-11, 1.0e-13 };

    private static void checkRealSolves() throws Exception {
        System.out.println("\n--- Check 2: finalGradientNorm consistency on real solves " + "-".repeat(19));

        for (String[] c : CASES) {
            Workspace workspace = new Workspace();
            CEWorkbenchContext context = new CEWorkbenchContext(workspace);
            ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                    .build(new SystemId(c[0], c[1], c[2]), EngineConfig.CVM, null);

            CvmGeometry geo = CvmGeometry.build(c[0], c[1], c[2], null);
            CvmNewtonSolver solver = new CvmNewtonSolver(new CVMGibbsModel(geo, session.cecEntry));
            int K = geo.numComponents;

            System.out.printf("%n  %s / %s / %s%n", c[0], c[1], c[2]);

            for (double t : TEMPERATURES) {
                for (double[] x : compositions(K)) {
                    for (double tol : TOLERANCES) {
                        CvmNewtonSolver.Result r = solver.solve(t, x, tol, null, null);

                        if (!r.converged()) continue;

                        // The returned norm must be the norm of the returned state.
                        double[] gu = r.state().gmu();
                        double recomputed = 0;
                        for (double g : gu) recomputed += Math.abs(g);

                        double diff = Math.abs(recomputed - r.finalGradientNorm());
                        double scale = Math.max(recomputed, 1.0e-30);
                        boolean consistent = diff <= 1.0e-9 * scale + 1.0e-30;

                        check(String.format("T=%.0f x=%s tol=%.0e  norm matches state", t, fmt(x), tol),
                                consistent,
                                String.format("reported=%.6e recomputed=%.6e", r.finalGradientNorm(), recomputed));

                        // A converged result must actually be within tolerance.
                        // (Slack factor for the early-exit degeneracy branch,
                        // which legitimately returns the current point's true
                        // gradient without a tolerance test -- that norm is
                        // still self-consistent, only possibly above tol.)
                        boolean withinTol = r.finalGradientNorm() <= tol
                                || isDegenerateEarlyExit(solver, t, x, r);

                        check(String.format("T=%.0f x=%s tol=%.0e  norm <= tol (or degenerate exit)", t, fmt(x), tol),
                                withinTol,
                                String.format("finalGradientNorm=%.6e tol=%.6e", r.finalGradientNorm(), tol));
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // Check 3 -- the NON-converged iteration-cap Result is self-consistent
    // ------------------------------------------------------------------------

    /**
     * Nb-Ti-V / BCC_A2 / T at x=[0.05, 0.05, 0.90], 1000 K, tol=1e-10: a
     * documented near-edge point where the solver stagnates in a limit-cycle
     * and runs out at {@code MAX_ITER}. Deterministic -- uses the stored
     * {@code Nb-Ti-V_BCC_A2_T_CVCF} Hamiltonian, no timing, no external file
     * that other work modifies.
     *
     * <p>Before this fix the iteration-cap return reported the loop variable
     * {@code errf} (the gradient at the iterate <em>before</em> the last
     * step-8 update) as {@code finalGradientNorm()}, while {@code state()}
     * described the iterate <em>after</em> it -- the two named different
     * points. This asserts they now agree, on a genuinely non-converged
     * Result.</p>
     */
    private static void checkIterationCapResultConsistency() throws Exception {
        System.out.println("\n--- Check 3: iteration-cap (non-converged) Result consistency " + "-".repeat(15));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Ti-V", "BCC_A2", "T"), EngineConfig.CVM, null);
        CvmGeometry geo = CvmGeometry.build("Nb-Ti-V", "BCC_A2", "T", null);
        CvmNewtonSolver solver = new CvmNewtonSolver(new CVMGibbsModel(geo, session.cecEntry));

        double[] x = { 0.05, 0.05, 0.90 };
        double tol = 1.0e-10;
        CvmNewtonSolver.Result r = solver.solve(1000.0, x, tol, null, null);

        System.out.printf("    solve -> converged=%s iterations=%d finalGradientNorm=%.6e%n",
                r.converged(), r.iterations(), r.finalGradientNorm());

        check("Check 3: run did NOT converge (hit the cap as intended)",
                !r.converged(), "converged=" + r.converged());
        check("Check 3: iterations() == MAX_ITER",
                r.iterations() == CvmNewtonSolver.MAX_ITER,
                r.iterations() + " vs " + CvmNewtonSolver.MAX_ITER);

        // Independently recompute the gradient L1 norm at the RETURNED state.
        double[] gu = r.state().gmu();
        double recomputed = 0;
        for (double g : gu) recomputed += Math.abs(g);

        double diff = Math.abs(recomputed - r.finalGradientNorm());
        double scale = Math.max(recomputed, 1.0e-30);
        boolean consistent = diff <= 1.0e-9 * scale + 1.0e-30;
        check(String.format("Check 3: finalGradientNorm == sum|state().gmu()|  (reported=%.6e recomputed=%.6e)",
                        r.finalGradientNorm(), recomputed),
                consistent, "diff=" + diff);

        // And it must be a genuinely non-stationary point (above tol) -- i.e.
        // this is the real cap path, not an accidental late convergence.
        check("Check 3: the returned point is genuinely non-stationary (norm > tol)",
                r.finalGradientNorm() > tol,
                "finalGradientNorm=" + r.finalGradientNorm() + " tol=" + tol);
    }

    /**
     * The step-3 degeneracy exit returns {@code converged = true} with the
     * current point's genuine gradient norm and no tolerance test -- that is
     * intentional and unchanged by this fix. Distinguish it from a step-5/step-9
     * exit by checking whether the minimum non-point cluster variable at the
     * returned point is non-positive.
     */
    private static boolean isDegenerateEarlyExit(
            CvmNewtonSolver solver, double t, double[] x, CvmNewtonSolver.Result r) {
        double[][][] cv = r.state().clusterVariables();
        CvmGeometry geo = solver.model().geometry();
        double minCv = Double.POSITIVE_INFINITY;
        for (int tt = 0; tt < geo.tcdis - 1; tt++) {
            double[][] a = cv[tt];
            if (a == null) continue;
            for (double[] b : a) {
                if (b == null) continue;
                for (double v : b) minCv = Math.min(minCv, v);
            }
        }
        return minCv <= 0;
    }

    private static double[][] compositions(int k) {
        if (k == 2) {
            return new double[][] { { 0.5, 0.5 }, { 0.2, 0.8 }, { 0.05, 0.95 } };
        }
        if (k == 3) {
            return new double[][] {
                    { 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 0.2, 0.3, 0.5 }, { 0.1, 0.1, 0.8 } };
        }
        return new double[][] { { 0.25, 0.25, 0.25, 0.25 }, { 0.1, 0.2, 0.3, 0.4 } };
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-58s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-58s [!] FAIL  %s%n", label, detail);
        }
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.3f", v[i]));
            if (i < v.length - 1) sb.append(",");
        }
        return sb.append(']').toString();
    }

    private CvmNewtonStepSizeCertification() {
    }
}
