package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;

/**
 * Validates both CVM minimisers against a Mathematica {@code phaseq} reference
 * point for a ternary system.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.TernaryReferenceValidation
 * </pre>
 *
 * <p><b>Reference point.</b> Mo-Nb-Ta / BCC_A2 / T-model, 1000 K,
 * x = [0.33, 0.33, 0.34], supplied by the user from the Mathematica
 * implementation. Scalars are quoted to 6 significant figures there, so
 * agreement is asserted at {@code 1e-5} relative, not to machine precision.</p>
 *
 * <p><b>What makes this more than a self-check.</b> The same point is solved
 * twice, by two solvers that share only the evaluator:</p>
 * <ul>
 *   <li>{@link CvmNewtonSolver} holds composition fixed and solves an
 *       {@code ncf}-dimensional system against {@code gmu}/{@code gmuu}.</li>
 *   <li>{@link HillertSolver} (np=1) treats composition as an unknown, solving
 *       the widened {@code ncf+K} system against {@code gmuFull}/{@code gmuuFull}
 *       with a Lagrange constraint and an outer chemical-potential loop.</li>
 * </ul>
 * <p>Agreement between them, and with an external reference, is therefore
 * evidence about {@code CVMGibbsModel} itself rather than about one loop.</p>
 *
 * <p><b>The CF index permutation is expected, not a defect.</b> Our CVCF basis
 * order and the reference's {@code u2List} order differ by four transpositions
 * (see {@link #REF_TO_OURS}): {@code v3ABC1}/{@code v3ABC3}, and the whole
 * {@code v21}/{@code v22} 1NN/2NN pair block. Both differences are long-known
 * and were deliberately left unfixed -- our ordering is internally consistent
 * and only the labels disagree. This gate therefore compares under the
 * permutation; comparing position-by-position would report eight false
 * mismatches. If the permutation is ever removed, this table is what must
 * change, and the mismatch will be loud rather than silent.</p>
 *
 * <p><b>Chemical potential is deliberately not compared entry-by-entry.</b>
 * For a single phase the only constraint on mu is one Gibbs-Duhem equation, so
 * with K unknowns there is a (K-1)-dimensional family of valid solutions; ours
 * and the reference's are different points in it. What is asserted is the
 * constraint itself, {@code sum(mu_i * x_i) == G}. Note also that
 * {@code CvmNewtonSolver} has no chemical potential at all -- composition is
 * fixed, so mu is not part of its problem.</p>
 *
 * <p><b>Not covered.</b> This validates one interior composition. The
 * near-edge band (e.g. x = [0.05, 0.05, 0.90] on this same system) is where
 * both solvers are known to fail, and is a separate open item.</p>
 */
public final class TernaryReferenceValidation {

    // ---- Reference point ------------------------------------------------
    private static final String ELEMENTS  = "Mo-Nb-Ta";
    private static final String STRUCTURE = "BCC_A2";
    private static final String MODEL     = "T";
    private static final double T         = 1000.0;
    private static final double[] REF_X   = { 0.33, 0.33, 0.34 };

    private static final double REF_G  = -69246.3;
    private static final double REF_GM = -20633.7;
    private static final double REF_HM = -11972.7;
    private static final double REF_SM = 8.66101;

    /** Reference mu -- one valid solution of the underdetermined np=1 system. */
    private static final double[] REF_MU = { -73827.8, -59613.7, -74148.9 };

    /** Converged CFs in the reference's own {@code u2List} order. */
    private static final double[] REF_CF = {
            0.01314,     0.0205594,   0.00865153,  0.016275,    0.0134986,   0.016106,
            0.00523559,  0.00647258,  0.000333726, 0.0453497,   0.0403362,   0.0437667,
            0.116077,    0.132253,    0.105641,    0.11387,     0.119071,    0.108481 };

    /**
     * {@code REF_TO_OURS[i]} is the index in <em>our</em> CVCF order holding the
     * quantity the reference stores at its index {@code i}. Identity except for
     * the four known transpositions described in the class documentation.
     */
    private static final int[] REF_TO_OURS = {
            0, 1, 2, 3, 4, 5, 6, 7, 8,
            11,           // ref v3ABC1 -> our v3ABC3 slot
            10,
            9,            // ref v3ABC3 -> our v3ABC1 slot
            15, 16, 17,   // ref v22{AB,AC,BC} -> our v21 block
            12, 13, 14 }; // ref v21{AB,AC,BC} -> our v22 block

    /** Scalars are quoted to 6 significant figures in the reference. */
    private static final double SCALAR_TOL = 1e-5;
    private static final double CF_TOL     = 1e-4;

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.printf("  %s / %s / %s   T = %.0f K   x = %s%n",
                ELEMENTS, STRUCTURE, MODEL, T, fmt(REF_X));
        System.out.println("  Both minimisers vs Mathematica phaseq reference");
        System.out.println("=".repeat(78));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);

        List<String> names = CvCfBasis.getNonPointCfNames(STRUCTURE, MODEL, REF_X.length);

        // Separate model instances: the two solvers must not share state.
        CVMGibbsModel nrModel = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, session.cecEntry, null);
        CVMGibbsModel hiModel = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, session.cecEntry, null);

        // ---- Newton-Raphson, fixed composition --------------------------
        CvmNewtonSolver.Result nr = new CvmNewtonSolver(nrModel).solve(T, REF_X, 1e-5, null, null);
        System.out.printf("%n--- Newton-Raphson: converged=%s in %d iterations, |grad|=%.3e ---%n",
                nr.converged(), nr.iterations(), nr.finalGradientNorm());
        check("NR converged", nr.converged());
        if (nr.converged()) {
            compareState("NR", nr.state(), null, names);
        }

        // ---- Hillert, single phase, same seed ---------------------------
        List<HillertSolver.Phase> phases = new ArrayList<>();
        phases.add(new HillertSolver.Phase(
                "bcc", session, hiModel, 1.0, hiModel.randomStateFull(REF_X)));
        HillertSolver.Result eq = HillertSolver.solve(phases, T, 50, 10, 1.0e-6, null);
        System.out.printf("%n--- Hillert (np=1): converged=%s in %d outer iterations, residual=%.3e ---%n",
                eq.overallConverged(), eq.outerIterations(), eq.finalResidualNorm());
        check("Hillert converged", eq.overallConverged());
        if (eq.overallConverged()) {
            HillertSolver.PhaseResult p = eq.phases().get(0);
            compareState("Hillert", p.state(), eq.mu(), names);

            // With one phase there is nothing to trade mass with, so the
            // composition must not have drifted off the input.
            check("Hillert holds x (np=1)", maxAbsDiff(p.composition(), REF_X) < 1e-9);
            check("Hillert amount == 1", Math.abs(p.amount() - 1.0) < 1e-9);
        }

        // ---- The two solvers against each other -------------------------
        if (nr.converged() && eq.overallConverged()) {
            CVMGibbsModel.State a = nr.state();
            CVMGibbsModel.State b = eq.phases().get(0).state();
            System.out.printf("%n--- NR vs Hillert (independent routes to the same point) ---%n");
            agree("Gm",  a.gm(),  b.gm());
            agree("Hm",  a.hm(),  b.hm());
            agree("Sm",  a.sm(),  b.sm());
            agree("G0m", a.g0m(), b.g0m());
        }

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " reference validation checks failed");
        }
    }

    /** Compares one converged state against the reference. */
    private static void compareState(String solver, CVMGibbsModel.State st,
            double[] mu, List<String> names) {

        System.out.printf("%n  %s vs reference:%n", solver);
        rel(solver + " G",  st.g(),  REF_G);
        rel(solver + " Gm", st.gm(), REF_GM);
        rel(solver + " Hm", st.hm(), REF_HM);
        rel(solver + " Sm", st.sm(), REF_SM);

        // G must decompose exactly -- an internal identity, not a reference check.
        check(solver + " G == G0m + Gm",
                Math.abs(st.g() - (st.g0m() + st.gm())) < 1e-9);

        // CFs, compared under the known index permutation.
        double[] u = st.u();
        int bad = 0;
        for (int i = 0; i < REF_CF.length; i++) {
            int j = REF_TO_OURS[i];
            double r = Math.abs((u[j] - REF_CF[i]) / REF_CF[i]);
            if (r >= CF_TOL) {
                bad++;
                System.out.printf("    [!] CF %-8s ours[%2d]=%.9f  ref[%2d]=%.7f  rel=%.2e%n",
                        names.get(j), j, u[j], i, REF_CF[i], r);
            }
        }
        check(solver + " all 18 CFs match reference (under permutation)", bad == 0);
        System.out.printf("    %d/%d CFs match to %.0e relative%n",
                REF_CF.length - bad, REF_CF.length, CF_TOL);

        if (mu != null) {
            // mu is underdetermined for np=1: assert the constraint, not the vector.
            double dot = 0.0;
            for (int i = 0; i < mu.length; i++) dot += mu[i] * REF_X[i];
            System.out.printf("    mu           = %s%n", fmt(mu));
            System.out.printf("    (reference)  = %s%n", fmt(REF_MU));
            System.out.printf("    sum(mu*x)    = %.6f   vs G = %.6f%n", dot, st.g());
            check(solver + " Gibbs-Duhem: sum(mu*x) == G",
                    Math.abs(dot - st.g()) / Math.abs(st.g()) < 1e-9);
        }
    }

    private static void rel(String what, double ours, double ref) {
        double r = Math.abs((ours - ref) / ref);
        System.out.printf("    %-14s ours=%16.6f  ref=%12.4f  rel=%.2e%n", what, ours, ref, r);
        check(what + " matches reference", r < SCALAR_TOL);
    }

    private static void agree(String what, double a, double b) {
        double d = Math.abs(a - b);
        System.out.printf("    |%-4s(NR) - %-4s(Hillert)| = %.3e%n", what, what, d);
        check("NR and Hillert agree on " + what, d < 1e-6);
    }

    private static double maxAbsDiff(double[] a, double[] b) {
        double m = 0.0;
        for (int i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
        return m;
    }

    private static void check(String what, boolean ok) {
        if (!ok) {
            failures++;
            System.out.println("    [!] FAIL  " + what);
        }
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.5f", v[i]));
            if (i < v.length - 1) sb.append(", ");
        }
        return sb.append(']').toString();
    }

    private TernaryReferenceValidation() {
    }
}
