package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.cvm.CvmNewtonSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

/**
 * Parity gate for Step B: proves {@link CvmNewtonSolver} reproduces the
 * fixed-composition Newton-Raphson loop that {@code CVMGibbsModel} used to own,
 * frozen here as {@link PreFacadeCVMGibbsModel}.
 *
 * <p>The loop was moved verbatim -- same early exit on a non-positive cluster
 * variable, same {@code sum |dGm/du|} convergence test, same {@code stpmx} step
 * clamp, same {@code errx} check on the raw Newton step -- so the converged
 * point, the iteration count and the final gradient norm should all match
 * exactly, not merely closely. A drift in any of them means the algorithm
 * changed during the move, not just its location.</p>
 *
 * <p>Iteration count is the sharpest of the three: two solvers can land on the
 * same minimum by different paths, but taking the same number of steps to get
 * there means every intermediate iterate agreed too.</p>
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.CvmNewtonSolverParity
 * </pre>
 */
public final class CvmNewtonSolverParity {

    private static final String[][] CASES = {
            { "Nb-Ti", "BCC_A2", "T" },
            { "Nb-Ti-V", "BCC_A2", "T" },
            { "Nb-Ti-V-Zr", "BCC_A2", "T" },
    };

    private static final double[] TEMPERATURES = { 800.0, 1000.0, 1273.0 };

    /** The port is verbatim, so require agreement to double rounding. */
    private static final double TOL = 1.0e-12;

    private static final double SOLVER_TOL = 1.0e-5;

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  CvmNewtonSolver parity vs the pre-refactor minimisation loop");
        System.out.println("=".repeat(78));

        for (String[] c : CASES) {
            runCase(c[0], c[1], c[2]);
        }

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " solver parity checks failed");
        }
    }

    private static void runCase(String elements, String structure, String model) throws Exception {
        System.out.printf("%n--- %s / %s / %s %s%n", elements, structure, model, "-".repeat(24));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId(elements, structure, model), EngineConfig.CVM, null);

        PreFacadeCVMGibbsModel legacy = new PreFacadeCVMGibbsModel();
        legacy.initialize(elements, structure, model, session.cecEntry, null);

        CvmGeometry geo = CvmGeometry.build(elements, structure, model, null);
        CvmNewtonSolver solver = new CvmNewtonSolver(new CVMGibbsModel(geo, session.cecEntry));

        int K = geo.numComponents;

        for (double t : TEMPERATURES) {
            for (double[] x : compositions(K)) {
                PreFacadeCVMGibbsModel.EquilibriumResult expected =
                        legacy.getEquilibriumState(t, x, SOLVER_TOL, null, null, null);
                CvmNewtonSolver.Result got = solver.solve(t, x, SOLVER_TOL, null, null);

                System.out.printf("  T=%.0f x=%s%n", t, fmt(x));
                check("converged", expected.converged == got.converged(),
                        expected.converged + " vs " + got.converged());
                check("iterations", expected.iterations == got.iterations(),
                        expected.iterations + " vs " + got.iterations());
                checkScalar("gradNorm", expected.finalGradientNorm, got.finalGradientNorm());
                checkScalar("Gm", expected.modelResult.G, got.state().gm());
                checkScalar("Hm", expected.modelResult.H, got.state().hm());
                checkScalar("Sm", expected.modelResult.S, got.state().sm());
                checkVector("u", expected.u, got.u());
            }
        }
    }

    /**
     * A spread of compositions per system: equiatomic, plus two asymmetric
     * points. The asymmetric ones matter because the early-exit and step-clamp
     * branches only trigger when some cluster variable approaches zero, which
     * an equiatomic point rarely does.
     */
    private static double[][] compositions(int k) {
        if (k == 2) {
            return new double[][] { { 0.5, 0.5 }, { 0.2, 0.8 }, { 0.05, 0.95 } };
        }
        if (k == 3) {
            return new double[][] { { 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 0.2, 0.3, 0.5 }, { 0.05, 0.05, 0.90 } };
        }
        return new double[][] { { 0.25, 0.25, 0.25, 0.25 }, { 0.1, 0.2, 0.3, 0.4 } };
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-11s OK  (%s)%n", label, detail);
        } else {
            failures++;
            System.out.printf("    %-11s [!] FAIL  %s%n", label, detail);
        }
    }

    private static void checkScalar(String label, double expected, double got) {
        double rel = Math.abs(expected - got) / Math.max(Math.abs(expected), 1.0);
        if (rel <= TOL) {
            System.out.printf("    %-11s OK  (%.10g)%n", label, got);
        } else {
            failures++;
            System.out.printf("    %-11s [!] FAIL  legacy=%.12g new=%.12g rel=%.3e%n",
                    label, expected, got, rel);
        }
    }

    private static void checkVector(String label, double[] expected, double[] got) {
        int n = Math.min(expected.length, got.length);
        double worst = 0;
        int at = -1;
        for (int i = 0; i < n; i++) {
            double rel = Math.abs(expected[i] - got[i]) / Math.max(Math.abs(expected[i]), 1.0);
            if (rel > worst) {
                worst = rel;
                at = i;
            }
        }
        if (worst <= TOL) {
            System.out.printf("    %-11s OK  (%d entries, worst rel %.1e)%n", label, n, worst);
        } else {
            failures++;
            System.out.printf("    %-11s [!] FAIL  worst rel=%.3e at [%d] legacy=%.12g new=%.12g%n",
                    label, worst, at, expected[at], got[at]);
        }
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.4f", v[i]));
            if (i < v.length - 1) sb.append(", ");
        }
        return sb.append(']').toString();
    }

    private CvmNewtonSolverParity() {
    }
}
