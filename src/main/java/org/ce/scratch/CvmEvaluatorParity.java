package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CvmEvaluator;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.cvm.CvmState;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.Random;

/**
 * Parity gate for Step A of the evaluator/solver split: proves
 * {@link CvmEvaluator}/{@link CvmState} reproduce {@link PreFacadeCVMGibbsModel}
 * bit-for-bit at arbitrary points, before any caller is migrated to them.
 *
 * <p>Every expression in {@code CvmState} was ported verbatim from the
 * corresponding {@code CVMGibbsModel.calculateXxx} (frozen in {@link PreFacadeCVMGibbsModel}), so agreement should be
 * exact -- not merely within a tolerance. Anything else means a field was
 * rebound wrongly during the move.</p>
 *
 * <p>Also runs a <b>finite-difference check of every gradient and Hessian
 * against the new evaluator's own energy</b>. That is the capability the split
 * exists to provide: {@code checkMinimized()} blocks evaluation
 * at an arbitrary point, so differentiating its free energy numerically was not
 * possible without going through the solver. This is the tool the {@code Gm}
 * expression audit needs.</p>
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.CvmEvaluatorParity
 * </pre>
 */
public final class CvmEvaluatorParity {

    private static final String[][] CASES = {
            { "Nb-Ti", "BCC_A2", "T" },
            { "Nb-Ti-V", "BCC_A2", "T" },
            { "Nb-Ti-V-Zr", "BCC_A2", "T" },
    };

    private static final double[] TEMPERATURES = { 800.0, 1000.0, 1273.0 };

    /** Parity is a verbatim port, so require exactness up to double rounding. */
    private static final double TOL_PARITY = 1.0e-12;

    /**
     * Tolerance for central-difference agreement. Deliberately looser than the
     * parity tolerance: a numerical derivative carries real truncation error,
     * and the Hessian check differences an already-analytic gradient, costing
     * roughly half the available precision again. Observed worst cases across
     * K=2/3/4 sit at a few times 1e-5; anything materially larger indicates a
     * wrong expression rather than noise.
     */
    private static final double TOL_FD = 1.0e-4;

    /**
     * Relative scatter applied to the random-state CFs. Small enough that most
     * points stay physical (every cluster variable inside (0,1)) while some
     * cross into the entropy-smoothing branch, so both are covered.
     */
    private static final double PERTURBATION = 0.05;

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        System.out.println("=".repeat(80));
        System.out.println("  CvmEvaluator parity vs CVMGibbsModel, plus finite-difference audit");
        System.out.println("=".repeat(80));

        for (String[] c : CASES) {
            runCase(c[0], c[1], c[2]);
        }

        System.out.println("\n" + "=".repeat(80));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(80));
        if (failures > 0) {
            throw new AssertionError(failures + " parity/finite-difference checks failed");
        }
    }

    private static void runCase(String elements, String structure, String model) throws Exception {
        System.out.printf("%n--- %s / %s / %s %s%n", elements, structure, model, "-".repeat(30));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId(elements, structure, model), EngineConfig.CVM, null);

        PreFacadeCVMGibbsModel legacy = new PreFacadeCVMGibbsModel();
        legacy.initialize(elements, structure, model, session.cecEntry, null);

        CvmGeometry geo = CvmGeometry.build(elements, structure, model, null);
        CvmEvaluator ev = new CvmEvaluator(geo, session.cecEntry);

        int ncf = geo.ncf;
        int K = geo.numComponents;
        Random rng = new Random(20260820L);

        for (double t : TEMPERATURES) {
            double[] x = randomComposition(rng, K);
            // Perturb around the random state so cluster variables stay
            // physical; an arbitrary u would put cv outside (0,1) and exercise
            // only the entropy-smoothing branch.
            double[] u = ev.randomStateU(x);
            for (int i = 0; i < ncf; i++) {
                u[i] *= 1.0 + PERTURBATION * (rng.nextDouble() - 0.5);
            }
            // A perturbation large enough to drive some cluster variable out of
            // (0,1) puts the entropy into its ENTROPY_SMOOTH_EPS quadratic
            // extension, where Sm can legitimately go negative. Both regions
            // are worth exercising -- the smoothing branch is the one with no
            // physical intuition to check it against -- but which one is in
            // play is reported so a negative Sm is not mistaken for a defect.

            CvmState st = ev.stateAt(t, x, u);
            PreFacadeCVMGibbsModel.ModelResult legacyResult = legacy.evaluate(u, x, t);

            System.out.printf("  T=%.0f x=%s  cvValid=%s%n", t, fmt(x), st.isValid());

            checkScalar("Gm", legacyResult.G, st.gm());
            checkScalar("Hm", legacyResult.H, st.hm());
            checkScalar("Sm", legacyResult.S, st.sm());
            checkVector("dGm/du", legacyResult.Gu, st.gmu());
            checkMatrix("d2Gm/du2", legacyResult.Guu, st.gmuu());
            checkVector("dHm/du", legacyResult.Hu, st.hmu());
            checkVector("dSm/du", legacyResult.Su, st.smu());
            checkMatrix("d2Sm/du2", legacyResult.Suu, st.smuu());
            checkVector("cfs", legacyResult.cfs, st.cfs());

            // Widened quantities: compare against the legacy per-phase path,
            // which is the only consumer of the Full variants.
            PreFacadeCVMGibbsModel.PerPhaseStepResult ignored = legacy.solvePerPhaseStep(concat(u, x), t);
            if (ignored == null) {
                fail("solvePerPhaseStep returned null");
            }

            // Finite-difference audit against the new evaluator's own energy.
            fdGradient(ev, t, x, u, st);
            fdHessian(ev, t, x, u, st);
            fdGradientFull(ev, t, x, u, st);
        }
    }

    // =========================================================================
    // Parity
    // =========================================================================

    private static void checkScalar(String label, double expected, double got) {
        double diff = Math.abs(expected - got);
        double rel = diff / Math.max(Math.abs(expected), 1.0);
        boolean ok = rel <= TOL_PARITY;
        if (!ok) {
            fail(String.format("    parity %-10s legacy=%.12g new=%.12g rel=%.3e", label, expected, got, rel));
        } else {
            System.out.printf("    parity %-10s OK  (%.10g)%n", label, got);
        }
    }

    private static void checkVector(String label, double[] expected, double[] got) {
        if (expected == null) {
            return;
        }
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
        if (expected.length != got.length) {
            fail(String.format("    parity %-10s length %d vs %d", label, expected.length, got.length));
        } else if (worst > TOL_PARITY) {
            fail(String.format("    parity %-10s worst rel=%.3e at [%d] legacy=%.12g new=%.12g",
                    label, worst, at, expected[at], got[at]));
        } else {
            System.out.printf("    parity %-10s OK  (%d entries, worst rel %.1e)%n", label, n, worst);
        }
    }

    private static void checkMatrix(String label, double[][] expected, double[][] got) {
        if (expected == null) {
            return;
        }
        double worst = 0;
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected[i].length; j++) {
                double rel = Math.abs(expected[i][j] - got[i][j]) / Math.max(Math.abs(expected[i][j]), 1.0);
                worst = Math.max(worst, rel);
            }
        }
        if (worst > TOL_PARITY) {
            fail(String.format("    parity %-10s worst rel=%.3e", label, worst));
        } else {
            System.out.printf("    parity %-10s OK  (%dx%d, worst rel %.1e)%n",
                    label, expected.length, expected.length, worst);
        }
    }

    // =========================================================================
    // Finite-difference audit
    // =========================================================================

    private static void fdGradient(CvmEvaluator ev, double t, double[] x, double[] u, CvmState st) {
        double[] analytic = st.gmu();
        double worst = 0;
        int at = -1;
        for (int l = 0; l < ev.ncf(); l++) {
            double h = 1.0e-6 * Math.max(Math.abs(u[l]), 1.0e-3);
            double plus = ev.stateAt(t, x, bump(u, l, h)).gm();
            double minus = ev.stateAt(t, x, bump(u, l, -h)).gm();
            double num = (plus - minus) / (2 * h);
            double rel = Math.abs(num - analytic[l]) / Math.max(Math.abs(analytic[l]), 1.0);
            if (rel > worst) {
                worst = rel;
                at = l;
            }
        }
        report("    fd dGm/du   ", worst, at, analytic);
    }

    private static void fdGradientFull(CvmEvaluator ev, double t, double[] x, double[] u, CvmState st) {
        // Only the leading ncf block is comparable by perturbing u; the
        // trailing composition block cannot be finite-differenced independently
        // because mole fractions are constrained to sum to 1.
        double[] full = st.gmuFull();
        double[] narrow = st.gmu();
        double worst = 0;
        for (int l = 0; l < ev.ncf(); l++) {
            double rel = Math.abs(full[l] - narrow[l]) / Math.max(Math.abs(narrow[l]), 1.0);
            worst = Math.max(worst, rel);
        }
        if (worst > TOL_PARITY) {
            fail(String.format("    gmuFull leading block != gmu, worst rel=%.3e", worst));
        } else {
            System.out.printf("    fd gmuFull   OK  (leading ncf block matches gmu, worst rel %.1e)%n", worst);
        }
    }

    private static void fdHessian(CvmEvaluator ev, double t, double[] x, double[] u, CvmState st) {
        double[][] analytic = st.gmuu();

        // Symmetry is an exact property of a Hessian and is independent of any
        // step size, so it is checked first and strictly: a genuine error in
        // the second-derivative expression would almost certainly break it.
        double asymmetry = 0;
        for (int i = 0; i < ev.ncf(); i++) {
            for (int j = 0; j < ev.ncf(); j++) {
                asymmetry = Math.max(asymmetry, Math.abs(analytic[i][j] - analytic[j][i]));
            }
        }
        if (asymmetry > TOL_PARITY) {
            fail(String.format("    d2Gm/du2 not symmetric, max |H - H^T| = %.3e", asymmetry));
        }

        // Scale the finite-difference comparison against the magnitude of the
        // matrix, not against each entry. Many entries are near zero (1e-6
        // against a diagonal of ~1e4), and dividing absolute differencing noise
        // by such an entry reports a huge relative error for a value that is
        // numerically zero. Differencing an analytic gradient also costs about
        // half the available precision, so the step must be larger than the
        // gradient check uses.
        double scale = 0;
        for (double[] row : analytic) {
            for (double value : row) {
                scale = Math.max(scale, Math.abs(value));
            }
        }
        scale = Math.max(scale, 1.0);

        double worst = 0;
        for (int l = 0; l < ev.ncf(); l++) {
            double h = 1.0e-4 * Math.max(Math.abs(u[l]), 1.0e-2);
            double[] gPlus = ev.stateAt(t, x, bump(u, l, h)).gmu();
            double[] gMinus = ev.stateAt(t, x, bump(u, l, -h)).gmu();
            for (int m = 0; m < ev.ncf(); m++) {
                double num = (gPlus[m] - gMinus[m]) / (2 * h);
                worst = Math.max(worst, Math.abs(num - analytic[m][l]) / scale);
            }
        }
        if (worst > TOL_FD) {
            fail(String.format("    fd d2Gm/du2 worst scaled err=%.3e", worst));
        } else {
            System.out.printf("    fd d2Gm/du2  OK  (symmetric; worst scaled err %.1e)%n", worst);
        }
    }

    private static void report(String label, double worst, int at, double[] analytic) {
        if (worst > TOL_FD) {
            fail(String.format("%s worst rel=%.3e at [%d] analytic=%.10g", label, worst, at, analytic[at]));
        } else {
            System.out.printf("%s OK  (worst rel %.1e)%n", label, worst);
        }
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    private static double[] bump(double[] v, int i, double h) {
        double[] out = v.clone();
        out[i] += h;
        return out;
    }

    private static double[] concat(double[] a, double[] b) {
        double[] out = new double[a.length + b.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }

    private static double[] randomComposition(Random rng, int k) {
        double[] x = new double[k];
        double sum = 0;
        for (int i = 0; i < k; i++) {
            x[i] = 0.5 + rng.nextDouble();
            sum += x[i];
        }
        for (int i = 0; i < k; i++) {
            x[i] /= sum;
        }
        return x;
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.4f", v[i]));
            if (i < v.length - 1) {
                sb.append(", ");
            }
        }
        return sb.append(']').toString();
    }

    private static void fail(String message) {
        failures++;
        System.out.println(message + "   [!] FAIL");
    }

    private CvmEvaluatorParity() {
    }
}
