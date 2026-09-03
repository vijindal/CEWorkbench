package org.ce.scratch;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.hamiltonian.CECEntry;

import java.util.Arrays;

/**
 * V2 STEP 2 regression gate for the generalised {@link HillertSolver.PhaseStep}
 * (M_A-based KKT system). The full {@link HillertSolver} is <b>not</b> run for
 * the ordered case -- it still rejects BCC_B2 at entry.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertGeneralizedPhaseStep
 * </pre>
 *
 * <h2>What is checked</h2>
 * <ul>
 *   <li><b>A -- V1 reduction.</b> For BCC_A2 (K=2,3,4) the generalised
 *       {@code PhaseStep} output must be <em>bit-identical</em> to the pre-V2
 *       system, which is re-assembled independently in this test (selector
 *       {@code J_M}, x-block simplex {@code c}, RHS {@code -Gu} / {@code e_{ncf+k}}).
 *       {@code M == x}, {@code J_M == selector} verified alongside.</li>
 *   <li><b>B -- BCC_B2 isolated KKT.</b> A valid ordered {@code State} is built
 *       with {@code atFullWide}; an independently assembled generalised KKT
 *       system {@code [[H, -c],[c^T, 0]]} with RHS {@code [-Gu ; 0]} (mu=0) and
 *       {@code [J_M[k,:] ; 0]} (unit mu) is solved and compared to production
 *       {@code PhaseStep.step}. The full linearised residual
 *       {@code H*deltaY - J_M^T*deltaMu - c*deltaLambda + Gu} is checked ~0.</li>
 *   <li><b>C -- simplex.</b> {@code c^T . deltaY == 0} for mu=0 and each unit mu,
 *       both phases.</li>
 *   <li><b>D -- M response.</b> {@code deltaM(mu) == J_M . deltaY(mu)} for mu=0
 *       and each unit mu (production {@code Step.deltaComponentAmountsAt} vs a
 *       direct contraction).</li>
 * </ul>
 */
public final class HillertGeneralizedPhaseStep {

    private static int failures = 0;
    private static final double T = 1000.0;
    private static final double EXACT = 0.0;
    private static final double KKT_TOL = 1.0e-7;   // independent-assembly vs production
    private static final double RESID_TOL = 1.0e-6; // linearised KKT residual

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(84));
        System.out.println("  V2 STEP 2 -- generalised PhaseStep (M_A-based KKT)");
        System.out.println("=".repeat(84));

        // A -- V1 reduction
        v1Reduction("Nb-Ti", "BCC_A2", "T", new double[] { 0.35, 0.65 });
        v1Reduction("Nb-Ti-V", "BCC_A2", "T", new double[] { 0.2, 0.3, 0.5 });
        v1Reduction("Nb-Ti-V-Zr", "BCC_A2", "T", new double[] { 0.1, 0.2, 0.3, 0.4 });

        // B/C/D -- BCC_B2 isolated
        bccB2Isolated();

        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(84));
        if (failures > 0) {
            throw new AssertionError(failures + " generalized-PhaseStep checks failed");
        }
    }

    // ==================================================================
    // A -- V1 reduction: production Step == independently re-assembled
    //      pre-V2 system, bit for bit
    // ==================================================================

    private static void v1Reduction(String elements, String structure, String model, double[] x) {
        System.out.printf("%n--- V1 reduction: %s / %s / %s  x=%s ---%n",
                elements, structure, model, Arrays.toString(x));
        CvmGeometry geo = CvmGeometry.build(elements, structure, model, null);
        CVMGibbsModel m = new CVMGibbsModel(geo, emptyEntry(elements, structure + "_" + model));
        int ncf = geo.ncf, K = geo.numComponents, w = ncf + K;

        double[] uFull = m.randomStateFull(x);
        // Build the reference state EXACTLY as PhaseStep.step does internally
        // (model.atFull(T, uFull)); the trailing-K slice of a random-state uFull
        // is the CVCF-transformed composition, ~1e-14 off the nominal x, so
        // re-slicing here is required for a bit-exact comparison.
        CVMGibbsModel.State st = m.atFull(T, uFull);

        // -- M == composition() (bit-exact selector), J_M == selector --
        double[] mA = st.componentAmountsPerFormulaUnit();
        double[] comp = st.composition();
        double worstMx = 0;
        for (int a = 0; a < K; a++) worstMx = Math.max(worstMx, Math.abs(mA[a] - comp[a]));
        check(structure + " K=" + K + ": M == composition() exactly", worstMx == EXACT,
                "worst=" + worstMx);

        double[][] jM = st.componentAmountsJacobian();
        boolean sel = jM.length == K;
        for (int a = 0; a < K && sel; a++)
            for (int j = 0; j < w; j++)
                if (jM[a][j] != ((j == ncf + a) ? 1.0 : 0.0)) { sel = false; break; }
        check(structure + " K=" + K + ": J_M == selector", sel, "jM=" + Arrays.deepToString(jM));

        // -- production step --
        HillertSolver.PhaseStep.Step prod = new HillertSolver.PhaseStep(m).step(uFull, T);

        // -- independently re-assemble the PRE-V2 system --
        double[] Gu = st.guFull();          // absolute dG/dY, width w
        double[][] H = st.gmuuFull();       // width w
        int n = w + 1;
        double[][] A = new double[n][n];
        for (int i = 0; i < w; i++) System.arraycopy(H[i], 0, A[i], 0, w);
        for (int i = ncf; i < w; i++) { A[i][w] = -1.0; A[w][i] = 1.0; }

        double[] b0 = new double[n];
        for (int i = 0; i < w; i++) b0[i] = -Gu[i];
        double[] sol0 = LinAlgSolve(A, b0);
        double[] refDeltaY0 = Arrays.copyOf(sol0, w);
        double[] refDeltaComp0 = Arrays.copyOfRange(sol0, ncf, w);

        double[][] refDeltaYSens = new double[K][];
        double[][] refDeltaCompSens = new double[K][];
        for (int k = 0; k < K; k++) {
            double[] ek = new double[n];
            ek[ncf + k] = 1.0;
            double[] solK = LinAlgSolve(A, ek);
            refDeltaYSens[k] = Arrays.copyOf(solK, w);
            refDeltaCompSens[k] = Arrays.copyOfRange(solK, ncf, w);
        }

        // -- compare bit-for-bit --
        check(structure + " K=" + K + ": deltaY0 identical to pre-V2",
                Arrays.equals(prod.deltaY0(), refDeltaY0),
                "max|d| = " + maxDiff(prod.deltaY0(), refDeltaY0));
        check(structure + " K=" + K + ": deltaComposition0 identical to pre-V2",
                Arrays.equals(prod.deltaComposition0(), refDeltaComp0),
                "max|d| = " + maxDiff(prod.deltaComposition0(), refDeltaComp0));
        boolean sensOk = true, compSensOk = true;
        for (int k = 0; k < K; k++) {
            if (!Arrays.equals(prod.deltaYSensitivity()[k], refDeltaYSens[k])) sensOk = false;
            if (!Arrays.equals(prod.deltaCompositionSensitivity()[k], refDeltaCompSens[k])) compSensOk = false;
        }
        check(structure + " K=" + K + ": every deltaYSensitivity column identical to pre-V2",
                sensOk, "column mismatch");
        check(structure + " K=" + K + ": every deltaCompositionSensitivity column identical to pre-V2",
                compSensOk, "column mismatch");
    }

    // ==================================================================
    // B / C / D -- BCC_B2 isolated
    // ==================================================================

    private static void bccB2Isolated() {
        System.out.printf("%n--- BCC_B2 isolated generalised KKT (HillertSolver NOT run) ---%n");
        CvmGeometry geo = CvmGeometry.build("Nb-Ti", "BCC_B2", "T", null);
        CVMGibbsModel m = new CVMGibbsModel(geo, emptyEntry("Nb-Ti", "BCC_B2_T"));
        int ncf = geo.ncf, K = geo.numComponents, w = geo.tcf;
        int idxEta = geo.basis.indexOfCf("eta"),
            idxXA = geo.basis.indexOfCf("xA"),
            idxXB = geo.basis.indexOfCf("xB");

        // A genuinely valid ordered state: random state at xA, then a small eta.
        // (The pair/triangle CFs are left at their random-state values, so only a
        // modest eta keeps every cluster variable inside (0,1) -- larger eta
        // needs a fully relaxed ordered state, which needs the solver this step
        // does not yet run. eta = 0.04 still fully exercises the J_M coupling and
        // the eta-column KKT row: M_A = xA - 0.02 != xA.)
        double xA = 0.45, xB = 0.55, eta = 0.04;
        double[] uFull = m.randomStateFull(new double[] { xA, xB });
        uFull[idxXA] = xA; uFull[idxXB] = xB; uFull[idxEta] = eta;

        CVMGibbsModel.State st = m.atFullWide(T, uFull);
        check("BCC_B2 seed state is valid (eta=" + eta + ")", st.isValidIncludingPoints(),
                "invalid seed");
        check("BCC_B2 M != x at this ordered state",
                Math.abs(st.componentAmountsPerFormulaUnit()[0] - xA) > 1e-6,
                "M_A - xA = " + (st.componentAmountsPerFormulaUnit()[0] - xA));

        // production
        HillertSolver.PhaseStep.Step prod = new HillertSolver.PhaseStep(m).step(uFull, T);

        // independent generalised KKT
        double[] Gu = st.guFull(w);
        double[][] H = st.gmuuFull(w);
        double[][] jM = st.componentAmountsJacobian();       // K x w
        double[] c = new double[w];
        for (int a = 0; a < K; a++) c[ncf + a] = 1.0;

        int n = w + 1;
        double[][] A = new double[n][n];
        for (int i = 0; i < w; i++) {
            System.arraycopy(H[i], 0, A[i], 0, w);
            A[i][w] = -c[i];
            A[w][i] = c[i];
        }
        double[] b0 = new double[n];
        for (int i = 0; i < w; i++) b0[i] = -Gu[i];
        double[] sol0 = LinAlgSolve(A, b0);
        double[] refDeltaY0 = Arrays.copyOf(sol0, w);
        double refLambda0 = sol0[w];

        double[][] refDeltaYSens = new double[K][];
        double[] refLambdaSens = new double[K];
        for (int k = 0; k < K; k++) {
            double[] ek = new double[n];
            for (int j = 0; j < w; j++) ek[j] = jM[k][j];
            double[] solK = LinAlgSolve(A, ek);
            refDeltaYSens[k] = Arrays.copyOf(solK, w);
            refLambdaSens[k] = solK[w];
        }

        // B: production == independent assembly
        check("BCC_B2: production deltaY0 matches independent generalised KKT",
                close(prod.deltaY0(), refDeltaY0, KKT_TOL),
                "max|d| = " + maxDiff(prod.deltaY0(), refDeltaY0));
        check("BCC_B2: production lambda0 matches",
                Math.abs(prod.lambda0() - refLambda0) < KKT_TOL,
                "|d| = " + Math.abs(prod.lambda0() - refLambda0));
        boolean sensOk = true;
        double worstSens = 0;
        for (int k = 0; k < K; k++) {
            worstSens = Math.max(worstSens, maxDiff(prod.deltaYSensitivity()[k], refDeltaYSens[k]));
            if (!close(prod.deltaYSensitivity()[k], refDeltaYSens[k], KKT_TOL)) sensOk = false;
        }
        check("BCC_B2: every deltaYSensitivity column matches independent KKT (worst "
                + String.format("%.2e", worstSens) + ")", sensOk, "mismatch");

        // B (residual): H*deltaY - J_M^T*deltaMu - c*deltaLambda + Gu ~ 0
        //   mu = 0:  deltaMu = 0
        residualCheck("BCC_B2 mu=0", H, jM, c, Gu, prod.deltaY0(), new double[K], prod.lambda0());
        for (int k = 0; k < K; k++) {
            double[] muUnit = new double[K];
            muUnit[k] = 1.0;
            double[] dY = prod.deltaYAt(muUnit);
            // deltaLambda at this mu = lambda0 + lambdaSensitivity[k]
            double dLam = prod.lambda0() + prod.lambdaSensitivity()[k];
            residualCheck("BCC_B2 mu=e" + k, H, jM, c, Gu, dY, muUnit, dLam);
        }

        // C: simplex c^T deltaY == 0
        check("BCC_B2: c^T deltaY0 == 0", Math.abs(dot(c, prod.deltaY0())) < 1e-9,
                "c^T dY0 = " + dot(c, prod.deltaY0()));
        for (int k = 0; k < K; k++) {
            double[] muUnit = new double[K];
            muUnit[k] = 1.0;
            double v = dot(c, prod.deltaYAt(muUnit));
            check("BCC_B2: c^T deltaY(mu=e" + k + ") == 0", Math.abs(v) < 1e-9, "= " + v);
        }
        // also for a generic mu
        double[] muGen = { 3000.0, -5000.0 };
        check("BCC_B2: c^T deltaY(mu=generic) == 0",
                Math.abs(dot(c, prod.deltaYAt(muGen))) < 1e-6,
                "= " + dot(c, prod.deltaYAt(muGen)));

        // D: deltaM(mu) == J_M . deltaY(mu)
        for (double[] mu : new double[][] { new double[K], { 1, 0 }, { 0, 1 }, { 3000, -5000 } }) {
            double[] dM = prod.deltaComponentAmountsAt(mu);
            double[] dY = prod.deltaYAt(mu);
            double[] jmdy = new double[K];
            for (int a = 0; a < K; a++)
                for (int j = 0; j < w; j++) jmdy[a] += jM[a][j] * dY[j];
            check("BCC_B2: deltaM(mu=" + Arrays.toString(mu) + ") == J_M . deltaY",
                    close(dM, jmdy, 1e-12), "max|d| = " + maxDiff(dM, jmdy));
        }

        // Affine property: deltaY(mu) linear in mu
        double[] mu1 = { 1000, 0 }, mu2 = { 0, 2000 };
        double[] sum = new double[K];
        for (int i = 0; i < K; i++) sum[i] = mu1[i] + mu2[i];
        double[] lhs = prod.deltaYAt(sum);
        double[] rhs = prod.deltaY0().clone();
        double[] a1 = prod.deltaYAt(mu1), a2 = prod.deltaYAt(mu2), z = prod.deltaY0();
        for (int i = 0; i < lhs.length; i++) rhs[i] = a1[i] + a2[i] - z[i];
        check("BCC_B2: deltaY is affine in mu (superposition)",
                close(lhs, rhs, 1e-6), "max|d| = " + maxDiff(lhs, rhs));
    }

    /** ||H*deltaY - J_M^T*deltaMu - c*deltaLambda + Gu||_inf ~ 0. */
    private static void residualCheck(String tag, double[][] H, double[][] jM, double[] c,
            double[] Gu, double[] deltaY, double[] mu, double deltaLambda) {
        int w = deltaY.length, K = mu.length;
        double worst = 0;
        for (int j = 0; j < w; j++) {
            double v = Gu[j] - c[j] * deltaLambda;
            for (int i = 0; i < w; i++) v += H[j][i] * deltaY[i];
            for (int a = 0; a < K; a++) v -= jM[a][j] * mu[a];
            worst = Math.max(worst, Math.abs(v));
        }
        check(tag + " : linearised KKT residual ~ 0 (worst " + String.format("%.2e", worst) + ")",
                worst < RESID_TOL, "worst=" + worst);
    }

    // ==================================================================
    // helpers
    // ==================================================================

    /**
     * Solve via the SAME production solver {@code PhaseStep} uses
     * ({@link LinearAlgebra#solveChecked}), so that for a V1 phase the only
     * possible difference between this independent re-assembly and production is
     * the matrix / RHS <em>construction</em> -- which must be identical, giving a
     * bit-exact match.
     */
    private static double[] LinAlgSolve(double[][] Ain, double[] bin) {
        return LinearAlgebra.solveChecked(Ain, bin).x();
    }

    private static double dot(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) s += a[i] * b[i];
        return s;
    }

    private static double maxDiff(double[] a, double[] b) {
        double m = 0;
        for (int i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
        return m;
    }

    private static boolean close(double[] a, double[] b, double tol) {
        return maxDiff(a, b) <= tol;
    }

    private static CECEntry emptyEntry(String elements, String structurePhase) {
        CECEntry e = new CECEntry();
        e.elements = elements;
        e.structurePhase = structurePhase;
        e.model = "T";
        e.cecTerms = new CECEntry.CECTerm[0];
        return e;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-76s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-76s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertGeneralizedPhaseStep() {
    }
}
