package org.ce.scratch;

import org.ce.model.equilibrium.HillertSolver.EquilibriumMatrix;
import org.ce.model.equilibrium.HillertSolver.EquilibriumMatrix.EquilibriumStepResult;
import org.ce.model.equilibrium.HillertSolver.EquilibriumMatrix.PhaseContribution;

import java.util.List;

/**
 * V2 STEP 2 regression gate for the generalised {@link EquilibriumMatrix}
 * (Gibbs-Duhem + mass-balance rows written on {@code M_A} and its affine
 * response {@code deltaM(mu)}, not on the raw composition coordinate).
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertGeneralizedEquilibriumMatrix
 * </pre>
 *
 * <p>No CVM model and no {@link org.ce.model.equilibrium.HillertSolver} run:
 * {@code PhaseContribution} records are supplied directly with hand-chosen
 * numbers, and the assembled {@code (K+np)} system is re-derived independently
 * and compared.</p>
 *
 * <h2>Cases</h2>
 * <ul>
 *   <li><b>E -- V1 reduction.</b> With {@code m == x} and
 *       {@code deltaM == deltaComposition} the generalised assembly must equal
 *       the pre-V2 assembly (re-implemented here verbatim from the old code).
 *       Checked for np=1 and np=2, with and without a target.</li>
 *   <li><b>F -- synthetic nontrivial M.</b> A two-phase, K=2 problem where
 *       {@code M != x} (an ordered-like phase: {@code M_A = xA - eta/2}). The
 *       Gibbs-Duhem rows, mass-balance rows (matrix + RHS), {@code deltaN}
 *       coefficients and target residual are each computed independently and
 *       compared to what {@code EquilibriumMatrix.solve} produced -- via the
 *       solved {@code (mu, deltaN)} satisfying the independently assembled
 *       system to tight precision.</li>
 * </ul>
 */
public final class HillertGeneralizedEquilibriumMatrix {

    private static int failures = 0;
    private static final double TOL = 1.0e-7;

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(84));
        System.out.println("  V2 STEP 2 -- generalised EquilibriumMatrix");
        System.out.println("=".repeat(84));

        v1ReductionNp1();
        v1ReductionNp2();
        syntheticNontrivialM();

        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(84));
        if (failures > 0) {
            throw new AssertionError(failures + " generalized-EquilibriumMatrix checks failed");
        }
    }

    // ==================================================================
    // E -- V1 reduction (M == x, deltaM == deltaComposition)
    // ==================================================================

    private static void v1ReductionNp1() {
        System.out.printf("%n--- E: V1 reduction, np=1, K=2 ---%n");
        int K = 2;
        // one phase; m == x
        double[] x = { 0.4, 0.6 };
        double[] dC0 = { 0.01, -0.01 };
        double[][] dCs = { { 0.002, -0.002 }, { -0.0015, 0.0015 } };
        PhaseContribution pc = new PhaseContribution(1.0, x, -5000.0, dC0, dCs);

        double[] target = { 0.4, 0.6 };
        compareToPreV2("np1 with target", List.of(pc), K, target);
        compareToPreV2("np1 no target", List.of(pc), K, null);
    }

    private static void v1ReductionNp2() {
        System.out.printf("%n--- E: V1 reduction, np=2, K=3 ---%n");
        int K = 3;
        PhaseContribution a = new PhaseContribution(
                0.6, new double[] { 0.5, 0.3, 0.2 }, -7000.0,
                new double[] { 0.01, -0.005, -0.005 },
                new double[][] {
                        { 0.003, -0.001, -0.002 },
                        { -0.001, 0.004, -0.003 },
                        { -0.002, -0.002, 0.004 } });
        PhaseContribution b = new PhaseContribution(
                0.4, new double[] { 0.2, 0.5, 0.3 }, -6500.0,
                new double[] { -0.008, 0.004, 0.004 },
                new double[][] {
                        { 0.002, -0.0005, -0.0015 },
                        { -0.0005, 0.0025, -0.002 },
                        { -0.0015, -0.002, 0.0035 } });
        double[] target = { 0.38, 0.38, 0.24 }; // = 0.6*a + 0.4*b
        compareToPreV2("np2 with target", List.of(a, b), K, target);
        compareToPreV2("np2 no target", List.of(a, b), K, null);
        // also perturb the target so r_i != 0
        compareToPreV2("np2 drifted target", List.of(a, b), K,
                new double[] { 0.40, 0.36, 0.24 });
    }

    /**
     * The pre-V2 EquilibriumMatrix.solve, re-implemented verbatim on the same
     * PhaseContribution fields (m plays the role of the old `composition`,
     * deltaM0/deltaMSensitivity the role of the old deltaComposition0/...). For
     * M == x these are the same numbers, so this must match the production
     * result exactly.
     */
    private static void compareToPreV2(String tag, List<PhaseContribution> phases,
            int K, double[] target) {
        int np = phases.size();
        int n = K + np;
        double[][] A = new double[n][n];
        double[] b = new double[n];
        for (int p = 0; p < np; p++) {
            PhaseContribution ph = phases.get(p);
            for (int i = 0; i < K; i++) A[p][i] = ph.m()[i];
            b[p] = ph.g();
        }
        for (int i = 0; i < K; i++) {
            int row = np + i;
            double represented = 0, rhs = 0;
            for (int p = 0; p < np; p++) {
                PhaseContribution ph = phases.get(p);
                A[row][K + p] = ph.m()[i];
                for (int k = 0; k < K; k++) A[row][k] += ph.amount() * ph.deltaMSensitivity()[k][i];
                rhs -= ph.amount() * ph.deltaM0()[i];
                represented += ph.amount() * ph.m()[i];
            }
            double ri = (target != null) ? target[i] - represented : 0.0;
            b[row] = ri + rhs;
        }
        double[] ref = gauss(A, b);
        double[] refMu = new double[K];
        double[] refDN = new double[np];
        System.arraycopy(ref, 0, refMu, 0, K);
        System.arraycopy(ref, K, refDN, 0, np);

        EquilibriumStepResult got = EquilibriumMatrix.solve(phases, K, target);

        check(tag + ": mu matches pre-V2 assembly", close(got.mu(), refMu, TOL),
                "max|d| = " + maxDiff(got.mu(), refMu));
        check(tag + ": deltaN matches pre-V2 assembly", close(got.deltaN(), refDN, TOL),
                "max|d| = " + maxDiff(got.deltaN(), refDN));
    }

    // ==================================================================
    // F -- synthetic nontrivial M(Y): M_A = xA - eta/2
    // ==================================================================

    private static void syntheticNontrivialM() {
        System.out.printf("%n--- F: synthetic nontrivial M (M_A = xA - eta/2), np=2, K=2 ---%n");
        int K = 2;

        // Two phases. We choose x, eta and thus M directly; deltaM0 / deltaMSens
        // are chosen arbitrarily but MUST sum to zero across components (M stays
        // on the simplex: M_A + M_B = 1 for this representation) so the assembly
        // is physically consistent -- the test does not require that, but it
        // makes the numbers realistic.
        double xA1 = 0.45, eta1 = 0.10;                 // M_A1 = 0.40, M_B1 = 0.60
        double xA2 = 0.70, eta2 = -0.20;                // M_A2 = 0.80, M_B2 = 0.20
        double[] m1 = { xA1 - eta1 / 2, (1 - xA1) + eta1 / 2 };
        double[] m2 = { xA2 - eta2 / 2, (1 - xA2) + eta2 / 2 };

        double[] dM0_1 = { 0.012, -0.012 };
        double[] dM0_2 = { -0.008, 0.008 };
        double[][] dMs_1 = { { 0.0030, -0.0030 }, { -0.0020, 0.0020 } };
        double[][] dMs_2 = { { 0.0025, -0.0025 }, { -0.0018, 0.0018 } };

        double N1 = 0.5, N2 = 0.5;
        PhaseContribution p1 = new PhaseContribution(N1, m1, -5200.0, dM0_1, dMs_1);
        PhaseContribution p2 = new PhaseContribution(N2, m2, -4800.0, dM0_2, dMs_2);
        List<PhaseContribution> phases = List.of(p1, p2);

        // target deliberately OFF the currently represented inventory
        double[] represented = {
                N1 * m1[0] + N2 * m2[0],
                N1 * m1[1] + N2 * m2[1] };
        double[] target = { represented[0] + 0.03, represented[1] - 0.03 };

        // ---- independently assemble the generalised (K+np) system ----
        int np = 2, n = K + np;
        double[][] A = new double[n][n];
        double[] b = new double[n];

        // Gibbs-Duhem rows: sum_A mu_A * M_A^p = G_p
        for (int p = 0; p < np; p++) {
            PhaseContribution ph = phases.get(p);
            for (int i = 0; i < K; i++) A[p][i] = ph.m()[i];
            b[p] = ph.g();
        }
        // Mass-balance rows: sum_p deltaN_p M_A^p
        //                    + Σ_k mu_k ( sum_p N_p dMSens^p[k][A] )
        //                    = r_A - sum_p N_p dM0^p[A]
        for (int i = 0; i < K; i++) {
            int row = np + i;
            double rep = 0, rhs = 0;
            for (int p = 0; p < np; p++) {
                PhaseContribution ph = phases.get(p);
                A[row][K + p] = ph.m()[i];
                for (int k = 0; k < K; k++) A[row][k] += ph.amount() * ph.deltaMSensitivity()[k][i];
                rhs -= ph.amount() * ph.deltaM0()[i];
                rep += ph.amount() * ph.m()[i];
            }
            double ri = target[i] - rep;
            b[row] = ri + rhs;
        }

        // check individual coefficients explicitly (not just the solved answer)
        check("F: Gibbs-Duhem row 0 coeffs == M^0", A[0][0] == m1[0] && A[0][1] == m1[1], "GD0");
        check("F: Gibbs-Duhem row 1 coeffs == M^1", A[1][0] == m2[0] && A[1][1] == m2[1], "GD1");
        check("F: Gibbs-Duhem RHS == G", b[0] == -5200.0 && b[1] == -4800.0, "GD rhs");
        // deltaN coefficient in mass-balance row A for phase p is M_A^p
        check("F: mass-balance deltaN coeff [A=0][p=0] == M_0^0",
                A[2][K + 0] == m1[0], "= " + A[2][K + 0] + " vs " + m1[0]);
        check("F: mass-balance deltaN coeff [A=1][p=1] == M_1^1",
                A[3][K + 1] == m2[1], "= " + A[3][K + 1] + " vs " + m2[1]);
        // mu coefficient in mass-balance row A=0 is sum_p N_p dMSens^p[k][0]
        double expMuCoef00 = N1 * dMs_1[0][0] + N2 * dMs_2[0][0];
        check("F: mass-balance mu coeff [row A=0][mu_0] == sum_p N_p dMSens^p[0][0]",
                Math.abs(A[2][0] - expMuCoef00) < 1e-15, "= " + A[2][0] + " vs " + expMuCoef00);
        // target residual r_A in the RHS
        double rA0 = target[0] - represented[0];
        double expRhs0 = rA0 - (N1 * dM0_1[0] + N2 * dM0_2[0]);
        check("F: mass-balance RHS row A=0 == r_0 - sum_p N_p dM0^p[0]",
                Math.abs(b[2] - expRhs0) < 1e-15, "= " + b[2] + " vs " + expRhs0);

        double[] ref = gauss(A, b);
        double[] refMu = { ref[0], ref[1] };
        double[] refDN = { ref[2], ref[3] };

        EquilibriumStepResult got = EquilibriumMatrix.solve(phases, K, target);
        check("F: production mu matches independent generalised assembly",
                close(got.mu(), refMu, TOL), "max|d| = " + maxDiff(got.mu(), refMu));
        check("F: production deltaN matches independent generalised assembly",
                close(got.deltaN(), refDN, TOL), "max|d| = " + maxDiff(got.deltaN(), refDN));

        // Gibbs-Duhem must hold at the solution: sum_A mu_A M_A^p == G_p
        for (int p = 0; p < np; p++) {
            PhaseContribution ph = phases.get(p);
            double gd = got.mu()[0] * ph.m()[0] + got.mu()[1] * ph.m()[1];
            check("F: Gibbs-Duhem residual phase " + p + " (sum mu*M == G)",
                    Math.abs(gd - ph.g()) < 1e-6, "sum mu*M = " + gd + " vs G = " + ph.g());
        }

        // Mass-balance Newton row must be satisfied by (mu, deltaN):
        //   sum_p deltaN_p M_A^p + Σ_k mu_k sum_p N_p dMSens^p[k][A]
        //     == r_A - sum_p N_p dM0^p[A]
        for (int i = 0; i < K; i++) {
            double lhs = 0;
            for (int p = 0; p < np; p++) {
                PhaseContribution ph = phases.get(p);
                lhs += got.deltaN()[p] * ph.m()[i];
                for (int k = 0; k < K; k++)
                    lhs += got.mu()[k] * ph.amount() * ph.deltaMSensitivity()[k][i];
            }
            double rA = target[i] - represented[i];
            double rhs = rA;
            for (int p = 0; p < np; p++)
                rhs -= phases.get(p).amount() * phases.get(p).deltaM0()[i];
            check("F: mass-balance Newton row A=" + i + " satisfied",
                    Math.abs(lhs - rhs) < 1e-6, "lhs=" + lhs + " rhs=" + rhs);
        }
    }

    // ==================================================================
    // helpers
    // ==================================================================

    private static double[] gauss(double[][] Ain, double[] bin) {
        int n = bin.length;
        double[][] A = new double[n][n];
        double[] b = bin.clone();
        for (int i = 0; i < n; i++) A[i] = Ain[i].clone();
        for (int col = 0; col < n; col++) {
            int piv = col;
            for (int r = col + 1; r < n; r++)
                if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) piv = r;
            double[] t = A[col]; A[col] = A[piv]; A[piv] = t;
            double tb = b[col]; b[col] = b[piv]; b[piv] = tb;
            for (int r = col + 1; r < n; r++) {
                double f = A[r][col] / A[col][col];
                for (int j = col; j < n; j++) A[r][j] -= f * A[col][j];
                b[r] -= f * b[col];
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double s = b[i];
            for (int j = i + 1; j < n; j++) s -= A[i][j] * x[j];
            x[i] = s / A[i][i];
        }
        return x;
    }

    private static double maxDiff(double[] a, double[] b) {
        double m = 0;
        for (int i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
        return m;
    }

    private static boolean close(double[] a, double[] b, double tol) {
        return maxDiff(a, b) <= tol;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-74s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-74s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertGeneralizedEquilibriumMatrix() {
    }
}
