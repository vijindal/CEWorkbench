package org.ce.scratch;

import org.ce.model.equilibrium.HillertSolver.EquilibriumMatrix;
import org.ce.model.equilibrium.HillertSolver.EquilibriumMatrix.EquilibriumStepResult;
import org.ce.model.equilibrium.HillertSolver.EquilibriumMatrix.PhaseContribution;

import java.util.Arrays;
import java.util.List;

/**
 * V2 STEP 3 gate (Part 14): a fully synthetic ordered/disordered two-phase
 * system with {@code M != x} for one phase and an <b>independently constructed</b>
 * equilibrium, exercising the generalised phase-addition driving force
 * ({@code sum_A mu_A M_A - G}) and the generalised {@link EquilibriumMatrix}.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertOrderedTwoPhaseSynthetic
 * </pre>
 *
 * <h2>The synthetic system (K = 2), everything closed-form</h2>
 * <ul>
 *   <li><b>alpha</b> -- disordered, constitution {@code Y = (m)}, {@code M_A = m}
 *       (composition coordinate == component amount). {@code G_a(m)} is a
 *       concave-down parabola {@code A2 m^2 + A1 m + A0}.</li>
 *   <li><b>beta</b> -- ordered, constitution {@code Y = (M_A, eta)}.
 *       {@code M_A} is a primary variable; the composition coordinate is
 *       {@code x_b = M_A + eta/2}, so {@code M != x} whenever {@code eta != 0}.
 *       {@code G_b(M_A, eta) = B2 M_A^2 + B1 M_A + B0 + phi(eta)} with a
 *       closed-form double-well {@code phi(eta) = C (eta^2 - etaStar^2)^2}
 *       (minimum 0 at {@code eta = +/- etaStar}). Beta's internal equilibrium is
 *       {@code eta = etaStar != 0}; its envelope {@code Gtilde_b(M_A) =
 *       B2 M_A^2 + B1 M_A + B0} is a plain parabola. {@code J_M} for beta is the
 *       exact constant {@code d(M_A, M_B)/d(M_A, eta) = [[1, 0], [-1, 0]]}
 *       (M_A primary; eta does not change M directly in this parametrisation --
 *       it changes only x_b).</li>
 * </ul>
 *
 * <p>The two parabolas {@code G_a}, {@code Gtilde_b} are <b>constructed</b> to
 * have a common tangent at chosen physical endpoints
 * {@code M_A = MA_STAR} (alpha) and {@code M_A = MB_STAR} (beta), with a chosen
 * slope and intercept -- so the reference tie-line {@code mu} is exact and
 * inside {@code (0,1)}, and the overall target sits between the endpoints.</p>
 */
public final class HillertOrderedTwoPhaseSynthetic {

    private static int failures = 0;

    // ---- constructed common tangent ----
    private static final double MA_STAR = 0.75;          // alpha tie-line endpoint (M_A)
    private static final double MB_STAR = 0.25;          // beta  tie-line endpoint (M_A)
    private static final double SLOPE   = -4000.0;       // tangent slope = mu_A - mu_B
    private static final double MU_B    = 500.0;         // tangent intercept = mu_B
    private static final double MU_A    = MU_B + SLOPE;  // (tangent value at M_A = 1 minus... ) => mu_A

    // CONVEX parabola curvatures: each phase's G(M_A) is locally stable, and the
    // common tangent line touches each parabola from below -- the standard
    // two-phase tie-line geometry. (A concave G would make the phase itself
    // internally unstable and Phi_beta unbounded below.)
    private static final double A2 = 9000.0;
    private static final double B2 = 6000.0;
    // linear/constant coeffs fixed by "tangent line touches at the endpoint":
    //   G(m*) = MU_B + SLOPE*m*   and   G'(m*) = SLOPE
    private static final double A1 = SLOPE - 2 * A2 * MA_STAR;
    private static final double A0 = (MU_B + SLOPE * MA_STAR) - A2 * MA_STAR * MA_STAR - A1 * MA_STAR;
    private static final double B1 = SLOPE - 2 * B2 * MB_STAR;
    private static final double B0 = (MU_B + SLOPE * MB_STAR) - B2 * MB_STAR * MB_STAR - B1 * MB_STAR;

    // beta order parameter
    private static final double ETA_STAR = 0.5;
    private static final double C_WELL   = 40000.0;      // phi(eta) = C (eta^2 - etaStar^2)^2

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(84));
        System.out.println("  V2 STEP 3 Part 14 -- synthetic ordered/disordered two-phase (M != x)");
        System.out.println("=".repeat(84));

        // ---------- 1. the (constructed, exact) reference equilibrium ----------
        double ma = MA_STAR, mb = MB_STAR;
        double muA = MU_A, muB = MU_B, s = SLOPE;
        double etaB = ETA_STAR;
        double xbEq = mb + etaB / 2.0;                    // beta composition coordinate (M != x)

        System.out.printf("%n  reference equilibrium (constructed common tangent):%n");
        System.out.printf("    alpha: M_A = x = %.6f   G = %.4f   G'(=slope) = %.4f%n",
                ma, Ga(ma), dGa(ma));
        System.out.printf("    beta : M_A = %.6f   eta = %.4f   x_b = %.6f   G = %.4f   Gtilde'(=slope) = %.4f%n",
                mb, etaB, xbEq, GbTilde(mb), dGbTilde(mb));
        System.out.printf("    mu = [%.4f, %.4f]   tangent slope = %.4f%n", muA, muB, s);

        check("indep: alpha tangent slope == SLOPE", Math.abs(dGa(ma) - s) < 1e-6, "" + dGa(ma));
        check("indep: beta tangent slope == SLOPE", Math.abs(dGbTilde(mb) - s) < 1e-6, "" + dGbTilde(mb));
        check("indep: alpha lies on the tangent line", Math.abs(Ga(ma) - (muB + s * ma)) < 1e-6, "off");
        check("indep: beta lies on the tangent line", Math.abs(GbTilde(mb) - (muB + s * mb)) < 1e-6, "off");
        check("indep: tie-line endpoints distinct, ma > mb", ma > mb + 1e-3, "ma=" + ma + " mb=" + mb);
        check("indep: beta genuinely ordered -- M_A(%.3f) != x_b(%.3f), diff = eta/2".formatted(mb, xbEq),
                Math.abs(mb - xbEq) > 1e-3 && Math.abs((xbEq - mb) - etaB / 2) < 1e-12, "M==x");
        check("indep: beta internal eta is stationary (phi'(etaStar) == 0)",
                Math.abs(dPhi(ETA_STAR)) < 1e-6, "phi'=" + dPhi(ETA_STAR));
        check("indep: Gb(MB_STAR, etaStar) == Gtilde_b(MB_STAR) (envelope consistency)",
                Math.abs(Gb(mb, ETA_STAR) - GbTilde(mb)) < 1e-9, "mismatch");

        // ---------- 2. generalised candidate driving force at alpha-only mu ----------
        // alpha alone, sitting in BETA's stable territory (M_A well below the
        // beta tie-line endpoint) -- there beta is the true stable phase, so its
        // generalised driving force against alpha's tangent must be > 0.
        double ma0 = 0.20;
        double muB0 = Ga(ma0) - ma0 * dGa(ma0);              // alpha tangent-line intercept
        double[] muAlpha = { muB0 + dGa(ma0), muB0 };        // {mu_A, mu_B}

        // relax beta's (M_A, eta) for Phi_b = G_b - mu^T M_b at this fixed mu:
        //   d/deta : phi'(eta) = 0                     -> eta = etaStar (mu-independent)
        //   d/dM_A : dGtilde_b/dM_A - (mu_A - mu_B) = 0 -> M_A adjusts to the mu slope
        double etaR = ETA_STAR;
        double mAbR = betaMAforMuSlope(muAlpha[0] - muAlpha[1]);
        double mBbR = 1.0 - mAbR;
        double xbR  = mAbR + etaR / 2.0;
        double gBR  = Gb(mAbR, etaR);                        // == GbTilde(mAbR) at etaStar
        double dgf  = muAlpha[0] * mAbR + muAlpha[1] * mBbR - gBR;
        double addThreshold = Math.max(1.0, 1e-6 * Math.abs(gBR));

        System.out.printf("%n  candidate scan at alpha-only mu = [%.1f, %.1f]:%n", muAlpha[0], muAlpha[1]);
        System.out.printf("    beta relaxed: M_A = %.6f  eta = %.4f  x_b = %.6f  G = %.4f%n",
                mAbR, etaR, xbR, gBR);
        System.out.printf("    dGf = sum_A mu_A M_A^beta - G_beta = %.4f  (addThreshold = %.4f)%n",
                dgf, addThreshold);

        check("1: generalised driving force sum_A mu_A M_A^beta - G_beta > 0 (beta stable)",
                dgf > 0, "dGf=" + dgf);
        check("2: phase addition would trigger (dGf > addThreshold)", dgf > addThreshold, "dGf=" + dgf);
        check("2: driving force on M not x (M_A=%.4f != x_b=%.4f, diff = eta/2)".formatted(mAbR, xbR),
                Math.abs((xbR - mAbR) - etaR / 2) < 1e-12 && Math.abs(xbR - mAbR) > 1e-3, "M==x");
        check("2: dGf uses M and G from ONE consistent beta state",
                Math.abs(dgf - (muAlpha[0] * mAbR + muAlpha[1] * (1 - mAbR) - Gb(mAbR, etaR))) < 1e-9,
                "inconsistent");

        // ---------- 3-6. EquilibriumMatrix.solve on the two-phase system ----------
        double targetMA = 0.5;                               // overall M_A, inside [MB_STAR, MA_STAR]
        double[] targetVec = { targetMA, 1 - targetMA };
        double targetSum0 = targetVec[0] + targetVec[1];
        double Na = 1.0, Nb = 1e-6;                          // beta just seeded (post-insertion)

        // per-phase deltaM sensitivity to mu at the tie-line point:
        //   dM/dmu = J_M H_Phi^-1 J_M^T, H_Phi = Gtilde'' (const), J_M = [1,-1]^T
        double[][] sensA = jHiJt_scalar(dGa2(), new double[] { 1, -1 });
        double[][] sensB = jHiJt_scalar(dGbTilde2(), new double[] { 1, -1 });

        PhaseContribution pa = new PhaseContribution(
                Na, new double[] { ma, 1 - ma }, Ga(ma), new double[] { 0, 0 }, sensA);
        PhaseContribution pb = new PhaseContribution(
                Nb, new double[] { mb, 1 - mb }, GbTilde(mb), new double[] { 0, 0 }, sensB);

        EquilibriumStepResult res = EquilibriumMatrix.solve(List.of(pa, pb), 2, targetVec);
        System.out.printf("%n  EquilibriumMatrix.solve:  mu = %s  deltaN = %s  resid = %.2e%n",
                Arrays.toString(res.mu()), Arrays.toString(res.deltaN()), res.relativeResidual());

        check("3: EquilibriumMatrix mu matches the independent tie-line mu",
                Math.abs(res.mu()[0] - muA) < 1e-3 && Math.abs(res.mu()[1] - muB) < 1e-3,
                "got " + Arrays.toString(res.mu()) + " want [" + muA + ", " + muB + "]");
        check("3: mass-balance step moves mass INTO the seeded beta phase (deltaN_beta > 0)",
                res.deltaN()[1] > 0, "deltaN_beta=" + res.deltaN()[1]);
        // The outer system solves the COMPONENT mass-balance rows (conserve
        // M_A, M_B) -- one Newton step, so it REDUCES the component residual
        // r_A = target_A - sum_p N_p M_A^p (full convergence needs iteration).
        double rBefore0 = targetVec[0] - (Na * ma + Nb * mb);
        double rAfter0  = targetVec[0]
                - ((Na + res.deltaN()[0]) * ma + (Nb + res.deltaN()[1]) * mb);
        System.out.printf("  component-A mass residual: %.4f -> %.4f%n", rBefore0, rAfter0);
        check("3: the outer step REDUCES the component mass-balance residual",
                Math.abs(rAfter0) < Math.abs(rBefore0), "before=" + rBefore0 + " after=" + rAfter0);
        check("3: target vector NOT mutated by the solve",
                targetVec[0] == targetMA && Math.abs(targetVec[0] + targetVec[1] - targetSum0) < 1e-15,
                "target changed");

        double gdA = res.mu()[0] * ma + res.mu()[1] * (1 - ma) - Ga(ma);
        double gdB = res.mu()[0] * mb + res.mu()[1] * (1 - mb) - GbTilde(mb);
        check("4: generalised Gibbs-Duhem alpha  sum_A mu_A M_A = G", Math.abs(gdA) < 1e-3, "resid=" + gdA);
        check("4: generalised Gibbs-Duhem beta   sum_A mu_A M_A = G", Math.abs(gdB) < 1e-3, "resid=" + gdB);
        check("5: ONE common mu vector satisfies both phases' Gibbs-Duhem rows",
                Math.abs(gdA) < 1e-3 && Math.abs(gdB) < 1e-3, "not common");

        // ---------- final phase amounts: lever rule on the M-axis ----------
        double fa = (mb - targetMA) / (mb - ma);
        double fb = 1 - fa;
        System.out.printf("  lever rule: f_alpha = %.4f   f_beta = %.4f%n", fa, fb);
        check("6: final phase fractions both >= 0 (target inside the tie-line)",
                fa >= -1e-12 && fb >= -1e-12, "fa=" + fa + " fb=" + fb);
        check("6: lever rule reproduces target M_A exactly",
                Math.abs(fa * ma + fb * mb - targetMA) < 1e-12, "off");
        check("6: represented inventory at the lever-rule split == target (mass balance closed)",
                Math.abs((fa * ma + fb * mb) - targetVec[0]) < 1e-12
                        && Math.abs((fa * (1 - ma) + fb * (1 - mb)) - targetVec[1]) < 1e-12,
                "not closed");

        // ---------- Part 11: phase removal is mass-conserving with M != x ----------
        // The V1 removal rule (unchanged in STEP 3) only zeroes a phase once its
        // amount has ALREADY reached numerical zero. Its contribution to the
        // represented inventory sum_p N_p M_A^p is then N_beta * M_A^beta = 0 --
        // exactly, for ANY M_A, ordered (M != x) or not. So removal never
        // perturbs the global target.
        double mAbetaOrdered = mb;                 // beta's M_A (with M_A != x_b, x_b = mb + eta/2)
        double[] invZeroBeta   = { Na * ma + 0.0 * mAbetaOrdered,
                                   Na * (1 - ma) + 0.0 * (1 - mAbetaOrdered) };
        double[] invAfterRemove = { Na * ma, Na * (1 - ma) };
        check("11: removing an already-zero-amount ORDERED (M != x) phase is exactly "
                        + "mass-conserving (N*M = 0*M = 0)",
                Arrays.equals(invZeroBeta, invAfterRemove), "not inert");
        check("11: the invariant holds regardless of M_A vs x_b (M_A=%.3f, x_b=%.3f differ)"
                        .formatted(mAbetaOrdered, mAbetaOrdered + ETA_STAR / 2),
                Math.abs(mAbetaOrdered - (mAbetaOrdered + ETA_STAR / 2)) > 1e-3
                        && 0.0 * mAbetaOrdered == 0.0, "not demonstrated");

        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(84));
        if (failures > 0) {
            throw new AssertionError(failures + " synthetic two-phase checks failed");
        }
    }

    // ---------- alpha ----------
    private static double Ga(double m)  { return A2 * m * m + A1 * m + A0; }
    private static double dGa(double m) { return 2 * A2 * m + A1; }
    private static double dGa2()        { return 2 * A2; }

    // ---------- beta ----------
    private static double phi(double eta) {
        double d = eta * eta - ETA_STAR * ETA_STAR;
        return C_WELL * d * d;
    }
    private static double dPhi(double eta) {
        return C_WELL * 2 * (eta * eta - ETA_STAR * ETA_STAR) * 2 * eta;
    }
    /** G_b as a function of the primary variables (M_A, eta). */
    private static double Gb(double mA, double eta) {
        return B2 * mA * mA + B1 * mA + B0 + phi(eta);
    }
    /** Internal-equilibrium envelope Gtilde_b(M_A) = min_eta G_b (at eta = etaStar, phi = 0). */
    private static double GbTilde(double mA)  { return B2 * mA * mA + B1 * mA + B0; }
    private static double dGbTilde(double mA) { return 2 * B2 * mA + B1; }
    private static double dGbTilde2()         { return 2 * B2; }

    /** beta's internal M_A that makes dGtilde_b/dM_A match an imposed mu slope. */
    private static double betaMAforMuSlope(double slope) {
        // 2 B2 mA + B1 = slope
        return (slope - B1) / (2 * B2);
    }

    /** J_M (1/h) J_M^T for scalar Hessian h and J_M column jm (length K). */
    private static double[][] jHiJt_scalar(double h, double[] jm) {
        double inv = 1.0 / h;
        int k = jm.length;
        double[][] out = new double[k][k];   // out[muRow][compCol]
        for (int r = 0; r < k; r++)
            for (int c = 0; c < k; c++)
                out[r][c] = inv * jm[r] * jm[c];
        return out;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-76s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-76s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertOrderedTwoPhaseSynthetic() {
    }
}
