package org.ce.model.equilibrium;

import org.ce.model.equilibrium.HillertSolver.CandidatePhaseState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * V2 STEP 4 gate (Part 10): a small independent synthetic test of the
 * deterministic multi-start + duplicate-minimum filtering + best-candidate
 * selection, on a {@code Phi(Y)} with TWO analytically-known local minima.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertOrderedBasinSynthetic
 * </pre>
 *
 * <p>In {@code org.ce.model.equilibrium} so it can drive the package-private
 * {@link HillertSolver#selectBestDistinctMinimum} /
 * {@link HillertSolver#countDistinctMinima} / {@link HillertSolver#minimumForTest}.
 * No CVM model, no production solver -- the expected minima are known in closed
 * form.</p>
 *
 * <h2>The synthetic objective</h2>
 * <p>A 1-D asymmetric double well
 * {@code Phi(y) = KW * (y^2 - 1)^2 - TILT * y}
 * (the "constitution vector" is the scalar {@code y}). For
 * {@code KW = 1000}, {@code TILT = 400} it has:</p>
 * <ul>
 *   <li>a shallow minimum near {@code y = -1} (call it A),</li>
 *   <li>a deeper minimum near {@code y = +1} (call it B),</li>
 *   <li>a barrier at {@code y ~ 0},</li>
 * </ul>
 * <p>so {@code Phi(A) > Phi(B)}. With the sign convention {@code dGf = -Phi}
 * (the driving force the production search maximises), {@code dGf(B) > dGf(A)}.
 * A single seed at {@code y = -0.5} descends only to A; a seed at
 * {@code y = +0.5} descends only to B -- two deterministic seeds cover both
 * basins.</p>
 *
 * <h2>Checks</h2>
 * <ol>
 *   <li>Newton descent from {@code y = -0.5} reaches basin A
 *       ({@code y_A ~ -1.07}, analytically the negative root of
 *       {@code 4 KW y (y^2 - 1) - TILT = 0}); from {@code y = +0.5} reaches B
 *       ({@code y_B ~ +1.06}).</li>
 *   <li>Both minima are found (the multi-start covers both basins).</li>
 *   <li>{@link HillertSolver#countDistinctMinima} sees exactly 2 after dedup.</li>
 *   <li>A near-duplicate of A (within {@code SEED_DEDUP_REL}) is merged -- still 2.</li>
 *   <li>A distinct minimum with nearly the same {@code Phi} as B is NOT merged
 *       (dedup is on {@code Y}, not energy).</li>
 *   <li>{@link HillertSolver#selectBestDistinctMinimum} returns the
 *       lower-{@code Phi} / higher-{@code dGf} minimum (B).</li>
 * </ol>
 */
public final class HillertOrderedBasinSynthetic {

    private static int failures = 0;

    private static final double KW = 1000.0;
    private static final double TILT = 400.0;

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(84));
        System.out.println("  V2 STEP 4 Part 10 -- synthetic two-basin multi-start / dedup / selection");
        System.out.println("=".repeat(84));

        // ---- analytic minima: 4 KW y (y^2 - 1) - TILT = 0 ----
        double yA = newtonRoot(-1.0);   // negative-side root  (shallow well, higher Phi)
        double yB = newtonRoot(+1.0);   // positive-side root  (deep well, lower Phi)
        System.out.printf("  analytic minima:  A y=%.6f Phi=%.4f   B y=%.6f Phi=%.4f%n",
                yA, phi(yA), yB, phi(yB));
        check("analytic: A and B are both stationary (Phi'(y) ~ 0)",
                Math.abs(dPhi(yA)) < 1e-6 && Math.abs(dPhi(yB)) < 1e-6,
                "dPhi(A)=" + dPhi(yA) + " dPhi(B)=" + dPhi(yB));
        check("analytic: A and B are both minima (Phi''(y) > 0)",
                d2Phi(yA) > 0 && d2Phi(yB) > 0, "curvature");
        check("analytic: Phi(A) > Phi(B) (A is the shallower basin)",
                phi(yA) > phi(yB) + 1.0, "PhiA=" + phi(yA) + " PhiB=" + phi(yB));
        check("analytic: the two minima are far apart (distinct basins)",
                Math.abs(yA - yB) > 1.0, "|yA-yB|=" + Math.abs(yA - yB));

        // ---- 1. two deterministic seeds, one per basin, local descent ----
        double[] seeds = { -0.5, +0.5 };   // fixed, no random
        double relaxedFromNeg = descend(seeds[0]);
        double relaxedFromPos = descend(seeds[1]);
        System.out.printf("  descent: seed %.2f -> %.6f ;  seed %.2f -> %.6f%n",
                seeds[0], relaxedFromNeg, seeds[1], relaxedFromPos);
        check("1: seed y=-0.5 descends to basin A", Math.abs(relaxedFromNeg - yA) < 1e-6,
                "got " + relaxedFromNeg);
        check("1: seed y=+0.5 descends to basin B", Math.abs(relaxedFromPos - yB) < 1e-6,
                "got " + relaxedFromPos);
        check("2: the multi-start found BOTH minima (single seed would miss one)",
                Math.abs(relaxedFromNeg - relaxedFromPos) > 1.0, "same basin");

        // ---- 2-3. feed the discovered minima to the production dedup/select ----
        List<CandidatePhaseState> found = new ArrayList<>();
        found.add(HillertSolver.minimumForTest(new double[] { relaxedFromNeg }, -phi(relaxedFromNeg)));
        found.add(HillertSolver.minimumForTest(new double[] { relaxedFromPos }, -phi(relaxedFromPos)));
        check("3: countDistinctMinima sees exactly 2 after dedup",
                HillertSolver.countDistinctMinima(found) == 2,
                "got " + HillertSolver.countDistinctMinima(found));

        // ---- 4. a near-duplicate of A is merged ----
        List<CandidatePhaseState> withDup = new ArrayList<>(found);
        withDup.add(HillertSolver.minimumForTest(
                new double[] { relaxedFromNeg + 1e-9 }, -phi(relaxedFromNeg)));   // ~A again
        check("4: a near-duplicate of A (within SEED_DEDUP_REL) is merged -- still 2",
                HillertSolver.countDistinctMinima(withDup) == 2,
                "got " + HillertSolver.countDistinctMinima(withDup));

        // ---- 5. a DISTINCT minimum with nearly the same Phi as B is NOT merged ----
        // place it well away in Y but hand it dGf ~ dGf(B)
        List<CandidatePhaseState> withEnergyTwin = new ArrayList<>(found);
        withEnergyTwin.add(HillertSolver.minimumForTest(
                new double[] { relaxedFromPos - 0.4 }, -phi(relaxedFromPos) - 0.01));
        check("5: a physically distinct minimum with ~equal Phi is NOT merged (dedup on Y, not energy)",
                HillertSolver.countDistinctMinima(withEnergyTwin) == 3,
                "got " + HillertSolver.countDistinctMinima(withEnergyTwin));

        // ---- 6. selection picks the lower-Phi / higher-dGf minimum ----
        CandidatePhaseState best = HillertSolver.selectBestDistinctMinimum(found);
        check("6: selectBestDistinctMinimum returns basin B (lowest Phi / highest dGf)",
                Math.abs(best.uFull()[0] - yB) < 1e-6, "picked y=" + best.uFull()[0]);
        check("6: selected dGf == -Phi(B) (max over the distinct minima)",
                Math.abs(best.drivingForce() - (-phi(yB))) < 1e-6,
                "dGf=" + best.drivingForce());
        // order-independence of selection
        List<CandidatePhaseState> reversed = new ArrayList<>();
        reversed.add(found.get(1));
        reversed.add(found.get(0));
        check("6: selection is order-independent (same winner if the list is reversed)",
                Math.abs(HillertSolver.selectBestDistinctMinimum(reversed).uFull()[0] - yB) < 1e-6,
                "different winner");

        // ---- determinism: repeat everything, identical results ----
        double d1 = descend(seeds[0]), d2 = descend(seeds[0]);
        check("8: local descent is deterministic (identical on repeat)", d1 == d2, d1 + " vs " + d2);

        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(84));
        if (failures > 0) {
            throw new AssertionError(failures + " basin-test checks failed");
        }
    }

    // Phi(y) = KW (y^2 - 1)^2 - TILT y
    private static double phi(double y)   { double a = y * y - 1; return KW * a * a - TILT * y; }
    private static double dPhi(double y)  { return 4 * KW * y * (y * y - 1) - TILT; }
    private static double d2Phi(double y) { return 4 * KW * (3 * y * y - 1); }

    /** Newton on dPhi = 0 from a starting guess -- the analytic minimum. */
    private static double newtonRoot(double y0) {
        double y = y0;
        for (int i = 0; i < 200; i++) {
            double f = dPhi(y), fp = d2Phi(y);
            double dy = -f / fp;
            y += dy;
            if (Math.abs(dy) < 1e-14) break;
        }
        return y;
    }

    /**
     * Bounded Newton descent on Phi from a seed -- the 1-D analogue of
     * {@link HillertSolver#relaxWide}'s inner loop (Newton step on the
     * stationarity condition, halve the step if it would increase Phi).
     */
    private static double descend(double seed) {
        double y = seed;
        for (int it = 0; it < 200; it++) {
            double g = dPhi(y);
            if (Math.abs(g) < 1e-12) break;
            double step = -g / Math.abs(d2Phi(y));   // guarded Newton (descent direction)
            double lambda = 1.0;
            for (int bt = 0; bt < 60; bt++) {
                if (phi(y + lambda * step) < phi(y)) break;
                lambda *= 0.5;
            }
            double dy = lambda * step;
            y += dy;
            if (Math.abs(dy) < 1e-14) break;
        }
        // polish to the exact stationary point
        return newtonRoot(y);
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-74s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-74s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertOrderedBasinSynthetic() {
    }
}
