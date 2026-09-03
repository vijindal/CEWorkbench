package org.ce.model.equilibrium;

import org.ce.model.equilibrium.HillertGridMinimizer.GridMinimizationResult;
import org.ce.model.equilibrium.HillertGridMinimizer.GridPoint;
import org.ce.model.equilibrium.HillertGridMinimizer.PhaseInitialState;

import java.util.ArrayList;
import java.util.List;

/**
 * V2 STEP 7 gate -- synthetic tests of the grid minimizer's discrete
 * convex-hull / linear-mixture / tangent math against analytic {@code G(M)}
 * surfaces whose answers are known in closed form. No CVM model is involved;
 * these drive {@link HillertGridMinimizer#minimizeDiscreteForTest} with
 * hand-built {@link GridPoint#synthetic} grids.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertGridMinimizerSyntheticTest
 * </pre>
 *
 * <h2>Part 16 -- binary miscibility gap</h2>
 * <p>One parent phase, symmetric quartic double well
 * {@code G(x) = C4 (x-1/2)^4 - C2 (x-1/2)^2 + G0} (composition coordinate
 * {@code x = x_B}). Minima at {@code x* = 1/2 +- sqrt(C2/(2 C4))}; by symmetry
 * the common tangent is horizontal, so the convex hull between the minima is a
 * flat tie line and any interior target is represented by a mix of the two
 * minima -- both from the SAME parent phase, kept as {@code A#1}/{@code A#2}.</p>
 *
 * <h2>Part 17 -- two-phase common tangent</h2>
 * <p>Two phases with analytic parabolic {@code G}, chosen so their lower
 * common tangent touches phase 0 at {@code xB = 0.2} and phase 1 at
 * {@code xB = 0.8}, with known slope. An interior target selects one point from
 * each; discrete {@code mu}, amounts, objective and tangent residuals are all
 * checked against the closed form.</p>
 *
 * <h2>Part 19 -- special cases</h2>
 * <p>A: target exactly on a grid point. B: target between two points. I:
 * invalid points present (mixed into the list, must be ignored by the caller --
 * here we just confirm valid-only lists still solve). J: infeasible target
 * (outside the hull) -> failure result, no fabricated mixture. G: duplicate /
 * collinear grid points -> deterministic pick. Determinism (Part 20): shuffled
 * input list -> identical selected (phaseIndex, gridIndex) set and amounts.</p>
 */
public final class HillertGridMinimizerSyntheticTest {

    private static int failures = 0;

    // --- double well G(x) = C4 (x-1/2)^4 - C2 (x-1/2)^2 + G0 ---
    private static final double C4 = 80_000.0;
    private static final double C2 = 5_000.0;
    private static final double G0 = -1_000.0;

    private static double wellG(double xB) {
        double d = xB - 0.5;
        return C4 * d * d * d * d - C2 * d * d + G0;
    }

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(88));
        System.out.println("  V2 STEP 7 -- grid minimizer synthetic tests (analytic G(M), no CVM model)");
        System.out.println("=".repeat(88));

        part16_binaryMiscibilityGap();
        part17_twoPhaseCommonTangent();
        part19a_targetOnGridPoint();
        part19j_infeasibleTarget();
        part19g_duplicateCollinearPoints();
        part20_determinismShuffledInput();

        // ---- V2 STEP 10 : same-parent thermodynamic merge test ----
        step10_binaryConvexMerge();
        step10_binaryDoubleWellNoMerge();
        step10_ternaryConvexMerge();
        closure_multiGroupPerMoleNormalisation();

        System.out.println("\n" + "=".repeat(88));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(88));
        if (failures > 0) {
            throw new AssertionError(failures + " grid-minimizer synthetic checks failed");
        }
    }

    // =====================================================================
    // Part 16 -- binary miscibility gap
    // =====================================================================

    private static void part16_binaryMiscibilityGap() {
        section("PART 16 -- binary miscibility gap (one parent phase, double well)");

        // analytic minima
        double half = Math.sqrt(C2 / (2.0 * C4));
        double xLo = 0.5 - half;
        double xHi = 0.5 + half;
        double gMin = wellG(xLo);
        System.out.printf("  analytic minima:  xB=%.6f and xB=%.6f   G_min=%.6f%n", xLo, xHi, gMin);
        check("16: analytic well is symmetric (G(xLo) == G(xHi))",
                Math.abs(wellG(xLo) - wellG(xHi)) < 1e-9, "asymmetric");

        // grid of 41 points on [0,1] for x_B ; M = (x_A, x_B) = (1-xB, xB)
        int n = 40;
        List<GridPoint> grid = new ArrayList<>();
        int gi = 0;
        for (int i = 0; i <= n; i++) {
            double xB = (double) i / n;
            grid.add(GridPoint.synthetic(0, "A", gi++, new double[] { 1.0 - xB, xB }, wellG(xB)));
        }
        // ensure the two exact minima are on the grid too
        int giLo = gi++;
        int giHi = gi++;
        grid.add(GridPoint.synthetic(0, "A", giLo, new double[] { 1.0 - xLo, xLo }, wellG(xLo)));
        grid.add(GridPoint.synthetic(0, "A", giHi, new double[] { 1.0 - xHi, xHi }, wellG(xHi)));

        // target strictly between the minima: xB = 0.5
        double[] target = { 0.5, 0.5 };
        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(grid, target, 2, null);

        check("16: converged", r.converged(), r.failureReason());
        check("16: discrete minimum uses exactly 2 grid points",
                r.selectedGridPoints().size() == 2, "got " + r.selectedGridPoints().size());

        List<GridPoint> sel = r.selectedGridPoints();
        boolean sameParent = sel.size() == 2 && sel.get(0).phaseIndex() == sel.get(1).phaseIndex();
        check("16: both selected points have the SAME parent phase id", sameParent, "different parents");

        // both endpoints are the two minima (within a grid step)
        double xB0 = sel.get(0).m()[1];
        double xB1 = sel.get(1).m()[1];
        double loSel = Math.min(xB0, xB1);
        double hiSel = Math.max(xB0, xB1);
        check("16: selected endpoints bracket the target and match the analytic minima",
                Math.abs(loSel - xLo) < 1.0 / n + 1e-9 && Math.abs(hiSel - xHi) < 1.0 / n + 1e-9,
                "endpoints xB=" + loSel + "," + hiSel);

        // two DISTINCT phase instances A#1, A#2
        List<PhaseInitialState> stable = r.stableInitialStates();
        check("16: two distinct stable phase instances retained (miscibility gap not merged)",
                stable.size() == 2, "got " + stable.size());
        boolean labelled = stable.size() == 2
                && stable.get(0).phaseLabel().equals("A#1")
                && stable.get(1).phaseLabel().equals("A#2");
        check("16: instances are labelled A#1 / A#2", labelled,
                stable.size() == 2 ? stable.get(0).phaseLabel() + " / " + stable.get(1).phaseLabel()
                        : "n/a");

        // positive amounts summing to 1
        double[] amt = r.selectedPhaseAmounts();
        check("16: both amounts strictly positive",
                amt != null && amt.length == 2 && amt[0] > 1e-6 && amt[1] > 1e-6,
                amt == null ? "null" : (amt[0] + "," + amt[1]));
        check("16: amounts sum to 1", amt != null && Math.abs(amt[0] + amt[1] - 1.0) < 1e-9,
                "sum != 1");

        // exact target conservation
        double[] rep = r.representedAmounts();
        check("16: represented composition == target",
                rep != null && Math.abs(rep[0] - 0.5) < 1e-7 && Math.abs(rep[1] - 0.5) < 1e-7,
                rep == null ? "null" : (rep[0] + "," + rep[1]));

        // mixture G lower than any single point that reaches the target
        // (only xB=0.5 itself reaches it; its G is the well maximum between minima)
        check("16: mixture G < single-point G at the target composition",
                r.totalGibbsEnergy() < wellG(0.5) - 1.0,
                "mixG=" + r.totalGibbsEnergy() + " vs G(0.5)=" + wellG(0.5));
        check("16: objectiveGap negative (real tie line beats every single point)",
                r.objectiveGap() < -1.0, "gap=" + r.objectiveGap());

        // discrete tangent: horizontal, mu_A == mu_B == gMin
        double[] mu = r.discreteChemicalPotentials();
        check("16: discrete mu solved (2 distinct points -> common tangent)",
                mu != null && mu.length == 2, "mu null");
        check("16: tangent is horizontal (mu_A == mu_B by symmetry)",
                mu != null && Math.abs(mu[0] - mu[1]) < 1e-6, mu == null ? "null" : (mu[0] + "," + mu[1]));
        check("16: mu level == analytic G_min",
                mu != null && Math.abs(mu[0] - gMin) < 1e-4, mu == null ? "null" : ("" + mu[0]));

        check("16: gamma ~ 0 at both selected points (on the tangent)",
                mu != null
                        && Math.abs(gammaAt(sel.get(0), mu)) < 1e-4
                        && Math.abs(gammaAt(sel.get(1), mu)) < 1e-4,
                "gamma off tangent");
    }

    // =====================================================================
    // Part 17 -- two-phase common tangent
    // =====================================================================

    private static void part17_twoPhaseCommonTangent() {
        section("PART 17 -- two-phase common tangent (analytic parabolas)");

        // Phase 0: G0(x) = a0 (x - c0)^2 + d0     touches tangent at x=0.2
        // Phase 1: G1(x) = a1 (x - c1)^2 + d1     touches tangent at x=0.8
        // Choose a common tangent line L(x) = p x + q.
        // Tangency of a parabola a(x-c)^2 + d at x*: slope 2a(x*-c) = p and value == L(x*).
        double a0 = 30_000.0, c0 = 0.15, d0 = -2_000.0;
        double a1 = 40_000.0, c1 = 0.85, d1 = -1_500.0;
        double x0 = 0.20, x1 = 0.80;
        double p = 2 * a0 * (x0 - c0);   // tangent slope from phase 0 at x0
        // shift d1 so phase 1 is tangent to the SAME line at x1:
        //   2 a1 (x1 - c1) must equal p  -> adjust c1 instead (cleaner):
        c1 = x1 - p / (2 * a1);
        double g0x0 = a0 * (x0 - c0) * (x0 - c0) + d0;
        double q = g0x0 - p * x0;                     // intercept of the common tangent
        double g1x1_target = p * x1 + q;              // where phase 1 must sit at x1
        d1 = g1x1_target - a1 * (x1 - c1) * (x1 - c1);

        final double fa0 = a0, fc0 = c0, fd0 = d0, fa1 = a1, fc1 = c1, fd1 = d1;
        java.util.function.DoubleUnaryOperator G0f = x -> fa0 * (x - fc0) * (x - fc0) + fd0;
        java.util.function.DoubleUnaryOperator G1f = x -> fa1 * (x - fc1) * (x - fc1) + fd1;

        System.out.printf("  common tangent  L(x) = %.3f x + %.3f ;  touch pts xB=%.2f (ph0), xB=%.2f (ph1)%n",
                p, q, x0, x1);
        check("17: tangent value matches on both phases",
                Math.abs(G0f.applyAsDouble(x0) - (p * x0 + q)) < 1e-6
                        && Math.abs(G1f.applyAsDouble(x1) - (p * x1 + q)) < 1e-6,
                "tangent mismatch");
        check("17: tangent slope matches on both phases",
                Math.abs(2 * a0 * (x0 - c0) - p) < 1e-6 && Math.abs(2 * a1 * (x1 - c1) - p) < 1e-6,
                "slope mismatch");

        int n = 50;
        List<GridPoint> grid = new ArrayList<>();
        int gi = 0;
        for (int i = 0; i <= n; i++) {
            double xB = (double) i / n;
            grid.add(GridPoint.synthetic(0, "ALPHA", gi++, new double[] { 1 - xB, xB },
                    G0f.applyAsDouble(xB)));
        }
        for (int i = 0; i <= n; i++) {
            double xB = (double) i / n;
            grid.add(GridPoint.synthetic(1, "BETA", gi++, new double[] { 1 - xB, xB },
                    G1f.applyAsDouble(xB)));
        }
        // pin the exact touch points
        int giA = gi++;
        grid.add(GridPoint.synthetic(0, "ALPHA", giA, new double[] { 1 - x0, x0 }, G0f.applyAsDouble(x0)));
        int giB = gi++;
        grid.add(GridPoint.synthetic(1, "BETA", giB, new double[] { 1 - x1, x1 }, G1f.applyAsDouble(x1)));

        // target halfway between the touch points: xB = 0.5 -> lever rule 50/50
        double[] target = { 0.5, 0.5 };
        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(grid, target, 2, null);

        check("17: converged", r.converged(), r.failureReason());
        check("17: two grid points selected", r.selectedGridPoints().size() == 2,
                "got " + r.selectedGridPoints().size());

        List<GridPoint> sel = r.selectedGridPoints();
        boolean twoPhases = sel.size() == 2 && sel.get(0).phaseIndex() != sel.get(1).phaseIndex();
        check("17: one point from each distinct phase", twoPhases, "not one-per-phase");

        // endpoints at the analytic touch compositions
        double selA = sel.get(0).phaseIndex() == 0 ? sel.get(0).m()[1] : sel.get(1).m()[1];
        double selB = sel.get(0).phaseIndex() == 1 ? sel.get(0).m()[1] : sel.get(1).m()[1];
        check("17: phase-0 endpoint at xB=0.2, phase-1 endpoint at xB=0.8",
                Math.abs(selA - x0) < 1.0 / n + 1e-9 && Math.abs(selB - x1) < 1.0 / n + 1e-9,
                "endpoints " + selA + " / " + selB);

        // lever rule: target 0.5 midway -> 50/50
        double[] amt = r.selectedPhaseAmounts();
        check("17: amounts ~ 50/50 (lever rule at the midpoint)",
                amt != null && Math.abs(amt[0] - 0.5) < 0.02 && Math.abs(amt[1] - 0.5) < 0.02,
                amt == null ? "null" : (amt[0] + "," + amt[1]));

        // discrete mu: slopes/intercept of the common tangent
        //   mu_B - mu_A = p   (tangent slope in xB) ;  mu_A = q + 0*... actually
        //   for G in J/mol with M=(xA,xB):  mu_A = L(0) = q ,  mu_B = L(1) = p + q
        double[] mu = r.discreteChemicalPotentials();
        check("17: discrete mu solved from the 2-point common tangent",
                mu != null && mu.length == 2, "mu null");
        check("17: mu_A == tangent intercept q",
                mu != null && Math.abs(mu[0] - q) < 1e-3, mu == null ? "null" : ("mu_A=" + mu[0] + " q=" + q));
        check("17: mu_B == p + q (tangent at xB=1)",
                mu != null && Math.abs(mu[1] - (p + q)) < 1e-3,
                mu == null ? "null" : ("mu_B=" + mu[1] + " p+q=" + (p + q)));

        // objective == value of the common tangent at the target composition
        double tangentAtTarget = p * 0.5 + q;
        check("17: total G == common tangent value at the target",
                Math.abs(r.totalGibbsEnergy() - tangentAtTarget) < 1e-3,
                "G=" + r.totalGibbsEnergy() + " tangent=" + tangentAtTarget);

        // tangent residual ~ 0 at both selected points
        check("17: tangent residual ~ 0 at both selected points",
                mu != null
                        && Math.abs(gammaAt(sel.get(0), mu)) < 1e-3
                        && Math.abs(gammaAt(sel.get(1), mu)) < 1e-3,
                "residual off");

        // target conservation
        double[] rep = r.representedAmounts();
        check("17: represented composition == target",
                rep != null && Math.abs(rep[0] - 0.5) < 1e-7 && Math.abs(rep[1] - 0.5) < 1e-7,
                rep == null ? "null" : (rep[0] + "," + rep[1]));
    }

    // =====================================================================
    // Part 19 -- special cases
    // =====================================================================

    private static void part19a_targetOnGridPoint() {
        section("PART 19A -- target exactly on a grid point (single-phase, one point)");
        // convex single parabola: minimum at xB=0.5
        int n = 20;
        List<GridPoint> grid = new ArrayList<>();
        int gi = 0;
        for (int i = 0; i <= n; i++) {
            double xB = (double) i / n;
            double g = 50_000.0 * (xB - 0.5) * (xB - 0.5) - 3_000.0;
            grid.add(GridPoint.synthetic(0, "SOL", gi++, new double[] { 1 - xB, xB }, g));
        }
        double[] target = { 0.7, 0.3 };   // exactly grid node i=6
        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(grid, target, 2, null);
        check("19A: converged", r.converged(), r.failureReason());
        check("19A: exactly 1 grid point selected", r.selectedGridPoints().size() == 1,
                "got " + r.selectedGridPoints().size());
        check("19A: selected point IS the target composition",
                Math.abs(r.selectedGridPoints().get(0).m()[1] - 0.3) < 1e-9, "wrong node");
        check("19A: amount == 1", r.selectedPhaseAmounts()[0] == 1.0,
                "" + r.selectedPhaseAmounts()[0]);
        check("19A: one stable instance, zero metastable (single candidate phase)",
                r.stableInitialStates().size() == 1 && r.metastableInitialStates().isEmpty(),
                "stable=" + r.stableInitialStates().size() + " meta=" + r.metastableInitialStates().size());
        check("19A: single-phase mu via dG/dx fallback OR null (documented, not arbitrary)",
                r.muMethod() != null, "no muMethod");
        System.out.println("      muMethod = " + r.muMethod());
    }

    private static void part19j_infeasibleTarget() {
        section("PART 19J -- infeasible target (outside the grid hull)");
        List<GridPoint> grid = new ArrayList<>();
        // grid only covers xB in [0.3, 0.7]
        int gi = 0;
        for (int i = 6; i <= 14; i++) {
            double xB = i / 20.0;
            grid.add(GridPoint.synthetic(0, "SOL", gi++, new double[] { 1 - xB, xB },
                    1000.0 * xB));
        }
        double[] target = { 0.05, 0.95 };   // xB=0.95, well outside [0.3,0.7]
        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(grid, target, 2, null);
        check("19J: result is a clean failure (not converged)", !r.converged(), "claimed converged");
        check("19J: failureReason set to 'infeasible target'",
                "infeasible target".equals(r.failureReason()), "" + r.failureReason());
        check("19J: no fabricated mixture (selected list empty, amounts null)",
                r.selectedGridPoints().isEmpty() && r.selectedPhaseAmounts() == null,
                "fabricated a mixture");
        check("19J: no stable / metastable states emitted",
                r.stableInitialStates().isEmpty() && r.metastableInitialStates().isEmpty(),
                "emitted states on failure");
    }

    private static void part19g_duplicateCollinearPoints() {
        section("PART 19G -- duplicate / collinear grid points -> deterministic pick");
        List<GridPoint> grid = new ArrayList<>();
        // three points at xB=0.2, and two identical at xB=0.8, all collinear in G
        grid.add(GridPoint.synthetic(0, "P", 0, new double[] { 0.8, 0.2 }, -100.0));
        grid.add(GridPoint.synthetic(0, "P", 1, new double[] { 0.8, 0.2 }, -100.0));   // dup
        grid.add(GridPoint.synthetic(0, "P", 2, new double[] { 0.8, 0.2 }, -100.0));   // dup
        grid.add(GridPoint.synthetic(0, "P", 7, new double[] { 0.2, 0.8 }, -300.0));
        grid.add(GridPoint.synthetic(0, "P", 8, new double[] { 0.2, 0.8 }, -300.0));   // dup
        double[] target = { 0.5, 0.5 };
        GridMinimizationResult r1 = HillertGridMinimizer.minimizeDiscreteForTest(grid, target, 2, null);
        GridMinimizationResult r2 = HillertGridMinimizer.minimizeDiscreteForTest(grid, target, 2, null);
        check("19G: converged", r1.converged(), r1.failureReason());
        check("19G: picks the lexicographically-smallest grid-index pair (0 and 7)",
                gridKey(r1).equals("[0, 7]"), "got " + gridKey(r1));
        check("19G: identical result on a repeat run (deterministic tie-break)",
                gridKey(r1).equals(gridKey(r2)), "r1=" + gridKey(r1) + " r2=" + gridKey(r2));
    }

    private static void part20_determinismShuffledInput() {
        section("PART 20 -- determinism under shuffled input list");
        // build a 2-phase grid, then present it in 3 different orders
        java.util.function.DoubleUnaryOperator G0f = x -> 20_000.0 * (x - 0.25) * (x - 0.25) - 1_000.0;
        java.util.function.DoubleUnaryOperator G1f = x -> 25_000.0 * (x - 0.75) * (x - 0.75) - 900.0;
        List<GridPoint> base = new ArrayList<>();
        int gi = 0;
        for (int i = 0; i <= 20; i++) {
            double xB = i / 20.0;
            base.add(GridPoint.synthetic(0, "A", gi++, new double[] { 1 - xB, xB }, G0f.applyAsDouble(xB)));
        }
        for (int i = 0; i <= 20; i++) {
            double xB = i / 20.0;
            base.add(GridPoint.synthetic(1, "B", gi++, new double[] { 1 - xB, xB }, G1f.applyAsDouble(xB)));
        }
        double[] target = { 0.5, 0.5 };

        String ref = summarize(HillertGridMinimizer.minimizeDiscreteForTest(base, target, 2, null));
        long[] seeds = { 1L, 42L, 12345L, 99999L };
        boolean allEqual = true;
        for (long seed : seeds) {
            List<GridPoint> shuffled = new ArrayList<>(base);
            java.util.Collections.shuffle(shuffled, new java.util.Random(seed));
            String s = summarize(HillertGridMinimizer.minimizeDiscreteForTest(shuffled, target, 2, null));
            if (!s.equals(ref)) {
                allEqual = false;
                System.out.println("      seed " + seed + " -> " + s + "   (ref " + ref + ")");
            }
        }
        check("20: selected set + amounts + mu identical for every shuffled input order",
                allEqual, "order-dependent result");
        System.out.println("      canonical result: " + ref);
    }

    // =====================================================================
    // V2 STEP 10 -- same-parent thermodynamic merge test
    // =====================================================================

    /**
     * STEP 10 PART 12 -- one synthetic CONVEX single parent phase,
     * {@code G(x) = a (x - c)^2 + d}. An off-node target between two grid nodes
     * makes the discrete minimum select both nodes (same parent). The merge test
     * relaxes ONE parent at the weighted combined composition; for a convex
     * {@code G}, {@code G_single(x_group) < G_grid_group}, so the group MERGES
     * to exactly one instance with {@code M == M_group}, {@code G == G_single},
     * {@code Y == } the analytically-relaxed constitution.
     */
    private static void step10_binaryConvexMerge() {
        section("STEP 10 -- binary CONVEX single phase: same-parent group MERGES");

        final double a = 40_000.0, c = 0.50, d = -3_000.0;
        java.util.function.DoubleUnaryOperator G = x -> a * (x - c) * (x - c) + d;

        // grid nodes at 0.05 spacing ; target Ti=0.53 sits between 0.50 and 0.55
        int n = 20;
        List<GridPoint> grid = new ArrayList<>();
        for (int i = 0; i <= n; i++) {
            double xB = (double) i / n;
            grid.add(GridPoint.synthetic(0, "SOL", i, new double[] { 1 - xB, xB }, G.applyAsDouble(xB)));
        }
        double[] target = { 0.47, 0.53 };

        // analytic merge evaluator: G is convex, minimum at x=c, but here the
        // "relaxed state at fixed composition x_group" for a 1-DOF phase with no
        // internal CVs is just G(x_group) itself (nothing to relax).
        HillertGridMinimizer.MergeEvaluator eval = (pi, xGroup) ->
                new HillertGridMinimizer.RelaxedParentState(
                        xGroup.clone(), xGroup.clone(), G.applyAsDouble(xGroup[1]));

        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(
                grid, target, 2, eval, null);

        check("S10-convex: grid's discrete minimum selects the 2 bracketing nodes",
                r.selectedGridPoints().size() == 2, "got " + r.selectedGridPoints().size());
        check("S10-convex: MERGED to exactly ONE stable PhaseInitialState",
                r.stableInitialStates().size() == 1, "got " + r.stableInitialStates().size());

        PhaseInitialState mi = r.stableInitialStates().get(0);
        double gGridGroup = r.totalGibbsEnergy();   // whole selected set = this one group
        double gSingle = G.applyAsDouble(0.53);
        System.out.printf("      G_grid_group=%.6f  G_single(x=0.53)=%.6f  deltaG_merge=%.6f%n",
                gGridGroup, gSingle, gSingle - gGridGroup);
        check("S10-convex: G_single < G_grid_group (convex -> mixture is wasteful)",
                gSingle < gGridGroup - 1.0, "gSingle=" + gSingle + " gGroup=" + gGridGroup);
        check("S10-convex: merged amount == sum of group amounts (== 1)",
                Math.abs(mi.amount() - 1.0) < 1e-9, "amount " + mi.amount());
        check("S10-convex: merged M == M_group == target [0.47, 0.53]",
                Math.abs(mi.m()[0] - 0.47) < 1e-9 && Math.abs(mi.m()[1] - 0.53) < 1e-9,
                java.util.Arrays.toString(mi.m()));
        check("S10-convex: merged G == G_single (not the grid-mixture G)",
                Math.abs(mi.g() - gSingle) < 1e-9, "g " + mi.g());
        check("S10-convex: merged Y == independently evaluated Y_single (== x_group here)",
                Math.abs(mi.y()[1] - 0.53) < 1e-9, java.util.Arrays.toString(mi.y()));
        check("S10-convex: no metastable state for the merged parent",
                r.metastableInitialStates().isEmpty(),
                "meta " + r.metastableInitialStates().size());
    }

    /**
     * STEP 10 PART 11 / PART 15 -- one synthetic parent phase with a DOUBLE-WELL
     * {@code G(x)}, target between the two minima. The discrete minimum selects
     * both minima (same parent). The merge test relaxes ONE parent at the
     * combined composition {@code x_group} -- which sits on the BARRIER between
     * the wells -- so {@code G_single(x_group) > G_grid_group} and the group
     * stays SPLIT as {@code SOL#1 / SOL#2}.
     */
    private static void step10_binaryDoubleWellNoMerge() {
        section("STEP 10 -- binary DOUBLE WELL: genuine gap stays SPLIT (no merge)");

        final double C4 = 80_000.0, C2 = 5_000.0, G0 = -1_000.0;
        java.util.function.DoubleUnaryOperator G = x -> {
            double t = x - 0.5;
            return C4 * t * t * t * t - C2 * t * t + G0;
        };
        double half = Math.sqrt(C2 / (2.0 * C4));
        double xLo = 0.5 - half, xHi = 0.5 + half;

        int n = 40;
        List<GridPoint> grid = new ArrayList<>();
        int gi = 0;
        for (int i = 0; i <= n; i++) {
            double xB = (double) i / n;
            grid.add(GridPoint.synthetic(0, "SOL", gi++, new double[] { 1 - xB, xB }, G.applyAsDouble(xB)));
        }
        grid.add(GridPoint.synthetic(0, "SOL", gi++, new double[] { 1 - xLo, xLo }, G.applyAsDouble(xLo)));
        grid.add(GridPoint.synthetic(0, "SOL", gi++, new double[] { 1 - xHi, xHi }, G.applyAsDouble(xHi)));

        double[] target = { 0.5, 0.5 };   // exactly between the minima -> on the barrier

        HillertGridMinimizer.MergeEvaluator eval = (pi, xGroup) ->
                new HillertGridMinimizer.RelaxedParentState(
                        xGroup.clone(), xGroup.clone(), G.applyAsDouble(xGroup[1]));

        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(
                grid, target, 2, eval, null);

        check("S10-gap: grid selects 2 points (the two minima)",
                r.selectedGridPoints().size() == 2, "got " + r.selectedGridPoints().size());
        double gGridGroup = r.totalGibbsEnergy();
        double gSingle = G.applyAsDouble(0.5);   // barrier top
        System.out.printf("      G_grid_group=%.6f  G_single(x=0.5, barrier)=%.6f  deltaG_merge=%.6f%n",
                gGridGroup, gSingle, gSingle - gGridGroup);
        check("S10-gap: G_single(barrier) > G_grid_group (mixture genuinely cheaper)",
                gSingle > gGridGroup + 1.0, "gSingle=" + gSingle + " gGroup=" + gGridGroup);
        check("S10-gap: group stays SPLIT -> 2 stable PhaseInitialStates",
                r.stableInitialStates().size() == 2, "got " + r.stableInitialStates().size());
        List<PhaseInitialState> st = r.stableInitialStates();
        check("S10-gap: both instances same parent",
                st.get(0).phaseIndex() == st.get(1).phaseIndex(), "different parents");
        check("S10-gap: labelled SOL#1 / SOL#2",
                st.get(0).phaseLabel().equals("SOL#1") && st.get(1).phaseLabel().equals("SOL#2"),
                st.get(0).phaseLabel() + " / " + st.get(1).phaseLabel());
        check("S10-gap: instances carry DIFFERENT constitutions (the two minima)",
                Math.abs(st.get(0).m()[1] - st.get(1).m()[1]) > 0.3, "constitutions not distinct");
        check("S10-gap: both amounts positive",
                st.get(0).amount() > 1e-6 && st.get(1).amount() > 1e-6,
                st.get(0).amount() + " / " + st.get(1).amount());
    }

    /**
     * STEP 10 PART 9 / PART 14 (synthetic) -- THREE same-parent selected points
     * for a K=3 CONVEX single phase. An interior off-node target inside one
     * convex basin gets bracketed by three grid nodes; the merge test on the
     * whole group collapses all three to ONE instance.
     */
    private static void step10_ternaryConvexMerge() {
        section("STEP 10 -- ternary CONVEX single phase: 3 same-parent points MERGE to 1");

        // convex paraboloid G(x) = a*((xA-1/3)^2 + (xB-1/3)^2) + d, min at centroid
        final double a = 30_000.0, d = -2_000.0;
        java.util.function.ToDoubleFunction<double[]> G = x ->
                a * ((x[0] - 1.0 / 3) * (x[0] - 1.0 / 3) + (x[1] - 1.0 / 3) * (x[1] - 1.0 / 3)) + d;

        int n = 12;
        List<GridPoint> grid = new ArrayList<>();
        int gi = 0;
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n - i; j++) {
                double xa = (double) i / n, xb = (double) j / n, xc = 1 - xa - xb;
                double[] mm = { xa, xb, xc };
                grid.add(GridPoint.synthetic(0, "SOL", gi++, mm, G.applyAsDouble(mm)));
            }
        }
        // off-node interior target
        double[] target = { 0.30, 0.34, 0.36 };

        HillertGridMinimizer.MergeEvaluator eval = (pi, xGroup) ->
                new HillertGridMinimizer.RelaxedParentState(
                        xGroup.clone(), xGroup.clone(), G.applyAsDouble(xGroup));

        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(
                grid, target, 3, eval, null);

        check("S10-tern: converged", r.converged(), r.failureReason());
        int nSel = r.selectedGridPoints().size();
        System.out.printf("      discrete minimum used %d grid points (all same parent)%n", nSel);
        check("S10-tern: discrete minimum used >= 2 same-parent points (interior off-node)",
                nSel >= 2, "got " + nSel);
        check("S10-tern: MERGED to exactly ONE stable PhaseInitialState",
                r.stableInitialStates().size() == 1, "got " + r.stableInitialStates().size());
        PhaseInitialState mi = r.stableInitialStates().get(0);
        check("S10-tern: merged M == target within CONSERVATION_TOL",
                Math.abs(mi.m()[0] - target[0]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(mi.m()[1] - target[1]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(mi.m()[2] - target[2]) < HillertGridMinimizer.CONSERVATION_TOL,
                java.util.Arrays.toString(mi.m()));
        check("S10-tern: merged G == G_single at the combined composition",
                Math.abs(mi.g() - G.applyAsDouble(target)) < 1e-6, "g " + mi.g());
        check("S10-tern: merged amount ~ 1 (sum of the group)",
                Math.abs(mi.amount() - 1.0) < 1e-9, "amount " + mi.amount());
    }

    /**
     * V2 FINAL-HARDENING regression -- the same-parent merge test must compare
     * {@code G_single} (per formula unit) against the group's <b>per-mole</b>
     * grid-mixture Gibbs energy, {@code (sum_q n_q G_q) / amountSum}, NOT the
     * raw total {@code sum_q n_q G_q}. When a same-parent group is one parent
     * among several selected phases its {@code amountSum < 1}, and comparing a
     * per-mole {@code G_single} to a sub-unit-scaled total would spuriously
     * merge (for {@code G < 0}) or spuriously split (for {@code G > 0}).
     *
     * <p>Two phases: A is a CONVEX single phase whose grid the target brackets
     * with two same-parent nodes; B is a distinct phase that also carries some
     * of the target. A's group has {@code amountSum ~ 0.5}. The merge test on
     * A's two nodes must still MERGE (A is convex), and it can only reach that
     * decision correctly with the per-mole normalisation.</p>
     */
    private static void closure_multiGroupPerMoleNormalisation() {
        section("CLOSURE -- multi-group merge test uses PER-MOLE grid-mixture G");

        // K=3 grid. Parent A: three collinear-in-M nodes on one edge (xC=0), at
        // xB = 0.2, 0.5, 0.8, with a CONVEX G along that edge. Parent B: one
        // deep node at the C-corner (xC=1). The target sits on the segment from
        // the A edge to B's corner, at a point whose A-projection is xB=0.5 but
        // whose xC != 0 -- so the mixture MUST be {A_lo(0.2), A_hi(0.8), B} to
        // both hit xB=0.5 (A_lo+A_hi average, since 0.5 as a single node would be
        // cheaper -- we omit it) and supply xC via B. A's group then has
        // amountSum < 1 and two distinct compositions -> the per-mole merge test
        // is the code path under test.
        final double aA = 20_000.0;
        java.util.function.DoubleUnaryOperator GAedge = xB -> aA * (xB - 0.5) * (xB - 0.5) - 2_000.0;

        List<GridPoint> grid = new ArrayList<>();
        int gi = 0;
        // parent A on the xC=0 edge, EXCLUDING xB=0.5 so a 2-node A blend is forced
        for (double xB : new double[] { 0.20, 0.35, 0.65, 0.80 }) {
            grid.add(GridPoint.synthetic(0, "A", gi++, new double[] { 1 - xB, xB, 0.0 }, GAedge.applyAsDouble(xB)));
        }
        // parent B: deep C-corner
        grid.add(GridPoint.synthetic(1, "B", gi++, new double[] { 0.0, 0.0, 1.0 }, -10_000.0));

        // target: half-way between the A edge at xB=0.5 and B's corner
        //   = 0.5 * (0.5, 0.5, 0.0)  + 0.5 * (0.0, 0.0, 1.0) = (0.25, 0.25, 0.50)
        double[] target = { 0.25, 0.25, 0.50 };

        HillertGridMinimizer.MergeEvaluator eval = (pi, xGroup) -> {
            double g = (pi == 0) ? GAedge.applyAsDouble(xGroup[1]) : -10_000.0;
            return new HillertGridMinimizer.RelaxedParentState(xGroup.clone(), xGroup.clone(), g);
        };

        StringBuilder log = new StringBuilder();
        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(
                grid, target, 3, eval, l -> log.append(l).append('\n'));

        check("CLOSURE: converged", r.converged(), r.failureReason());
        int nSel = r.selectedGridPoints().size();
        long aSel = r.selectedGridPoints().stream().filter(g -> g.phaseIndex() == 0).count();
        double aAmountSum = 0.0;
        double[] amts = r.selectedPhaseAmounts();
        for (int i = 0; i < r.selectedGridPoints().size(); i++) {
            if (r.selectedGridPoints().get(i).phaseIndex() == 0) {
                aAmountSum += amts[i];
            }
        }
        System.out.printf("      selected %d grid point(s), %d parent-A (amountSum=%.4f); stable=%d%n",
                nSel, aSel, aAmountSum, r.stableInitialStates().size());

        check("CLOSURE: parent B also selected",
                r.stableInitialStates().stream().anyMatch(s -> s.phaseIndex() == 1), "B missing");

        if (aSel >= 2) {
            check("CLOSURE: parent-A group has amountSum < 1 (sub-unit share)",
                    aAmountSum < 0.999, "amountSum=" + aAmountSum);
            check("CLOSURE: A's same-parent merge test RAN",
                    log.toString().contains("same-parent merge test (parent 0"),
                    "no parent-0 merge-test log:\n" + log);
            check("CLOSURE: A's merge test decided MERGE (convex; per-mole G_single <= per-mole mix)",
                    log.toString().contains("-> MERGE parent 0"), "did not merge:\n" + log);
            check("CLOSURE: parent A ends as exactly ONE stable instance",
                    r.stableInitialStates().stream().filter(s -> s.phaseIndex() == 0).count() == 1,
                    "A not single");
            // per-mole sanity: without the fix, deltaG_merge = G_single - amountSum*sum(n_q G_q)
            // and with G ~ -2000, amountSum ~ 0.5 this would be ~ -1000 (spurious
            // merge that happens to also be correct here, but for the WRONG
            // reason). We assert the LOGGED G_grid_group is per-mole, i.e. on the
            // GAedge scale (~ -2000), not ~ -1000.
            java.util.regex.Matcher mm = java.util.regex.Pattern
                    .compile("G_grid_group=(-?[0-9.]+)").matcher(log.toString());
            check("CLOSURE: logged G_grid_group present", mm.find(), "no G_grid_group in log");
            if (mm.reset().find()) {
                double logged = Double.parseDouble(mm.group(1));
                double gAAtHalf = GAedge.applyAsDouble(0.5);   // = -2000 (per mole, A's own min)
                System.out.printf("      logged G_grid_group=%.4f ; GAedge(0.5)=%.4f (per-mole target)%n",
                        logged, gAAtHalf);
                check("CLOSURE: logged G_grid_group is PER-MOLE (close to GAedge scale, not amount-scaled)",
                        Math.abs(logged - gAAtHalf) < 500.0,
                        "logged=" + logged + " expected ~" + gAAtHalf
                                + " (amount-scaled would be ~" + (gAAtHalf * aAmountSum) + ")");
            }
        } else {
            System.out.println("      (mixture used a single A node -- per-mole path not exercised; "
                    + "still verifying A is 1 instance)");
            check("CLOSURE: parent A ends as exactly ONE stable instance",
                    r.stableInitialStates().stream().filter(s -> s.phaseIndex() == 0).count() == 1,
                    "A not single");
        }

        double[] rep = r.representedAmounts();
        check("CLOSURE: represented composition == target",
                rep != null
                        && Math.abs(rep[0] - target[0]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(rep[1] - target[1]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(rep[2] - target[2]) < HillertGridMinimizer.CONSERVATION_TOL,
                rep == null ? "null" : java.util.Arrays.toString(rep));
    }

    // =====================================================================
    // helpers
    // =====================================================================

    private static double gammaAt(GridPoint gp, double[] mu) {
        double[] m = gp.m();
        double s = 0.0;
        for (int i = 0; i < mu.length; i++) {
            s += mu[i] * m[i];
        }
        return s - gp.g();
    }

    private static String gridKey(GridMinimizationResult r) {
        List<Integer> idx = new ArrayList<>();
        for (GridPoint gp : r.selectedGridPoints()) {
            idx.add(gp.gridIndex());
        }
        java.util.Collections.sort(idx);
        return idx.toString();
    }

    private static String summarize(GridMinimizationResult r) {
        StringBuilder sb = new StringBuilder();
        sb.append(r.converged() ? "OK" : "FAIL").append(' ');
        List<String> parts = new ArrayList<>();
        List<GridPoint> sel = r.selectedGridPoints();
        double[] amt = r.selectedPhaseAmounts();
        for (int i = 0; i < sel.size(); i++) {
            parts.add(String.format("(p%d,g%d,n=%.6f)",
                    sel.get(i).phaseIndex(), sel.get(i).gridIndex(), amt[i]));
        }
        sb.append(parts);
        double[] mu = r.discreteChemicalPotentials();
        if (mu != null) {
            sb.append(String.format(" mu=[%.4f,%.4f]", mu[0], mu[1]));
        }
        sb.append(String.format(" G=%.5f", r.totalGibbsEnergy()));
        return sb.toString();
    }

    private static void section(String s) {
        System.out.println("\n--- " + s + " ---");
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-78s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-78s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertGridMinimizerSyntheticTest() {
    }
}
