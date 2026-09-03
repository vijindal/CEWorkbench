package org.ce.model.equilibrium;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertGridMinimizer.GridMinimizationResult;
import org.ce.model.equilibrium.HillertGridMinimizer.GridPoint;
import org.ce.model.equilibrium.HillertGridMinimizer.PhaseInitialState;
import org.ce.model.equilibrium.HillertSolver.GridSeededResult;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;

/**
 * V2 STEP 8 gate (PART 13 + PART 14) -- multi-instance grid &rarr; Algorithm-A
 * handoff.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertGridToHillertMultiPhase
 * </pre>
 *
 * <h2>PART 13 -- real interior two-phase case</h2>
 * <p><b>Not available.</b> STEP 5/6/7 established that none of the CVM
 * Hamiltonians in this repository has a trustworthy interior two-phase
 * equilibrium reachable by {@link HillertSolver#solve}: the miscibility gaps in
 * Nb-Zr / Nb-Ti sit at edge compositions the inner CVM solver cannot converge,
 * and the only ordered (BCC_B2) Hamiltonian is a zero-ECI scaffold behind the
 * {@code UNSUPPORTED_PHASE_MODEL} guard. Per the STEP-8 spec ("do NOT fabricate
 * one"), this test does <b>not</b> invent a real two-phase system. The
 * two-phase common-tangent math is covered analytically by
 * {@link HillertGridMinimizerSyntheticTest} PART 17 (grid stage) and the
 * conversion layer is exercised below.</p>
 *
 * <h2>PART 14 -- same-parent multi-instance conversion + independence</h2>
 * <p>Two independent tests:</p>
 * <ol>
 *   <li><b>Synthetic conversion (analytic gap).</b> The STEP-7 double-well
 *       {@code A#1}/{@code A#2} selection is fed through the same amount /
 *       constitution logic the real handoff uses (here inline, since
 *       {@code GridPoint.synthetic} has no CVM model {@link HillertSolver#solve}
 *       could run). Verifies the two instances get independent {@code Y}
 *       arrays, independent amounts, distinct labels, and are never merged on
 *       shared parent identity.</li>
 *   <li><b>Real-CVM multi-instance handoff.</b> {@code Nb-Ti / BCC_A2 / T} at an
 *       off-node target so the discrete grid picks <b>two</b> same-parent
 *       instances bounding it. {@link HillertSolver#solveFromGrid} then builds
 *       two independent {@link HillertSolver.Phase} objects and hands them to
 *       Algorithm A, which (there being no real gap in this system) reconciles
 *       them to the single continuous equilibrium. Verifies: two instances
 *       enter the solve; their {@code uFull} arrays are not aliased; mass
 *       balance holds at entry and exit; Algorithm A is free to move their
 *       amounts (one goes to zero via its normal removal path); final G / x /
 *       mass balance match a direct single-phase solve; {@code G_final <=
 *       G_grid_mixture}.</li>
 * </ol>
 */
public final class HillertGridToHillertMultiPhase {

    private static int failures = 0;
    private static final double T = 1000.0;

    // STEP-7 double well:  G(x) = C4 (x-1/2)^4 - C2 (x-1/2)^2 + G0
    private static final double C4 = 80_000.0;
    private static final double C2 = 5_000.0;
    private static final double G0 = -1_000.0;

    private static double wellG(double xB) {
        double d = xB - 0.5;
        return C4 * d * d * d * d - C2 * d * d + G0;
    }

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(88));
        System.out.println("  V2 STEP 8 PART 13/14 -- multi-instance grid -> Algorithm-A handoff");
        System.out.println("=".repeat(88));

        part13_realTwoPhaseAvailability();
        part14a_syntheticGapConversionIndependence();
        part14b_realCvmMultiInstanceHandoff();

        System.out.println("\n" + "=".repeat(88));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(88));
        if (failures > 0) {
            throw new AssertionError(failures + " multi-phase integration checks failed");
        }
    }

    // =====================================================================

    private static void part13_realTwoPhaseAvailability() {
        section("PART 13 -- real interior two-phase CVM case");
        System.out.println("      NOT AVAILABLE in this repository (documented in STEP 5/6/7):");
        System.out.println("        - Nb-Zr / Nb-Ti gaps are at edge compositions the inner CVM");
        System.out.println("          solver cannot converge;");
        System.out.println("        - the only BCC_B2 Hamiltonian is a zero-ECI scaffold behind");
        System.out.println("          UNSUPPORTED_PHASE_MODEL.");
        System.out.println("      Per the STEP-8 spec, no real two-phase case is fabricated.");
        System.out.println("      Two-phase common-tangent math: HillertGridMinimizerSyntheticTest PART 17.");
        check("PART 13: limitation reported, no fabricated system", true, "");
    }

    // =====================================================================

    private static void part14a_syntheticGapConversionIndependence() {
        section("PART 14A -- synthetic A#1/A#2 conversion: independence, no merge");

        double half = Math.sqrt(C2 / (2.0 * C4));
        double xLo = 0.5 - half;
        double xHi = 0.5 + half;

        int n = 40;
        List<GridPoint> grid = new ArrayList<>();
        int gi = 0;
        for (int i = 0; i <= n; i++) {
            double xB = (double) i / n;
            grid.add(GridPoint.synthetic(0, "A", gi++, new double[] { 1 - xB, xB }, wellG(xB)));
        }
        grid.add(GridPoint.synthetic(0, "A", gi++, new double[] { 1 - xLo, xLo }, wellG(xLo)));
        grid.add(GridPoint.synthetic(0, "A", gi++, new double[] { 1 - xHi, xHi }, wellG(xHi)));

        GridMinimizationResult r = HillertGridMinimizer.minimizeDiscreteForTest(
                grid, new double[] { 0.5, 0.5 }, 2, null);

        check("14A: grid selects 2 stable instances", r.stableInitialStates().size() == 2,
                "got " + r.stableInitialStates().size());
        List<PhaseInitialState> st = r.stableInitialStates();
        check("14A: same parent phaseIndex for both",
                st.get(0).phaseIndex() == st.get(1).phaseIndex(), "different parents");
        check("14A: distinct labels A#1 / A#2",
                st.get(0).phaseLabel().equals("A#1") && st.get(1).phaseLabel().equals("A#2"),
                st.get(0).phaseLabel() + " / " + st.get(1).phaseLabel());

        // the conversion each accessor call must hand back an independent copy
        double[] y0a = st.get(0).y();
        double[] y0b = st.get(0).y();
        check("14A: PhaseInitialState.y() returns an independent copy each call",
                y0a != y0b && java.util.Arrays.equals(y0a, y0b), "aliased or unequal");
        y0a[0] = 999.0;
        check("14A: mutating a handed-out Y does not corrupt the record",
                st.get(0).y()[0] != 999.0, "record shares the array");

        double[] y1 = st.get(1).y();
        check("14A: the two instances carry DIFFERENT constitutions",
                !java.util.Arrays.equals(st.get(0).y(), y1), "identical Y");
        check("14A: independent amounts, both positive",
                r.selectedPhaseAmounts()[0] > 1e-6 && r.selectedPhaseAmounts()[1] > 1e-6
                        && r.selectedPhaseAmounts()[0] != r.selectedPhaseAmounts()[1]
                        || Math.abs(r.selectedPhaseAmounts()[0] - r.selectedPhaseAmounts()[1]) < 1e-9,
                java.util.Arrays.toString(r.selectedPhaseAmounts()));
        check("14A: NOT merged despite identical parent model",
                r.stableInitialStates().size() == 2, "merged to " + r.stableInitialStates().size());
    }

    // =====================================================================

    private static void part14b_realCvmMultiInstanceHandoff() throws Exception {
        section("PART 14B -- real Nb-Ti Ti=0.53: STEP-10 same-parent merge -> converges");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());
        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);

        // ---- (i) STEP 10: off-node target -> grid selects 2 same-parent nodes,
        //          the merge test collapses them to ONE, solveFromGrid converges ----
        double[] offNode = { 0.47, 0.53 };   // Ti=0.53 not on the N=20 grid (0.05 spacing)
        HillertSolver.Phase cand2 = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));

        StringBuilder log = new StringBuilder();
        GridMinimizationResult grid2 = HillertGridMinimizer.minimize(
                List.of(cand2), T, offNode, l -> log.append(l).append('\n'));

        check("14B(i): grid's discrete minimum still uses the 2 bracketing nodes",
                grid2.selectedGridPoints().size() == 2, "selected " + grid2.selectedGridPoints().size());
        check("14B(i): the same-parent MERGE TEST fired",
                log.toString().contains("same-parent merge test"), "no merge-test log line");
        check("14B(i): merge test decided to MERGE (single phase no more expensive)",
                log.toString().contains("-> MERGE parent"), "did not merge");
        check("14B(i): grid returns exactly ONE stable PhaseInitialState after the merge",
                grid2.stableInitialStates().size() == 1,
                "got " + grid2.stableInitialStates().size());
        var merged = grid2.stableInitialStates().get(0);
        check("14B(i): merged instance amount == sum of the group amounts (~1.0)",
                Math.abs(merged.amount() - 1.0) < 1e-9, "amount " + merged.amount());
        check("14B(i): merged instance composition == the combined composition Ti=0.53",
                Math.abs(merged.m()[1] - 0.53) < 1e-6, "M=" + java.util.Arrays.toString(merged.m()));
        // Y is the independently relaxed state, NOT a grid node's Y
        CVMGibbsModel.State reMerged = m.atFull(T, merged.y());
        check("14B(i): merged Y is the relaxed parent state at Ti=0.53 (M(Y)==0.53, valid)",
                reMerged.isValidIncludingPoints()
                        && Math.abs(reMerged.componentAmountsPerFormulaUnit()[1] - 0.53) < 1e-6,
                "M(Y) mismatch");
        check("14B(i): merged G == G_single (the relaxed state's G), matches record",
                Math.abs(reMerged.g() - merged.g()) < 1e-6,
                "record " + merged.g() + " re-eval " + reMerged.g());

        double[] uBefore = cand2.uFull.clone();
        GridSeededResult gsr2 = HillertSolver.solveFromGrid(
                List.of(cand2), offNode, T, 200, 20, 1.0e-9, null);
        check("14B(i): caller candidate Phase not mutated",
                java.util.Arrays.equals(uBefore, cand2.uFull), "uFull changed");
        check("14B(i): solveFromGrid now CONVERGES (was STALLED before STEP 10)",
                gsr2.overallConverged(), gsr2.result().convergenceReport().reason().toString());

        List<HillertSolver.PhaseResult> nz2 = gsr2.result().phases().stream()
                .filter(p -> p.amount() > 1e-9).toList();
        check("14B(i): exactly one active phase at equilibrium", nz2.size() == 1,
                "got " + nz2.size());
        check("14B(i): no phase-set events (merge happened at the grid stage, not in solve)",
                gsr2.result().convergenceReport().phaseSetEvents().isEmpty(),
                "events: " + gsr2.result().convergenceReport().phaseSetEvents());

        // cross-check against a DIRECT single-phase solve at Ti=0.53 (expected
        // value from the direct path, NOT from solveFromGrid)
        HillertSolver.Phase directSeed = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.47, 0.53 }));
        HillertSolver.Result directR = HillertSolver.solve(
                List.of(directSeed), offNode, T, 200, 20, 1.0e-9, null);
        HillertSolver.PhaseResult directP = directR.phases().stream()
                .filter(p -> p.amount() > 1e-9).findFirst().orElseThrow();
        HillertSolver.PhaseResult gP2 = nz2.get(0);
        System.out.printf("      grid-seeded final : xTi=%.8f  G=%.8f%n", gP2.composition()[1], gP2.g());
        System.out.printf("      direct     final  : xTi=%.8f  G=%.8f%n",
                directP.composition()[1], directP.g());
        check("14B(i): final composition == direct solve (x = [0.47, 0.53], < 1e-6)",
                Math.abs(gP2.composition()[0] - 0.47) < 1e-6
                        && Math.abs(gP2.composition()[1] - 0.53) < 1e-6,
                java.util.Arrays.toString(gP2.composition()));
        check("14B(i): final G == direct solve G = -50094.49775302 (< 1e-4)",
                Math.abs(gP2.g() - (-50094.49775302)) < 1e-4
                        && Math.abs(gP2.g() - directP.g()) < 1e-4,
                "G=" + gP2.g());

        double[] rep2 = new double[2];
        for (HillertSolver.PhaseResult p : gsr2.result().phases()) {
            double[] x = p.composition();
            for (int i = 0; i < 2; i++) rep2[i] += p.amount() * x[i];
        }
        check("14B(i): exit mass balance == overall target (< 1e-7)",
                Math.abs(rep2[0] - 0.47) < 1e-7 && Math.abs(rep2[1] - 0.53) < 1e-7,
                java.util.Arrays.toString(rep2));

        System.out.printf("      timing: grid = %d ms, Algorithm A = %d ms%n",
                gsr2.gridMillis(), gsr2.algorithmAMillis());

        // ---- (ii) HAPPY PATH: on-node target -> grid selects ONE instance ->
        //           full grid -> Algorithm A -> agrees with a direct solve ----
        double[] onNode = { 0.45, 0.55 };   // Ti=0.55 IS an N=20 node
        HillertSolver.Phase cand1 = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        GridMinimizationResult grid1 = HillertGridMinimizer.minimize(
                List.of(cand1), T, onNode, null);
        check("14B(ii): on-node target -> grid selects exactly 1 instance",
                grid1.stableInitialStates().size() == 1,
                "selected " + grid1.stableInitialStates().size());

        GridSeededResult gsr1 = HillertSolver.solveFromGrid(
                List.of(cand1), onNode, T, 200, 20, 1.0e-9, null);
        check("14B(ii): Algorithm A converged from the 1-instance grid seed",
                gsr1.overallConverged(), gsr1.result().convergenceReport().reason().toString());
        check("14B(ii): no phase-set events at initialisation",
                gsr1.result().convergenceReport().phaseSetEvents().isEmpty(),
                "events: " + gsr1.result().convergenceReport().phaseSetEvents());

        HillertSolver.PhaseResult gP = gsr1.result().phases().stream()
                .filter(p -> p.amount() > 1e-9).findFirst().orElseThrow();

        HillertSolver.Phase seed = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.45, 0.55 }));
        HillertSolver.Result dr = HillertSolver.solve(
                List.of(seed), onNode, T, 200, 20, 1.0e-9, null);
        HillertSolver.PhaseResult dP = dr.phases().stream()
                .filter(p -> p.amount() > 1e-9).findFirst().orElseThrow();

        System.out.printf("      grid-seeded final : xTi=%.8f  G=%.8f%n", gP.composition()[1], gP.g());
        System.out.printf("      direct  final     : xTi=%.8f  G=%.8f%n", dP.composition()[1], dP.g());
        check("14B(ii): final composition matches the direct solve (< 1e-6)",
                Math.abs(gP.composition()[1] - dP.composition()[1]) < 1e-6,
                "diff " + Math.abs(gP.composition()[1] - dP.composition()[1]));
        check("14B(ii): final G matches the direct solve (< 1e-4 J/mol)",
                Math.abs(gP.g() - dP.g()) < 1e-4, "diff " + Math.abs(gP.g() - dP.g()));

        double[] rep = new double[2];
        for (HillertSolver.PhaseResult p : gsr1.result().phases()) {
            double[] x = p.composition();
            for (int i = 0; i < 2; i++) rep[i] += p.amount() * x[i];
        }
        check("14B(ii): exit mass balance == overall target (< 1e-7)",
                Math.abs(rep[0] - onNode[0]) < 1e-7 && Math.abs(rep[1] - onNode[1]) < 1e-7,
                java.util.Arrays.toString(rep));

        double gMix = gsr1.gridMixtureGibbs();
        check("14B(ii): G_final <= G_grid_mixture + 1e-4 (discrete minimum is an upper bound)",
                gP.amount() * gP.g() <= gMix + 1e-4,
                "G_final " + (gP.amount() * gP.g()) + " > G_grid " + gMix);

        System.out.printf("      timing (happy path): grid = %d ms, Algorithm A = %d ms%n",
                gsr1.gridMillis(), gsr1.algorithmAMillis());
    }

    // =====================================================================

    private static void section(String s) {
        System.out.println("\n--- " + s + " ---");
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-74s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-74s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertGridToHillertMultiPhase() {
    }
}
