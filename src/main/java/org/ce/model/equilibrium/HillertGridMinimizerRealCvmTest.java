package org.ce.model.equilibrium;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertGridMinimizer.GridMinimizationResult;
import org.ce.model.equilibrium.HillertGridMinimizer.GridPoint;
import org.ce.model.equilibrium.HillertGridMinimizer.PhaseInitialState;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;

/**
 * V2 STEP 7 gate -- the grid minimizer against a real V1 CVM phase
 * ({@code Nb-Ti / BCC_A2 / T}, and a ternary {@code Mo-Nb-Ta} smoke), end to
 * end through {@link HillertGridMinimizer#minimize}.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertGridMinimizerRealCvmTest
 * </pre>
 *
 * <p>Checks (STEP 7 PART 18):</p>
 * <ul>
 *   <li>grid generation is deterministic (two runs -> identical valid /
 *       skipped counts, identical selected {@code (phaseIndex,gridIndex)} set
 *       and amounts);</li>
 *   <li>every emitted grid point has finite {@code M} and {@code G} and a
 *       valid CVM state; invalid points are skipped, not repaired;</li>
 *   <li>the discrete minimum's represented composition equals the (normalised)
 *       target to {@link HillertGridMinimizer#CONSERVATION_TOL};</li>
 *   <li>a selected grid point's {@code Y} builds a {@link HillertSolver.Phase}
 *       whose state {@link HillertSolver} would consume <b>unchanged</b> --
 *       {@code phase.uFull} equals the grid point's {@code Y} bit-for-bit and
 *       {@code M(Y)} / {@code G(Y)} re-evaluate to the reported values;</li>
 *   <li>{@code M == x} for this V1 phase (STEP 7 PART 2 / PART 24).</li>
 * </ul>
 *
 * <p>This does NOT require a real miscibility gap and does NOT wire the result
 * into {@code HillertSolver.solve} -- that is a later step.</p>
 */
public final class HillertGridMinimizerRealCvmTest {

    private static int failures = 0;
    private static final double T = 1000.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(88));
        System.out.println("  V2 STEP 7 -- grid minimizer on a real V1 CVM phase");
        System.out.println("=".repeat(88));

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        binaryNbTi(builder);
        ternaryMoNbTaSmoke(builder);

        System.out.println("\n" + "=".repeat(88));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(88));
        if (failures > 0) {
            throw new AssertionError(failures + " grid-minimizer real-CVM checks failed");
        }
    }

    // =====================================================================

    private static void binaryNbTi(ModelSession.Builder builder) throws Exception {
        section("Nb-Ti / BCC_A2 / T  (K=2, single disordered phase)");

        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);

        HillertSolver.Phase phase = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        List<HillertSolver.Phase> phases = List.of(phase);

        double[] target = { 0.5, 0.5 };   // Ti = 0.5

        long t0 = System.nanoTime();
        GridMinimizationResult r1 = HillertGridMinimizer.minimize(phases, T, target, null);
        long t1 = System.nanoTime();
        GridMinimizationResult r2 = HillertGridMinimizer.minimize(phases, T, target, null);

        check("converged", r1.converged(), r1.failureReason());
        System.out.printf("      %d valid grid points, %d skipped, %.1f ms%n",
                r1.gridPointCount(), r1.skippedGridPointCount(), (t1 - t0) / 1e6);
        check("expected point count for K=2, N=" + HillertGridMinimizer.GRID_DIVISIONS
                        + " is N+1=" + (HillertGridMinimizer.GRID_DIVISIONS + 1)
                        + " minus skips",
                r1.gridPointCount() + r1.skippedGridPointCount()
                        == HillertGridMinimizer.GRID_DIVISIONS + 1,
                "got " + r1.gridPointCount() + " + " + r1.skippedGridPointCount());

        // determinism
        check("deterministic: identical valid/skipped counts on a repeat run",
                r1.gridPointCount() == r2.gridPointCount()
                        && r1.skippedGridPointCount() == r2.skippedGridPointCount(),
                "counts differ");
        check("deterministic: identical selected grid set + amounts on a repeat run",
                sameSelection(r1, r2), "selection differs");

        // represented composition == target
        double[] rep = r1.representedAmounts();
        check("represented composition == target (within CONSERVATION_TOL)",
                rep != null
                        && Math.abs(rep[0] - 0.5) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(rep[1] - 0.5) < HillertGridMinimizer.CONSERVATION_TOL,
                rep == null ? "null" : (rep[0] + "," + rep[1]));

        // every selected grid point: finite, valid, M == x, and Hillert-consumable unchanged
        for (int i = 0; i < r1.selectedGridPoints().size(); i++) {
            GridPoint gp = r1.selectedGridPoints().get(i);
            double[] y = gp.y();
            double[] mm = gp.m();
            check("selected[" + i + "]: Y width == ncf+K", y.length == phase.ncf + 2,
                    "width " + y.length);
            check("selected[" + i + "]: M and G finite",
                    Double.isFinite(gp.g()) && Double.isFinite(mm[0]) && Double.isFinite(mm[1]),
                    "non-finite");
            check("selected[" + i + "]: M sums to 1", Math.abs(mm[0] + mm[1] - 1.0) < 1e-9,
                    "sum " + (mm[0] + mm[1]));

            // M == x for a V1 disordered phase (STEP 7 PART 2 / 24)
            double[] x = gp.composition();
            check("selected[" + i + "]: M == x bit-for-bit (V1 disordered phase)",
                    mm[0] == x[0] && mm[1] == x[1], mm[0] + "," + mm[1] + " vs " + x[0] + "," + x[1]);

            // build a Phase from the grid point's Y and confirm Hillert would
            // consume it unchanged
            HillertSolver.Phase seeded = new HillertSolver.Phase(
                    "BCC_A2_seeded", s, m, Math.max(1e-6, r1.selectedPhaseAmounts()[i]), y);
            check("selected[" + i + "]: Phase.uFull == grid point Y bit-for-bit",
                    java.util.Arrays.equals(seeded.uFull, y), "uFull changed on construction");

            CVMGibbsModel.State reeval = m.atFull(T, y);
            double gRe = reeval.g();
            double[] mRe = reeval.componentAmountsPerFormulaUnit();
            check("selected[" + i + "]: G(Y) re-evaluates to the reported value",
                    Math.abs(gRe - gp.g()) < 1e-6, "reported " + gp.g() + " re-eval " + gRe);
            check("selected[" + i + "]: M(Y) re-evaluates to the reported value",
                    Math.abs(mRe[0] - mm[0]) < 1e-9 && Math.abs(mRe[1] - mm[1]) < 1e-9,
                    "M mismatch");
        }

        // metastable initial states: none here (only one candidate phase, and it is selected)
        check("no metastable states (single candidate phase, selected)",
                r1.metastableInitialStates().isEmpty(),
                "got " + r1.metastableInitialStates().size());

        // single-phase mu method is documented (not arbitrary)
        System.out.println("      muMethod = " + r1.muMethod()
                + "   mu = " + java.util.Arrays.toString(r1.discreteChemicalPotentials()));
        check("muMethod recorded", r1.muMethod() != null && !r1.muMethod().isBlank(), "blank");

        // an off-eutectic target still solves and conserves
        GridMinimizationResult r3 = HillertGridMinimizer.minimize(
                phases, T, new double[] { 0.7, 0.3 }, null);
        check("off-centre target (Ti=0.3) converges", r3.converged(), r3.failureReason());
        double[] rep3 = r3.representedAmounts();
        check("off-centre target conserved",
                rep3 != null && Math.abs(rep3[1] - 0.3) < HillertGridMinimizer.CONSERVATION_TOL,
                rep3 == null ? "null" : ("" + rep3[1]));
    }

    private static void ternaryMoNbTaSmoke(ModelSession.Builder builder) throws Exception {
        section("Mo-Nb-Ta / BCC_A2 / T  (K=3 smoke: grid builds, minimum reproducible)");

        ModelSession s = builder.build(new SystemId("Mo-Nb-Ta", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Mo-Nb-Ta", "BCC_A2", "T", s.cecEntry, null);

        HillertSolver.Phase phase = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.34, 0.33, 0.33 }));
        List<HillertSolver.Phase> phases = List.of(phase);

        double[] target = { 0.34, 0.33, 0.33 };

        long t0 = System.nanoTime();
        GridMinimizationResult r1 = HillertGridMinimizer.minimize(phases, T, target, null);
        long t1 = System.nanoTime();
        GridMinimizationResult r2 = HillertGridMinimizer.minimize(phases, T, target, null);

        check("K=3 converged", r1.converged(), r1.failureReason());
        System.out.printf("      %d valid grid points, %d skipped, %.0f ms  (expected ~%d nodes)%n",
                r1.gridPointCount(), r1.skippedGridPointCount(), (t1 - t0) / 1e6,
                (HillertGridMinimizer.GRID_DIVISIONS + 1) * (HillertGridMinimizer.GRID_DIVISIONS + 2) / 2);
        check("K=3 node budget: valid + skipped == (N+1)(N+2)/2",
                r1.gridPointCount() + r1.skippedGridPointCount()
                        == (HillertGridMinimizer.GRID_DIVISIONS + 1)
                        * (HillertGridMinimizer.GRID_DIVISIONS + 2) / 2,
                "got " + (r1.gridPointCount() + r1.skippedGridPointCount()));
        check("K=3 minimum reproducible (same selection + amounts)",
                sameSelection(r1, r2), "selection differs between runs");

        double[] rep = r1.representedAmounts();
        check("K=3 represented composition == target",
                rep != null
                        && Math.abs(rep[0] - target[0]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(rep[1] - target[1]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(rep[2] - target[2]) < HillertGridMinimizer.CONSERVATION_TOL,
                rep == null ? "null" : java.util.Arrays.toString(rep));

        for (PhaseInitialState st : r1.allInitialStates()) {
            check("K=3 state '" + st.phaseLabel() + "': Y finite, width ncf+3",
                    st.y().length == phase.ncf + 3 && allFinite(st.y()), "bad Y");
        }

        // ---- V2 STEP 10: the interior equiatomic target is one convex BCC_A2
        // basin, so the same-parent merge test collapses the (formerly 3)
        // bracketing grid instances to ONE. ----
        check("K=3 STEP-10: exactly ONE stable instance after the same-parent merge",
                r1.stableInitialStates().size() == 1,
                "got " + r1.stableInitialStates().size() + " instances (merge did not collapse them)");
        check("K=3 STEP-10: the single instance is unsuffixed 'BCC_A2' (not BCC_A2#1)",
                r1.stableInitialStates().get(0).phaseLabel().equals("BCC_A2"),
                "label " + r1.stableInitialStates().get(0).phaseLabel());
        PhaseInitialState merged = r1.stableInitialStates().get(0);
        check("K=3 STEP-10: merged M == target",
                Math.abs(merged.m()[0] - target[0]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(merged.m()[1] - target[1]) < HillertGridMinimizer.CONSERVATION_TOL
                        && Math.abs(merged.m()[2] - target[2]) < HillertGridMinimizer.CONSERVATION_TOL,
                java.util.Arrays.toString(merged.m()));
        // re-evaluate the merged Y: it must be the relaxed parent state at the target
        CVMGibbsModel.State reMerged = m.atFull(T, merged.y());
        check("K=3 STEP-10: merged Y is a valid relaxed state with M(Y)==target",
                reMerged.isValidIncludingPoints()
                        && Math.abs(reMerged.componentAmountsPerFormulaUnit()[0] - target[0]) < 1e-6
                        && Math.abs(reMerged.componentAmountsPerFormulaUnit()[1] - target[1]) < 1e-6,
                "M(Y) mismatch");
        check("K=3 STEP-10: merged G == the relaxed state's G (record consistent)",
                Math.abs(reMerged.g() - merged.g()) < 1e-6,
                "record " + merged.g() + " re-eval " + reMerged.g());

        // ---- full pipeline: solveFromGrid must now converge for this interior target ----
        HillertSolver.Phase cand = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.34, 0.33, 0.33 }));
        HillertSolver.GridSeededResult gsr = HillertSolver.solveFromGrid(
                List.of(cand), target, T, 200, 20, 1.0e-9, null);
        check("K=3 STEP-10: solveFromGrid CONVERGES for the interior equiatomic target",
                gsr.overallConverged(), gsr.result().convergenceReport().reason().toString());
        long nzTern = gsr.result().phases().stream().filter(p -> p.amount() > 1e-9).count();
        check("K=3 STEP-10: exactly one active phase at equilibrium", nzTern == 1,
                "got " + nzTern);
        // cross-check final G against a direct single-phase solve at the target
        HillertSolver.Phase directSeed = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(target));
        HillertSolver.Result directR = HillertSolver.solve(
                List.of(directSeed), target, T, 200, 20, 1.0e-9, null);
        var gPt = gsr.result().phases().stream().filter(p -> p.amount() > 1e-9).findFirst().orElseThrow();
        var dPt = directR.phases().stream().filter(p -> p.amount() > 1e-9).findFirst().orElseThrow();
        System.out.printf("      grid-seeded final G=%.8f ;  direct final G=%.8f%n", gPt.g(), dPt.g());
        check("K=3 STEP-10: grid-seeded final G matches the direct solve (< 1e-4)",
                Math.abs(gPt.g() - dPt.g()) < 1e-4, "diff " + Math.abs(gPt.g() - dPt.g()));
    }

    // =====================================================================

    private static boolean sameSelection(GridMinimizationResult a, GridMinimizationResult b) {
        if (a.selectedGridPoints().size() != b.selectedGridPoints().size()) {
            return false;
        }
        for (int i = 0; i < a.selectedGridPoints().size(); i++) {
            GridPoint ga = a.selectedGridPoints().get(i);
            GridPoint gb = b.selectedGridPoints().get(i);
            if (ga.phaseIndex() != gb.phaseIndex() || ga.gridIndex() != gb.gridIndex()) {
                return false;
            }
            if (Math.abs(a.selectedPhaseAmounts()[i] - b.selectedPhaseAmounts()[i]) > 1e-12) {
                return false;
            }
        }
        return true;
    }

    private static boolean allFinite(double[] v) {
        for (double d : v) {
            if (!Double.isFinite(d)) {
                return false;
            }
        }
        return true;
    }

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

    private HillertGridMinimizerRealCvmTest() {
    }
}
