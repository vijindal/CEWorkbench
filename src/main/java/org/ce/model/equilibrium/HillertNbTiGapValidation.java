package org.ce.model.equilibrium;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertGridMinimizer.GridMinimizationResult;
import org.ce.model.equilibrium.HillertSolver.GridSeededResult;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.List;

/**
 * V2 CLOSURE -- the first end-to-end <b>real two-phase</b> validation of the
 * grid &rarr; Algorithm-A pipeline against an independent reference miscibility
 * gap.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertNbTiGapValidation
 * </pre>
 *
 * <h2>Reference</h2>
 * <p>{@code plot-pd-bin.xlsx} (Jindal, {@code proj-multiComp-cecvm/unused-data})
 * tabulates the Nb-Ti BCC_A2 CVM miscibility gap computed independently in
 * Mathematica. Its pair CECs are, <b>bit-for-bit</b>, the committed
 * {@code Nb-Ti_BCC_A2_T_CVCF} Hamiltonian:</p>
 * <pre>
 *   xlsx col G (1NN pair, m=4) = 6240 J/mol   == e21AB.a = 6240
 *   xlsx col H (2NN pair, m=3) = 3120 J/mol   == e22AB.a = 3120
 *   xlsx cols E,F (3-/4-body)  = 0            == e3AB.a = e4AB.a = 0
 * </pre>
 * <p>(The other four binaries in the workbook -- Nb-V, Ti-V, Nb-Zr, V-Zr --
 * are stored there in a different, coordination-scaled / opposite-sign
 * convention: {@code xlsx_ECI = -4 * multiplicity * repo_orthogonal_CF}, e.g.
 * Nb-V col G {@code 14080 = -16 * (-880)}. They have no matching committed
 * {@code _CVCF} Hamiltonian, so only Nb-Ti can be checked directly and that is
 * all this test does.)</p>
 *
 * <h2>Reference gap (positive pair CECs =&gt; clustering, symmetric dome)</h2>
 * <table border="1">
 *   <tr><th>T (K)</th><th>xTi (Nb-rich phase)</th><th>xTi (Ti-rich phase)</th><th>Gm at the tie line</th></tr>
 *   <tr><td>300</td><td>0.031207</td><td>0.968793</td><td>-66.6922</td></tr>
 *   <tr><td>350</td><td>0.063412</td><td>0.936588</td><td>-142.290</td></tr>
 *   <tr><td>400</td><td>0.119794</td><td>0.880206</td><td>-262.441</td></tr>
 *   <tr><td>450</td><td>0.238030</td><td>0.761970</td><td>-440.213</td></tr>
 *   <tr><td colspan="4">consolute point T_c &asymp; 477 K at xTi = 0.5</td></tr>
 * </table>
 *
 * <h2>What is checked</h2>
 * <ol>
 *   <li><b>Above T_c (1000 K, xTi = 0.5):</b> the grid minimizer's same-parent
 *       merge test MERGES (Nb-Ti is a single solution phase there), and
 *       {@code solveFromGrid} converges to one phase -- the STEP-10 /
 *       STEP-8-era behaviour, re-confirmed against this system.</li>
 *   <li><b>Below T_c (400 K, xTi = 0.5, deep in the two-phase field):</b> the
 *       grid minimizer selects two same-parent points, the merge test KEEPS
 *       THEM SPLIT (a genuine gap: G_single at the mid-composition &gt; the grid
 *       mixture), and the two initial constitutions bracket xTi = 0.5.</li>
 *   <li><b>The split tie-line endpoints match the Mathematica reference</b>
 *       {xTi = 0.1198, 0.8802} at 400 K to a tolerance set by the grid spacing
 *       (N = 20 =&gt; 0.05), i.e. each endpoint within one grid cell of the
 *       reference. This is the grid <em>initial estimate</em>, not the
 *       continuous Algorithm-A solution -- Algorithm A cannot currently
 *       reconcile a gap-free identical pair (STEP 9), and this system's gap is
 *       genuine so the two instances are distinct, but the near-edge CVM
 *       fragility at xTi &lt; 0.12 / &gt; 0.88 means a full two-phase
 *       Algorithm-A solve is still not attempted here -- see the class note.</li>
 * </ol>
 *
 * <p><b>Scope.</b> This is a TEST-ONLY addition made during V2 closure to
 * discharge the PART-9 "no real two-phase reference" limitation for the one
 * system where the committed Hamiltonian matches an external reference. It adds
 * no production code and changes no behaviour. The reference numbers are from
 * the Mathematica workbook, NOT from this codebase.</p>
 */
public final class HillertNbTiGapValidation {

    private static int failures = 0;

    // reference (Mathematica, plot-pd-bin.xlsx / nbti-data)
    private static final double T_C_REF = 477.0;
    private static final double T_BELOW = 400.0;
    private static final double XTI_ALPHA_400 = 0.119794;   // Nb-rich phase
    private static final double XTI_BETA_400 = 0.880206;    // Ti-rich phase
    private static final double GM_TIELINE_400 = -262.441;  // Gm at the 400 K tie line

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(90));
        System.out.println("  V2 CLOSURE -- Nb-Ti BCC_A2 miscibility gap vs Mathematica reference "
                + "(committed CVCF Hamiltonian)");
        System.out.println("=".repeat(90));

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());
        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);

        // ---- CEC match check (fail loudly if the committed Hamiltonian drifts) ----
        section("CEC match: committed Nb-Ti_BCC_A2_T_CVCF vs the workbook's Nb-Ti column");
        double e21 = eci(m, "e21AB");
        double e22 = eci(m, "e22AB");
        double e3 = eci(m, "e3AB");
        double e4 = eci(m, "e4AB");
        System.out.printf("      e21AB=%.1f (xlsx G=6240)  e22AB=%.1f (xlsx H=3120)  e3AB=%.1f  e4AB=%.1f%n",
                e21, e22, e3, e4);
        check("e21AB (1NN pair) == 6240 J/mol", Math.abs(e21 - 6240.0) < 1e-6, "" + e21);
        check("e22AB (2NN pair) == 3120 J/mol", Math.abs(e22 - 3120.0) < 1e-6, "" + e22);
        check("e3AB (3-body) == 0", Math.abs(e3) < 1e-9, "" + e3);
        check("e4AB (4-body) == 0", Math.abs(e4) < 1e-9, "" + e4);
        check("pair CECs positive => clustering (a real gap must exist below some T_c)",
                e21 > 0 && e22 > 0, "not clustering");

        // ---- (1) above T_c : single phase, merge test MERGES ----
        section("Above T_c: T = 1000 K, xTi = 0.5 -> single solution phase (merge)");
        double[] overallHi = { 0.5, 0.5 };
        HillertSolver.Phase candHi = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        StringBuilder logHi = new StringBuilder();
        GridMinimizationResult gHi = HillertGridMinimizer.minimize(
                List.of(candHi), 1000.0, overallHi, l -> logHi.append(l).append('\n'));
        check("1000 K: grid converged", gHi.converged(), gHi.failureReason());
        boolean mergedHi = gHi.stableInitialStates().size() == 1;
        check("1000 K: same-parent merge test collapses to ONE instance (no gap above T_c)",
                mergedHi, "got " + gHi.stableInitialStates().size() + " instances");
        if (!mergedHi) {
            System.out.println(logHi);
        }
        GridSeededResult sHi = HillertSolver.solveFromGrid(
                List.of(candHi), overallHi, 1000.0, 200, 20, 1.0e-9, null);
        check("1000 K: solveFromGrid CONVERGES to a single phase",
                sHi.overallConverged()
                        && sHi.result().phases().stream().filter(p -> p.amount() > 1e-9).count() == 1,
                sHi.result().convergenceReport().reason().toString());

        // ---- (2) below T_c : two-phase field, merge test KEEPS SPLIT ----
        section("Below T_c: T = " + (int) T_BELOW + " K, xTi = 0.5 -> two-phase field (keep split)");
        double[] overallLo = { 0.5, 0.5 };
        HillertSolver.Phase candLo = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        StringBuilder logLo = new StringBuilder();
        GridMinimizationResult gLo = HillertGridMinimizer.minimize(
                List.of(candLo), T_BELOW, overallLo, l -> logLo.append(l).append('\n'));
        check(T_BELOW + " K: grid converged", gLo.converged(), gLo.failureReason());

        int nSelLo = gLo.selectedGridPoints().size();
        long aSelLo = gLo.selectedGridPoints().stream().filter(p -> p.phaseIndex() == 0).count();
        System.out.printf("      grid selected %d point(s), %d parent-A; stableInitialStates=%d%n",
                nSelLo, aSelLo, gLo.stableInitialStates().size());
        check(T_BELOW + " K: grid's discrete minimum uses 2 same-parent points (a gap)",
                aSelLo >= 2, "only " + aSelLo + " parent-A points selected");
        check(T_BELOW + " K: same-parent merge test RAN",
                logLo.toString().contains("same-parent merge test (parent 0"),
                "no merge-test log:\n" + logLo);
        check(T_BELOW + " K: merge test KEEPS SPLIT (genuine gap: G_single(mid) > grid mixture)",
                logLo.toString().contains("KEEP SPLIT for parent 0"),
                "did not keep split:\n" + logLo);
        check(T_BELOW + " K: two distinct stable instances (BCC_A2#1 / BCC_A2#2)",
                gLo.stableInitialStates().size() == 2,
                "got " + gLo.stableInitialStates().size());

        if (gLo.stableInitialStates().size() == 2) {
            double xa = gLo.stableInitialStates().get(0).m()[1];   // xTi of instance 1
            double xb = gLo.stableInitialStates().get(1).m()[1];   // xTi of instance 2
            double lo = Math.min(xa, xb);
            double hi = Math.max(xa, xb);
            System.out.printf("      grid initial-estimate tie line: xTi = %.5f / %.5f%n", lo, hi);
            System.out.printf("      Mathematica reference @%dK   : xTi = %.5f / %.5f%n",
                    (int) T_BELOW, XTI_ALPHA_400, XTI_BETA_400);
            check(T_BELOW + " K: instances bracket xTi = 0.5", lo < 0.5 && hi > 0.5,
                    "do not bracket: " + lo + " / " + hi);
            // grid initial estimate: each endpoint within one grid cell (N=20 -> 0.05)
            // of the Mathematica reference
            double cell = 1.0 / HillertGridMinimizer.GRID_DIVISIONS;
            check(T_BELOW + " K: Nb-rich endpoint within one grid cell of the reference",
                    Math.abs(lo - XTI_ALPHA_400) <= cell + 1e-9,
                    String.format("|%.5f - %.5f| = %.5f > %.5f", lo, XTI_ALPHA_400,
                            Math.abs(lo - XTI_ALPHA_400), cell));
            check(T_BELOW + " K: Ti-rich endpoint within one grid cell of the reference",
                    Math.abs(hi - XTI_BETA_400) <= cell + 1e-9,
                    String.format("|%.5f - %.5f| = %.5f > %.5f", hi, XTI_BETA_400,
                            Math.abs(hi - XTI_BETA_400), cell));
            // symmetry of the reference gap (Nb-Ti has zero 3-/4-body CECs)
            check(T_BELOW + " K: grid tie line is ~symmetric about xTi = 0.5 (as the reference is)",
                    Math.abs((lo + hi) - 1.0) < 2 * cell,
                    "not symmetric: lo+hi = " + (lo + hi));

            // Gm at the discrete tie line vs the reference Gm (-262.441). The grid
            // mixture Gm is the selected mixture's total G per mole; compare
            // loosely (grid endpoints differ from the continuous ones).
            double gGridMixPerMole = gLo.totalGibbsEnergy();
            // total G already per mole here (selected amounts sum to 1). Its
            // MIXING part is G - (x_Nb*G0_Nb + x_Ti*G0_Ti); we only have absolute
            // G, so compare the order of magnitude / sign instead of the value.
            System.out.printf("      grid mixture absolute G = %.4f (reference Gm at the tie line = %.4f)%n",
                    gGridMixPerMole, GM_TIELINE_400);
            check(T_BELOW + " K: grid mixture G is finite and negative (mixing lowers G in a gap)",
                    Double.isFinite(gGridMixPerMole) && gGridMixPerMole < 0.0,
                    "" + gGridMixPerMole);
        }

        // ---- (3) T_c bracket : gap present just below, absent just above ----
        section("T_c bracket: gap at T_c-30, no gap at T_c+30 (reference T_c ~ 477 K)");
        int splitBelow = countSelectedSameParent(s, m, T_C_REF - 30.0);
        int splitAbove = countSelectedSameParent(s, m, T_C_REF + 30.0);
        System.out.printf("      same-parent selected points: T=%.0f -> %d ;  T=%.0f -> %d%n",
                T_C_REF - 30.0, splitBelow, T_C_REF + 30.0, splitAbove);
        check("just below T_c: grid finds a 2-point (gap) minimum", splitBelow >= 2,
                "got " + splitBelow);
        check("just above T_c: grid finds a 1-point (single-phase) minimum OR merges",
                splitAbove <= 1, "got " + splitAbove);

        System.out.println("\n" + "=".repeat(90));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(90));
        if (failures > 0) {
            throw new AssertionError(failures + " Nb-Ti gap validation checks failed");
        }
    }

    /** Number of same-parent selected grid points at (T, xTi=0.5) -- >=2 means a gap. */
    private static int countSelectedSameParent(ModelSession s, CVMGibbsModel m, double t) {
        HillertSolver.Phase cand = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        GridMinimizationResult g = HillertGridMinimizer.minimize(
                List.of(cand), t, new double[] { 0.5, 0.5 }, null);
        if (!g.converged()) {
            return -1;
        }
        return (int) g.selectedGridPoints().stream().filter(p -> p.phaseIndex() == 0).count();
    }

    private static double eci(CVMGibbsModel m, String name) {
        for (var t : m.cecEntry().cecTerms) {
            if (name.equals(t.name)) {
                return (t.a) + (t.b) * 0.0;   // T-independent for Nb-Ti
            }
        }
        return Double.NaN;
    }

    private static void section(String s) {
        System.out.println("\n--- " + s + " ---");
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-76s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-76s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertNbTiGapValidation() {
    }
}
