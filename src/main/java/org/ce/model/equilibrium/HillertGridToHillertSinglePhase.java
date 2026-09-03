package org.ce.model.equilibrium;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertSolver.GridSeededResult;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.List;

/**
 * V2 STEP 8 gate (PART 12) -- the grid &rarr; Algorithm-A pipeline end to end
 * on a real single-phase V1 CVM equilibrium.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertGridToHillertSinglePhase
 * </pre>
 *
 * <p>{@code Nb-Ti / BCC_A2 / T} at {@code Ti = 0.5}, 1000 K, has a single
 * stable BCC_A2 phase. This test runs
 * {@link HillertSolver#solveFromGrid} (grid minimizer picks the initial
 * constitution, then the STEP-6 Algorithm-A loop) and cross-checks the final
 * state against a <b>direct</b> {@link HillertSolver#solve} started from an
 * equivalent good seed. The two must agree within the solver's own
 * tolerances.</p>
 *
 * <p>Also asserted (STEP 8 PART 6/15/16/17):</p>
 * <ul>
 *   <li>no {@code PHASE_ADDED} / {@code PHASE_REMOVED} events at initialisation
 *       (the seed is not a phase-set change);</li>
 *   <li>Algorithm A computes its own {@code mu} -- the grid's discrete
 *       {@code mu} is returned only as {@code gridMu} diagnostics, and the two
 *       are "reasonably close" for a well-resolved grid (no hard tolerance);</li>
 *   <li>{@code G_final <= G_grid_mixture} within tolerance (the discrete grid
 *       minimum is an upper bound);</li>
 *   <li>the caller's candidate {@link HillertSolver.Phase} objects are not
 *       mutated by {@code solveFromGrid}.</li>
 * </ul>
 */
public final class HillertGridToHillertSinglePhase {

    private static int failures = 0;
    private static final double T = 1000.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(88));
        System.out.println("  V2 STEP 8 PART 12 -- grid -> Algorithm A, single-phase real CVM");
        System.out.println("=".repeat(88));

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);

        double[] overall = { 0.5, 0.5 };   // Ti = 0.5

        // ---------- grid-seeded path ----------
        double[] candidateUFullBefore;
        double candidateAmountBefore;
        boolean candidateActiveBefore;
        HillertSolver.Phase candidate = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        candidateUFullBefore = candidate.uFull.clone();
        candidateAmountBefore = candidate.amount;
        candidateActiveBefore = candidate.active;

        GridSeededResult gsr = HillertSolver.solveFromGrid(
                List.of(candidate), overall, T, 120, 20, 1.0e-9, null);

        check("grid stage succeeded", !gsr.gridFailed(), "grid failed: "
                + (gsr.gridResult() == null ? "null" : gsr.gridResult().failureReason()));
        check("Algorithm A converged from the grid seed", gsr.overallConverged(),
                gsr.result().convergenceReport().reason().toString());

        HillertSolver.Result gr = gsr.result();

        // caller's candidate Phase not mutated
        check("caller candidate Phase.uFull not mutated",
                java.util.Arrays.equals(candidateUFullBefore, candidate.uFull), "uFull changed");
        check("caller candidate Phase.amount not mutated",
                candidate.amount == candidateAmountBefore, "amount changed");
        check("caller candidate Phase.active not mutated",
                candidate.active == candidateActiveBefore, "active changed");

        // no phase-set events at initialisation (single phase, no add/remove)
        check("no PHASE_ADDED / PHASE_REMOVED events (init is not a phase-set change)",
                gr.convergenceReport().phaseSetEvents().isEmpty(),
                "events: " + gr.convergenceReport().phaseSetEvents());

        // exactly one active phase at equilibrium
        long activeCount = gr.phases().stream().filter(p -> p.amount() > 1e-9).count();
        check("exactly one phase carries material at equilibrium", activeCount == 1,
                "got " + activeCount);

        HillertSolver.PhaseResult gPhase = gr.phases().stream()
                .filter(p -> p.amount() > 1e-9).findFirst().orElseThrow();

        // ---------- direct solve() from an equivalent good seed ----------
        HillertSolver.Phase seed = new HillertSolver.Phase(
                "BCC_A2", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        HillertSolver.Result dr = HillertSolver.solve(
                List.of(seed), overall, T, 120, 20, 1.0e-9, null);
        check("direct solve() also converged", dr.overallConverged(),
                dr.convergenceReport().reason().toString());
        HillertSolver.PhaseResult dPhase = dr.phases().stream()
                .filter(p -> p.amount() > 1e-9).findFirst().orElseThrow();

        // ---------- the two paths must agree ----------
        double gG = gPhase.g();
        double dG = dPhase.g();
        System.out.printf("      grid-seeded  : x=%s  N=%.8f  G=%.8f  mu=%s%n",
                java.util.Arrays.toString(round(gPhase.composition())), gPhase.amount(), gG,
                java.util.Arrays.toString(round(gr.mu())));
        System.out.printf("      direct solve : x=%s  N=%.8f  G=%.8f  mu=%s%n",
                java.util.Arrays.toString(round(dPhase.composition())), dPhase.amount(), dG,
                java.util.Arrays.toString(round(dr.mu())));

        check("final phase composition agrees (< 1e-6)",
                maxAbsDiff(gPhase.composition(), dPhase.composition()) < 1e-6,
                "diff " + maxAbsDiff(gPhase.composition(), dPhase.composition()));
        check("final phase amount agrees (< 1e-6)",
                Math.abs(gPhase.amount() - dPhase.amount()) < 1e-6,
                "diff " + Math.abs(gPhase.amount() - dPhase.amount()));
        check("final G agrees (< 1e-4 J/mol)", Math.abs(gG - dG) < 1e-4,
                "diff " + Math.abs(gG - dG));
        check("final mu agrees (< 1e-3 J/mol)",
                maxAbsDiff(gr.mu(), dr.mu()) < 1e-3, "diff " + maxAbsDiff(gr.mu(), dr.mu()));

        // mass balance satisfied on the grid-seeded result
        double[] rep = new double[2];
        for (HillertSolver.PhaseResult p : gr.phases()) {
            double[] x = p.composition();
            for (int i = 0; i < 2; i++) rep[i] += p.amount() * x[i];
        }
        check("grid-seeded result satisfies overall mass balance (< 1e-7)",
                Math.abs(rep[0] - 0.5) < 1e-7 && Math.abs(rep[1] - 0.5) < 1e-7,
                java.util.Arrays.toString(rep));

        // ---------- diagnostics: grid mu vs Algorithm-A mu (PART 15) ----------
        // For a SINGLE stable phase mu is underdetermined -- one Gibbs-Duhem
        // equation (sum_i mu_i x_i = G), K unknowns. The grid's dG/dx-tangent
        // fallback and Algorithm A's equilibrium-matrix mu are BOTH valid points
        // on that line and need not be equal entry-by-entry. The spec (PART 15)
        // asks only for a diagnostic, "no tolerance unless scientifically
        // justified" -- so the check is the physically meaningful invariant,
        // sum_i mu_i x_i == G, for each mu, plus the raw gap reported.
        double[] gridMu = gsr.gridMu();
        System.out.printf("      gridMu (diagnostic) = %s   Algorithm-A mu = %s%n",
                java.util.Arrays.toString(round(gridMu)), java.util.Arrays.toString(round(gr.mu())));
        check("gridMu is available as a diagnostic (not injected into solve)",
                gridMu != null && gridMu.length == 2, "gridMu null");
        double[] xFinal = gPhase.composition();
        double gridMuDotX = gridMu[0] * xFinal[0] + gridMu[1] * xFinal[1];
        double algMuDotX = gr.mu()[0] * xFinal[0] + gr.mu()[1] * xFinal[1];
        check("grid mu lies on the Gibbs-Duhem line (sum mu_i x_i == G, < 1e-3)",
                Math.abs(gridMuDotX - gG) < 1e-3, "sum=" + gridMuDotX + " G=" + gG);
        check("Algorithm-A mu lies on the same Gibbs-Duhem line (sum mu_i x_i == G, < 1e-3)",
                Math.abs(algMuDotX - gG) < 1e-3, "sum=" + algMuDotX + " G=" + gG);
        double muGap = maxAbsDiff(gridMu, gr.mu());
        System.out.printf("      |gridMu - AlgorithmA mu|_inf = %.4f J/mol  "
                + "(underdetermined for a single phase; both satisfy sum mu.x = G)%n", muGap);

        // ---------- grid mixture G is an upper bound (PART 16) ----------
        double gGridMix = gsr.gridMixtureGibbs();
        System.out.printf("      G_grid_mixture = %.8f   G_final = %.8f%n", gGridMix, gG);
        check("G_final <= G_grid_mixture + 1e-4 (discrete grid minimum is an upper bound)",
                gG <= gGridMix + 1e-4, "G_final " + gG + " > G_grid " + gGridMix);

        // ---------- performance (PART 19) ----------
        System.out.printf("      timing: grid = %d ms, Algorithm A = %d ms%n",
                gsr.gridMillis(), gsr.algorithmAMillis());

        System.out.println("\n" + "=".repeat(88));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(88));
        if (failures > 0) {
            throw new AssertionError(failures + " single-phase integration checks failed");
        }
    }

    private static double maxAbsDiff(double[] a, double[] b) {
        double m = 0.0;
        for (int i = 0; i < a.length; i++) {
            m = Math.max(m, Math.abs(a[i] - b[i]));
        }
        return m;
    }

    private static double[] round(double[] v) {
        if (v == null) {
            return null;
        }
        double[] r = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            r[i] = Math.round(v[i] * 1e6) / 1e6;
        }
        return r;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-72s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-72s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertGridToHillertSinglePhase() {
    }
}
