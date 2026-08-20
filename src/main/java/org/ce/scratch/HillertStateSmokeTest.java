package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.equilibrium.PhaseEquilibriumResult;
import org.ce.model.equilibrium.PhaseState;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Exercises the Hillert multi-phase path end to end and checks that each
 * phase's retained {@link CVMGibbsModel.State} behaves the same way the
 * single-phase Newton solver's does -- any property, on demand, with no
 * re-solve.
 *
 * <p>This path had no end-to-end coverage: its pieces were each verified, but
 * nothing ran the full outer loop. Convergence on a real two-phase system is
 * not asserted here (that is a physics question this test cannot settle);
 * what is asserted is that the machinery runs, that every reported quantity
 * comes from one consistent state, and that the state stays usable after the
 * solve returns.</p>
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertStateSmokeTest
 * </pre>
 */
public final class HillertStateSmokeTest {

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(76));
        System.out.println("  Hillert multi-phase path: state reuse after convergence");
        System.out.println("=".repeat(76));

        double T = 1000.0;
        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());

        ModelSession sessionA2 = builder.build(
                new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);

        CVMGibbsModel a2 = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", sessionA2.cecEntry, null);

        // Two phases of the same structure at different starting compositions:
        // enough to drive the outer loop without depending on a second
        // structure's Hamiltonian being physically meaningful.
        List<PhaseState> phases = new ArrayList<>();
        phases.add(phase("alpha", sessionA2, a2, 0.5, new double[] { 0.35, 0.65 }));
        phases.add(phase("beta", sessionA2, a2, 0.5, new double[] { 0.65, 0.35 }));

        System.out.printf("%n  T = %.1f K, %d phases%n", T, phases.size());
        for (PhaseState p : phases) {
            System.out.printf("    %-6s start x=%s amount=%.3f%n",
                    p.label, fmt(p.composition()), p.amount);
        }

        PhaseEquilibriumResult eq = HillertSolver.solve(phases, T, 50, 20, 1.0e-6, null);

        System.out.printf("%n  overallConverged = %s   outerIterations = %d   residual = %.4e%n",
                eq.overallConverged(), eq.outerIterations(), eq.finalResidualNorm());
        System.out.println("  mu = " + fmt(eq.mu()));

        for (PhaseEquilibriumResult.PhaseResultEntry p : eq.phases()) {
            System.out.printf("%n  --- %s ---%n", p.label());
            System.out.printf("    amount      = %.6f%n", p.amount());
            System.out.printf("    composition = %s%n", fmt(p.composition()));
            System.out.printf("    G (abs)     = %.6f%n", p.g());

            // The point of retaining the state: everything else, on demand.
            CVMGibbsModel.State st = p.state();
            System.out.printf("    Gm          = %.6f%n", st.gm());
            System.out.printf("    G0m         = %.6f%n", st.g0m());
            System.out.printf("    Hm          = %.6f%n", st.hm());
            System.out.printf("    Sm          = %.6f%n", st.sm());
            System.out.printf("    SRO 1NN     = %s%n",
                    st.unlikePairSro(CVMGibbsModel.Shell.FIRST));

            check("g() == state().g()", p.g() == st.g());
            check("G == G0m + Gm", Math.abs(st.g() - (st.g0m() + st.gm())) < 1e-9);
            check("state names its model", p.model() == a2);
            check("composition matches state", Arrays.equals(p.composition(), st.composition()));

            // Widened gradient must be available -- it is what the per-phase
            // step solved against, and is meaningless on a single-phase state
            // that never had composition as an unknown.
            double[] guFull = st.gmuFull();
            check("gmuFull width = ncf+K", guFull.length == a2.ncf() + a2.numComponents());
        }

        System.out.println("\n" + "=".repeat(76));
        System.out.printf("RESULT: %s   (%d assertion failures)%n",
                failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(76));
        if (failures > 0) {
            throw new AssertionError(failures + " assertions failed");
        }
    }

    /** Builds a phase whose joint vector starts at the random state for {@code x}. */
    private static PhaseState phase(String label, ModelSession session,
            CVMGibbsModel model, double amount, double[] x) {
        double[] uFull = model.randomStateFull(x);
        return new PhaseState(label, session, model, amount, uFull);
    }

    private static void check(String what, boolean ok) {
        if (!ok) {
            failures++;
            System.out.println("    [!] FAIL  " + what);
        }
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.5f", v[i]));
            if (i < v.length - 1) sb.append(", ");
        }
        return sb.append(']').toString();
    }

    private HillertStateSmokeTest() {
    }
}
