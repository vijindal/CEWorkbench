package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Regression gate for the defensive-copy contract of
 * {@link HillertSolver.Result#mu()} and
 * {@link HillertSolver.PhaseResult#composition()} (STEP 6, PARTS 3-4).
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertResultImmutability
 * </pre>
 *
 * <p>A Java record auto-generates by-reference accessors and stores its
 * arguments by reference. Before this fix, a caller could mutate a
 * {@code Result}'s chemical-potential vector, or a {@code PhaseResult}'s
 * composition, either through the array handed back by the accessor or (for
 * {@code Result.mu}) through the array passed to the constructor. Both records
 * now clone on the way in and on the way out.</p>
 *
 * <p>The {@code Result} used here comes from a real (tiny) two-phase
 * {@link HillertSolver#solve} call, so the test also confirms the copying does
 * not disturb a normal solve.</p>
 */
public final class HillertResultImmutability {

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  HillertSolver.Result / PhaseResult defensive-copy contract");
        System.out.println("=".repeat(78));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel model = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", session.cecEntry, null);

        List<HillertSolver.Phase> phases = new ArrayList<>();
        phases.add(new HillertSolver.Phase("alpha", session, model, 0.5,
                model.randomStateFull(new double[] { 0.35, 0.65 })));
        phases.add(new HillertSolver.Phase("beta", session, model, 0.5,
                model.randomStateFull(new double[] { 0.65, 0.35 })));
        double[] target = { 0.5 * 0.35 + 0.5 * 0.65, 0.5 * 0.65 + 0.5 * 0.35 };

        HillertSolver.Result result =
                HillertSolver.solve(phases, target, 1000.0, 50, 20, 1.0e-6, null);

        checkResultMu(result);
        checkPhaseResultComposition(result);
        checkConstructorInputIsolation(model);

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " immutability checks failed");
        }
    }

    // ---- PART 3: Result.mu() -----------------------------------------

    private static void checkResultMu(HillertSolver.Result result) {
        System.out.println("\n--- Result.mu() ---");

        double[] pristine = result.mu().clone();
        check("mu() is non-null and non-empty", pristine != null && pristine.length > 0,
                "mu=" + Arrays.toString(pristine));

        double[] handed = result.mu();
        check("two mu() calls return distinct arrays", handed != result.mu(), "same reference");

        Arrays.fill(handed, 123456.0);
        check("mu() unchanged after scribbling on a prior return value",
                Arrays.equals(result.mu(), pristine),
                Arrays.toString(result.mu()));
    }

    // ---- PART 4: PhaseResult.composition() --------------------------

    private static void checkPhaseResultComposition(HillertSolver.Result result) {
        System.out.println("\n--- PhaseResult.composition() ---");

        check("result has phase entries", !result.phases().isEmpty(), "empty");

        for (HillertSolver.PhaseResult pr : result.phases()) {
            double[] pristine = pr.composition().clone();
            double[] handed = pr.composition();
            check("[" + pr.label() + "] two composition() calls return distinct arrays",
                    handed != pr.composition(), "same reference");

            Arrays.fill(handed, -1.0);
            check("[" + pr.label() + "] composition() unchanged after scribbling on a prior return",
                    Arrays.equals(pr.composition(), pristine),
                    Arrays.toString(pr.composition()));
        }
    }

    // ---- Constructor-input isolation (direct record construction) --

    private static void checkConstructorInputIsolation(CVMGibbsModel model) {
        System.out.println("\n--- Direct construction: input arrays isolated ---");

        double[] muIn = { 1.0, 2.0 };
        double[] compIn = { 0.4, 0.6 };
        double[] muSnap = muIn.clone();
        double[] compSnap = compIn.clone();

        CVMGibbsModel.State st = model.at(1000.0, new double[] { 0.4, 0.6 },
                model.randomStateU(new double[] { 0.4, 0.6 }));
        HillertSolver.PhaseResult pr =
                new HillertSolver.PhaseResult("x", 1.0, compIn, st.g(), st, false);
        HillertSolver.ConvergenceReport rep = new HillertSolver.ConvergenceReport(
                HillertSolver.ConvergenceReason.MAX_ITERATIONS,
                Double.POSITIVE_INFINITY, Double.NaN, Double.NaN, Double.NaN,
                false, false, 0, List.of(),
                new HillertSolver.MassBalanceReport(null, null, Double.NaN, Double.NaN, Double.NaN));
        HillertSolver.Result res =
                new HillertSolver.Result(List.of(pr), muIn, rep);

        Arrays.fill(muIn, -9.0);
        Arrays.fill(compIn, -9.0);

        check("Result.mu() unaffected by mutating the constructor input",
                Arrays.equals(res.mu(), muSnap), Arrays.toString(res.mu()));
        check("PhaseResult.composition() unaffected by mutating the constructor input",
                Arrays.equals(res.phases().get(0).composition(), compSnap),
                Arrays.toString(res.phases().get(0).composition()));

        check("phaseSetEvents is unmodifiable",
                isUnmodifiable(rep.phaseSetEvents()), "modifiable list");
    }

    private static boolean isUnmodifiable(List<?> list) {
        try {
            list.add(null);
            return false;   // add succeeded -> modifiable
        } catch (UnsupportedOperationException e) {
            return true;
        }
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-70s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-70s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertResultImmutability() {
    }
}
