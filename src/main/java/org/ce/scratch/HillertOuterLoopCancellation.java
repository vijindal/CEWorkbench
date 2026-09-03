package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;

/**
 * Regression gate for {@link HillertSolver}'s outer-loop cooperative
 * cancellation (STEP 4).
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertOuterLoopCancellation
 * </pre>
 *
 * <p>Before STEP 4, {@code HillertSolver.solve} polled nothing itself --
 * cancellation was observed only if the run happened to enter the
 * phase-addition candidate scan, whose inner {@link CvmNewtonSolver} calls
 * check the interrupt flag. A run with no inactive candidates (or one that
 * never settles enough to scan) could execute every {@code maxOuterIterations}
 * with zero responsiveness. STEP 4 adds one check at the top of the outer loop,
 * before any phase state is touched.</p>
 *
 * <p>Deterministic -- no sleeps, no timing windows:</p>
 * <ul>
 *   <li><b>A</b> -- interrupt flag set before {@code solve()} is entered:
 *       {@link CancellationException} on the very first outer-iteration
 *       boundary, no iteration executed, no phase mutation.</li>
 *   <li><b>C</b> -- interrupt set from a progress-sink callback the instant the
 *       first outer iteration reports: the <em>next</em> iteration's
 *       top-of-loop check throws. This exercises the outer-loop poll
 *       specifically, not the transitive {@code CvmNewtonSolver} path. Phases
 *       reflect only the one accepted iteration; the interrupt flag is left
 *       set.</li>
 * </ul>
 *
 * <p>The existing candidate-scan cancellation coverage (R2 in
 * {@code HillertAbsoluteMuValidation}) is not duplicated here; this suite adds
 * only the outer-loop-specific cases.</p>
 */
public final class HillertOuterLoopCancellation {

    private static int failures = 0;
    private static final double T = 1000.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  HillertSolver outer-loop cooperative cancellation");
        System.out.println("=".repeat(78));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel model = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", session.cecEntry, null);

        testA(session, model);
        testC(session, model);

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " cancellation checks failed");
        }
    }

    // ---- Test A: interrupt before solve() -------------------------------

    private static void testA(ModelSession s, CVMGibbsModel m) throws Exception {
        System.out.println("\n--- Test A: interrupt set BEFORE solve() ---");

        final Throwable[] caught = new Throwable[1];
        final boolean[] normalReturn = { false };
        final boolean[] flagStillSetInsideCatch = { false };
        final double[][] uFullBefore = new double[2][];
        final double[] amountBefore = new double[2];

        Thread worker = new Thread(() -> {
            List<HillertSolver.Phase> ph = twoPhase(s, m);
            uFullBefore[0] = ph.get(0).uFull.clone();
            uFullBefore[1] = ph.get(1).uFull.clone();
            amountBefore[0] = ph.get(0).amount;
            amountBefore[1] = ph.get(1).amount;

            Thread.currentThread().interrupt();               // set BEFORE entering solve
            try {
                HillertSolver.solve(ph, targetOf(ph), T, 80, 20, 1.0e-8, null);
                normalReturn[0] = true;
            } catch (Throwable t) {
                caught[0] = t;
                flagStillSetInsideCatch[0] = Thread.currentThread().isInterrupted();
            }
            // Compare after: solve must not have mutated the phases at all.
            check("A: phase 0 uFull unchanged", Arrays.equals(ph.get(0).uFull, uFullBefore[0]),
                    "mutated");
            check("A: phase 1 uFull unchanged", Arrays.equals(ph.get(1).uFull, uFullBefore[1]),
                    "mutated");
            check("A: phase 0 amount unchanged", ph.get(0).amount == amountBefore[0],
                    amountBefore[0] + " -> " + ph.get(0).amount);
            check("A: phase 1 amount unchanged", ph.get(1).amount == amountBefore[1],
                    amountBefore[1] + " -> " + ph.get(1).amount);
        }, "hillert-outer-cancel-A");

        worker.start();
        worker.join(30_000);

        System.out.printf("    A: alive=%s normalReturn=%s caught=%s%n",
                worker.isAlive(), normalReturn[0],
                caught[0] == null ? "(none)" : caught[0].getClass().getName());
        check("A: solve did not return normally", !normalReturn[0], "returned normally");
        check("A: a CancellationException propagated",
                caught[0] instanceof CancellationException,
                caught[0] == null ? "nothing thrown" : caught[0].getClass().getName());
        check("A: interrupt flag left set (not cleared)", flagStillSetInsideCatch[0], "flag cleared");
        check("A: worker finished (no deadlock)", !worker.isAlive(), "still alive");
    }

    // ---- Test C: interrupt from the progress sink, mid-run --------------

    private static void testC(ModelSession s, CVMGibbsModel m) throws Exception {
        System.out.println("\n--- Test C: interrupt from progress sink after outer iter 1 ---");

        final Throwable[] caught = new Throwable[1];
        final boolean[] normalReturn = { false };
        final AtomicInteger iterLinesSeen = new AtomicInteger();
        final int[] interruptedAtLine = { -1 };

        Thread worker = new Thread(() -> {
            List<HillertSolver.Phase> ph = twoPhase(s, m);

            Consumer<String> sink = line -> {
                // The per-iteration summary line is emitted once near the end of
                // each outer iteration: "Hillert outer iter N: lambda=...".
                if (line.startsWith("Hillert outer iter ") && line.contains("accepted=")) {
                    int n = iterLinesSeen.incrementAndGet();
                    if (n == 1) {
                        // We are still on the solver thread, at the end of
                        // iteration 1. Setting the flag now must make
                        // iteration 2's top-of-loop check throw.
                        Thread.currentThread().interrupt();
                        interruptedAtLine[0] = n;
                    }
                }
            };

            try {
                HillertSolver.solve(ph, targetOf(ph), T, 80, 20, 1.0e-8, sink);
                normalReturn[0] = true;
            } catch (Throwable t) {
                caught[0] = t;
            }

            check("C: interrupt flag left set after solve threw",
                    Thread.currentThread().isInterrupted(), "flag cleared");
        }, "hillert-outer-cancel-C");

        worker.start();
        worker.join(60_000);

        System.out.printf("    C: alive=%s normalReturn=%s caught=%s iterLinesSeen=%d%n",
                worker.isAlive(), normalReturn[0],
                caught[0] == null ? "(none)" : caught[0].getClass().getName(),
                iterLinesSeen.get());

        check("C: the run actually started iterating (>=1 iteration ran)",
                iterLinesSeen.get() >= 1, "no iteration ran");
        check("C: solve did not return normally after mid-run interrupt",
                !normalReturn[0], "returned normally");
        check("C: a CancellationException propagated",
                caught[0] instanceof CancellationException,
                caught[0] == null ? "nothing thrown" : caught[0].getClass().getName());
        // If the outer-loop check were absent, the run would either finish
        // normally (normalReturn) or only stop if a candidate scan happened to
        // run -- neither is what we asserted above.
        check("C: cancellation was observed within a bounded number of iterations",
                iterLinesSeen.get() <= 5,
                "kept iterating for " + iterLinesSeen.get() + " lines after interrupt");
        check("C: worker finished (no deadlock)", !worker.isAlive(), "still alive");
    }

    // ---- helpers ------------------------------------------------------

    private static List<HillertSolver.Phase> twoPhase(ModelSession s, CVMGibbsModel m) {
        List<HillertSolver.Phase> ph = new ArrayList<>();
        ph.add(new HillertSolver.Phase("alpha", s, m, 0.5,
                m.randomStateFull(new double[] { 0.35, 0.65 })));
        ph.add(new HillertSolver.Phase("beta", s, m, 0.5,
                m.randomStateFull(new double[] { 0.65, 0.35 })));
        return ph;
    }

    private static double[] targetOf(List<HillertSolver.Phase> ph) {
        int k = ph.get(0).numComponents;
        double[] t = new double[k];
        for (HillertSolver.Phase p : ph) {
            double[] x = p.composition();
            for (int i = 0; i < k; i++) t[i] += p.amount * x[i];
        }
        return t;
    }

    private static synchronized void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-58s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-58s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertOuterLoopCancellation() {
    }
}
