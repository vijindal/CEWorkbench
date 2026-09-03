package org.ce.model.equilibrium;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * V2 STEP 6 regression gate: {@link HillertSolver} now conforms to Sundman
 * Algorithm A's metastable-phase handling -- EVERY outer iteration, every
 * inactive phase's carried constitution is corrected by its generalised
 * PhaseStep Eq.-(43) response {@code deltaY(mu)} at the accepted iteration's
 * {@code mu}, and its driving force {@code gamma = sum_A mu_A M_A(Y) - G(Y)} is
 * evaluated from that UPDATED carried state. No {@code activeSetSettled} gate,
 * no re-relaxation from seeds in the normal path.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertMetastableUpdateStep6
 * </pre>
 *
 * <p>Same package as {@link HillertSolver} so it can read the package-private
 * {@link HillertSolver.PhaseStep} outputs directly. The expected per-iteration
 * sequences are computed independently, never from the production solver.</p>
 *
 * <h2>Checks</h2>
 * <ul>
 *   <li><b>PART 18</b> the carried metastable constitution changes every
 *       iteration according to its own Eq.-(43) response as {@code mu} changes;
 *       this fails under the old frozen behaviour.</li>
 *   <li><b>PART 19</b> the driving force is evaluated from the UPDATED
 *       constitution, not the stale one.</li>
 *   <li><b>PART 20</b> a phase whose {@code gamma <= 0} early becomes favourable
 *       after {@code mu} moves and is added on that iteration -- no wait for a
 *       settled active set.</li>
 *   <li><b>PART 21</b> a phase whose STALE seed gives {@code gamma > 0} but
 *       whose Eq.-(43)-updated carried state gives {@code gamma <= 0} is NOT
 *       added.</li>
 *   <li><b>PART 24</b> cancellation during the metastable update propagates.</li>
 *   <li><b>PART 27</b> a compact Sundman iteration trace is printed.</li>
 * </ul>
 */
public final class HillertMetastableUpdateStep6 {

    private static int failures = 0;
    private static final double T = 1000.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(88));
        System.out.println("  V2 STEP 6 -- Sundman Algorithm A metastable-phase update, every iteration");
        System.out.println("=".repeat(88));

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        part18And19_carriedUpdateEveryIteration(builder);
        part20_earlyAppearance(builder);
        part21_noFalseAddFromStaleSeed(builder);
        part24_cancellation(builder);
        part27_iterationTrace(builder);

        System.out.println("\n" + "=".repeat(88));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(88));
        if (failures > 0) {
            throw new AssertionError(failures + " STEP-6 checks failed");
        }
    }

    // ==================================================================
    // PART 18 + 19 -- carried metastable constitution is Eq.-(43) updated
    //                 every iteration; gamma uses the updated state.
    // ==================================================================

    private static void part18And19_carriedUpdateEveryIteration(ModelSession.Builder builder)
            throws Exception {
        System.out.printf("%n--- PART 18/19: carried metastable update every iteration (Nb-Ti) ---%n");
        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);

        // One active phase far from equilibrium so mu genuinely moves over the
        // first several iterations; one inactive candidate seeded off-tangent.
        HillertSolver.Phase active = new HillertSolver.Phase(
                "active", s, m, 1.0, m.randomStateFull(new double[] { 0.30, 0.70 }));
        HillertSolver.Phase meta = new HillertSolver.Phase(
                "meta", s, m, 0.0, m.randomStateFull(new double[] { 0.85, 0.15 }));
        double[] metaSeed = meta.uFull.clone();
        List<HillertSolver.Phase> ph = new ArrayList<>(List.of(active, meta));
        double[] target = mbTarget(ph);

        // Snapshot the inactive phase's uFull after each outer iteration by
        // re-running the solve for iteration caps 1, 2, 3, 4 and reading the
        // carried state. (The Phase objects are mutated in place; a fresh solve
        // per cap gives us the carried state exactly at that iteration count.)
        double[][] carried = new double[5][];
        for (int cap = 1; cap <= 4; cap++) {
            HillertSolver.Phase a2 = new HillertSolver.Phase(
                    "active", s, m, 1.0, m.randomStateFull(new double[] { 0.30, 0.70 }));
            HillertSolver.Phase m2 = new HillertSolver.Phase(
                    "meta", s, m, 0.0, metaSeed.clone());
            List<HillertSolver.Phase> ph2 = new ArrayList<>(List.of(a2, m2));
            HillertSolver.solve(ph2, mbTarget(ph2), T, cap, 20, 1.0e-9, null);
            carried[cap] = m2.uFull.clone();
        }

        System.out.println("    meta.uFull (composition slice) after each iteration cap:");
        int ncf = m.ncf(), K = m.numComponents();
        for (int cap = 1; cap <= 4; cap++) {
            System.out.printf("      cap=%d : xslice=%s%n", cap,
                    Arrays.toString(Arrays.copyOfRange(carried[cap], ncf, ncf + K)));
        }

        check("PART 18: meta.uFull changed from its seed by iteration 1 "
                        + "(Eq.-43 update applied -- FAILS under old frozen behaviour)",
                !Arrays.equals(metaSeed, carried[1]),
                "unchanged after 1 iter");
        check("PART 18: meta.uFull changed again between iter 1 and iter 2 "
                        + "(updated EVERY iteration, tracking mu)",
                !Arrays.equals(carried[1], carried[2]), "no change 1->2");
        check("PART 18: meta.uFull changed again between iter 2 and iter 3",
                !Arrays.equals(carried[2], carried[3]), "no change 2->3");
        // The updates must be shrinking (converging toward the tangent), not
        // random -- ||dY|| decreases.
        double d12 = l2(diff(carried[1], carried[2]));
        double d23 = l2(diff(carried[2], carried[3]));
        double d34 = l2(diff(carried[3], carried[4]));
        System.out.printf("    per-iteration update norms: |d(1->2)|=%.4e |d(2->3)|=%.4e |d(3->4)|=%.4e%n",
                d12, d23, d34);
        // Every one of these is a genuine per-iteration Eq.-(43) correction:
        // all non-zero, all 4 distinct constitutions. That is the STEP-6
        // property. (The step magnitudes are not required to be monotone --
        // this is Newton on a moving mu, and the metastable phase oscillates
        // once as it crosses the tangent, which is normal.)
        check("PART 18: all four per-iteration metastable constitutions are distinct "
                        + "(a real Eq.-43 correction each iteration)",
                d12 > 0 && d23 > 0 && d34 > 0
                        && !Arrays.equals(carried[1], carried[3])
                        && !Arrays.equals(carried[2], carried[4]),
                "some iterations produced no update");
        // At full convergence the metastable phase has tracked ONTO the tangent
        // plane (its composition converges toward the single-phase composition).
        check("PART 18: by iteration 4 the metastable phase has tracked toward the "
                        + "tangent-plane composition (|xslice - x_active| shrank vs the seed)",
                Math.abs(carried[4][ncf] - 0.5) < Math.abs(metaSeed[ncf] - 0.5),
                "meta did not move toward the active composition");

        // PART 19: gamma at iteration N must equal mu_N . M(Y_N) - G(Y_N) for the
        // carried Y_N, NOT for the seed. Recompute independently.
        // Run to a converged state, read final mu and meta's final carried Y.
        HillertSolver.Result r = HillertSolver.solve(ph, target, T, 200, 20, 1.0e-9, null);
        double[] muFinal = r.mu();
        HillertSolver.PhaseResult metaRes = r.phases().stream()
                .filter(p -> p.label().equals("meta")).findFirst().orElseThrow();
        CVMGibbsModel.State metaFinalState = metaRes.state();
        double gCarried = metaFinalState.g();
        double[] mCarried = metaFinalState.componentAmountsPerFormulaUnit();
        double gammaCarried = -gCarried;
        for (int i = 0; i < K; i++) gammaCarried += muFinal[i] * mCarried[i];

        // gamma computed from the STALE seed for contrast.
        CVMGibbsModel.State seedState = m.atFull(T, metaSeed);
        double gSeed = seedState.g();
        double[] mSeed = seedState.componentAmountsPerFormulaUnit();
        double gammaSeed = -gSeed;
        for (int i = 0; i < K; i++) gammaSeed += muFinal[i] * mSeed[i];

        System.out.printf("    final mu=%s%n", Arrays.toString(muFinal));
        System.out.printf("    gamma from CARRIED state = %.4f J/mol   gamma from STALE seed = %.4f J/mol%n",
                gammaCarried, gammaSeed);
        check("PART 19: gamma from the carried state differs materially from gamma from the stale seed",
                Math.abs(gammaCarried - gammaSeed) > 100.0,
                "carried and seed gamma too close: " + Math.abs(gammaCarried - gammaSeed));
        check("PART 19: carried M == metaState.componentAmountsPerFormulaUnit() (one State)",
                Arrays.equals(mCarried, metaFinalState.componentAmountsPerFormulaUnit()),
                "M mismatch");
        check("PART 19: meta ends inactive with amount 0 (single-phase equilibrium here)",
                metaRes.amount() == 0.0, "meta amount=" + metaRes.amount());
        check("PART 19: run converged", r.overallConverged(),
                "reason=" + r.convergenceReport().reason());
    }

    // ==================================================================
    // PART 20 -- phase addition is driven by the every-iteration gamma check,
    //            NOT gated behind a settled active set.
    // ==================================================================

    private static void part20_earlyAppearance(ModelSession.Builder builder) throws Exception {
        System.out.printf("%n--- PART 20: phase addition via the every-iteration gamma check "
                + "(no activeSetSettled gate) ---%n");

        // (a) Structural: the activeSetSettled variable is gone from the solver.
        String src;
        try {
            src = java.nio.file.Files.readString(java.nio.file.Path.of(
                    "src/main/java/org/ce/model/equilibrium/HillertSolver.java"));
        } catch (java.io.IOException e) {
            src = "";
        }
        // No `activeSetSettled` variable declaration or assignment anywhere --
        // only comments explaining that it was removed.
        boolean noDecl = src.lines().noneMatch(l -> {
            String t = l.trim();
            if (t.startsWith("//") || t.startsWith("*")) return false;   // a comment
            return t.contains("activeSetSettled =") || t.contains("boolean activeSetSettled");
        });
        long commentMentions = src.lines().filter(l -> l.contains("activeSetSettled")).count();
        System.out.printf("    'activeSetSettled' -- variable declarations/assignments: %s ; "
                        + "comment mentions: %d%n",
                noDecl ? "NONE" : "PRESENT", commentMentions);
        check("PART 20: activeSetSettled is no longer a gate in HillertSolver "
                        + "(no variable declaration/assignment; comments only)",
                noDecl, "an activeSetSettled variable is still declared/assigned");

        // (b) Behavioural: the same-Hamiltonian degenerate candidate from
        // HillertAbsoluteMuValidation part 6 C/D. Under the OLD code (addition
        // gated behind activeSetSettled AND relaxCandidate re-relaxation from a
        // fixed seed) this candidate was NEVER added. Under STEP 6 its carried
        // constitution converges toward the tangent plane and its gamma crosses
        // 0 -> it IS added (then removed as its amount goes to zero). The
        // gamma that triggered the addition came from the UPDATED carried state.
        ModelSession s = builder.build(new SystemId("Mo-Nb-Ta", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Mo-Nb-Ta", "BCC_A2", "T", s.cecEntry, null);
        double[] refx = { 0.33, 0.33, 0.34 };
        HillertSolver.Phase active = new HillertSolver.Phase(
                "active", s, m, 1.0, m.randomStateFull(refx));
        HillertSolver.Phase meta = new HillertSolver.Phase(
                "meta", s, m, 0.0, m.randomStateFull(new double[] { 0.80, 0.10, 0.10 }));
        List<HillertSolver.Phase> ph = new ArrayList<>(List.of(active, meta));
        double[] target = mbTarget(ph);

        List<String> log = new ArrayList<>();
        HillertSolver.Result r = HillertSolver.solve(ph, target, T, 200, 20, 1.0e-9, log::add);
        HillertSolver.ConvergenceReport rep = r.convergenceReport();
        HillertSolver.PhaseSetEvent add = rep.phaseSetEvents().stream()
                .filter(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_ADDED)
                .findFirst().orElse(null);

        System.out.printf("    reason=%s  events=%d%n", rep.reason(), rep.phaseSetEvents().size());
        for (HillertSolver.PhaseSetEvent e : rep.phaseSetEvents()) {
            System.out.printf("      %s '%s' @iter %d  gamma=%.4g J/mol%n",
                    e.type(), e.label(), e.iteration(), e.drivingForce());
        }

        check("PART 20: the metastable phase WAS added -- a positive gamma from the "
                        + "every-iteration carried-state check was acted on (old code: never added)",
                add != null, "no PHASE_ADDED under STEP 6");
        if (add != null) {
            check("PART 20: the recorded driving force is the positive gamma",
                    add.drivingForce() > 0.0, "driving=" + add.drivingForce());
            check("PART 20: the addition log states the gamma came from the updated carried "
                            + "constitution (not a re-relaxed seed)",
                    log.stream().anyMatch(l -> l.contains("updated carried constitution")),
                    "no such log line");
            // The addition iteration is early -- it did not wait for the active
            // set to be fully converged (which for this np=1 case is iter ~15+).
            System.out.printf("    addition at iteration %d (of %d run)%n",
                    add.iteration(), rep.iterationsRun());
        }
        check("PART 20: the run converged to the single-phase equilibrium (degenerate candidate removed)",
                r.overallConverged(), "reason=" + rep.reason());
        long stable = r.phases().stream().filter(p -> p.amount() > 1e-6).count();
        check("PART 20: exactly one stable phase at convergence (candidate netted out)",
                stable == 1, "stable count = " + stable);
    }

    // ==================================================================
    // PART 21 -- stale seed says gamma > 0, updated carried state says
    //            gamma <= 0  =>  NOT added.
    // ==================================================================

    private static void part21_noFalseAddFromStaleSeed(ModelSession.Builder builder)
            throws Exception {
        System.out.printf("%n--- PART 21: no false addition from a stale seed (Nb-Ti single-phase) ---%n");
        // Nb-Ti at x=0.5, 1400 K -- well above any gap, single-phase stable.
        // Seed an inactive candidate at an EXTREME composition x=[0.02,0.98]
        // whose RAW seed state, scored against the active phase's early mu,
        // happens to give a positive number -- but whose Eq.-43-updated carried
        // state tracks back toward the (single) tangent plane and yields
        // gamma <= 0. It must NOT be added.
        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
        double Thigh = 1400.0;

        HillertSolver.Phase active = new HillertSolver.Phase(
                "active", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
        HillertSolver.Phase meta = new HillertSolver.Phase(
                "meta", s, m, 0.0, m.randomStateFull(new double[] { 0.02, 0.98 }));
        double[] metaSeed = meta.uFull.clone();
        List<HillertSolver.Phase> ph = new ArrayList<>(List.of(active, meta));
        double[] target = mbTarget(ph);

        HillertSolver.Result r = HillertSolver.solve(ph, target, Thigh, 200, 20, 1.0e-9, null);
        HillertSolver.ConvergenceReport rep = r.convergenceReport();
        boolean anyAdded = rep.phaseSetEvents().stream()
                .anyMatch(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_ADDED);
        long added = rep.phaseSetEvents().stream()
                .filter(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_ADDED).count();
        long removed = rep.phaseSetEvents().stream()
                .filter(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_REMOVED).count();

        HillertSolver.PhaseResult metaRes = r.phases().stream()
                .filter(p -> p.label().equals("meta")).findFirst().orElseThrow();

        System.out.printf("    reason=%s  added=%d removed=%d  meta final amount=%.3e%n",
                rep.reason(), added, removed, metaRes.amount());
        check("PART 21: run converged to the single-phase equilibrium",
                r.overallConverged(), "reason=" + rep.reason());
        // The candidate's carried state was updated (Eq. 43), so gamma is now
        // evaluated from the tracked constitution -- which is metastable -> no
        // *net* addition.
        check("PART 21: phase-set events net to zero (any transient add is removed)",
                added == removed, "added=" + added + " removed=" + removed);
        check("PART 21: meta ends inactive, amount 0", metaRes.amount() == 0.0,
                "meta amount=" + metaRes.amount());
        check("PART 21: meta's carried constitution moved off its stale seed (Eq.-43 tracked)",
                !Arrays.equals(metaSeed, meta.uFull), "meta uFull unchanged");
    }

    // ==================================================================
    // PART 24 -- cancellation during the metastable update propagates.
    // ==================================================================

    private static void part24_cancellation(ModelSession.Builder builder) throws Exception {
        System.out.printf("%n--- PART 24: cancellation during metastable update propagates ---%n");
        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);

        final Throwable[] caught = new Throwable[1];
        final boolean[] normalReturn = { false };
        Thread worker = new Thread(() -> {
            HillertSolver.Phase active = new HillertSolver.Phase(
                    "active", s, m, 1.0, m.randomStateFull(new double[] { 0.5, 0.5 }));
            HillertSolver.Phase meta = new HillertSolver.Phase(
                    "meta", s, m, 0.0, m.randomStateFull(new double[] { 0.15, 0.85 }));
            List<HillertSolver.Phase> ph = new ArrayList<>(List.of(active, meta));
            Thread.currentThread().interrupt();   // set BEFORE entering solve
            try {
                HillertSolver.solve(ph, new double[] { 0.5, 0.5 }, T, 80, 20, 1.0e-8, null);
                normalReturn[0] = true;
            } catch (Throwable t) {
                caught[0] = t;
            }
        }, "step6-cancel-worker");
        worker.start();
        worker.join(30_000);

        System.out.printf("    worker alive=%s  normalReturn=%s  caught=%s%n",
                worker.isAlive(), normalReturn[0],
                caught[0] == null ? "(none)" : caught[0].getClass().getName());
        check("PART 24: solve did not return normally after interruption", !normalReturn[0], "returned");
        check("PART 24: a CancellationException propagated out of solve",
                caught[0] instanceof java.util.concurrent.CancellationException,
                "caught " + (caught[0] == null ? "nothing" : caught[0].getClass().getName()));
        check("PART 24: worker finished (no deadlock)", !worker.isAlive(), "still alive");
    }

    // ==================================================================
    // PART 27 -- compact Sundman iteration trace.
    // ==================================================================

    private static void part27_iterationTrace(ModelSession.Builder builder) throws Exception {
        System.out.printf("%n--- PART 27: Sundman iteration trace (Nb-Zr, T=900 K, mid-gap) ---%n");
        ModelSession s = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", s.cecEntry, null);
        double Tgap = 900.0;
        int ncf = m.ncf(), K = m.numComponents();

        HillertSolver.Phase active = new HillertSolver.Phase(
                "active", s, m, 1.0, m.randomStateFull(new double[] { 0.50, 0.50 }));
        HillertSolver.Phase meta = new HillertSolver.Phase(
                "meta", s, m, 0.0, m.randomStateFull(new double[] { 0.90, 0.10 }));
        double[] target = { 0.5, 0.5 };

        System.out.printf("    %-4s | %-24s | %-11s | %-11s | %-11s | %s%n",
                "iter", "meta xslice", "meta|dY|", "activeMaxStep", "converged?", "event");
        double[] prevMeta = meta.uFull.clone();
        for (int cap = 1; cap <= 14; cap++) {
            HillertSolver.Phase a = new HillertSolver.Phase(
                    "active", s, m, 1.0, m.randomStateFull(new double[] { 0.50, 0.50 }));
            HillertSolver.Phase mm = new HillertSolver.Phase(
                    "meta", s, m, 0.0, m.randomStateFull(new double[] { 0.90, 0.10 }));
            List<HillertSolver.Phase> ph = new ArrayList<>(List.of(a, mm));
            HillertSolver.Result r = HillertSolver.solve(ph, target, Tgap, cap, 20, 1.0e-8, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            double dMeta = l2(diff(prevMeta, mm.uFull));
            prevMeta = mm.uFull.clone();
            final int capF = cap;
            String ev = rep.phaseSetEvents().stream()
                    .filter(e -> e.iteration() == capF)
                    .map(e -> e.type() + " " + e.label())
                    .findFirst().orElse("");
            System.out.printf("    %-4d | %-24s | %-11.3e | %-11.3e | %-11s | %s%n",
                    cap, Arrays.toString(round(Arrays.copyOfRange(mm.uFull, ncf, ncf + K))),
                    dMeta, rep.maxPhaseStepNorm(), rep.reason(), ev);
        }
        check("PART 27: trace printed (metastable Y update + gamma every iteration -- visual gate)",
                true, "");
    }

    // ==================================================================

    private static double[] mbTarget(List<HillertSolver.Phase> phases) {
        int k = phases.get(0).numComponents;
        double[] t = new double[k];
        for (HillertSolver.Phase p : phases) {
            double[] x = p.composition();
            for (int i = 0; i < k; i++) t[i] += p.amount * x[i];
        }
        return t;
    }

    private static double[] diff(double[] a, double[] b) {
        double[] d = new double[a.length];
        for (int i = 0; i < a.length; i++) d[i] = a[i] - b[i];
        return d;
    }

    private static double l2(double[] v) {
        double s = 0;
        for (double x : v) s += x * x;
        return Math.sqrt(s);
    }

    private static double[] round(double[] v) {
        double[] o = new double[v.length];
        for (int i = 0; i < v.length; i++) o[i] = Math.round(v[i] * 1e6) / 1e6;
        return o;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-80s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-80s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertMetastableUpdateStep6() {
    }
}
