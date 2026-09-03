package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Regression gate for {@link HillertSolver}'s V1 phase-model scope guard.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertUnsupportedPhaseGuard
 * </pre>
 *
 * <p>V1 mass balance is {@code target_i = sum_p N_p x^p_i}, valid only when a
 * phase's mole fraction equals its moles of component {@code i} per formula unit
 * ({@code M_A = x_A}). That holds for single-site, one-atom-per-site
 * <b>disordered</b> phases, whose CVCF point block is exactly the K mole
 * fractions: {@code geometry.tcf - geometry.ncf == K}. An ordered /
 * multi-sublattice phase (BCC_B2 binary registers an extra long-range-order
 * point CF {@code eta}, so {@code tcf - ncf == 3} at {@code K == 2}) breaks that
 * assumption -- the solver must reject the problem before iterating rather than
 * silently treat {@code x} as {@code M_A}.</p>
 *
 * <p>Checks:</p>
 * <ol>
 *   <li>a lone ordered (BCC_B2) phase -> {@code UNSUPPORTED_PHASE_MODEL}, zero
 *       iterations, no phase-state mutation;</li>
 *   <li>a supported (BCC_A2) active phase together with an <em>inactive</em>
 *       ordered candidate -> still {@code UNSUPPORTED_PHASE_MODEL} (Part 4: the
 *       guard inspects candidates too, since phase addition could activate one);</li>
 *   <li>negative control: two supported BCC_A2 phases are NOT rejected by the
 *       guard (the run proceeds and its reason is anything but
 *       {@code UNSUPPORTED_PHASE_MODEL}).</li>
 * </ol>
 */
public final class HillertUnsupportedPhaseGuard {

    private static int failures = 0;
    private static final double T = 1000.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  HillertSolver V1 phase-model scope guard");
        System.out.println("=".repeat(78));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());

        // --- Supported reference: BCC_A2 binary ---
        ModelSession sessionA2 = builder.build(
                new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel a2 = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", sessionA2.cecEntry, null);

        // --- Unsupported: BCC_B2 binary (ordered, carries the eta LRO point CF) ---
        CvmGeometry b2geo = CvmGeometry.build("Nb-Ti", "BCC_B2", "T", null);
        System.out.printf("%n  BCC_B2 geometry: K=%d ncf=%d tcf=%d  (tcf-ncf=%d, expected > K for an ordered phase)%n",
                b2geo.numComponents, b2geo.ncf, b2geo.tcf, b2geo.tcf - b2geo.ncf);
        check("BCC_B2 really is out of V1 scope (tcf-ncf != K)",
                (b2geo.tcf - b2geo.ncf) != b2geo.numComponents,
                "tcf-ncf=" + (b2geo.tcf - b2geo.ncf) + " K=" + b2geo.numComponents);

        // A B2 model needs a CECEntry; the A2 one has the wrong basis names, so
        // build a permissive empty entry -- the guard fires before any ECI is
        // evaluated, so its contents never matter here.
        CECEntry b2entry = emptyEntry("Nb-Ti", "BCC_B2_T");
        CVMGibbsModel b2 = new CVMGibbsModel(b2geo, b2entry);

        checkLoneOrderedPhaseRejected(b2);
        checkInactiveOrderedCandidateRejected(a2, sessionA2, b2, b2geo);
        checkSupportedNotRejected(a2, sessionA2);

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " scope-guard checks failed");
        }
    }

    // -- Check 1: a lone ordered phase is rejected, untouched --------------

    private static void checkLoneOrderedPhaseRejected(CVMGibbsModel b2) {
        System.out.println("\n--- Check 1: lone BCC_B2 phase ---");

        double[] uFull0 = orderedSeed(b2);
        HillertSolver.Phase p = new HillertSolver.Phase("b2", null, b2, 1.0, uFull0);
        double[] uFullBefore = p.uFull.clone();
        double amountBefore = p.amount;
        boolean activeBefore = p.active;

        HillertSolver.Result r = HillertSolver.solve(
                List.of(p), new double[] { 0.5, 0.5 }, T, 50, 20, 1.0e-6, null);

        check("reason == UNSUPPORTED_PHASE_MODEL",
                r.convergenceReport().reason() == HillertSolver.ConvergenceReason.UNSUPPORTED_PHASE_MODEL,
                String.valueOf(r.convergenceReport().reason()));
        check("not overallConverged", !r.overallConverged(), "converged");
        check("zero outer iterations", r.convergenceReport().iterationsRun() == 0,
                "iterationsRun=" + r.convergenceReport().iterationsRun());
        check("no phase-set events", r.convergenceReport().phaseSetEvents().isEmpty(),
                "events=" + r.convergenceReport().phaseSetEvents().size());
        check("phase uFull not mutated", Arrays.equals(p.uFull, uFullBefore), "uFull changed");
        check("phase amount not mutated", p.amount == amountBefore,
                amountBefore + " -> " + p.amount);
        check("phase active flag not mutated", p.active == activeBefore,
                activeBefore + " -> " + p.active);
        check("Result still structurally valid (1 phase entry, finite mu vector)",
                r.phases().size() == 1 && isFinite(r.mu()), "malformed Result");
    }

    // -- Check 2: an inactive ordered candidate alongside a supported phase --

    private static void checkInactiveOrderedCandidateRejected(
            CVMGibbsModel a2, ModelSession sessionA2, CVMGibbsModel b2, CvmGeometry b2geo) {
        System.out.println("\n--- Check 2: supported active phase + inactive ordered candidate ---");

        HillertSolver.Phase active = new HillertSolver.Phase(
                "a2-active", sessionA2, a2, 1.0, a2.randomStateFull(new double[] { 0.5, 0.5 }));
        // Inactive candidate: amount 0 so Phase's constructor leaves active=false.
        HillertSolver.Phase candidate = new HillertSolver.Phase(
                "b2-candidate", null, b2, 0.0, orderedSeed(b2));
        check("candidate really starts inactive", !candidate.active, "active");

        double[] activeUFullBefore = active.uFull.clone();
        double activeAmountBefore = active.amount;

        HillertSolver.Result r = HillertSolver.solve(
                List.of(active, candidate), new double[] { 0.5, 0.5 }, T, 50, 20, 1.0e-6, null);

        check("reason == UNSUPPORTED_PHASE_MODEL (candidate counts)",
                r.convergenceReport().reason() == HillertSolver.ConvergenceReason.UNSUPPORTED_PHASE_MODEL,
                String.valueOf(r.convergenceReport().reason()));
        check("zero outer iterations", r.convergenceReport().iterationsRun() == 0,
                "iterationsRun=" + r.convergenceReport().iterationsRun());
        check("supported phase not mutated (uFull)",
                Arrays.equals(active.uFull, activeUFullBefore), "uFull changed");
        check("supported phase not mutated (amount)", active.amount == activeAmountBefore,
                activeAmountBefore + " -> " + active.amount);
        check("both phase entries echoed back", r.phases().size() == 2, "size=" + r.phases().size());
    }

    // -- Check 3: negative control -- supported phases are not rejected ----

    private static void checkSupportedNotRejected(CVMGibbsModel a2, ModelSession sessionA2) {
        System.out.println("\n--- Check 3: negative control -- two BCC_A2 phases ---");

        List<HillertSolver.Phase> phases = new ArrayList<>();
        phases.add(new HillertSolver.Phase(
                "alpha", sessionA2, a2, 0.5, a2.randomStateFull(new double[] { 0.35, 0.65 })));
        phases.add(new HillertSolver.Phase(
                "beta", sessionA2, a2, 0.5, a2.randomStateFull(new double[] { 0.65, 0.35 })));

        double[] target = new double[2];
        for (HillertSolver.Phase p : phases) {
            double[] x = p.composition();
            for (int i = 0; i < 2; i++) target[i] += p.amount * x[i];
        }

        HillertSolver.Result r = HillertSolver.solve(phases, target, T, 50, 20, 1.0e-6, null);

        check("supported problem NOT rejected as UNSUPPORTED_PHASE_MODEL",
                r.convergenceReport().reason() != HillertSolver.ConvergenceReason.UNSUPPORTED_PHASE_MODEL,
                "reason=" + r.convergenceReport().reason());
        check("supported problem actually iterated",
                r.convergenceReport().iterationsRun() > 0,
                "iterationsRun=" + r.convergenceReport().iterationsRun());
    }

    // -- helpers ---------------------------------------------------------

    /**
     * A physical-ish joint vector for an ordered phase, length {@code tcf}: the
     * random state at x=[0.5,0.5] built by the basis (which does fill the wider
     * ordered point block), so the pre-iteration {@code model.atFull} call
     * inside the rejection path evaluates cleanly.
     */
    private static double[] orderedSeed(CVMGibbsModel b2) {
        try {
            return b2.randomStateFull(new double[] { 0.5, 0.5 });
        } catch (RuntimeException e) {
            // Fallback: a plain zero point block of the right width. The guard
            // fires before this vector is used for anything numeric.
            CvmGeometry g = b2.geometry();
            double[] v = new double[g.tcf];
            for (int i = g.ncf; i < g.tcf; i++) v[i] = 0.5;
            return v;
        }
    }

    private static CECEntry emptyEntry(String elements, String structurePhase) {
        CECEntry e = new CECEntry();
        e.elements = elements;
        e.structurePhase = structurePhase;
        e.model = "T";
        e.cecTerms = new CECEntry.CECTerm[0];
        return e;
    }

    private static boolean isFinite(double[] v) {
        for (double x : v) if (!Double.isFinite(x)) return false;
        return true;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-62s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-62s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertUnsupportedPhaseGuard() {
    }
}
