package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;

/**
 * Computes the converged bcc-bcc two-phase equilibrium for the binary Nb-Zr
 * system at T = 1073 K, seeded at (xNb) = 0.98 and 0.01, per user
 * instruction.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.NbZrMiscibilityGapTieLine
 * </pre>
 */
public final class NbZrMiscibilityGapTieLine {

    private static final double T = 1073.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(76));
        System.out.println("  Nb-Zr binary bcc-bcc two-phase equilibrium at T = 1073 K");
        System.out.println("=".repeat(76));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);

        CVMGibbsModel model = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", session.cecEntry, null);
        System.out.printf("%n  ncf = %d, numComponents = %d%n", model.ncf(), model.numComponents());

        List<HillertSolver.Phase> phases = new ArrayList<>();
        phases.add(phase("alpha", session, model, 0.5, new double[] { 0.98, 0.02 }));
        phases.add(phase("beta",  session, model, 0.5, new double[] { 0.01, 0.99 }));

        System.out.printf("%n  T = %.1f K, %d phases%n", T, phases.size());
        for (HillertSolver.Phase p : phases) {
            System.out.printf("    %-6s start x(Nb,Zr)=%s amount=%.3f%n",
                    p.label, fmt(p.composition()), p.amount);
        }

        double[] target = new double[model.numComponents()];
        for (HillertSolver.Phase p : phases) {
            double[] x = p.composition();
            for (int i = 0; i < target.length; i++) target[i] += p.amount * x[i];
        }
        HillertSolver.Result eq = HillertSolver.solve(phases, target, T, 200, 30, 1.0e-8, null);

        System.out.printf("%n  overallConverged = %s   outerIterations = %d   residual = %.4e%n",
                eq.overallConverged(), eq.outerIterations(), eq.finalResidualNorm());
        System.out.println("  mu(Nb,Zr) = " + fmt(eq.mu()));

        boolean bothStable = true;
        for (HillertSolver.PhaseResult p : eq.phases()) {
            System.out.printf("%n  --- %s ---%n", p.label());
            System.out.printf("    converged   = %s%n", p.phaseConverged());
            System.out.printf("    amount      = %.6f%n", p.amount());
            System.out.printf("    x(Nb,Zr)    = %s%n", fmt(p.composition()));
            System.out.printf("    G (abs)     = %.6f%n", p.g());
            System.out.printf("    Gm          = %.6f%n", p.state().gm());
            System.out.printf("    Hm          = %.6f%n", p.state().hm());
            System.out.printf("    Sm          = %.6f%n", p.state().sm());
            if (p.amount() <= 0) bothStable = false;
        }

        double sep = compositionDistance(eq.phases().get(0).composition(), eq.phases().get(1).composition());
        System.out.printf("%n  composition separation = %.6f   bothPhasesStable = %s%n", sep, bothStable);

        System.out.println("\n" + "=".repeat(76));
        if (eq.overallConverged() && bothStable && sep > 1.0e-3) {
            System.out.println("RESULT: CONVERGED, GENUINE TWO-PHASE TIE-LINE -- ready to use.");
        } else if (!eq.overallConverged()) {
            System.out.println("RESULT: DID NOT CONVERGE -- do not use these numbers.");
        } else {
            System.out.println("RESULT: Converged but NOT a genuine two-phase split (one phase unstable, "
                    + "or compositions collapsed together) -- confirms single-phase stability, consistent "
                    + "with the convexity check finding no gap on this binary edge at 1073 K.");
        }
        System.out.println("=".repeat(76));
    }

    private static HillertSolver.Phase phase(String label, ModelSession session,
            CVMGibbsModel model, double amount, double[] x) {
        double[] uFull = model.randomStateFull(x);
        return new HillertSolver.Phase(label, session, model, amount, uFull);
    }

    private static double compositionDistance(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            double d = a[i] - b[i];
            s += d * d;
        }
        return Math.sqrt(s);
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.6f", v[i]));
            if (i < v.length - 1) sb.append(", ");
        }
        return sb.append(']').toString();
    }

    private NbZrMiscibilityGapTieLine() {
    }
}
