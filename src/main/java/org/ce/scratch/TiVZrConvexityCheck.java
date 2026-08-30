package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.workflow.TernarySubsystemExtractor;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.List;

/**
 * Directly tests whether single-phase G_m(x) is convex along a straight
 * composition line through the region flagged as non-convergent by
 * {@link TiVZrGapLocator} -- the unambiguous, solver-independent signature of
 * a two-phase (miscibility-gap) region: if the chord between two endpoints
 * ever lies below the sampled G_m curve, the system phase-separates somewhere
 * on that line, regardless of whether {@link org.ce.model.equilibrium.HillertSolver}
 * can be seeded to find it directly.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.TiVZrConvexityCheck
 * </pre>
 */
public final class TiVZrConvexityCheck {

    private static final double T = 1073.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession quatSession = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Ti-V-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);

        List<String> parentElements = List.of("Nb", "Ti", "V", "Zr");
        CECEntry ternaryEntry = TernarySubsystemExtractor.extractTernary(
                quatSession.cecEntry, parentElements, "Ti", "V", "Zr");
        CVMGibbsModel model = CVMGibbsModel.of("Ti-V-Zr", "BCC_A2", "T", ternaryEntry, null);

        // Endpoints: a converged Ti-rich point and a converged V-rich point,
        // with the suspect non-convergent band lying between them.
        double[] xA = { 0.55, 0.25, 0.20 };
        double[] xB = { 0.10, 0.85, 0.05 };

        System.out.println("=".repeat(90));
        System.out.printf("  Convexity check along the line x(t) = (1-t)*A + t*B, t in [0,1], T=%.1f K%n", T);
        System.out.printf("  A = %s   B = %s%n", fmt(xA), fmt(xB));
        System.out.println("=".repeat(90));

        int n = 40;
        double[] gValues = new double[n + 1];
        boolean[] conv = new boolean[n + 1];
        for (int i = 0; i <= n; i++) {
            double t = i / (double) n;
            double[] x = interp(xA, xB, t);
            CvmNewtonSolver.Result r = new CvmNewtonSolver(model).solve(T, x, 1.0e-9, null, null);
            gValues[i] = r.state().gm();
            conv[i] = r.converged();
            System.out.printf("  t=%.3f  x=%s  converged=%-5s  Gm=%.4f%n",
                    t, fmt(x), r.converged(), gValues[i]);
        }

        System.out.println();
        System.out.println("  Convex-hull check: does any interior point lie above the chord");
        System.out.println("  between two other (possibly non-adjacent) sampled points?");

        double worstViolation = 0;
        int worstI = -1, worstJ = -1, worstK = -1;
        for (int i = 0; i <= n; i++) {
            for (int j = i + 2; j <= n; j++) {
                for (int k = i + 1; k < j; k++) {
                    double tk = (k - i) / (double) (j - i);
                    double chord = gValues[i] * (1 - tk) + gValues[j] * tk;
                    double violation = gValues[k] - chord; // positive => above chord => non-convex signature is chord BELOW curve, i.e. curve dips below chord means convex; we want gValues[k] < chord for a stable two-phase split (lower G by splitting)
                    double dip = chord - gValues[k]; // positive means splitting into i,j lowers G below single-phase k -- i.e. i,j more stable combined
                    if (dip > worstViolation) {
                        worstViolation = dip;
                        worstI = i; worstJ = j; worstK = k;
                    }
                }
            }
        }

        System.out.printf("%n  Largest (chord - Gm) found: %.4f J/mol at i=%d (t=%.3f), k=%d (t=%.3f), j=%d (t=%.3f)%n",
                worstViolation, worstI, worstI / (double) n, worstK, worstK / (double) n, worstJ, worstJ / (double) n);

        if (worstViolation > 1.0) {
            System.out.println();
            System.out.println("  ==> NON-CONVEX: a two-phase split between the endpoints at those t values");
            System.out.println("      lowers G below the single-phase curve -- a genuine miscibility gap");
            System.out.println("      exists somewhere in this composition range at T=1073 K.");
            System.out.printf("  ==> Candidate tie-line endpoints: x(t=%.3f)=%s   x(t=%.3f)=%s%n",
                    worstI / (double) n, fmt(interp(xA, xB, worstI / (double) n)),
                    worstJ / (double) n, fmt(interp(xA, xB, worstJ / (double) n)));
        } else {
            System.out.println();
            System.out.println("  ==> CONVEX (within tolerance): no evidence of a two-phase split along");
            System.out.println("      this line at T=1073 K. The non-convergence seen in TiVZrGapLocator");
            System.out.println("      is most likely a solver-robustness artifact, not a real gap, OR the");
            System.out.println("      gap (if any) lies off this particular line.");
        }
        System.out.println("=".repeat(90));
    }

    private static double[] interp(double[] a, double[] b, double t) {
        double[] x = new double[a.length];
        for (int i = 0; i < a.length; i++) x[i] = a[i] * (1 - t) + b[i] * t;
        return x;
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.4f", v[i]));
            if (i < v.length - 1) sb.append(", ");
        }
        return sb.append(']').toString();
    }

    private TiVZrConvexityCheck() {
    }
}
