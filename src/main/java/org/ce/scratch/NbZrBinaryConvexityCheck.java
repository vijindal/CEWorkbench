package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

/**
 * Convexity check on the single-phase Nb-Zr binary G_m(x) curve at T = 1073 K,
 * to confirm/locate a bcc-bcc miscibility gap near this edge before seeding
 * the (ternary or binary) Hillert two-phase solver.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.NbZrBinaryConvexityCheck
 * </pre>
 */
public final class NbZrBinaryConvexityCheck {

    private static final double T = 1073.0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);

        CVMGibbsModel model = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", session.cecEntry, null);

        System.out.println("=".repeat(90));
        System.out.printf("  Nb-Zr binary G_m(x) scan and convexity check at T = %.1f K%n", T);
        System.out.println("=".repeat(90));

        int n = 100;
        double[] xNbArr = new double[n + 1];
        double[] gVals = new double[n + 1];
        boolean[] conv = new boolean[n + 1];
        for (int i = 0; i <= n; i++) {
            double xNb = i / (double) n;
            if (xNb == 0.0 || xNb == 1.0) {
                xNbArr[i] = xNb;
                gVals[i] = 0.0;
                conv[i] = true;
                continue;
            }
            double[] x = { xNb, 1.0 - xNb };
            CvmNewtonSolver.Result r = new CvmNewtonSolver(model).solve(T, x, 1.0e-9, null, null);
            xNbArr[i] = xNb;
            gVals[i] = r.state().gm();
            conv[i] = r.converged();
            System.out.printf("  xNb=%.3f  converged=%-5s  Gm=%.4f%n", xNb, r.converged(), gVals[i]);
        }

        System.out.println();
        System.out.println("  Common-tangent scan: for every pair (i,j), does the chord dip below any k between them?");
        double worstDip = 0;
        int wi = -1, wj = -1, wk = -1;
        for (int i = 0; i <= n; i++) {
            for (int j = i + 2; j <= n; j++) {
                if (!conv[i] || !conv[j]) continue;
                for (int k = i + 1; k < j; k++) {
                    if (!conv[k]) continue;
                    double tk = (k - i) / (double) (j - i);
                    double chord = gVals[i] * (1 - tk) + gVals[j] * tk;
                    double dip = chord - gVals[k];
                    if (dip > worstDip) { worstDip = dip; wi = i; wj = j; wk = k; }
                }
            }
        }

        System.out.printf("%n  Largest (chord - Gm) = %.4f J/mol at xNb(i)=%.3f, xNb(k)=%.3f, xNb(j)=%.3f%n",
                worstDip, xNbArr[wi], xNbArr[wk], xNbArr[wj]);

        if (worstDip > 1.0) {
            System.out.println();
            System.out.println("  ==> NON-CONVEX: a genuine two-phase (miscibility gap) region exists here.");
            System.out.printf("  ==> Candidate tie-line endpoints: xNb=%.3f  and  xNb=%.3f%n",
                    xNbArr[wi], xNbArr[wj]);
        } else {
            System.out.println();
            System.out.println("  ==> CONVEX (within tolerance): no evidence of a binary Nb-Zr gap at this T.");
        }
        System.out.println("=".repeat(90));
    }

    private NbZrBinaryConvexityCheck() {
    }
}
