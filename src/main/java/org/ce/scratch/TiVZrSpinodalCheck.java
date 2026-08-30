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
 * HillertSolver keeps collapsing both phases to the SAME Ti-rich composition
 * (~0.40,0.30,0.29) regardless of starting separation, OR flinging one phase
 * to an unphysical V-rich corner with negative amount as it goes unstable.
 * Neither outcome demonstrates a genuine two-phase split. This does a
 * completely independent check: scan single-phase Gm along a very fine line
 * through the region and look at the SECOND DERIVATIVE (curvature) directly,
 * to see if there is a genuine spinodal (concave-down) region at all, as
 * opposed to isolated non-convergence noise creating spurious "dips" in the
 * coarser scans.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.TiVZrSpinodalCheck
 * </pre>
 */
public final class TiVZrSpinodalCheck {

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

        // Fine line through the region where HillertSeed4 landed:
        // (0.658,0.027,0.314) -> (0.42,0.30,0.28)-ish -> further to V-rich
        double[] a = { 0.66, 0.02, 0.32 };
        double[] b = { 0.10, 0.88, 0.02 };

        System.out.println("=".repeat(100));
        System.out.printf("  Fine curvature scan: A=%s -> B=%s, T=%.1fK%n", fmt(a), fmt(b), T);
        System.out.println("=".repeat(100));

        int n = 100;
        double[] g = new double[n + 1];
        boolean[] conv = new boolean[n + 1];
        double[][] xs = new double[n + 1][];
        for (int i = 0; i <= n; i++) {
            double t = i / (double) n;
            double[] x = interp(a, b, t);
            xs[i] = x;
            CvmNewtonSolver.Result r = new CvmNewtonSolver(model).solve(T, x, 1.0e-9, null, null);
            g[i] = r.state().gm();
            conv[i] = r.converged();
        }

        // Print with local second-difference curvature (only for consecutive
        // converged triples), and flag chord violations among ALL converged
        // pairs/triples regardless of adjacency.
        System.out.printf("  %-6s %-28s %-10s %14s %14s%n", "t", "x", "converged", "Gm", "curvature");
        for (int i = 0; i <= n; i++) {
            String curv = "";
            if (i > 0 && i < n && conv[i - 1] && conv[i] && conv[i + 1]) {
                double dt = 1.0 / n;
                double d2 = (g[i - 1] - 2 * g[i] + g[i + 1]) / (dt * dt);
                curv = String.format("%14.1f", d2);
            }
            System.out.printf("  %.3f  %-28s %-10s %14.4f %s%n",
                    i / (double) n, fmt(xs[i]), conv[i], g[i], curv);
        }

        double worstDip = 0;
        int wi = -1, wj = -1, wk = -1;
        for (int i = 0; i <= n; i++) {
            if (!conv[i]) continue;
            for (int j = i + 2; j <= n; j++) {
                if (!conv[j]) continue;
                for (int k = i + 1; k < j; k++) {
                    if (!conv[k]) continue;
                    double tk = (k - i) / (double) (j - i);
                    double chord = g[i] * (1 - tk) + g[j] * tk;
                    double dip = chord - g[k];
                    if (dip > worstDip) { worstDip = dip; wi = i; wj = j; wk = k; }
                }
            }
        }
        System.out.println();
        System.out.printf("  worst dip (converged only) = %.4f J/mol", worstDip);
        if (wi >= 0) {
            System.out.printf(" at i=%d(t=%.3f,x=%s) k=%d(t=%.3f,x=%s) j=%d(t=%.3f,x=%s)%n",
                    wi, wi/(double)n, fmt(xs[wi]), wk, wk/(double)n, fmt(xs[wk]), wj, wj/(double)n, fmt(xs[wj]));
        }
        System.out.println("=".repeat(100));
    }

    private static double[] interp(double[] p, double[] q, double t) {
        double[] x = new double[p.length];
        for (int i = 0; i < p.length; i++) x[i] = p[i] * (1 - t) + q[i] * t;
        return x;
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.4f", v[i]));
            if (i < v.length - 1) sb.append(",");
        }
        return sb.append(']').toString();
    }

    private TiVZrSpinodalCheck() {
    }
}
