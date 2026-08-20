package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.equilibrium.LatticeStability;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Smoke test for the G = G0m + Gm split and its derivatives, cross-checked
 * against a Mathematica reference gradient at the same state.
 *
 * <p><b>No minimization is performed.</b> Everything is evaluated at one
 * fixed, caller-supplied CF set via {@link CVMGibbsModel#evaluate}, which
 * only calls {@code setT}/{@code setX}/{@code setU} and then evaluates — the
 * Newton-Raphson loop in {@code minimize()} is never entered.</p>
 *
 * <h2>CF ordering</h2>
 *
 * <p>The Mathematica reference and this codebase order the two pair blocks
 * <em>oppositely</em>:</p>
 * <ul>
 *   <li>Mathematica: {@code v2AB1 v2AC1 v2BC1} (I-pair / 1NN) then
 *       {@code v2AB2 v2AC2 v2BC2} (II-pair / 2NN)</li>
 *   <li>Java (see {@code hamiltonian.json}): {@code v22AB v22AC v22BC}
 *       (II-pair, multiplicity 3) then {@code v21AB v21AC v21BC}
 *       (I-pair, multiplicity 4)</li>
 * </ul>
 *
 * <p>So {@code v2XY1 -> v21XY} and {@code v2XY2 -> v22XY}, and the blocks
 * swap position. This test therefore keys both the input CFs and the
 * reference gradient <em>by name</em> and reorders into Java's basis order
 * via {@link CvCfBasis#getNonPointCfNames}, rather than assuming positional
 * agreement — feeding a Mathematica-ordered array in positionally evaluates
 * a different physical state and silently produces a mismatched gradient.</p>
 */
public class SmokeTestG0mDerivatives {

    /** Relative tolerance for comparing against the 6-significant-digit reference. */
    static final double REL_TOL = 1e-5;
    static boolean allPass = true;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        String elements = "Mo-Nb-Ta";
        String structure = "BCC_A2";
        String model = "T";
        double T = 1273.0;
        double[] x = {0.05, 0.05, 0.90};

        // --- "Updated CFs" at this state, keyed by Mathematica's names -----
        //
        // Full double precision, taken from the Java trace of the same state
        // (the two agreed to all printed digits at this point). Do NOT
        // substitute the 6-significant-digit values from the Mathematica
        // printout: several cluster variables here are ~1e-6, and the entropy
        // gradient goes through log(cv) and 1/cv, so 6-digit inputs swing
        // ||dGm/du|| from 3.1e5 to 4.5e7. The comparison below is only
        // meaningful at full input precision.
        Map<String, Double> mmaCfs = new LinkedHashMap<>();
        mmaCfs.put("v4AB",   4.573970112248618E-6);
        mmaCfs.put("v4AC",   0.0018683795937251065);
        mmaCfs.put("v4BC",   0.0020350050982310517);
        mmaCfs.put("v4ABC1", 9.317425280014712E-5);
        mmaCfs.put("v4ABC2", 1.0517157239269508E-4);
        mmaCfs.put("v4ABC3", 0.00204714842775464);
        mmaCfs.put("v3AB",   1.50101984131667E-5);
        mmaCfs.put("v3AC",   0.03947912283817138);
        mmaCfs.put("v3BC",   0.03841912939954517);
        mmaCfs.put("v3ABC1", 0.0022423259201821617);
        mmaCfs.put("v3ABC2", 0.002187982543764565);
        mmaCfs.put("v3ABC3", 0.002195784220487935);
        mmaCfs.put("v2AB2",  0.0024544555163575417);
        mmaCfs.put("v2AC2",  0.04536121025512943);
        mmaCfs.put("v2BC2",  0.045054643001294054);
        mmaCfs.put("v2AB1",  0.0024057983970115168);
        mmaCfs.put("v2AC1",  0.04556667212648746);
        mmaCfs.put("v2BC1",  0.045106142173209865);

        // --- Mathematica GcuN at those CFs, same ordering ------------------
        Map<String, Double> mmaGrad = new LinkedHashMap<>();
        mmaGrad.put("v4AB",   -173571.0);
        mmaGrad.put("v4AC",   -191653.0);
        mmaGrad.put("v4BC",   -8.64047);
        mmaGrad.put("v4ABC1", -363200.0);
        mmaGrad.put("v4ABC2", -908.299);
        mmaGrad.put("v4ABC3", -1185.23);
        mmaGrad.put("v3AB",   -185503.0);
        mmaGrad.put("v3AC",   -200494.0);
        mmaGrad.put("v3BC",   -502.871);
        mmaGrad.put("v3ABC1", -385914.0);
        mmaGrad.put("v3ABC2",  8760.75);
        mmaGrad.put("v3ABC3", -6156.8);
        mmaGrad.put("v2AB1",   358379.0);
        mmaGrad.put("v2AC1",   366305.0);
        mmaGrad.put("v2BC1",  -777.666);
        mmaGrad.put("v2AB2",   5329.17);
        mmaGrad.put("v2AC2",  -8504.83);
        mmaGrad.put("v2BC2",  -396.594);

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());
        ModelSession session = builder.build(
                new SystemId(elements, structure, model), EngineConfig.CVM, null);

        CVMGibbsModel gm = new CVMGibbsModel();
        gm.initialize(elements, structure, model, session.cecEntry, null);

        int ncf = gm.getNcf();
        String[] els = elements.split("-");
        int K = els.length;
        List<String> javaNames = CvCfBasis.getNonPointCfNames(structure, model, K);

        System.out.println("=".repeat(78));
        System.out.printf("Smoke test (no minimize): %s %s %s at T=%.1f K, x=%s%n",
                elements, structure, model, T, Arrays.toString(x));
        System.out.println("  ncf = " + ncf + ", K = " + K);
        System.out.println("=".repeat(78));

        // --- Reorder Mathematica CFs into Java's basis order ---------------
        System.out.println("\n-- CF name mapping (Mathematica -> Java) -------------------------");
        System.out.printf("  %-4s %-10s %-10s %18s%n", "idx", "Java CF", "Mma CF", "value");
        System.out.println("  " + "-".repeat(48));

        double[] u = new double[ncf];
        for (int i = 0; i < ncf; i++) {
            String javaName = javaNames.get(i);
            String mmaName = toMathematicaName(javaName);
            Double v = mmaCfs.get(mmaName);
            if (v == null) {
                throw new IllegalStateException(
                        "No Mathematica CF for Java name " + javaName + " (mapped to " + mmaName + ")");
            }
            u[i] = v;
            System.out.printf("  %-4d %-10s %-10s %18.10g%n", i, javaName, mmaName, v);
        }

        // --- Evaluate at that state (pure evaluation, no solver) -----------
        org.ce.model.cvm.CvmState mixing = gm.getEvaluator().stateAt(T, x, u);
        double[] gmu = mixing.gmu();

        double g0m = LatticeStability.g0m(List.of(els), structure, x, T);

        // --- Scalars -------------------------------------------------------
        System.out.println("\n-- Scalars --------------------------------------------------------");
        System.out.printf("  G0m            = %22.10f%n", g0m);
        System.out.printf("  Gm             = %22.10f%n", mixing.gm());
        System.out.printf("  G  = G0m + Gm  = %22.10f%n", mixing.g());
        System.out.printf("  Hm             = %22.10f%n", mixing.hm());
        System.out.printf("  H  = H0m + Hm  = %22.10f%n", mixing.h());
        System.out.printf("  Sm  = S        = %22.10f   (S0m = 0)%n", mixing.sm());

        // --- Gradient comparison, component by component -------------------
        System.out.println("\n-- dGm/du vs Mathematica GcuN (matched by CF name) ----------------");
        System.out.printf("  %-4s %-10s %18s %18s %12s%n",
                "idx", "CF", "Java dGm/du", "Mma GcuN", "rel diff");
        System.out.println("  " + "-".repeat(68));

        int mismatches = 0;
        for (int i = 0; i < ncf; i++) {
            String javaName = javaNames.get(i);
            double ref = mmaGrad.get(toMathematicaName(javaName));
            double got = gmu[i];
            double denom = Math.max(Math.abs(ref), 1.0);
            double rel = Math.abs(got - ref) / denom;
            boolean ok = rel < REL_TOL;
            if (!ok) mismatches++;
            System.out.printf("  %-4d %-10s %18.6f %18.6f %12.2e%s%n",
                    i, javaName, got, ref, rel, ok ? "" : "   <-- MISMATCH");
        }
        System.out.println("  " + "-".repeat(68));
        System.out.printf("  ||dGm/du||     = %22.10f%n", norm(gmu));

        // --- G0m contribution ----------------------------------------------
        System.out.println("\n-- G0m derivatives -------------------------------------------------");
        System.out.println("  dG0m/du         = 0 (no u-dependence)  => dG/du == dGm/du");
        System.out.println("  d2G0m/du2       = 0                    => d2G/du2 == d2Gm/du2");
        System.out.println("  d2G0m/d(uFull)2 = 0                    => GuuFull == GmuuFull");
        System.out.println("  dG0m/dx_i (composition block of the widened gradient):");
        double dot = 0.0;
        for (int i = 0; i < K; i++) {
            double g0i = LatticeStability.g0(els[i], structure, T);
            dot += x[i] * g0i;
            System.out.printf("    G0(%-3s)      = %22.10f%n", els[i], g0i);
        }

        // --- Checks --------------------------------------------------------
        System.out.println("\n-- Checks ---------------------------------------------------------");

        check("All 18 gradient components match Mathematica (rel < " + REL_TOL + ")",
                mismatches == 0,
                mismatches == 0 ? "all components agree"
                                : mismatches + " component(s) differ — see MISMATCH rows above");

        check("Sum x_i*G0_i reproduces G0m",
                Math.abs(dot - g0m) < 1e-9,
                String.format("dot = %.10f, G0m = %.10f, diff = %.3e", dot, g0m, dot - g0m));

        check("G0m is nonzero (LatticeStability covers Mo/Nb/Ta in BCC_A2)",
                Math.abs(g0m) > 1.0,
                String.format("G0m = %.10f", g0m));

        check("G - Gm == G0m exactly",
                Math.abs((mixing.g() - mixing.gm()) - g0m) < 1e-9,
                String.format("G - Gm = %.10f", mixing.g() - mixing.gm()));

        check("H - Hm == G0m exactly (H0m == G0m)",
                Math.abs((mixing.h() - mixing.hm()) - g0m) < 1e-9,
                String.format("H - Hm = %.10f", mixing.h() - mixing.hm()));

        System.out.println("\n" + "=".repeat(78));
        System.out.println(allPass ? "RESULT: PASS" : "RESULT: FAIL");
        System.out.println("=".repeat(78));
    }

    /**
     * Java CVCF name -> Mathematica name. Only the pair block is spelled
     * differently; the correspondence is by suffix digit, and both
     * conventions list the blocks in the same order:
     * {@code v22XY <-> v2XY2} and {@code v21XY <-> v2XY1}.
     *
     * <p>Note this is <em>not</em> a block swap. An earlier version of this
     * test mapped {@code v21XY -> v2XY2} on the theory that the two
     * conventions ordered I-pair and II-pair blocks oppositely; that made
     * all 18 gradient components mismatch and inflated the gradient norm
     * from 3.1e5 to 4.5e7, confirming the straight suffix mapping below is
     * the correct one.</p>
     */
    static String toMathematicaName(String javaName) {
        if (javaName.startsWith("v22") && javaName.length() == 5) {
            return "v2" + javaName.substring(3) + "2";
        }
        if (javaName.startsWith("v21") && javaName.length() == 5) {
            return "v2" + javaName.substring(3) + "1";
        }
        return javaName;
    }

    static double norm(double[] v) {
        double s = 0;
        for (double a : v) s += a * a;
        return Math.sqrt(s);
    }

    static void check(String label, boolean ok, String detail) {
        System.out.printf("  [%s] %s%n", ok ? "PASS" : "FAIL", label);
        if (detail != null && !detail.isEmpty()) System.out.println("         " + detail);
        if (!ok) allPass = false;
    }
}
