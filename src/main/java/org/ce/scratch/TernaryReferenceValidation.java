package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Validates both CVM minimisers against a Mathematica {@code phaseq} reference
 * point for a ternary system.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.TernaryReferenceValidation
 * </pre>
 *
 * <p><b>Reference point.</b> Mo-Nb-Ta / BCC_A2 / T-model, 1000 K,
 * x = [0.33, 0.33, 0.34], supplied by the user from the Mathematica
 * implementation. Scalars are quoted to 6 significant figures there, so
 * agreement is asserted at {@code 1e-5} relative, not to machine precision.</p>
 *
 * <p><b>What makes this more than a self-check.</b> The same point is solved
 * twice, by two solvers that share only the evaluator:</p>
 * <ul>
 *   <li>{@link CvmNewtonSolver} holds composition fixed and solves an
 *       {@code ncf}-dimensional system against {@code gmu}/{@code gmuu}.</li>
 *   <li>{@link HillertSolver} (np=1) treats composition as an unknown, solving
 *       the widened {@code ncf+K} system against {@code gmuFull}/{@code gmuuFull}
 *       with a Lagrange constraint and an outer chemical-potential loop.</li>
 * </ul>
 * <p>Agreement between them, and with an external reference, is therefore
 * evidence about {@code CVMGibbsModel} itself rather than about one loop.</p>
 *
 * <p><b>The reference is stored as (ECI, CF) pairs</b>, in {@link #REFERENCE}.
 * The two codes differ in cluster-algebra conventions -- labels may be
 * exchanged, and block order differs -- so neither a position nor a label
 * identifies a cluster across them on its own. What must correspond is the
 * cluster: the ECI and the equilibrium CF have to be the same one. Checking the
 * ECI at the slot a name resolves to is what establishes that; a CF checked
 * alone could pass while attached to the wrong cluster, if a label were
 * exchanged consistently in both input and output. See {@link #REFERENCE} for
 * the full reasoning.</p>
 *
 * <p>Two conventions are reconciled to make the match, and they are different
 * in kind:</p>
 * <ol>
 *   <li><b>A spelling difference, not a reordering.</b> The reference writes
 *       pair CFs shell-last ({@code v2AB1} = pair AB, 1st shell); we write them
 *       shell-first ({@code v21AB}). {@link #canonical} decodes both. Every pair
 *       CF and ECI agrees once decoded -- the {@code v21}/{@code v22} blocks are
 *       <em>not</em> transposed, though the reference emits 1NN before 2NN while
 *       we emit 2NN first, which makes a positional read of an unlabelled vector
 *       look as though they are.</li>
 *   <li><b>One real labelling difference</b>, {@link #REF_TO_OURS_NAME}:
 *       {@code v3ABC1} and {@code v3ABC3} are exchanged. Deliberately unfixed.</li>
 * </ol>
 *
 * <p><b>Chemical potential is deliberately not compared entry-by-entry.</b>
 * For a single phase the only constraint on mu is one Gibbs-Duhem equation, so
 * with K unknowns there is a (K-1)-dimensional family of valid solutions; ours
 * and the reference's are different points in it. What is asserted is the
 * constraint itself, {@code sum(mu_i * x_i) == G}. Note also that
 * {@code CvmNewtonSolver} has no chemical potential at all -- composition is
 * fixed, so mu is not part of its problem.</p>
 *
 * <p><b>Not covered.</b> This validates one interior composition. The
 * near-edge band (e.g. x = [0.05, 0.05, 0.90] on this same system) is where
 * both solvers are known to fail, and is a separate open item.</p>
 */
public final class TernaryReferenceValidation {

    // ---- Reference point ------------------------------------------------
    private static final String ELEMENTS  = "Mo-Nb-Ta";
    private static final String STRUCTURE = "BCC_A2";
    private static final String MODEL     = "T";
    private static final double T         = 1000.0;
    private static final double[] REF_X   = { 0.33, 0.33, 0.34 };

    private static final double REF_G  = -69246.3;
    private static final double REF_GM = -20633.7;
    private static final double REF_HM = -11972.7;
    private static final double REF_SM = 8.66101;

    /** Reference mu -- one valid solution of the underdetermined np=1 system. */
    private static final double[] REF_MU = { -73827.8, -59613.7, -74148.9 };

    /**
     * The reference point, stored as {@code (name, ECI, converged CF)} triples
     * in the reference's own emission order.
     *
     * <p><b>Why ECI and CF are stored together.</b> The two codes differ in
     * cluster-algebra conventions: labels may be exchanged, and block order
     * differs (the reference emits the 1NN pair block before 2NN, we emit 2NN
     * first). Each code is internally consistent -- within it, ECI slot
     * {@code i} and CF slot {@code i} describe the same cluster -- but neither
     * a position nor a label is meaningful across the two on its own.</p>
     *
     * <p>The ECI is what identifies the cluster. A CF value checked alone could
     * pass while attached to the wrong cluster, if a label were exchanged
     * consistently in both the ECI input and the CF output. Requiring the ECI to
     * match at the same slot removes that possibility: the pair is the unit of
     * meaning, so the gate validates pairs.</p>
     *
     * <p>ECIs are from the Mo-Nb-Ta Hamiltonian (J/mol, temperature-independent
     * here); CFs are the converged values at the reference point, cross-checked
     * against both a labelled Mathematica NR dump and the positional
     * {@code phaseq} trace, which agree exactly.</p>
     *
     * <p><b>Provenance of the ECIs.</b> They are the three binary CEC sets
     * inherited into the ternary, verified cell by cell against the published
     * binary table (element order A=Mo, B=Nb, C=Ta, so AB=Mo-Nb, AC=Mo-Ta,
     * BC=Nb-Ta):</p>
     * <pre>
     *   binary   e21      e22      e3     e4
     *   Mo-Nb   -23774   -11887   +894    0
     *   Mo-Ta   -40107   -20054   -1802   0
     *   Nb-Ta    -1002     -501   -630    0
     * </pre>
     * <p>All ternary-specific terms ({@code e3ABC*}, {@code e4ABC*}) are zero:
     * binary-cluster CECs are inherited unchanged into any higher-order system
     * containing that pair, and there is no Mo-Nb-Ta ternary DFT data. Because
     * this gate reads the ECIs back from the evaluator rather than from the
     * file, it also confirms every name resolved -- an unmatched ECI name is
     * left silently at 0.0 by {@code CECEvaluator}.</p>
     */
    private static final RefTerm[] REFERENCE = {
            new RefTerm("v4AB",        0.0, 0.01314),
            new RefTerm("v4AC",        0.0, 0.0205594),
            new RefTerm("v4BC",        0.0, 0.00865153),
            new RefTerm("v4ABC1",      0.0, 0.016275),
            new RefTerm("v4ABC2",      0.0, 0.0134986),
            new RefTerm("v4ABC3",      0.0, 0.016106),
            new RefTerm("v3AB",      894.0, 0.00523559),
            new RefTerm("v3AC",    -1802.0, 0.00647258),
            new RefTerm("v3BC",     -630.0, 0.000333726),
            new RefTerm("v3ABC1",      0.0, 0.0453497),
            new RefTerm("v3ABC2",      0.0, 0.0403362),
            new RefTerm("v3ABC3",      0.0, 0.0437667),
            new RefTerm("v2AB1",  -23774.0, 0.116077),
            new RefTerm("v2AC1",  -40107.0, 0.132253),
            new RefTerm("v2BC1",   -1002.0, 0.105641),
            new RefTerm("v2AB2",  -11887.0, 0.11387),
            new RefTerm("v2AC2",  -20054.0, 0.119071),
            new RefTerm("v2BC2",    -501.0, 0.108481) };

    /** One reference cluster: its name, its ECI, and its converged CF. */
    private record RefTerm(String name, double eci, double cf) {
    }

    /**
     * The one genuine labelling difference: our {@code v3ABC1} holds what the
     * reference calls {@code v3ABC3}, and vice versa. Same three numbers, two
     * labels transposed -- confirmed by cross-comparing the triple, where each
     * of our values matches a reference value to ~1e-6 under exactly this
     * exchange and no other.
     *
     * <p>These are {@code diff}-type CVCF correlation functions (signed
     * combinations of occupation probabilities, per {@code CvCfBasis.VSpec}),
     * so which of the three symmetry-distinct ABC triangle arrangements gets
     * index 1 versus 3 is a convention, not physics. Both carry ECI 0 in this
     * Hamiltonian, so the ECI cannot discriminate between them -- this mapping
     * rests on the CF cross-comparison alone, and is noted as the weaker of the
     * two identifications.</p>
     */
    private static final Map<String, String> REF_TO_OURS_NAME = Map.of(
            "v3ABC1", "v3ABC3",
            "v3ABC3", "v3ABC1");

    /**
     * The published binary CEC table: {@code {binary, pair suffix, e21, e22,
     * e3, e4}}. Element order A=Mo, B=Nb, C=Ta fixes the suffixes.
     */
    private static final String[][] BINARY_TABLE = {
            { "Mo-Nb", "AB", "-23774", "-11887",   "894", "0" },
            { "Mo-Ta", "AC", "-40107", "-20054", "-1802", "0" },
            { "Nb-Ta", "BC",  "-1002",   "-501",  "-630", "0" } };

    /**
     * Reference {@code Gm} at the near-edge point x=[0.05, 0.05, 0.90], from a
     * Mathematica {@code phaseq} trace that converges there to residual
     * 4.88e-16. Both of our solvers stalled at this point until the entropy
     * smoothing was removed -- the true solution has a cluster variable at
     * 1.9e-08, well inside the old 1e-6 threshold, so the smoothed function
     * being minimised had its minimum somewhere else entirely.
     */
    private static final double EDGE_REF_GM = -6300.14;

    /** Scalars are quoted to 6 significant figures in the reference. */
    private static final double SCALAR_TOL = 1e-5;
    private static final double CF_TOL     = 1e-4;
    /** ECIs are exact integers in the Hamiltonian; allow only rounding slack. */
    private static final double ECI_TOL    = 0.5;

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.printf("  %s / %s / %s   T = %.0f K   x = %s%n",
                ELEMENTS, STRUCTURE, MODEL, T, fmt(REF_X));
        System.out.println("  Both minimisers vs Mathematica phaseq reference");
        System.out.println("=".repeat(78));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);

        List<String> names = CvCfBasis.getNonPointCfNames(STRUCTURE, MODEL, REF_X.length);

        // Separate model instances: the two solvers must not share state.
        CVMGibbsModel nrModel = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, session.cecEntry, null);
        CVMGibbsModel hiModel = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, session.cecEntry, null);

        // ---- Newton-Raphson, fixed composition --------------------------
        CvmNewtonSolver.Result nr = new CvmNewtonSolver(nrModel).solve(T, REF_X, 1e-5, null, null);
        System.out.printf("%n--- Newton-Raphson: converged=%s in %d iterations, |grad|=%.3e ---%n",
                nr.converged(), nr.iterations(), nr.finalGradientNorm());
        check("NR converged", nr.converged());
        if (nr.converged()) {
            compareState("NR", nr.state(), null, names);
        }

        // ---- Hillert, single phase, same seed ---------------------------
        List<HillertSolver.Phase> phases = new ArrayList<>();
        phases.add(new HillertSolver.Phase(
                "bcc", session, hiModel, 1.0, hiModel.randomStateFull(REF_X)));
        HillertSolver.Result eq = HillertSolver.solve(phases, T, 50, 10, 1.0e-6, null);
        System.out.printf("%n--- Hillert (np=1): converged=%s in %d outer iterations, residual=%.3e ---%n",
                eq.overallConverged(), eq.outerIterations(), eq.finalResidualNorm());
        check("Hillert converged", eq.overallConverged());
        if (eq.overallConverged()) {
            HillertSolver.PhaseResult p = eq.phases().get(0);
            compareState("Hillert", p.state(), eq.mu(), names);

            // With one phase there is nothing to trade mass with, so the
            // composition must not have drifted off the input.
            check("Hillert holds x (np=1)", maxAbsDiff(p.composition(), REF_X) < 1e-9);
            check("Hillert amount == 1", Math.abs(p.amount() - 1.0) < 1e-9);
        }

        // ---- The published binary table ---------------------------------
        // Checked against the evaluator's resolved ECIs, so this covers the
        // whole path from stored Hamiltonian through name matching, not just
        // the file's contents.
        if (nr.converged()) {
            System.out.printf("%n--- Published binary CEC table ---%n");
            double[] eci = nr.state().eci();
            for (String[] row : BINARY_TABLE) {
                String pair = row[1];
                for (int k = 0; k < 4; k++) {
                    String term = new String[] { "v21", "v22", "v3", "v4" }[k] + pair;
                    double expect = Double.parseDouble(row[2 + k]);
                    int j = indexOfCanonical(names, term);
                    check(String.format("%s %s = %.0f", row[0], term, expect),
                            j >= 0 && Math.abs(eci[j] - expect) <= ECI_TOL);
                }
                System.out.printf("    %-7s e21=%-8s e22=%-8s e3=%-6s e4=%s  ok%n",
                        row[0], row[2], row[3], row[4], row[5]);
            }
        }

        // ---- The two solvers against each other -------------------------
        if (nr.converged() && eq.overallConverged()) {
            CVMGibbsModel.State a = nr.state();
            CVMGibbsModel.State b = eq.phases().get(0).state();
            System.out.printf("%n--- NR vs Hillert (independent routes to the same point) ---%n");
            agree("Gm",  a.gm(),  b.gm());
            agree("Hm",  a.hm(),  b.hm());
            agree("Sm",  a.sm(),  b.sm());
            agree("G0m", a.g0m(), b.g0m());
        }

        // ---- Near-edge point: convergence regression --------------------
        // Both solvers stalled here until the entropy smoothing was removed;
        // the reference converges to residual 4.88e-16. Checked as its own
        // case because it is the failure mode, not a second sample.
        System.out.printf("%n--- Near-edge point x=[0.05, 0.05, 0.90] ---%n");
        double[] xEdge = { 0.05, 0.05, 0.90 };
        CVMGibbsModel em = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, session.cecEntry, null);
        CvmNewtonSolver.Result enr = new CvmNewtonSolver(em).solve(T, xEdge, 1e-5, null, null);
        System.out.printf("    NR      converged=%s iters=%d Gm=%.5f  (ref %.2f)%n",
                enr.converged(), enr.iterations(),
                enr.converged() ? enr.state().gm() : Double.NaN, EDGE_REF_GM);
        check("NR converges at the near-edge point", enr.converged());
        if (enr.converged()) {
            check("NR near-edge Gm matches reference",
                    Math.abs((enr.state().gm() - EDGE_REF_GM) / EDGE_REF_GM) < SCALAR_TOL);
        }

        CVMGibbsModel ehm = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, session.cecEntry, null);
        List<HillertSolver.Phase> ep = new ArrayList<>();
        ep.add(new HillertSolver.Phase("bcc", session, ehm, 1.0, ehm.randomStateFull(xEdge)));
        HillertSolver.Result eeq = HillertSolver.solve(ep, T, 50, 10, 1.0e-6, null);
        System.out.printf("    Hillert converged=%s outer=%d Gm=%.5f  (ref %.2f)%n",
                eeq.overallConverged(), eeq.outerIterations(),
                eeq.overallConverged() ? eeq.phases().get(0).state().gm() : Double.NaN,
                EDGE_REF_GM);
        check("Hillert converges at the near-edge point", eeq.overallConverged());
        if (eeq.overallConverged()) {
            check("Hillert near-edge Gm matches reference",
                    Math.abs((eeq.phases().get(0).state().gm() - EDGE_REF_GM) / EDGE_REF_GM)
                            < SCALAR_TOL);
        }

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " reference validation checks failed");
        }
    }

    /** Compares one converged state against the reference. */
    private static void compareState(String solver, CVMGibbsModel.State st,
            double[] mu, List<String> names) {

        System.out.printf("%n  %s vs reference:%n", solver);
        rel(solver + " G",  st.g(),  REF_G);
        rel(solver + " Gm", st.gm(), REF_GM);
        rel(solver + " Hm", st.hm(), REF_HM);
        rel(solver + " Sm", st.sm(), REF_SM);

        // G must decompose exactly -- an internal identity, not a reference check.
        check(solver + " G == G0m + Gm",
                Math.abs(st.g() - (st.g0m() + st.gm())) < 1e-9);

        // Each reference cluster is validated as an (ECI, CF) pair. Matching
        // is by name -- the two naming conventions are decoded by canonical()
        // -- and the ECI must agree at the slot the name resolves to, which is
        // what pins the CF to the right cluster rather than merely the right
        // label. An unresolvable name is a failure, never a skip.
        double[] u = st.u();
        double[] eci = st.eci();
        int bad = 0;
        for (RefTerm t : REFERENCE) {
            String ourName = REF_TO_OURS_NAME.getOrDefault(t.name(), t.name());
            int j = indexOfCanonical(names, ourName);
            if (j < 0) {
                bad++;
                System.out.printf("    [!] %-8s has no counterpart in our basis%n", t.name());
                continue;
            }
            // ECI first: it identifies the cluster the CF belongs to.
            double dEci = Math.abs(eci[j] - t.eci());
            if (dEci > ECI_TOL) {
                bad++;
                // Same resolved name means the cluster is agreed and the
                // value disagrees; different names mean the name resolved to
                // the wrong cluster. Say which, and keep each value with the
                // side it came from.
                boolean sameCluster = canonical(names.get(j)).equals(canonical(t.name()));
                System.out.printf("    [!] ECI %-7s: ref=%.4f  ours[%s]=%.4f  (%s)%n",
                        t.name(), t.eci(), names.get(j), eci[j],
                        sameCluster ? "value differs" : "resolved to wrong cluster");
                continue;   // CF is not meaningful if the cluster disagrees
            }
            double r = Math.abs((u[j] - t.cf()) / t.cf());
            if (r >= CF_TOL) {
                bad++;
                System.out.printf("    [!] CF  %-7s: ref=%.7f  ours[%s]=%.9f  rel=%.2e%n",
                        t.name(), t.cf(), names.get(j), u[j], r);
            }
        }
        check(solver + " all 18 (ECI, CF) pairs match reference", bad == 0);

        // The nonzero ECIs, as the evaluator resolved them. Printed because a
        // name that fails to match leaves its interaction silently at 0.0 --
        // seeing the values is what distinguishes "matched" from "defaulted".
        StringBuilder sb = new StringBuilder("    nonzero ECI: ");
        for (RefTerm t : REFERENCE) {
            if (t.eci() == 0.0) continue;
            String on = REF_TO_OURS_NAME.getOrDefault(t.name(), t.name());
            int j = indexOfCanonical(names, on);
            sb.append(String.format("%s=%.0f ", names.get(j), eci[j]));
        }
        System.out.println(sb.toString().trim());
        System.out.printf("    %d/%d (ECI, CF) pairs match  [ECI exact to %.0f J/mol, CF to %.0e rel]%n",
                REFERENCE.length - bad, REFERENCE.length, ECI_TOL, CF_TOL);

        if (mu != null) {
            // mu is underdetermined for np=1: assert the constraint, not the vector.
            double dot = 0.0;
            for (int i = 0; i < mu.length; i++) dot += mu[i] * REF_X[i];
            System.out.printf("    mu           = %s%n", fmt(mu));
            System.out.printf("    (reference)  = %s%n", fmt(REF_MU));
            System.out.printf("    sum(mu*x)    = %.6f   vs G = %.6f%n", dot, st.g());
            check(solver + " Gibbs-Duhem: sum(mu*x) == G",
                    Math.abs(dot - st.g()) / Math.abs(st.g()) < 1e-9);
        }
    }

    /**
     * Reduces either naming convention to {@code (order, pair, shell)} so the
     * two can be compared: the reference spells a pair CF {@code v2<pair><shell>}
     * and we spell it {@code v2<shell><pair>}.
     */
    private static String canonical(String name) {
        if (name.startsWith("v2") && name.length() > 3) {
            String rest = name.substring(2);
            return Character.isDigit(rest.charAt(0))
                    ? "2|" + rest.substring(1) + "|" + rest.charAt(0)   // ours
                    : "2|" + rest.substring(0, rest.length() - 1)
                            + "|" + rest.charAt(rest.length() - 1);      // reference
        }
        return name;
    }

    /** Index in {@code names} whose canonical form equals that of {@code wanted}. */
    private static int indexOfCanonical(List<String> names, String wanted) {
        String target = canonical(wanted);
        for (int i = 0; i < names.size(); i++) {
            if (canonical(names.get(i)).equals(target)) return i;
        }
        return -1;
    }

    private static void rel(String what, double ours, double ref) {
        double r = Math.abs((ours - ref) / ref);
        System.out.printf("    %-14s ours=%16.6f  ref=%12.4f  rel=%.2e%n", what, ours, ref, r);
        check(what + " matches reference", r < SCALAR_TOL);
    }

    private static void agree(String what, double a, double b) {
        double d = Math.abs(a - b);
        System.out.printf("    |%-4s(NR) - %-4s(Hillert)| = %.3e%n", what, what, d);
        check("NR and Hillert agree on " + what, d < 1e-6);
    }

    private static double maxAbsDiff(double[] a, double[] b) {
        double m = 0.0;
        for (int i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
        return m;
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

    private TernaryReferenceValidation() {
    }
}
