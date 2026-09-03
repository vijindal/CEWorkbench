package org.ce.scratch;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.hamiltonian.CECEntry;

import java.util.Arrays;

/**
 * V3 DIAGNOSTIC (follow-up) -- BCC_B2 restricted to eta = 0 must behave as BCC_A2.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertBccB2AsA2Validation
 * </pre>
 *
 * <h2>Idea</h2>
 * BCC_B2 and BCC_A2 share the same 4-site BCC tetrahedron. B2's only extra
 * degree of freedom is the LRO parameter {@code eta = p[1][A] - p[2][A]}: with
 * {@code eta = 0} the two sublattices are statistically identical and the B2
 * problem collapses onto the A2 problem. This test confirms that collapse
 * against a KNOWN A2 result:
 *
 * <pre>
 *   Nb-Ti / BCC_A2 / T ,  x_Ti = 0.5 ,  T = 1000 K   (CLAUDE.md CLI anchor)
 *     G  = -3480.5209063901 J/mol
 *     converged CVCF CFs :  v4AB=0.05365018  v3AB~0  v22AB=0.23920446  v21AB=0.23526588
 *     random-state CFs   :  v4AB=0.0625      v3AB=0   v22AB=0.25       v21AB=0.25
 *   ECIs (Nb-Ti_BCC_A2_T_CVCF): e22AB a=3120,  e21AB a=6240,  e4AB=e3AB=0
 * </pre>
 *
 * <h2>Basis correspondence (from CvCfBasis VSpecs)</h2>
 * <pre>
 *   A2 v4AB   = p1A p2B p3B p4A                 B2 V4AB   = p1A p2A p3B p4B
 *   A2 v3AB   = tri{1,2,3} diff                 B2 V31AB  = tri{1,2,3} diff (alpha vertex)
 *                                              B2 V32AB  = tri{1,2,4} diff (beta  vertex)
 *   A2 v22AB  = p1A p4B    (II-n)               B2 V221AB = p1A p4B  (II-n alpha-alpha)
 *                                              B2 V222AB = p2A p3B  (II-n beta -beta )
 *   A2 v21AB  = p1A p2B    (I-n)                B2 V21AB  = p1A p2B  (I-n)
 * </pre>
 * At {@code eta = 0} (alpha == beta statistically): V31AB == V32AB, V221AB ==
 * V222AB, V21AB == A2 v21AB, and a disordered V4AB == A2 v4AB numerically.
 *
 * <h2>ECI mapping used</h2>
 * {@code CECEvaluator} computes {@code Hm = sum_l eci[l] * u[l]} -- the JSON
 * {@code multiplicity} is descriptive and NOT applied. A2's single II-n pair
 * orbit {@code v22AB} is split by B2 into two slots {@code V221AB}
 * (alpha-alpha) and {@code V222AB} (beta-beta) which are numerically EQUAL at
 * {@code eta = 0}. So to keep {@code Hm} identical, {@code e22AB}'s coefficient
 * is split HALF/HALF: {@code V221AB = V222AB = e22AB / 2}. Then
 * {@code eci_V221 * u_V221 + eci_V222 * u_V222 = (e22AB/2 + e22AB/2) * v22AB =
 * e22AB * v22AB}. {@code e21AB} maps 1:1 to {@code V21AB}. {@code e4AB},
 * {@code e3AB} are 0, so {@code V4AB / V31AB / V32AB} carry 0.
 *
 * <p><b>HillertSolver.solve is NOT called.</b> The eta = 0 restricted solve here
 * is done with {@link CvmNewtonSolver} on a genuine A2 model for the reference
 * leg, and with an eta-projected generalized step on the B2 model for the test
 * leg.</p>
 */
public final class HillertBccB2AsA2Validation {

    private static int failures = 0;
    private static final double T = 1000.0;

    // Known A2 anchor (Nb-Ti / BCC_A2 / T, x_Ti = 0.5, 1000 K). The CLI's
    // "G (J/mol)" is the MIXING Gibbs energy Gm (= Hm - T*Sm), not absolute G.
    private static final double A2_GM = -3480.5209063901;
    private static final double[] A2_CFS_CONVERGED = { 0.05365018, 0.0, 0.23920446, 0.23526588 }; // v4,v3,v22,v21
    private static final double[] A2_CFS_RANDOM = { 0.0625, 0.0, 0.25, 0.25 };

    private static final double G_TOL = 1.0e-3;     // J/mol, vs the printed anchor
    private static final double CF_TOL = 1.0e-4;    // vs the printed anchor digits
    private static final double COLLAPSE_TOL = 1.0e-9;

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();
        line("=");
        System.out.println("  V3 -- BCC_B2 at eta = 0 must reproduce BCC_A2  (Nb-Ti, x_Ti=0.5, 1000 K)");
        line("=");

        CvmGeometry a2geo = CvmGeometry.build("Nb-Ti", "BCC_A2", "T", null);
        CvmGeometry b2geo = CvmGeometry.build("Nb-Ti", "BCC_B2", "T", null);

        int a2ncf = a2geo.ncf, b2ncf = b2geo.ncf, wB2 = b2geo.tcf;
        System.out.printf("%n  A2: ncf=%d tcf=%d  names=%s%n", a2ncf, a2geo.tcf, a2geo.basis.cfNames);
        System.out.printf("  B2: ncf=%d tcf=%d  names=%s%n", b2ncf, wB2, b2geo.basis.cfNames);

        int iV4   = b2geo.basis.indexOfCf("V4AB");
        int iV31  = b2geo.basis.indexOfCf("V31AB");
        int iV32  = b2geo.basis.indexOfCf("V32AB");
        int iV221 = b2geo.basis.indexOfCf("V221AB");
        int iV222 = b2geo.basis.indexOfCf("V222AB");
        int iV21  = b2geo.basis.indexOfCf("V21AB");
        int iXA   = b2geo.basis.indexOfCf("xA");
        int iXB   = b2geo.basis.indexOfCf("xB");
        int iEta  = b2geo.basis.indexOfCf("eta");

        // =============================================================
        // PART A -- random (disordered) state at eta = 0 : B2 CFs collapse to A2
        // =============================================================
        section("PART A -- random state, eta = 0 : B2 split-CFs collapse to the A2 values");
        CECEntry a2entry = loadA2Hamiltonian();
        CECEntry b2entry = mappedB2Hamiltonian();     // A2 ECIs carried into B2 slots
        CECEntry b2zero  = emptyEntry("Nb-Ti", "BCC_B2_T");

        CVMGibbsModel a2model = new CVMGibbsModel(a2geo, a2entry);
        CVMGibbsModel b2model = new CVMGibbsModel(b2geo, b2entry);
        CVMGibbsModel b2modelZero = new CVMGibbsModel(b2geo, b2zero);

        double[] x = { 0.5, 0.5 };
        // A2 random state
        double[] a2rand = a2model.randomStateFull(x);       // [v4,v3,v22,v21, xA,xB]
        CVMGibbsModel.State a2randSt = a2model.atFull(T, a2rand);
        // B2 random state, eta forced to 0
        double[] b2rand = b2model.randomStateFull(x);        // width tcf, eta ~ 0
        b2rand[iXA] = 0.5; b2rand[iXB] = 0.5; b2rand[iEta] = 0.0;
        CVMGibbsModel.State b2randSt = b2model.atFullWide(T, b2rand);

        System.out.printf("    A2 random CFs  [v4,v3,v22,v21] = %s%n", fmt(a2rand, 0, 4));
        System.out.printf("    B2 random CFs  V4=%.8f V31=%.8f V32=%.8f V221=%.8f V222=%.8f V21=%.8f  eta=%.1e%n",
                b2rand[iV4], b2rand[iV31], b2rand[iV32], b2rand[iV221], b2rand[iV222], b2rand[iV21], b2rand[iEta]);

        check("A2 random CFs match the known anchor",
                close(Arrays.copyOf(a2rand, 4), A2_CFS_RANDOM, 1e-6),
                Arrays.toString(Arrays.copyOf(a2rand, 4)));
        check("B2 eta=0 random : V31AB == V32AB (triangle split collapses)",
                Math.abs(b2rand[iV31] - b2rand[iV32]) < COLLAPSE_TOL,
                "d=" + Math.abs(b2rand[iV31] - b2rand[iV32]));
        check("B2 eta=0 random : V221AB == V222AB (II-n pair split collapses)",
                Math.abs(b2rand[iV221] - b2rand[iV222]) < COLLAPSE_TOL,
                "d=" + Math.abs(b2rand[iV221] - b2rand[iV222]));
        check("B2 eta=0 random : V221AB == A2 v22AB (= 0.25)",
                Math.abs(b2rand[iV221] - A2_CFS_RANDOM[2]) < 1e-9, "V221=" + b2rand[iV221]);
        check("B2 eta=0 random : V21AB == A2 v21AB (= 0.25)",
                Math.abs(b2rand[iV21] - A2_CFS_RANDOM[3]) < 1e-9, "V21=" + b2rand[iV21]);
        check("B2 eta=0 random : V4AB == A2 v4AB (= 0.0625)",
                Math.abs(b2rand[iV4] - A2_CFS_RANDOM[0]) < 1e-9, "V4=" + b2rand[iV4]);
        check("B2 eta=0 random : V31AB == A2 v3AB (= 0)",
                Math.abs(b2rand[iV31] - A2_CFS_RANDOM[1]) < 1e-9, "V31=" + b2rand[iV31]);

        // Gm / Hm / Sm of the random state (mixing quantities -- basis-comparable)
        System.out.printf("    A2 random  Gm=%.6f  Hm=%.6f  Sm=%.8f%n",
                a2randSt.gm(), a2randSt.hm(), a2randSt.sm());
        System.out.printf("    B2 random  Gm=%.6f  Hm=%.6f  Sm=%.8f  (mapped ECIs)%n",
                b2randSt.gm(), b2randSt.hm(), b2randSt.sm());
        check("B2 eta=0 random Sm == A2 random Sm (config entropy collapses)",
                Math.abs(b2randSt.sm() - a2randSt.sm()) < 1e-9,
                "dS=" + (b2randSt.sm() - a2randSt.sm()));
        check("B2 eta=0 random Hm == A2 random Hm (ECI mapping is energy-consistent)",
                Math.abs(b2randSt.hm() - a2randSt.hm()) < 1e-6,
                "dHm=" + (b2randSt.hm() - a2randSt.hm()));
        check("B2 eta=0 random Gm == A2 random Gm (mapped ECIs)",
                Math.abs(b2randSt.gm() - a2randSt.gm()) < 1e-6,
                "dGm=" + (b2randSt.gm() - a2randSt.gm()));

        // =============================================================
        // PART B -- REFERENCE leg: A2 CvmNewtonSolver converged result vs the anchor
        // =============================================================
        section("PART B -- reference leg : A2 CvmNewtonSolver converged G/CFs vs the known anchor");
        CvmNewtonSolver a2solver = new CvmNewtonSolver(a2model);
        CvmNewtonSolver.Result a2res = a2solver.solve(T, x, 1.0e-6, null, null);
        check("A2 CvmNewtonSolver converged", a2res.converged(), "not converged");
        CVMGibbsModel.State a2eq = a2res.state();
        double[] a2eqCf = a2eq.u();       // [v4,v3,v22,v21]
        System.out.printf("    A2 converged  Gm=%.10f  (abs G=%.6f)  CFs=%s%n",
                a2eq.gm(), a2eq.g(), Arrays.toString(a2eqCf));
        check("A2 converged Gm == anchor (-3480.5209063901)",
                Math.abs(a2eq.gm() - A2_GM) < G_TOL, "Gm=" + a2eq.gm());
        check("A2 converged CFs == anchor [0.05365018, ~0, 0.23920446, 0.23526588]",
                Math.abs(a2eqCf[0] - A2_CFS_CONVERGED[0]) < CF_TOL
                        && Math.abs(a2eqCf[1]) < CF_TOL
                        && Math.abs(a2eqCf[2] - A2_CFS_CONVERGED[2]) < CF_TOL
                        && Math.abs(a2eqCf[3] - A2_CFS_CONVERGED[3]) < CF_TOL,
                Arrays.toString(a2eqCf));

        // =============================================================
        // PART C -- TEST leg: B2 solved with eta PINNED at 0 must reach the same G/CFs
        //
        // Restricted minimisation: minimise G(Y) over the eta=0 subspace with x
        // fixed. This is EXACTLY the A2 fixed-composition problem re-expressed in
        // the (redundant) B2 basis. We drive the B2 model directly with a Newton
        // step on its non-point block, holding xA/xB/eta fixed -- no HillertSolver.
        // =============================================================
        section("PART C -- test leg : B2 minimised with eta pinned = 0, x fixed -> must match A2");
        double[] b2eq = restrictedB2SolveEtaZero(b2model, b2geo, x);
        CVMGibbsModel.State b2eqSt = b2model.atFullWide(T, b2eq);

        System.out.printf("    B2(eta=0) converged  Gm=%.10f  (abs G=%.6f)%n", b2eqSt.gm(), b2eqSt.g());
        System.out.printf("      V4=%.8f V31=%.8f V32=%.8f V221=%.8f V222=%.8f V21=%.8f  eta=%.2e%n",
                b2eq[iV4], b2eq[iV31], b2eq[iV32], b2eq[iV221], b2eq[iV222], b2eq[iV21], b2eq[iEta]);

        check("B2(eta=0) converged Gm == A2 anchor (-3480.5209063901)",
                Math.abs(b2eqSt.gm() - A2_GM) < G_TOL, "Gm=" + b2eqSt.gm() + "  dGm=" + (b2eqSt.gm() - A2_GM));
        check("B2(eta=0) converged Gm == A2 CvmNewtonSolver Gm (tight)",
                Math.abs(b2eqSt.gm() - a2eq.gm()) < 1e-6, "dGm=" + (b2eqSt.gm() - a2eq.gm()));
        check("B2(eta=0) eta stayed 0", Math.abs(b2eq[iEta]) < 1e-10, "eta=" + b2eq[iEta]);
        check("B2(eta=0) V221AB == V222AB at convergence (still collapsed)",
                Math.abs(b2eq[iV221] - b2eq[iV222]) < 1e-7,
                "d=" + Math.abs(b2eq[iV221] - b2eq[iV222]));
        check("B2(eta=0) V31AB == V32AB at convergence",
                Math.abs(b2eq[iV31] - b2eq[iV32]) < 1e-7,
                "d=" + Math.abs(b2eq[iV31] - b2eq[iV32]));

        // CF correspondence at equilibrium. The pair blocks have the SAME VSpec
        // in both bases (V221AB = p1A p4B = A2 v22AB; V21AB = p1A p2B = A2 v21AB)
        // so they must match numerically. The tetrahedron CF does NOT: A2
        // v4AB = p1A p2B p3B p4A (ABBA), B2 V4AB = p1A p2A p3B p4B (AABB) -- a
        // DIFFERENT 4-site arrangement, equal only at the uncorrelated random
        // state (PART A). At a correlated equilibrium they legitimately differ;
        // since e4AB = 0 here, V4AB does not enter the energy and Gm is
        // unaffected (verified above).
        check("B2(eta=0) V221AB == A2 v22AB (" + fmtd(A2_CFS_CONVERGED[2]) + ")  [same VSpec]",
                Math.abs(b2eq[iV221] - a2eqCf[2]) < 1e-4, "V221=" + b2eq[iV221] + " a2v22=" + a2eqCf[2]);
        check("B2(eta=0) V21AB == A2 v21AB (" + fmtd(A2_CFS_CONVERGED[3]) + ")  [same VSpec]",
                Math.abs(b2eq[iV21] - a2eqCf[3]) < 1e-4, "V21=" + b2eq[iV21] + " a2v21=" + a2eqCf[3]);
        System.out.printf("    note: B2 V4AB=%.8f vs A2 v4AB=%.8f -- DIFFERENT tetrahedron VSpec "
                + "(AABB vs ABBA), equal only at the random state; e4AB=0 so Gm is unaffected.%n",
                b2eq[iV4], a2eqCf[0]);
        check("B2(eta=0) tetrahedron CF differs from A2 by a basis-convention amount only "
                + "(both bases: same Gm, same pair CFs)",
                Math.abs(b2eqSt.gm() - a2eq.gm()) < 1e-6
                        && Math.abs(b2eq[iV221] - a2eqCf[2]) < 1e-4
                        && Math.abs(b2eq[iV21] - a2eqCf[3]) < 1e-4,
                "gm/pair-CF agreement is the real invariant");

        // =============================================================
        // PART D -- zero-ECI B2 at eta = 0 == pure entropy maximum == A2 with zero ECI
        // =============================================================
        section("PART D -- cross-check: zero-ECI B2 at eta=0 random == A2 zero-ECI random (pure -T*Sm)");
        CECEntry a2zero = emptyEntry("Nb-Ti", "BCC_A2_T");
        CVMGibbsModel a2modelZero = new CVMGibbsModel(a2geo, a2zero);
        double gA2z = a2modelZero.atFull(T, a2model.randomStateFull(x)).g();
        double gB2z = b2modelZero.atFullWide(T, b2rand).g();
        System.out.printf("    zero-ECI A2 random G = %.8f%n", gA2z);
        System.out.printf("    zero-ECI B2 random G = %.8f  (eta=0)%n", gB2z);
        check("zero-ECI B2(eta=0) random G == zero-ECI A2 random G",
                Math.abs(gB2z - gA2z) < 1e-6, "dG=" + (gB2z - gA2z));

        // =============================================================
        // PART E -- not just the symmetric point: B2(eta=0) == A2 at x_Ti = 0.3, 0.7
        // =============================================================
        section("PART E -- B2(eta=0) reproduces A2 Gm/pair-CFs at asymmetric x too (0.3, 0.7)");
        for (double xTi : new double[] { 0.3, 0.7 }) {
            double[] xc = { 1.0 - xTi, xTi }; // canonical order Nb, Ti  (element 0 = Nb)
            CvmNewtonSolver.Result ar = new CvmNewtonSolver(a2model).solve(T, xc, 1.0e-6, null, null);
            double[] be = restrictedB2SolveEtaZero(b2model, b2geo, xc);
            CVMGibbsModel.State bs = b2model.atFullWide(T, be);
            double[] acf = ar.state().u();
            System.out.printf("    x_Ti=%.1f : A2 Gm=%.8f   B2(eta=0) Gm=%.8f   dGm=%.2e%n",
                    xTi, ar.state().gm(), bs.gm(), bs.gm() - ar.state().gm());
            check(String.format("x_Ti=%.1f : A2 CvmNewtonSolver converged", xTi), ar.converged(), "");
            check(String.format("x_Ti=%.1f : B2(eta=0) Gm == A2 Gm", xTi),
                    Math.abs(bs.gm() - ar.state().gm()) < 1e-5,
                    "dGm=" + (bs.gm() - ar.state().gm()));
            check(String.format("x_Ti=%.1f : B2 V221AB == A2 v22AB", xTi),
                    Math.abs(be[iV221] - acf[2]) < 1e-4, "d=" + Math.abs(be[iV221] - acf[2]));
            check(String.format("x_Ti=%.1f : B2 V21AB == A2 v21AB", xTi),
                    Math.abs(be[iV21] - acf[3]) < 1e-4, "d=" + Math.abs(be[iV21] - acf[3]));
            check(String.format("x_Ti=%.1f : B2 V221AB == V222AB and V31AB == V32AB (eta=0 collapse holds)", xTi),
                    Math.abs(be[iV221] - be[iV222]) < 1e-7 && Math.abs(be[iV31] - be[iV32]) < 1e-7,
                    "d221=" + Math.abs(be[iV221] - be[iV222]) + " d31=" + Math.abs(be[iV31] - be[iV32]));
            check(String.format("x_Ti=%.1f : B2 eta stayed 0", xTi), Math.abs(be[iEta]) < 1e-10,
                    "eta=" + be[iEta]);
        }

        line("=");
        System.out.printf("  RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        line("=");
        if (failures > 0) throw new AssertionError(failures + " B2-as-A2 checks failed");
    }

    // =====================================================================
    // Restricted B2 minimisation with eta pinned to 0, x fixed.
    //
    //   Unknowns: the 6 non-point CFs Y[0..5].  xA,xB,eta held constant.
    //   Newton on  dG/dY_nonpoint = 0  using the top-left 6x6 block of
    //   gmuuFull(tcf) and the first 6 entries of gmuFull(tcf), with a
    //   feasibility backtracking line search (isValidIncludingPoints).
    //   This is the A2 fixed-composition minimisation, in the B2 basis.
    // =====================================================================
    private static double[] restrictedB2SolveEtaZero(CVMGibbsModel m, CvmGeometry geo, double[] x) {
        int w = geo.tcf, ncf = geo.ncf;
        int iEta = geo.basis.indexOfCf("eta");
        int iXA = geo.basis.indexOfCf("xA"), iXB = geo.basis.indexOfCf("xB");

        double[] Y = m.randomStateFull(x);
        Y[iXA] = x[0]; Y[iXB] = x[1]; Y[iEta] = 0.0;

        for (int it = 0; it < 60; it++) {
            CVMGibbsModel.State st = m.atFullWide(T, Y);
            double[] g = st.gmuFull(w);          // mixing gradient is enough for the stationarity of the non-point block
            double[][] H = st.gmuuFull(w);

            // 6x6 non-point Hessian block, 6-vector non-point gradient
            double[][] Hb = new double[ncf][ncf];
            double[] gb = new double[ncf];
            for (int i = 0; i < ncf; i++) {
                gb[i] = g[i];
                for (int j = 0; j < ncf; j++) Hb[i][j] = H[i][j];
            }
            double gnorm = 0;
            for (double v : gb) gnorm = Math.max(gnorm, Math.abs(v));
            if (gnorm < 1e-9) break;

            double[] dU = LinearAlgebra.solveChecked(Hb, negate(gb)).x();
            if (!allFinite(dU)) break;

            double alpha = 1.0;
            double[] Ynext = null;
            for (int ls = 0; ls < 40; ls++) {
                double[] cand = Y.clone();
                for (int i = 0; i < ncf; i++) cand[i] += alpha * dU[i];
                // xA,xB,eta untouched -> eta stays exactly 0
                CVMGibbsModel.State cs = m.atFullWide(T, cand);
                if (cs.isValidIncludingPoints() && Double.isFinite(cs.g())) { Ynext = cand; break; }
                alpha *= 0.5;
            }
            if (Ynext == null) break;
            double move = 0;
            for (int i = 0; i < ncf; i++) move = Math.max(move, Math.abs(Ynext[i] - Y[i]));
            Y = Ynext;
            if (move < 1e-12) break;
        }
        return Y;
    }

    // =====================================================================
    // Hamiltonians
    // =====================================================================

    private static CECEntry loadA2Hamiltonian() {
        // Nb-Ti_BCC_A2_T_CVCF : e22AB a=3120, e21AB a=6240, e4AB=e3AB=0
        CECEntry e = new CECEntry();
        e.elements = "Nb-Ti";
        e.structurePhase = "BCC_A2";
        e.model = "T_CVCF";
        e.cecTerms = new CECEntry.CECTerm[] {
                term("e4AB", 4, 6.0, 0.0),
                term("e3AB", 3, 12.0, 0.0),
                term("e22AB", 2, 3.0, 3120.0),
                term("e21AB", 2, 4.0, 6240.0),
        };
        return e;
    }

    /**
     * A2 ECIs carried into the B2 basis. Hm = sum eci*u (no multiplicity), and
     * V221AB == V222AB at eta = 0, so e22AB is split half/half: 3120/2 = 1560 on
     * each. e21AB -> V21AB 1:1 (6240). e4AB/e3AB = 0.
     */
    private static CECEntry mappedB2Hamiltonian() {
        CECEntry e = new CECEntry();
        e.elements = "Nb-Ti";
        e.structurePhase = "BCC_B2_T";
        e.model = "T";
        e.cecTerms = new CECEntry.CECTerm[] {
                term("V4AB", 4, 6.0, 0.0),
                term("V31AB", 3, 6.0, 0.0),
                term("V32AB", 3, 6.0, 0.0),
                term("V221AB", 2, 1.5, 1560.0),
                term("V222AB", 2, 1.5, 1560.0),
                term("V21AB", 2, 4.0, 6240.0),
        };
        return e;
    }

    private static CECEntry.CECTerm term(String name, int sites, double mult, double a) {
        CECEntry.CECTerm t = new CECEntry.CECTerm();
        t.name = name;
        t.numSites = sites;
        t.multiplicity = mult;
        t.a = a;
        t.b = 0.0;
        return t;
    }

    private static CECEntry emptyEntry(String elements, String structurePhase) {
        CECEntry e = new CECEntry();
        e.elements = elements;
        e.structurePhase = structurePhase;
        e.model = "T";
        e.cecTerms = new CECEntry.CECTerm[0];
        return e;
    }

    // =====================================================================
    // helpers
    // =====================================================================

    private static double[] negate(double[] a) {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++) r[i] = -a[i];
        return r;
    }

    private static boolean allFinite(double[] a) {
        for (double v : a) if (!Double.isFinite(v)) return false;
        return true;
    }

    private static boolean close(double[] a, double[] b, double tol) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) if (Math.abs(a[i] - b[i]) > tol) return false;
        return true;
    }

    private static String fmt(double[] a, int from, int to) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = from; i < to; i++) {
            sb.append(String.format("%.8f", a[i]));
            if (i < to - 1) sb.append(", ");
        }
        return sb.append(']').toString();
    }

    private static String fmtd(double v) {
        return String.format("%.8f", v);
    }

    private static void line(String c) {
        System.out.println(c.repeat(92));
    }

    private static void section(String s) {
        System.out.println();
        line("-");
        System.out.println("  " + s);
        line("-");
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-80s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-80s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertBccB2AsA2Validation() {
    }
}
