package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.HillertSolver;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;

/**
 * Regression gate for the P1 correction: {@code HillertSolver.PhaseStep} solves
 * the stationarity of the <b>absolute</b> Gibbs energy {@code G = G0m + Gm}
 * (x-block gradient {@code dG/dx_i}, including the pure-element reference
 * {@code G0_i(T)}), so the {@code mu} it produces -- together with
 * {@code EquilibriumMatrix}'s absolute {@code G} in its Gibbs-Duhem rows -- is
 * the absolute (SER-referenced) component chemical potential.
 *
 * <p>Before the fix, {@code PhaseStep} used the mixing-only gradient
 * {@code dGm/dx_i} while {@code EquilibriumMatrix} used absolute {@code G}: an
 * inconsistent hybrid whose {@code mu} matched neither the {@code phaseq}
 * reference nor the absolute stationarity condition. This gate would have failed
 * on the hybrid (the entry-by-entry {@code mu} check misses by ~5400 J/mol RMS,
 * and the absolute-stationarity check by ~7000 J/mol).</p>
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertAbsoluteMuValidation
 * </pre>
 *
 * <p><b>Part 1 -- np=1 CVM against the {@code phaseq} reference.</b> Mo-Nb-Ta /
 * BCC_A2 / T, 1000 K, x = [0.33, 0.33, 0.34]. Asserts the absolute {@code mu}
 * vector <em>entry-by-entry</em> against the stored {@code phaseq} values
 * (not just {@code sum(mu*x)}), plus G / Gm / G0m / Hm / Sm / final x / the
 * converged non-point CF block, and the two independent thermodynamic
 * conditions: absolute Gibbs-Duhem {@code sum(mu*x) == G} and constant
 * {@code dG/dx_i - mu_i} across i (the simplex/Lagrange condition).</p>
 *
 * <p><b>Part 2 -- synthetic analytic two-phase.</b> Two binary regular-solution
 * phases with deliberately different linear {@code G0} coefficients. The
 * absolute-gradient + absolute-{@code G} formulation must recover the exact
 * common tangent to {@code G^r(x) = G0^r(x) + Gm^r(x)}. Independent of any CVM
 * model; guards the mathematical distinction itself.</p>
 */
public final class HillertAbsoluteMuValidation {

    // ---- np=1 reference point (same as TernaryReferenceValidation) ----------
    private static final String ELEMENTS = "Mo-Nb-Ta";
    private static final String STRUCTURE = "BCC_A2";
    private static final String MODEL = "T";
    private static final double T = 1000.0;
    private static final double[] REF_X = { 0.33, 0.33, 0.34 };

    /**
     * Absolute (SER-referenced) chemical potentials from the Mathematica
     * {@code phaseq} run, canonical order Mo, Nb, Ta, J/mol. Quoted to 6
     * significant figures; the diagnostic reproduced them with {@code guFull()}
     * to RMS 0.033 J/mol.
     */
    private static final double[] REF_MU = { -73827.8, -59613.7, -74148.9 };

    private static final double REF_G = -69246.3;
    private static final double REF_GM = -20633.7;
    private static final double REF_HM = -11972.7;
    private static final double REF_SM = 8.66101;

    /**
     * Converged non-point CF block (CVCF basis, our internal order) at the
     * reference point -- from CvmEvaluatorParity / TernaryReferenceValidation's
     * own converged run, used here only as a stability anchor on {@code u}.
     */
    private static final double[] REF_U_HEAD = {
            0.01314001, 0.02055940, 0.00865153, 0.01627499, 0.01349865, 0.01610601
    };

    // ---- tolerances -------------------------------------------------------
    /** mu is quoted to 6 sig figs (~0.1 J/mol on a ~1e5 value); allow 2 J/mol. */
    private static final double MU_ABS_TOL = 2.0;
    /** Scalars quoted to 6 sig figs. */
    private static final double SCALAR_REL_TOL = 1e-5;
    /** Gibbs-Duhem is an exact identity given (mu, x, G); allow only round-off. */
    private static final double GIBBS_DUHEM_REL_TOL = 1e-9;
    /** dG/dx_i - mu_i constant across i: solver residual scale. */
    private static final double STATIONARITY_ABS_TOL = 1e-3;
    /** Converged u block: 6-sig-fig anchor. */
    private static final double U_REL_TOL = 1e-3;

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(84));
        System.out.println("  Hillert absolute-mu regression  (P1 correction: PhaseStep uses guFull())");
        System.out.println("=".repeat(84));

        part1Np1CvmReference();
        part2SyntheticTwoPhase();
        part3OuterSystemRobustness();
        part4ConvergenceLogic();
        part5LineSearchFeasibility();
        part6PhaseSetManagement();
        part7MassBalance();
        part8MassBalanceNewtonRhs();
        part9PhaseAddition();
        part10ReleaseBlockerFixes();

        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(84));
        if (failures > 0) {
            throw new AssertionError(failures + " P1 regression checks failed");
        }
    }

    // ====================================================================
    // PART 1 : np=1 CVM against the phaseq reference
    // ====================================================================

    private static void part1Np1CvmReference() throws Exception {
        System.out.printf("%n--- PART 1 : np=1  %s / %s / %s   T=%.0f  x=%s ---%n",
                ELEMENTS, STRUCTURE, MODEL, T, fmt(REF_X));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
        CVMGibbsModel model = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, session.cecEntry, null);

        int ncf = model.ncf();
        int K = model.numComponents();

        List<HillertSolver.Phase> phases = new ArrayList<>();
        phases.add(new HillertSolver.Phase(
                "bcc", session, model, 1.0, model.randomStateFull(REF_X)));
        HillertSolver.Result eq = HillertSolver.solve(phases, mbTarget(phases), T, 80, 20, 1.0e-9, null);

        check("np=1 converged", eq.overallConverged());
        HillertSolver.PhaseResult p = eq.phases().get(0);
        CVMGibbsModel.State st = p.state();
        double[] mu = eq.mu();

        System.out.printf("    outerIterations=%d  residual=%.3e%n",
                eq.outerIterations(), eq.finalResidualNorm());
        System.out.printf("    mu   = %s%n", fmt(mu));
        System.out.printf("    ref  = %s%n", fmt(REF_MU));

        // ---- PART 3 (of the STEP-4 spec): mu entry-by-entry ----
        for (int i = 0; i < K; i++) {
            double d = Math.abs(mu[i] - REF_MU[i]);
            checkMsg(String.format("mu[%d] (%s) == phaseq ref  (|%.5f - %.1f| = %.4f <= %.1f)",
                            i, elementName(i), mu[i], REF_MU[i], d, MU_ABS_TOL),
                    d <= MU_ABS_TOL);
        }

        // ---- scalars ----
        relCheck("G  == phaseq ref", st.g(), REF_G, SCALAR_REL_TOL);
        relCheck("Gm == phaseq ref", st.gm(), REF_GM, SCALAR_REL_TOL);
        relCheck("Hm == phaseq ref", st.hm(), REF_HM, SCALAR_REL_TOL);
        relCheck("Sm == phaseq ref", st.sm(), REF_SM, SCALAR_REL_TOL);

        // G0m: identity G = G0m + Gm, and G0m must be the nonzero absolute
        // reference (guards against a silent revert to a mixing-only pipeline).
        double g0m = st.g0m();
        checkMsg(String.format("G == G0m + Gm  (|%.5f - (%.5f + %.5f)| small)",
                        st.g(), g0m, st.gm()),
                Math.abs(st.g() - (g0m + st.gm())) < 1e-6);
        checkMsg(String.format("G0m is a real absolute reference (%.3f, |G0m| > 1e4)", g0m),
                Math.abs(g0m) > 1e4);

        // ---- final composition ----
        checkMsg("final x == input x (np=1 cannot move composition)",
                maxAbsDiff(p.composition(), REF_X) < 1e-9);
        checkMsg("phase amount == 1", Math.abs(p.amount() - 1.0) < 1e-9);

        // ---- converged non-point CF block ----
        double[] u = st.u();
        double worstU = 0;
        for (int i = 0; i < REF_U_HEAD.length; i++) {
            worstU = Math.max(worstU, Math.abs((u[i] - REF_U_HEAD[i]) / REF_U_HEAD[i]));
        }
        checkMsg(String.format("converged u[0..%d] stable (worst rel %.2e <= %.0e)",
                        REF_U_HEAD.length - 1, worstU, U_REL_TOL),
                worstU <= U_REL_TOL);

        // ---- PART 5 (of the STEP-4 spec): absolute Gibbs-Duhem ----
        double dot = 0;
        for (int i = 0; i < K; i++) dot += mu[i] * st.composition()[i];
        double gdRel = Math.abs(dot - st.g()) / Math.abs(st.g());
        checkMsg(String.format("Gibbs-Duhem: sum(mu*x) == G  (%.5f vs %.5f, rel %.2e <= %.0e)",
                        dot, st.g(), gdRel, GIBBS_DUHEM_REL_TOL),
                gdRel <= GIBBS_DUHEM_REL_TOL);
        // And explicitly NOT the mixing identity -- protects the reference convention.
        double dotMinusGm = Math.abs(dot - st.gm());
        checkMsg(String.format("sum(mu*x) != Gm  (differ by %.1f ~ |G0m|, NOT ~0)", dotMinusGm),
                dotMinusGm > 1e3);

        // ---- PART 4 (of the STEP-4 spec): direct absolute-stationarity ----
        // Thermodynamic condition, expressed without assuming guFull == gmuFull + G0:
        //   at equilibrium  dG/dx_i - mu_i  is the same for every i  (= -pi, the
        //   simplex Lagrange multiplier). Equivalently its traceless part -> 0.
        double[] dGdx = new double[K];
        {
            // dG/dx_i via the model's absolute widened gradient, x-block.
            double[] guFull = st.guFull();
            for (int i = 0; i < K; i++) dGdx[i] = guFull[ncf + i];
        }
        double[] absResRaw = new double[K];
        double absMean = 0;
        for (int i = 0; i < K; i++) { absResRaw[i] = dGdx[i] - mu[i]; absMean += absResRaw[i]; }
        absMean /= K;
        double absTracelessMax = 0;
        for (int i = 0; i < K; i++) {
            absTracelessMax = Math.max(absTracelessMax, Math.abs(absResRaw[i] - absMean));
        }
        System.out.printf("    dG/dx_i - mu_i          = %s   (should be constant across i)%n",
                fmt(absResRaw));
        checkMsg(String.format("ABS stationarity: dG/dx_i - mu_i constant across i "
                        + "(traceless max %.3e <= %.0e)", absTracelessMax, STATIONARITY_ABS_TOL),
                absTracelessMax <= STATIONARITY_ABS_TOL);

        // Guard the OTHER direction: the mixing gradient minus mu is NOT constant
        // (i.e. mu is not being treated as a mixing potential against dGm/dx).
        double[] dGmdx = new double[K];
        {
            double[] gmuFull = st.gmuFull();
            for (int i = 0; i < K; i++) dGmdx[i] = gmuFull[ncf + i];
        }
        double[] mixResRaw = new double[K];
        double mixMean = 0;
        for (int i = 0; i < K; i++) { mixResRaw[i] = dGmdx[i] - mu[i]; mixMean += mixResRaw[i]; }
        mixMean /= K;
        double mixTracelessMax = 0;
        for (int i = 0; i < K; i++) {
            mixTracelessMax = Math.max(mixTracelessMax, Math.abs(mixResRaw[i] - mixMean));
        }
        System.out.printf("    dGm/dx_i - mu_i         = %s   (must NOT be constant: mu is absolute)%n",
                fmt(mixResRaw));
        checkMsg(String.format("mu is NOT a mixing potential: dGm/dx_i - mu_i varies across i "
                        + "(traceless max %.1f, well above tol)", mixTracelessMax),
                mixTracelessMax > 1.0);
    }

    // ====================================================================
    // PART 2 : synthetic analytic two-phase (independent of CVM)
    // ====================================================================

    private static void part2SyntheticTwoPhase() {
        System.out.printf("%n--- PART 2 : synthetic analytic two-phase (absolute grad + absolute G) ---%n");

        final double R = 8.3144598;
        final double Tt = 1000.0;
        final double RT = R * Tt;

        // Deliberately different linear reference energies, g0A - g0B traceless
        // and genuinely per-component (no global mu shift can absorb it).
        final double aA = -30000.0, bA = -20000.0;   // g0A = (aA, bA)
        final double aB = -24000.0, bB = -26000.0;   // g0B = (aB, bB)
        final double omA = 2.6 * RT, omB = 2.6 * RT;  // Omega > 2RT => gap

        // exact common tangent to the ABSOLUTE surfaces
        double[] tie = commonTangentAbs(aA, bA, omA, aB, bB, omB, RT);
        double xA = tie[0], xB = tie[1];
        double slope = dGabs(xA, aA, bA, omA, RT);
        double muAbs1 = Gabs(xA, aA, bA, omA, RT) + (1 - xA) * slope;
        double muAbs2 = Gabs(xA, aA, bA, omA, RT) - xA * slope;

        // benchmark values quoted in the STEP-3 diagnostic
        final double BM_XA = 0.95675999, BM_XB = 0.04324001;
        final double BM_MU1 = -30327.1033, BM_MU2 = -26327.1033;
        final double TIE_TOL = 1e-5;   // 2x2 Newton converges to ~1e-12; benchmark quoted to 8 dp
        final double MU_TOL = 0.1;

        System.out.printf("    analytic common tangent: xA=%.8f  xB=%.8f%n", xA, xB);
        System.out.printf("    analytic ABS mu        : [%.4f, %.4f]%n", muAbs1, muAbs2);

        // internal self-consistency of the analytic tangent
        checkMsg(String.format("tangent parallel: dG^A/dx(xA) == dG^B/dx(xB)  (%.2e)",
                        dGabs(xA, aA, bA, omA, RT) - dGabs(xB, aB, bB, omB, RT)),
                Math.abs(dGabs(xA, aA, bA, omA, RT) - dGabs(xB, aB, bB, omB, RT)) < 1e-6);
        checkMsg(String.format("tangent colinear: G^B(xB)-G^A(xA) == slope*(xB-xA)  (%.2e)",
                        Gabs(xB, aB, bB, omB, RT) - Gabs(xA, aA, bA, omA, RT) - slope * (xB - xA)),
                Math.abs(Gabs(xB, aB, bB, omB, RT) - Gabs(xA, aA, bA, omA, RT)
                        - slope * (xB - xA)) < 1e-6);

        // matches the diagnostic benchmark
        checkMsg(String.format("xA == benchmark %.8f  (%.2e <= %.0e)", BM_XA,
                        Math.abs(xA - BM_XA), TIE_TOL), Math.abs(xA - BM_XA) <= TIE_TOL);
        checkMsg(String.format("xB == benchmark %.8f  (%.2e <= %.0e)", BM_XB,
                        Math.abs(xB - BM_XB), TIE_TOL), Math.abs(xB - BM_XB) <= TIE_TOL);
        checkMsg(String.format("mu1 == benchmark %.4f  (%.4f <= %.1f)", BM_MU1,
                        Math.abs(muAbs1 - BM_MU1), MU_TOL), Math.abs(muAbs1 - BM_MU1) <= MU_TOL);
        checkMsg(String.format("mu2 == benchmark %.4f  (%.4f <= %.1f)", BM_MU2,
                        Math.abs(muAbs2 - BM_MU2), MU_TOL), Math.abs(muAbs2 - BM_MU2) <= MU_TOL);

        // The distinction this test guards: solving the SAME construction with
        // the mixing surfaces only (formulation C) gives a DIFFERENT tie-line,
        // because the two phases' linear G0 slopes differ. If someone reverts
        // PhaseStep to gmuFull() while EquilibriumMatrix keeps absolute G, the
        // np=1 CVM part above catches it; this documents why C is also wrong.
        double[] tieMix = commonTangentMix(omA, omB, RT);
        double dispA = Math.abs(tieMix[0] - xA);
        System.out.printf("    (mixing-surface common tangent xA=%.6f, displaced from physical by %.4f)%n",
                tieMix[0], dispA);
        checkMsg("mixing-surface tie-line is displaced from the absolute one "
                        + "(confirms G0 slope difference matters)",
                dispA > 1e-3);
    }

    // ====================================================================
    // PART 3 : outer equilibrium linear system robustness (STEP 5)
    //
    // The outer system  A z = b  ( z = [mu ; deltaN] ) is assembled by
    // HillertSolver.EquilibriumMatrix.solve. STEP 5 switched its inner solve
    // from LinearAlgebra.solve to LinearAlgebra.solveChecked (same algorithm,
    // same exact answer -- solve() is literally solveChecked().x()) and made
    // HillertSolver.solve catch a singular outer matrix instead of aborting.
    // These checks exercise the ACTUAL numerical residual, not just "did it
    // return".
    // ====================================================================

    private static void part3OuterSystemRobustness() throws Exception {
        System.out.printf("%n--- PART 3 : outer equilibrium system robustness ---%n");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        // ---- 3a: well-conditioned outer system -> backward-exact solve ----
        // Re-assemble the exact np=2 outer matrix (verbatim structure from
        // EquilibriumMatrix.solve) for Nb-Ti with two distinct compositions and
        // check ||A z - b|| / ||b|| is at machine epsilon.
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
            int ncf = m.ncf(), K = m.numComponents();

            OuterSystem os = assembleOuter(m, T, ncf, K,
                    new double[][] { { 0.35, 0.65 }, { 0.65, 0.35 } },
                    new double[] { 0.5, 0.5 });
            LinearAlgebra.Solution sol = LinearAlgebra.solveChecked(os.A, os.b);
            double resid = l2(matVecMinus(os.A, sol.x(), os.b));
            double relResid = resid / l2(os.b);
            System.out.printf("    3a well-conditioned np=2: ||b||=%.3e  ||Az-b||=%.3e  relResid=%.3e%n",
                    l2(os.b), resid, relResid);
            checkMsg(String.format("3a: outer solve is backward-exact (relResid %.2e <= 1e-10)", relResid),
                    relResid <= 1e-10);
            checkMsg(String.format("3a: solveChecked.relativeResidual agrees with recomputed (%.2e vs %.2e)",
                            sol.relativeResidual(), relResid),
                    Math.abs(sol.relativeResidual() - relResid) < 1e-12 + 1e-6 * relResid);
        }

        // ---- 3b: singular outer system -> detected, HillertSolver reports
        //           non-converged instead of throwing an unchecked exception ----
        // Two stable phases with IDENTICAL composition: the two Gibbs-Duhem rows
        // become identical -> rank-deficient outer matrix.
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
            int ncf = m.ncf(), K = m.numComponents();

            OuterSystem os = assembleOuter(m, T, ncf, K,
                    new double[][] { { 0.5, 0.5 }, { 0.5, 0.5 } },
                    new double[] { 0.5, 0.5 });
            boolean threw = false;
            try {
                LinearAlgebra.solveChecked(os.A, os.b);
            } catch (RuntimeException ex) {
                threw = true;
                System.out.printf("    3b raw solveChecked on identical-composition np=2: THROW (%s)%n",
                        ex.getMessage());
            }
            checkMsg("3b: raw outer solve on identical-composition phases throws (genuine rank deficiency)",
                    threw);

            // Now the whole solver: it must NOT propagate the unchecked
            // exception -- it returns a Result with overallConverged()==false.
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 0.5, m.randomStateFull(new double[] { 0.5, 0.5 })));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.5, m.randomStateFull(new double[] { 0.5, 0.5 })));
            HillertSolver.Result r;
            boolean solverThrew = false;
            try {
                r = HillertSolver.solve(ph, mbTarget(ph), T, 50, 20, 1.0e-9, null);
                System.out.printf("    3b HillertSolver on identical phases: overallConverged=%s  outerIterations=%d%n",
                        r.overallConverged(), r.outerIterations());
                checkMsg("3b: HillertSolver returns a Result (no unchecked abort) on a singular phase set",
                        r != null);
                checkMsg("3b: that Result is reported non-converged",
                        !r.overallConverged());
            } catch (RuntimeException ex) {
                solverThrew = true;
                System.out.printf("    3b HillertSolver THREW: %s%n", ex);
            }
            checkMsg("3b: HillertSolver does NOT throw an unchecked exception on a singular outer matrix",
                    !solverThrew);
        }

        // ---- 3c: nearly-singular outer system -> large but finite relResid /
        //           high condition, still solvable, diagnostic is meaningful ----
        // Two phases 2e-6 apart in composition: near-parallel Gibbs-Duhem rows.
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
            int ncf = m.ncf(), K = m.numComponents();
            double d = 1e-6;

            OuterSystem os = assembleOuter(m, T, ncf, K,
                    new double[][] { { 0.5 - d, 0.5 + d }, { 0.5 + d, 0.5 - d } },
                    new double[] { 0.5, 0.5 });
            LinearAlgebra.Solution sol = LinearAlgebra.solveChecked(os.A, os.b);
            double relResid = l2(matVecMinus(os.A, sol.x(), os.b)) / l2(os.b);
            double cond2 = condEstimate(os.A);
            System.out.printf("    3c near-identical (|dx|=%.0e): cond2~%.2e  relResid=%.2e%n",
                    2 * d, cond2, relResid);
            // The Jacobi scaling inside solveChecked keeps this backward-exact
            // even at cond2 ~ 1e7; the point of the check is that it does NOT
            // silently degrade -- relResid stays small, cond is high.
            checkMsg(String.format("3c: near-singular system still solves backward-exact (relResid %.2e <= 1e-9)",
                            relResid), relResid <= 1e-9);
            checkMsg(String.format("3c: condition estimate flags the near-degeneracy (cond2 %.2e >= 1e4)", cond2),
                    cond2 >= 1e4);
        }

        // ---- 3d: the STEP-4 np=1 Mo-Nb-Ta regression is unchanged by STEP 5 ----
        // (mu vector, G, Gm -- must be identical to Part 1 / STEP 4.)
        {
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0, m.randomStateFull(REF_X)));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 80, 20, 1.0e-9, null);
            double[] mu = r.mu();
            double worst = 0;
            for (int i = 0; i < mu.length; i++) worst = Math.max(worst, Math.abs(mu[i] - REF_MU[i]));
            System.out.printf("    3d Mo-Nb-Ta np=1 after STEP 5: mu=%s  worst|mu-ref|=%.4f%n", fmt(mu), worst);
            checkMsg(String.format("3d: mu still matches phaseq ref entry-by-entry (worst %.4f <= %.1f)",
                            worst, MU_ABS_TOL), worst <= MU_ABS_TOL);
            checkMsg("3d: still reported converged", r.overallConverged());
        }
    }

    /** Small holder for a re-assembled outer system. */
    private static final class OuterSystem {
        final double[][] A;
        final double[] b;
        OuterSystem(double[][] A, double[] b) { this.A = A; this.b = b; }
    }

    /**
     * Re-assembles the exact outer matrix A and RHS b that
     * HillertSolver.EquilibriumMatrix.solve would build, using the real
     * PhaseStep contributions, for the given per-phase compositions and amounts.
     * Structure copied verbatim from EquilibriumMatrix.solve.
     */
    private static OuterSystem assembleOuter(CVMGibbsModel m, double T, int ncf, int K,
            double[][] comps, double[] amounts) {
        int np = comps.length;
        int n = K + np;
        double[][] A = new double[n][n];
        double[] b = new double[n];

        double[][] c0 = new double[np][];
        double[][][] sens = new double[np][][];
        double[] g = new double[np];
        double[][] x = new double[np][];
        for (int p = 0; p < np; p++) {
            double[] uFull = m.randomStateFull(comps[p]);
            HillertSolver.PhaseStep.Step step = new HillertSolver.PhaseStep(m).step(uFull, T);
            c0[p] = step.deltaComposition0();
            sens[p] = step.deltaCompositionSensitivity();
            g[p] = m.atFull(T, uFull).g();
            x[p] = new double[K];
            System.arraycopy(uFull, ncf, x[p], 0, K);
        }
        for (int p = 0; p < np; p++) {
            for (int i = 0; i < K; i++) A[p][i] = x[p][i];
            b[p] = g[p];
        }
        for (int i = 0; i < K; i++) {
            int row = np + i;
            double rhs = 0.0;
            for (int p = 0; p < np; p++) {
                A[row][K + p] = x[p][i];
                for (int k = 0; k < K; k++) A[row][k] += amounts[p] * sens[p][k][i];
                rhs -= amounts[p] * c0[p][i];
            }
            b[row] = rhs;
        }
        return new OuterSystem(A, b);
    }

    private static double[] matVecMinus(double[][] A, double[] x, double[] b) {
        int n = b.length;
        double[] r = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) s += A[i][j] * x[j];
            r[i] = s - b[i];
        }
        return r;
    }

    /** Crude 2-norm condition estimate (power / inverse-power iteration on A^T A). */
    private static double condEstimate(double[][] A) {
        int n = A.length;
        double[][] ata = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double s = 0;
                for (int k = 0; k < n; k++) s += A[k][i] * A[k][j];
                ata[i][j] = s;
            }
        double lmax = powerIter(ata);
        double lmin;
        try {
            lmin = 1.0 / powerIterInv(ata);
        } catch (RuntimeException ex) {
            return Double.POSITIVE_INFINITY;
        }
        if (lmin <= 0) return Double.POSITIVE_INFINITY;
        return Math.sqrt(lmax / lmin);
    }

    private static double powerIter(double[][] M) {
        int n = M.length;
        double[] v = new double[n];
        java.util.Arrays.fill(v, 1.0 / Math.sqrt(n));
        double lam = 0;
        for (int it = 0; it < 300; it++) {
            double[] w = new double[n];
            for (int i = 0; i < n; i++) {
                double s = 0;
                for (int j = 0; j < n; j++) s += M[i][j] * v[j];
                w[i] = s;
            }
            double nw = l2(w);
            if (nw == 0) return 0;
            for (int i = 0; i < n; i++) w[i] /= nw;
            lam = nw;
            v = w;
        }
        return lam;
    }

    private static double powerIterInv(double[][] M) {
        int n = M.length;
        double[] v = new double[n];
        java.util.Arrays.fill(v, 1.0 / Math.sqrt(n));
        double lam = 0;
        for (int it = 0; it < 300; it++) {
            double[] w = LinearAlgebra.solve(M, v);
            double nw = l2(w);
            if (nw == 0 || Double.isNaN(nw) || Double.isInfinite(nw)) return 0;
            for (int i = 0; i < n; i++) w[i] /= nw;
            lam = nw;
            v = w;
        }
        return lam;
    }

    private static double l2(double[] v) {
        double s = 0;
        for (double x : v) s += x * x;
        return Math.sqrt(s);
    }

    // ====================================================================
    // PART 4 : convergence logic (STEP 6)
    //
    // The convergence rule changed from  min_p ||deltaY_p(mu)|| <= tol  to
    //   ( a step was accepted )  AND  ( max_{p stable} ||deltaY_p(mu)|| <= tol )
    //   AND  ( outer-solve relative residual <= max(1e-8, tol) ).
    // Plus: a failed line search can no longer report CONVERGED, and a
    // two-iteration numerical stall stops with reason STALLED.
    // These checks would FAIL under the old min()-based logic.
    // ====================================================================

    private static void part4ConvergenceLogic() throws Exception {
        System.out.printf("%n--- PART 4 : convergence logic ---%n");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        // ---- TEST A: two phases, one tiny Newton step, one large. NOT converged. ----
        // Nb-Ti np=2: phase alpha seeded AT its own single-phase CVM equilibrium
        // for x=[0.5,0.5] (so its first Newton step is ~0), phase beta seeded at
        // the raw random state for x=[0.2,0.8] and far from any tie-line with
        // alpha. Run a single outer iteration. Old min()-logic: min norm is
        // alpha's ~0 -> declares CONVERGED. New max()-logic: beta's large norm
        // -> NOT converged.
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);

            // alpha's converged single-phase u at x=[0.5,0.5]
            org.ce.model.equilibrium.CvmNewtonSolver nr =
                    new org.ce.model.equilibrium.CvmNewtonSolver(m);
            org.ce.model.equilibrium.CvmNewtonSolver.Result alphaEq =
                    nr.solve(T, new double[] { 0.5, 0.5 }, 1e-9, null, null);
            double[] alphaUFull = new double[m.ncf() + m.numComponents()];
            System.arraycopy(alphaEq.u(), 0, alphaUFull, 0, m.ncf());
            alphaUFull[m.ncf()] = 0.5;
            alphaUFull[m.ncf() + 1] = 0.5;

            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 0.5, alphaUFull));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.2, 0.8 })));

            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 1, 20, 1.0e-6, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    A: reason=%s  maxPhaseStepNorm=%.3e  overallConverged=%s%n",
                    rep.reason(), rep.maxPhaseStepNorm(), r.overallConverged());
            checkMsg("A: two-phase run with one tiny + one large step is NOT reported converged",
                    !r.overallConverged());
            checkMsg("A: the convergence metric is the MAX step norm (>> tol), not the MIN",
                    rep.maxPhaseStepNorm() > 1.0e-6);
        }

        // ---- TEST C: all corrections below tolerance -> converged. ----
        // Real Mo-Nb-Ta np=1: every metric is tiny at the fixed point.
        {
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0, m.randomStateFull(REF_X)));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 80, 20, 1.0e-9, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    C: reason=%s  maxStepNorm=%.3e  outerResid=%.3e  lambda=%.3f  valid=%s%n",
                    rep.reason(), rep.maxPhaseStepNorm(), rep.linearSolveResidual(),
                    rep.lastLambda(), rep.allStatesValid());
            checkMsg("C: converges when all metrics are below tolerance",
                    r.overallConverged()
                            && rep.reason() == HillertSolver.ConvergenceReason.CONVERGED);
            checkMsg("C: report metrics are consistent with convergence "
                            + "(maxStepNorm <= 1e-9, outerResid tiny, step accepted, all valid)",
                    rep.maxPhaseStepNorm() <= 1.0e-9
                            && rep.linearSolveResidual() <= 1.0e-8
                            && rep.lastStepAccepted()
                            && rep.allStatesValid());
        }

        // ---- TEST D: line search fails on the terminal iteration -> NOT
        //      converged, reason LINE_SEARCH_FAILED, no false convergence. ----
        // The near-edge Mo-Nb-Ta point that HILLERT_SOLVER_PLAN.md sec 6e
        // documents as where the per-phase Newton step is fragile, run with
        // innerBacktrackTries = 1 (only lambda = 1 tried) and a 1-iteration cap:
        // the single line-search attempt from the raw random seed leaves the
        // physical region, so nothing is accepted and the run ends with
        // lastAccepted == false.
        {
            double Tedge = 1273.0;
            ModelSession s = builder.build(new SystemId("Mo-Nb-Ta", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Mo-Nb-Ta", "BCC_A2", "T", s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0,
                    m.randomStateFull(new double[] { 0.05, 0.475, 0.475 })));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), Tedge, 1, 1, 1.0e-9, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    D: reason=%s  lastStepAccepted=%s  overallConverged=%s  "
                            + "lastLambda=%s  iters=%d%n",
                    rep.reason(), rep.lastStepAccepted(), r.overallConverged(),
                    Double.toString(rep.lastLambda()), rep.iterationsRun());
            checkMsg("D: a completely failed line search is NOT reported converged",
                    !r.overallConverged());
            if (rep.lastStepAccepted()) {
                // The seed happened to admit a valid lambda=1 step; the test's
                // intent (no false convergence from a failed line search) still
                // holds, but this instance didn't exercise LINE_SEARCH_FAILED.
                System.out.println("      (note: lambda=1 was valid for this seed; "
                        + "LINE_SEARCH_FAILED path not exercised by this instance)");
                checkMsg("D: reason is a non-converged condition",
                        rep.reason() != HillertSolver.ConvergenceReason.CONVERGED);
            } else {
                checkMsg("D: reason == LINE_SEARCH_FAILED (distinct from MAX_ITERATIONS)",
                        rep.reason() == HillertSolver.ConvergenceReason.LINE_SEARCH_FAILED);
                checkMsg("D: lastLambda is NaN when the line search failed",
                        Double.isNaN(rep.lastLambda()));
            }
        }

        // ---- TEST E: numerically negligible accepted step -> STALLED,
        //      not a false convergence, and it stops early (not at the cap). ----
        // Run a real np=1 case to convergence at a normal tolerance, then feed
        // that already-converged state back in with tol = 1e-16 (below the
        // step-norm floor of a well-converged CVM state, ~1e-13). Every further
        // "step" is then genuinely negligible (below stallRel), so the line
        // search finds only sub-threshold feasible steps -> STEP 7 routes it to
        // STEP 6's stall path, which stops with STALLED.
        {
            double Te = 1000.0;
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);

            List<HillertSolver.Phase> warm = new ArrayList<>();
            warm.add(new HillertSolver.Phase("bcc", s, m, 1.0, m.randomStateFull(REF_X)));
            HillertSolver.Result warmRes = HillertSolver.solve(warm, mbTarget(warm), Te, 80, 20, 1.0e-9, null);
            checkMsg("E(setup): warm-up run converged", warmRes.overallConverged());

            // Re-seed a fresh phase at the converged uFull.
            double[] convergedUFull = new double[m.ncf() + m.numComponents()];
            System.arraycopy(warmRes.phases().get(0).state().u(), 0, convergedUFull, 0, m.ncf());
            System.arraycopy(REF_X, 0, convergedUFull, m.ncf(), m.numComponents());

            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0, convergedUFull));
            int cap = 40;
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), Te, cap, 30, 1.0e-16, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    E: reason=%s  overallConverged=%s  iters=%d/%d  maxStepNorm=%.3e%n",
                    rep.reason(), r.overallConverged(), rep.iterationsRun(), cap, rep.maxPhaseStepNorm());
            checkMsg("E: a negligible-progress run is NOT reported converged",
                    !r.overallConverged());
            checkMsg("E: reason == STALLED (identified the numerical stall)",
                    rep.reason() == HillertSolver.ConvergenceReason.STALLED);
            checkMsg("E: it stopped early on the stall, not at the iteration cap",
                    rep.iterationsRun() < cap);
        }

        // ---- TEST F: Mo-Nb-Ta np=1 anchor -- still converges, mu byte-identical. ----
        {
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0, m.randomStateFull(REF_X)));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 80, 20, 1.0e-9, null);
            double[] mu = r.mu();
            double worst = 0;
            for (int i = 0; i < mu.length; i++) worst = Math.max(worst, Math.abs(mu[i] - REF_MU[i]));
            System.out.printf("    F: reason=%s  mu=%s  worst|mu-ref|=%.4f%n",
                    r.convergenceReport().reason(), fmt(mu), worst);
            checkMsg("F: Mo-Nb-Ta np=1 still converges after the convergence-logic change",
                    r.overallConverged());
            checkMsg(String.format("F: mu still matches STEP-4 phaseq ref (worst %.4f <= %.1f)",
                            worst, MU_ABS_TOL), worst <= MU_ABS_TOL);
        }

        // ---- TEST G: Nb-Ti np=2 smoke -- no false convergence, state unchanged. ----
        // Mirrors HillertStateSmokeTest's config (two BCC_A2 phases at 0.35/0.65
        // and 0.65/0.35, amount 0.5 each, T=1000).
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.35, 0.65 })));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.65, 0.35 })));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 50, 20, 1.0e-6, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    G: reason=%s  converged=%s  iters=%d  maxStepNorm=%.3e%n",
                    rep.reason(), r.overallConverged(), rep.iterationsRun(), rep.maxPhaseStepNorm());
            // Expected: this collapses to a single composition (documented
            // behavior of this seed), so it may converge OR report a
            // non-converged reason -- but it must NOT be a false positive: if
            // reason==CONVERGED then maxPhaseStepNorm must actually be <= tol.
            checkMsg("G: if reported converged, the MAX phase step norm is genuinely <= tol "
                            + "(no false positive)",
                    rep.reason() != HillertSolver.ConvergenceReason.CONVERGED
                            || rep.maxPhaseStepNorm() <= 1.0e-6 + 1e-12);
            checkMsg("G: final phase G values are finite and consistent (G == state.g())",
                    Double.isFinite(r.phases().get(0).g())
                            && Double.isFinite(r.phases().get(1).g())
                            && r.phases().get(0).g() == r.phases().get(0).state().g());
        }

        // ---- TEST H (iteration-count semantics): no off-by-one. ----
        // A run capped at N iterations that does NOT converge must report
        // iterationsRun == N; a run that converges reports the iteration it
        // converged on.
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.2, 0.8 })));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.8, 0.2 })));
            HillertSolver.Result capped = HillertSolver.solve(ph, mbTarget(ph), T, 3, 20, 1.0e-14, null);
            System.out.printf("    H: capped-at-3 run: reason=%s  iterationsRun=%d%n",
                    capped.convergenceReport().reason(), capped.outerIterations());
            checkMsg("H: a run capped at 3 non-converging iterations reports iterationsRun == 3 (no off-by-one)",
                    !capped.overallConverged() && capped.outerIterations() == 3);
        }

        // ---- TEST B (deferred): all corrections small but the GLOBAL mass
        //      balance still violated. ----
        // This case needs an overall-composition target passed to solve(), which
        // it does not take today (adding it is an API change, out of scope for
        // STEP 6 -- see the step's PART 4/PART 12 notes). The convergence rule
        // here is step-norm + accepted + outer-residual only; a hypothetical
        // state with tiny per-phase steps but a broken global balance is not
        // reachable without that target, and would be caught by the deferred
        // massBalanceResidual gate. Documented, not asserted.
        System.out.println("    B: (deferred -- needs an overall-composition target input to solve(); "
                + "out of scope for STEP 6)");
    }

    // ====================================================================
    // PART 5 : line-search / step-acceptance feasibility (STEP 7)
    //
    // Changes under test:
    //   - accepted trial never has a negative or non-finite phase amount;
    //   - the halving starts at min(1, 0.5*lambda_amount) so the first probe is
    //     already amount-feasible (lambda_amount = min over stepped p with
    //     deltaN<0 of -N/deltaN);
    //   - a physically negligible feasible step is not accepted while short of
    //     convergence (routes to STALLED, not silent no-progress);
    //   - a failed line search leaves EVERY phase's amount and uFull unchanged
    //     bit-for-bit;
    //   - amount<=0 phases are frozen during the line search (no-op today since
    //     all callers seed amount>0).
    // ====================================================================

    private static void part5LineSearchFeasibility() throws Exception {
        System.out.printf("%n--- PART 5 : line-search feasibility ---%n");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        // ---- TEST A + B: a phase driven toward N < 0 must never be COMMITTED
        //      with a negative amount; N ~ 0 is allowed (phase removal is NOT
        //      done in this step). ----
        // Nb-Zr, T=1073, seeded far apart -- the known config where one phase
        // goes unstable. Every iteration's committed amount must stay >= 0.
        {
            double T = 1073.0;
            ModelSession s = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", s.cecEntry, null);

            List<HillertSolver.Phase> ph = new ArrayList<>();
            HillertSolver.Phase alpha = new HillertSolver.Phase(
                    "alpha", s, m, 0.5, m.randomStateFull(new double[] { 0.98, 0.02 }));
            HillertSolver.Phase beta = new HillertSolver.Phase(
                    "beta", s, m, 0.5, m.randomStateFull(new double[] { 0.01, 0.99 }));
            ph.add(alpha);
            ph.add(beta);

            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 200, 30, 1.0e-8, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    A/B: reason=%s  iters=%d%n", rep.reason(), rep.iterationsRun());
            for (HillertSolver.PhaseResult pr : r.phases()) {
                System.out.printf("        %-6s amount=%.10f  x=%s%n",
                        pr.label(), pr.amount(), fmt(pr.composition()));
                checkMsg(String.format("A/B: committed amount for '%s' is >= 0 (%.3e), never negative",
                                pr.label(), pr.amount()),
                        pr.amount() >= 0.0);
                checkMsg(String.format("A/B: committed amount for '%s' is finite", pr.label()),
                        Double.isFinite(pr.amount()));
            }
            // The old code committed a negative amount here (~-0.07). Assert the
            // fix actually engaged: at least one phase should be at or near zero
            // OR both clearly positive -- but NONE negative (checked above).
        }

        // ---- TEST C: a cluster variable would cross a bound at lambda=1 ->
        //      the step is reduced and the accepted trial stays valid. ----
        // Mo-Nb-Ta near-edge np=1 (x with a 0.05 minority): the full Newton step
        // routinely pushes a cluster variable out of (0,1) at lambda=1, so the
        // line search must back off. Every committed state must be valid.
        {
            double T = 1273.0;
            ModelSession s = builder.build(new SystemId("Mo-Nb-Ta", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Mo-Nb-Ta", "BCC_A2", "T", s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0,
                    m.randomStateFull(new double[] { 0.05, 0.475, 0.475 })));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 60, 40, 1.0e-6, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    C: reason=%s  iters=%d  lastLambda=%s  allStatesValid=%s%n",
                    rep.reason(), rep.iterationsRun(),
                    Double.toString(rep.lastLambda()), rep.allStatesValid());
            checkMsg("C: the final committed state is physically valid (every CV in (0,1))",
                    rep.allStatesValid());
            checkMsg("C: the phase's final G is finite",
                    Double.isFinite(r.phases().get(0).g()));
        }

        // ---- TEST D: a valid step needs an extremely small lambda -> the
        //      solver detects the numerical stall rather than accepting an
        //      effectively-zero step indefinitely. ----
        // (Covered by Part 4 TEST E, which drives a converged state with
        //  tol=1e-16; re-checked here for the STEP-7 framing.)
        {
            double Te = 1000.0;
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);
            List<HillertSolver.Phase> warm = new ArrayList<>();
            warm.add(new HillertSolver.Phase("bcc", s, m, 1.0, m.randomStateFull(REF_X)));
            HillertSolver.Result warmRes = HillertSolver.solve(warm, mbTarget(warm), Te, 80, 20, 1.0e-9, null);
            double[] cu = new double[m.ncf() + m.numComponents()];
            System.arraycopy(warmRes.phases().get(0).state().u(), 0, cu, 0, m.ncf());
            System.arraycopy(REF_X, 0, cu, m.ncf(), m.numComponents());
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0, cu));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), Te, 40, 30, 1.0e-16, null);
            System.out.printf("    D: reason=%s  iters=%d%n",
                    r.convergenceReport().reason(), r.outerIterations());
            checkMsg("D: an unreachably-tight tol on a converged state -> STALLED, not a spin to the cap",
                    r.convergenceReport().reason() == HillertSolver.ConvergenceReason.STALLED
                            && r.outerIterations() < 40);
        }

        // ---- TEST E: all trial states invalid -> LINE_SEARCH_FAILED, no state
        //      mutation, no false convergence. Plus PART 11: state immutability
        //      on failure, checked bit-for-bit. ----
        {
            double T = 1273.0;
            ModelSession s = builder.build(new SystemId("Mo-Nb-Ta", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Mo-Nb-Ta", "BCC_A2", "T", s.cecEntry, null);
            double[] seed = m.randomStateFull(new double[] { 0.05, 0.475, 0.475 });
            double[] seedCopy = seed.clone();
            HillertSolver.Phase phase = new HillertSolver.Phase("bcc", s, m, 1.0, seed);
            double amountBefore = phase.amount;
            double[] uFullBefore = phase.uFull.clone();

            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(phase);
            // innerBacktrackTries = 1: only lambda = 1 is tried; the raw Newton
            // step from this near-edge seed leaves (0,1), so the single attempt
            // fails and nothing is accepted, every iteration.
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 1, 1, 1.0e-9, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    E: reason=%s  overallConverged=%s  lastStepAccepted=%s%n",
                    rep.reason(), r.overallConverged(), rep.lastStepAccepted());
            checkMsg("E: all-invalid line search -> LINE_SEARCH_FAILED",
                    rep.reason() == HillertSolver.ConvergenceReason.LINE_SEARCH_FAILED);
            checkMsg("E: not reported converged", !r.overallConverged());

            // PART 11: the phase must be bit-for-bit unchanged.
            checkMsg("E/imm: phase.amount unchanged after failed line search (bitwise)",
                    Double.doubleToRawLongBits(phase.amount)
                            == Double.doubleToRawLongBits(amountBefore));
            checkMsg("E/imm: phase.uFull unchanged after failed line search (bitwise, all entries)",
                    bitwiseEqual(phase.uFull, uFullBefore));
            checkMsg("E/imm: the seed array the caller passed is also untouched",
                    bitwiseEqual(seed, seedCopy));
        }

        // ---- TEST F: Mo-Nb-Ta np=1 anchor -- unchanged thermodynamic result. ----
        {
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0, m.randomStateFull(REF_X)));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 80, 20, 1.0e-9, null);
            double[] mu = r.mu();
            double worst = 0;
            for (int i = 0; i < mu.length; i++) worst = Math.max(worst, Math.abs(mu[i] - REF_MU[i]));
            System.out.printf("    F: reason=%s  worst|mu-ref|=%.4f  G=%.6f%n",
                    r.convergenceReport().reason(), worst, r.phases().get(0).g());
            checkMsg("F: still converges", r.overallConverged());
            checkMsg(String.format("F: mu unchanged vs STEP-4 (worst %.4f <= %.1f)", worst, MU_ABS_TOL),
                    worst <= MU_ABS_TOL);
            relCheck("F: G unchanged", r.phases().get(0).g(), REF_G, SCALAR_REL_TOL);
        }

        // ---- TEST G: Nb-Ti np=2 smoke -- result unchanged (the old behavior did
        //      NOT rely on a negative amount for this seed; it collapses to a
        //      single composition with both amounts positive). ----
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.35, 0.65 })));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.65, 0.35 })));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 50, 20, 1.0e-6, null);
            System.out.printf("    G: reason=%s  converged=%s  amounts=[%.6f, %.6f]  G=[%.6f, %.6f]%n",
                    r.convergenceReport().reason(), r.overallConverged(),
                    r.phases().get(0).amount(), r.phases().get(1).amount(),
                    r.phases().get(0).g(), r.phases().get(1).g());
            checkMsg("G: converges (unchanged)", r.overallConverged());
            checkMsg("G: both amounts >= 0 and finite",
                    r.phases().get(0).amount() >= 0.0 && r.phases().get(1).amount() >= 0.0
                            && Double.isFinite(r.phases().get(0).amount())
                            && Double.isFinite(r.phases().get(1).amount()));
            checkMsg("G: both phases at G = -50257.822809 +/- 1e-3 (STEP-6 baseline)",
                    Math.abs(r.phases().get(0).g() - (-50257.822809)) < 1e-3
                            && Math.abs(r.phases().get(1).g() - (-50257.822809)) < 1e-3);
        }
    }

    // ====================================================================
    // PART 6 : phase-set management (STEP 8)
    //
    // Phase REMOVAL is implemented: an active phase whose amount reaches
    // numerical zero (< 1e-9 of the represented total, still shrinking, after
    // >= 3 accepted iterations) is removed -- amount set to exact 0, frozen,
    // marked inactive, a PHASE_REMOVED event recorded -- and the iteration
    // continues on the reduced active set (this is NOT a STALLED outcome).
    // Removal is exactly mass-conserving (N*x = 0).
    //
    // Phase ADDITION is NOT implemented -- it cannot be made mass-conserving
    // without an overall-composition target in the API. Tests C/D below assert
    // that a candidate seeded inactive stays inactive.
    // ====================================================================

    private static void part6PhaseSetManagement() throws Exception {
        System.out.printf("%n--- PART 6 : phase-set management ---%n");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        // ---- TEST A/B/F/G: two active phases, one driven to zero -> removed,
        //      single-phase continuation converges, PHASE_REMOVED event, no
        //      negative amount, not STALLED. ----
        // Nb-Zr T=1073, seeded far apart -- there is no real two-phase split
        // here (the convexity check confirms), so beta must go to zero.
        {
            double T = 1073.0;
            ModelSession s = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.98, 0.02 })));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.01, 0.99 })));
            List<String> log = new ArrayList<>();
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 200, 30, 1.0e-8, log::add);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();

            System.out.printf("    A/B: reason=%s  iters=%d  events=%d%n",
                    rep.reason(), rep.iterationsRun(), rep.phaseSetEvents().size());
            for (HillertSolver.PhaseSetEvent e : rep.phaseSetEvents()) {
                System.out.printf("        event: %s '%s' @iter %d  amount %.3e -> %.3e%n",
                        e.type(), e.label(), e.iteration(), e.oldAmount(), e.newAmount());
            }
            for (HillertSolver.PhaseResult pr : r.phases()) {
                System.out.printf("        %-6s amount=%.10f  x=%s  G=%.6f%n",
                        pr.label(), pr.amount(), fmt(pr.composition()), pr.g());
            }

            // A: a PHASE_REMOVED event was recorded, for 'beta'.
            boolean removedBeta = rep.phaseSetEvents().stream().anyMatch(e ->
                    e.type() == HillertSolver.PhaseSetEventType.PHASE_REMOVED
                            && e.label().equals("beta"));
            checkMsg("A: a PHASE_REMOVED event for 'beta' was recorded", removedBeta);
            checkMsg("A: exactly one phase-set event (no chatter)",
                    rep.phaseSetEvents().size() == 1);
            checkMsg("A: the removal event's newAmount is exactly 0",
                    rep.phaseSetEvents().get(0).newAmount() == 0.0);

            // F: exact zero, deterministic removal, NOT STALLED.
            checkMsg("F: outcome is not STALLED (removal is progress, not a stall)",
                    rep.reason() != HillertSolver.ConvergenceReason.STALLED);

            HillertSolver.PhaseResult alpha = r.phases().stream()
                    .filter(p -> p.label().equals("alpha")).findFirst().orElseThrow();
            HillertSolver.PhaseResult beta = r.phases().stream()
                    .filter(p -> p.label().equals("beta")).findFirst().orElseThrow();

            // B (STEP-10): the target here is 0.5*[.98,.02] + 0.5*[.01,.99]
            // = [0.495, 0.505]. After beta is removed the reduced np=1 system
            // has r_i = [0.495,0.505] - N_alpha*x_alpha != 0. STEP 10 put the
            // mass-balance residual r_i into the Newton RHS, so the reduced
            // iteration drives alpha's composition BACK to [0.495, 0.505] -- the
            // single-phase equilibrium AT the overall composition (which exists
            // here: the Nb-Zr convexity check finds no gap at 1073 K). So it now
            // legitimately CONVERGES with an essentially zero mass-balance
            // residual, instead of STEP 9's MASS_BALANCE_DRIFT.
            checkMsg("B: reduced single-phase system CONVERGES to the target "
                            + "(STEP-10 Newton mass-balance correction recovers it)",
                    r.overallConverged()
                            && rep.reason() == HillertSolver.ConvergenceReason.CONVERGED);
            checkMsg("B: the recovered state matches the target inventory (maxAbsResidual < 1e-6)",
                    rep.massBalance().maxAbsResidual() < 1e-6);
            checkMsg("B: alpha ends at the overall composition [0.495, 0.505] +/- 1e-4",
                    Math.abs(alpha.composition()[0] - 0.495) < 1e-4
                            && Math.abs(alpha.composition()[1] - 0.505) < 1e-4);
            checkMsg("B: alpha's final state is physically valid",
                    alpha.state().isValidIncludingPoints());
            checkMsg("G: beta ends at exactly 0, never negative",
                    beta.amount() == 0.0);
            checkMsg("A/B: no phase has a negative amount",
                    r.phases().stream().allMatch(p -> p.amount() >= 0.0));
            checkMsg("B: alpha holds the full represented total amount (~1.0)",
                    Math.abs(alpha.amount() - 1.0) < 1e-6);
            checkMsg("B: total phase amount is conserved (alpha + beta == 1.0)",
                    Math.abs((alpha.amount() + beta.amount()) - 1.0) < 1e-6);
        }

        // ---- TEST C/D: a candidate seeded inactive stays inactive (phase
        //      ADDITION is not implemented). ----
        // Seed one active Mo-Nb-Ta phase and one INACTIVE candidate (amount 0).
        {
            double T = 1000.0;
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            HillertSolver.Phase active = new HillertSolver.Phase(
                    "active", s, m, 1.0, m.randomStateFull(REF_X));
            HillertSolver.Phase candidate = new HillertSolver.Phase(
                    "candidate", s, m, 0.0, m.randomStateFull(new double[] { 0.8, 0.1, 0.1 }));
            ph.add(active);
            ph.add(candidate);
            System.out.printf("    C/D: seeded active=%s (active) candidate=%s (active=%s, amount=%.1f)%n",
                    active.label, candidate.label, candidate.active, candidate.amount);

            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 80, 20, 1.0e-9, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            HillertSolver.PhaseResult cand = r.phases().stream()
                    .filter(p -> p.label().equals("candidate")).findFirst().orElseThrow();
            System.out.printf("    C/D: reason=%s  events=%d  candidate final amount=%.3e%n",
                    rep.reason(), rep.phaseSetEvents().size(), cand.amount());
            checkMsg("C/D: the inactive candidate is never activated (no PHASE_ADDED)",
                    rep.phaseSetEvents().stream().noneMatch(e ->
                            e.type() == HillertSolver.PhaseSetEventType.PHASE_ADDED));
            checkMsg("C/D: candidate final amount stays 0",
                    cand.amount() == 0.0);
            checkMsg("C/D: the active phase still converges normally (np=1 path)",
                    r.overallConverged());
            // the active phase's mu must still match the STEP-4 phaseq reference
            double worst = 0;
            for (int i = 0; i < r.mu().length; i++) {
                worst = Math.max(worst, Math.abs(r.mu()[i] - REF_MU[i]));
            }
            checkMsg(String.format("C/D: active-phase mu unchanged vs STEP-4 (worst %.4f <= %.1f)",
                            worst, MU_ABS_TOL), worst <= MU_ABS_TOL);
        }

        // ---- TEST E: phase-set oscillation cannot occur in this version. ----
        // With no phase addition, a removed phase never returns -- the removal
        // is one-way, so add/remove/add/remove is structurally impossible. This
        // test documents that by running the Nb-Zr case again and asserting at
        // most ONE removal event total.
        {
            double T = 1073.0;
            ModelSession s = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.985, 0.015 })));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.5,
                    m.randomStateFull(new double[] { 0.015, 0.985 })));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 200, 30, 1.0e-8, null);
            long removals = r.convergenceReport().phaseSetEvents().stream()
                    .filter(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_REMOVED)
                    .count();
            System.out.printf("    E: total PHASE_REMOVED events = %d (one-way removal, no oscillation)%n",
                    removals);
            checkMsg("E: at most one removal per phase -- no add/remove oscillation",
                    removals <= 1);
        }

        // ---- TEST H: Mo-Nb-Ta np=1 -- no regression, no phase-set events. ----
        {
            ModelSession s = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
            CVMGibbsModel m = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s, m, 1.0, m.randomStateFull(REF_X)));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), T, 80, 20, 1.0e-9, null);
            System.out.printf("    H: reason=%s  events=%d  mu=%s%n",
                    r.convergenceReport().reason(), r.convergenceReport().phaseSetEvents().size(),
                    fmt(r.mu()));
            checkMsg("H: np=1 still converges", r.overallConverged());
            checkMsg("H: np=1 produces no phase-set events",
                    r.convergenceReport().phaseSetEvents().isEmpty());
            double worst = 0;
            for (int i = 0; i < r.mu().length; i++) {
                worst = Math.max(worst, Math.abs(r.mu()[i] - REF_MU[i]));
            }
            checkMsg(String.format("H: mu byte-stable vs STEP-4 (worst %.4f <= %.1f)", worst, MU_ABS_TOL),
                    worst <= MU_ABS_TOL);
        }
    }

    // ====================================================================
    // PART 7 : explicit mass-balance target (STEP 9)
    //
    // solve() now takes an explicit overallAmounts target. It:
    //   - validates the target (finite, >= 0, sum > 0) and REJECTS otherwise;
    //   - validates the seed represents the target (INITIAL_MASS_BALANCE);
    //   - reports an INDEPENDENT mass-balance residual (reconstructed from the
    //     accepted phase states, not the Newton RHS) in ConvergenceReport;
    //   - downgrades CONVERGED -> MASS_BALANCE_DRIFT if the represented overall
    //     composition drifted past massRelExit.
    // The Newton mass-balance equations are unchanged.
    // ====================================================================

    private static void part7MassBalance() throws Exception {
        System.out.printf("%n--- PART 7 : explicit mass-balance target ---%n");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());

        double T = 1000.0;
        ModelSession s2 = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m2 = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s2.cecEntry, null);
        ModelSession s3 = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
        CVMGibbsModel m3 = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s3.cecEntry, null);

        // ---- TEST A: exact initial mass balance -> accepted. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("a", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.35, 0.65 })));
            ph.add(new HillertSolver.Phase("b", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.65, 0.35 })));
            double[] tgt = mbTarget(ph); // exactly [0.5, 0.5]
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, T, 50, 20, 1.0e-6, null);
            HillertSolver.MassBalanceReport mb = r.convergenceReport().massBalance();
            System.out.printf("    A: reason=%s  maxAbsResid=%.3e  maxRelResid=%.3e%n",
                    r.convergenceReport().reason(), mb.maxAbsResidual(), mb.maxRelResidual());
            checkMsg("A: an exact-initial-mass-balance problem is not rejected",
                    r.convergenceReport().reason() != HillertSolver.ConvergenceReason.INITIAL_MASS_BALANCE);
            checkMsg("A: it converges", r.overallConverged());
            checkMsg("A: final mass-balance residual is tiny (<= 1e-6 relative)",
                    mb.maxRelResidual() <= 1.0e-6);
        }

        // ---- TEST B: inconsistent initial mass balance -> rejected, no
        //      silent correction. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            HillertSolver.Phase a = new HillertSolver.Phase(
                    "a", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.35, 0.65 }));
            HillertSolver.Phase b = new HillertSolver.Phase(
                    "b", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.65, 0.35 }));
            ph.add(a);
            ph.add(b);
            double amtA0 = a.amount, amtB0 = b.amount;
            double[] xA0 = a.composition(), xB0 = b.composition();
            double[] wrongTarget = { 0.30, 0.70 }; // seed represents [0.5, 0.5]
            HillertSolver.Result r = HillertSolver.solve(ph, wrongTarget, T, 50, 20, 1.0e-6, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    B: reason=%s  iters=%d  maxAbsResid=%.3e%n",
                    rep.reason(), rep.iterationsRun(), rep.massBalance().maxAbsResidual());
            checkMsg("B: an inconsistent initial state is rejected (INITIAL_MASS_BALANCE)",
                    rep.reason() == HillertSolver.ConvergenceReason.INITIAL_MASS_BALANCE);
            checkMsg("B: no iterations were run", rep.iterationsRun() == 0);
            checkMsg("B: not reported converged", !r.overallConverged());
            checkMsg("B: the residual quantifies how far off the seed is (~0.2)",
                    rep.massBalance().maxAbsResidual() > 0.1);
            checkMsg("B: NO silent correction -- initial amounts unchanged",
                    a.amount == amtA0 && b.amount == amtB0);
            checkMsg("B: NO silent correction -- initial compositions unchanged",
                    bitwiseEqual(a.composition(), xA0) && bitwiseEqual(b.composition(), xB0));
        }

        // ---- TEST C: two-phase mass conservation at the final state. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("a", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.35, 0.65 })));
            ph.add(new HillertSolver.Phase("b", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.65, 0.35 })));
            double[] tgt = mbTarget(ph);
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, T, 50, 20, 1.0e-6, null);
            // independent reconstruction from the final PhaseResults
            double[] rep = new double[2];
            for (HillertSolver.PhaseResult pr : r.phases()) {
                double[] x = pr.composition();
                for (int i = 0; i < 2; i++) rep[i] += pr.amount() * x[i];
            }
            double maxAbs = Math.max(Math.abs(rep[0] - tgt[0]), Math.abs(rep[1] - tgt[1]));
            System.out.printf("    C: target=%s  reconstructed sum N x=%s  maxAbs=%.3e%n",
                    fmt(tgt), fmt(rep), maxAbs);
            checkMsg("C: final represented inventory matches target to <= 1e-6",
                    maxAbs <= 1.0e-6);
            checkMsg("C: the report's own residual agrees with the independent one",
                    Math.abs(r.convergenceReport().massBalance().maxAbsResidual() - maxAbs) < 1e-9);
        }

        // ---- TEST D: phase removal -- target is unchanged input, and STEP 9
        //      makes the removal's mass-balance effect explicit. ----
        {
            ModelSession sZ = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel mZ = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", sZ.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.98, 0.02 })));
            ph.add(new HillertSolver.Phase("beta", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.01, 0.99 })));
            double[] tgt = mbTarget(ph); // [0.495, 0.505]
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, 1073.0, 200, 30, 1.0e-8, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            System.out.printf("    D: reason=%s  target=%s  calc=%s  maxAbsResid=%.3e%n",
                    rep.reason(), fmt(rep.massBalance().targetOverall()),
                    fmt(rep.massBalance().calculatedOverall()), rep.massBalance().maxAbsResidual());
            // the target passed in is not mutated by the solver
            checkMsg("D: the target passed in is echoed back normalized, unchanged",
                    Math.abs(rep.massBalance().targetOverall()[0] - tgt[0] / (tgt[0] + tgt[1])) < 1e-12);
            checkMsg("D: beta was removed",
                    rep.phaseSetEvents().stream().anyMatch(e ->
                            e.type() == HillertSolver.PhaseSetEventType.PHASE_REMOVED));
            // STEP 10: the mass-balance residual is in the Newton RHS, so after
            // removal the reduced set is driven BACK to the target (a single-
            // phase equilibrium at [0.495, 0.505] exists -- no gap here). It now
            // converges with a machine-zero mass-balance residual, NOT a silent
            // target reset and NOT MASS_BALANCE_DRIFT.
            checkMsg("D: the solver does NOT silently reset the target -- targetOverall "
                            + "still normalizes to the input [0.495, 0.505]",
                    Math.abs(rep.massBalance().targetOverall()[0] - 0.495) < 1e-9);
            checkMsg("D: after removal the reduced set RECOVERS the target and converges",
                    rep.reason() == HillertSolver.ConvergenceReason.CONVERGED
                            && rep.massBalance().maxAbsResidual() < 1e-6);
        }

        // ---- TEST E: single phase, target == phase composition -> converges,
        //      zero mass-balance residual. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s3, m3, 1.0, m3.randomStateFull(REF_X)));
            double[] tgt = REF_X.clone();
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, T, 80, 20, 1.0e-9, null);
            HillertSolver.MassBalanceReport mb = r.convergenceReport().massBalance();
            System.out.printf("    E: reason=%s  maxAbsResid=%.3e%n",
                    r.convergenceReport().reason(), mb.maxAbsResidual());
            checkMsg("E: np=1 with target == composition converges", r.overallConverged());
            checkMsg("E: mass-balance residual is essentially zero (np=1 x is fixed at target)",
                    mb.maxAbsResidual() < 1.0e-12);
        }

        // ---- TEST F: a component near zero but positive -> residual scaling
        //      still works (no divide-by-zero, sensible relative residual). ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s3, m3, 1.0,
                    m3.randomStateFull(new double[] { 0.001, 0.4995, 0.4995 })));
            double[] tgt = { 0.001, 0.4995, 0.4995 };
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, T, 80, 20, 1.0e-8, null);
            HillertSolver.MassBalanceReport mb = r.convergenceReport().massBalance();
            System.out.printf("    F: reason=%s  maxAbsResid=%.3e  maxRelResid=%.3e  targetOverall=%s%n",
                    r.convergenceReport().reason(), mb.maxAbsResidual(), mb.maxRelResidual(),
                    fmt(mb.targetOverall()));
            checkMsg("F: residuals are finite with a near-zero component",
                    Double.isFinite(mb.maxAbsResidual()) && Double.isFinite(mb.maxRelResidual()));
            checkMsg("F: relative residual uses the system scale, not the tiny component "
                            + "(so it is not blown up)",
                    mb.maxRelResidual() <= Math.max(mb.maxAbsResidual(), 1e-30) + 1e-15);
        }

        // ---- TEST G: NaN / Inf target -> rejected. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s3, m3, 1.0, m3.randomStateFull(REF_X)));
            for (double[] bad : new double[][] {
                    { Double.NaN, 0.33, 0.34 },
                    { 0.33, Double.POSITIVE_INFINITY, 0.34 },
                    { 0.33, -0.01, 0.68 } }) {
                HillertSolver.Result r = HillertSolver.solve(ph, bad, T, 80, 20, 1.0e-9, null);
                System.out.printf("    G: target=%s -> reason=%s%n", fmt(bad), r.convergenceReport().reason());
                checkMsg("G: invalid target " + fmt(bad) + " is rejected (INITIAL_MASS_BALANCE)",
                        r.convergenceReport().reason() == HillertSolver.ConvergenceReason.INITIAL_MASS_BALANCE);
                checkMsg("G: no iterations run for " + fmt(bad),
                        r.convergenceReport().iterationsRun() == 0);
            }
        }

        // ---- TEST H: wrong-length / zero-sum target -> rejected with a
        //      useful reason. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s3, m3, 1.0, m3.randomStateFull(REF_X)));
            for (double[] bad : new double[][] {
                    { 0.5, 0.5 },              // wrong length (K=3)
                    { 0.0, 0.0, 0.0 } }) {     // zero sum
                HillertSolver.Result r = HillertSolver.solve(ph, bad, T, 80, 20, 1.0e-9, null);
                System.out.printf("    H: target=%s -> reason=%s%n", fmt(bad), r.convergenceReport().reason());
                checkMsg("H: malformed target " + fmt(bad) + " is rejected",
                        r.convergenceReport().reason() == HillertSolver.ConvergenceReason.INITIAL_MASS_BALANCE);
            }
            HillertSolver.Result rn = HillertSolver.solve(ph, null, T, 80, 20, 1.0e-9, null);
            checkMsg("H: a null target is rejected",
                    rn.convergenceReport().reason() == HillertSolver.ConvergenceReason.INITIAL_MASS_BALANCE);
        }

        // ---- TEST I: non-normalized target is ALLOWED (contract: amounts,
        //      not fractions) as long as the seed represents it. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            // two phases, amount 1.0 each -> represented inventory sums to 2.0
            ph.add(new HillertSolver.Phase("a", s2, m2, 1.0, m2.randomStateFull(new double[] { 0.35, 0.65 })));
            ph.add(new HillertSolver.Phase("b", s2, m2, 1.0, m2.randomStateFull(new double[] { 0.65, 0.35 })));
            double[] tgt = mbTarget(ph); // [1.0, 1.0], sum 2.0 -- NOT normalized
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, T, 50, 20, 1.0e-6, null);
            HillertSolver.MassBalanceReport mb = r.convergenceReport().massBalance();
            System.out.printf("    I: non-normalized target sum=%.1f  reason=%s  maxAbsResid=%.3e  "
                            + "targetOverall=%s%n",
                    tgt[0] + tgt[1], r.convergenceReport().reason(), mb.maxAbsResidual(),
                    fmt(mb.targetOverall()));
            checkMsg("I: a non-normalized target (sum != 1) is accepted (amounts contract)",
                    r.convergenceReport().reason() != HillertSolver.ConvergenceReason.INITIAL_MASS_BALANCE);
            checkMsg("I: targetOverall is normalized in the report (sums to 1)",
                    Math.abs(mb.targetOverall()[0] + mb.targetOverall()[1] - 1.0) < 1e-12);
        }
    }

    // ====================================================================
    // PART 8 : mass-balance Newton RHS (STEP 10)
    //
    // The outer mass-balance rows are the Newton linearization of the nonlinear
    // constraint  F_i = sum_p N_p x^p_i - b_i = 0.  STEP 10 put the current
    // residual  r_i = b_i - sum_p N_p x^p_i  into the RHS, so each step drives
    // r_i -> 0 (quadratically near a solution) instead of merely holding the
    // drifted value. When r_i == 0 the outer system is bit-identical to the
    // pre-STEP-10 formulation.
    // ====================================================================

    private static void part8MassBalanceNewtonRhs() throws Exception {
        System.out.printf("%n--- PART 8 : mass-balance Newton RHS ---%n");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());
        double T = 1000.0;
        ModelSession s2 = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m2 = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s2.cecEntry, null);
        ModelSession s3 = builder.build(new SystemId(ELEMENTS, STRUCTURE, MODEL), EngineConfig.CVM, null);
        CVMGibbsModel m3 = CVMGibbsModel.of(ELEMENTS, STRUCTURE, MODEL, s3.cecEntry, null);

        // ---- TEST A: initial state exactly balanced -> r_i = 0 at iter 1, so
        //      the first Newton system is identical to the pre-STEP-10 one, and
        //      the final thermodynamic state is unchanged. ----
        // Verified by: Mo-Nb-Ta np=1 with target == REF_X still converges to the
        // exact STEP-4 mu and G (the r_i term is 0 every iteration for np=1 with
        // x fixed at the target).
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s3, m3, 1.0, m3.randomStateFull(REF_X)));
            HillertSolver.Result r = HillertSolver.solve(ph, REF_X.clone(), T, 80, 20, 1.0e-9, null);
            double worst = 0;
            for (int i = 0; i < r.mu().length; i++) worst = Math.max(worst, Math.abs(r.mu()[i] - REF_MU[i]));
            HillertSolver.MassBalanceReport mb = r.convergenceReport().massBalance();
            System.out.printf("    A: reason=%s  worst|mu-ref|=%.4f  maxAbsMassResid=%.3e  "
                            + "residBeforeLastStep=%.3e%n",
                    r.convergenceReport().reason(), worst, mb.maxAbsResidual(), mb.residualBeforeLastStep());
            checkMsg("A: exactly-balanced np=1 still converges", r.overallConverged());
            checkMsg(String.format("A: mu byte-stable vs STEP-4 (worst %.4f <= %.1f) -- r_i=0 => "
                            + "first Newton system unchanged", worst, MU_ABS_TOL),
                    worst <= MU_ABS_TOL);
            relCheck("A: G byte-stable vs STEP-4", r.phases().get(0).g(), REF_G, SCALAR_REL_TOL);
            checkMsg("A: mass residual is ~0 the whole run (r_i term never bit-changes the system)",
                    mb.maxAbsResidual() < 1e-12
                            && (Double.isNaN(mb.residualBeforeLastStep())
                                || mb.residualBeforeLastStep() < 1e-9));
        }

        // ---- TEST B: small deliberate perturbation -> the corrected RHS drives
        //      the residual back toward zero. ----
        // Two Nb-Ti phases whose SEED is exactly balanced to [0.5,0.5], but we
        // pass a target perturbed by a small amount the entry gate still admits?
        // No -- the entry gate (1e-9) would reject that. Instead: exploit phase
        // removal. Seed Nb-Zr far apart (target [0.495,0.505]); STEP 8 removes
        // beta mid-run, injecting r_i ~ 0.028; STEP 10's RHS must then recover it.
        {
            ModelSession sZ = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel mZ = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", sZ.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.98, 0.02 })));
            ph.add(new HillertSolver.Phase("beta", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.01, 0.99 })));
            double[] tgt = mbTarget(ph); // [0.495, 0.505]
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, 1073.0, 200, 30, 1.0e-8, null);
            HillertSolver.MassBalanceReport mb = r.convergenceReport().massBalance();
            System.out.printf("    B: reason=%s  iters=%d  residBeforeLastStep=%.3e  finalMaxAbs=%.3e%n",
                    r.convergenceReport().reason(), r.outerIterations(),
                    mb.residualBeforeLastStep(), mb.maxAbsResidual());
            // during the two-phase phase the residual is O(1e-3); after beta is
            // removed the corrected RHS collapses it to machine zero.
            checkMsg("B: the corrected RHS drives the perturbed residual back to ~0",
                    mb.maxAbsResidual() < 1e-6);
            checkMsg("B: it converges (the target IS representable by the reduced set here)",
                    r.overallConverged());
        }

        // ---- TEST C: Nb-Ti two-phase -- final mass-balance residual materially
        //      smaller than / no worse than before (the +r_i term pulls the
        //      first-order drift back). ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("a", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.35, 0.65 })));
            ph.add(new HillertSolver.Phase("b", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.65, 0.35 })));
            double[] tgt = mbTarget(ph);
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, T, 50, 20, 1.0e-6, null);
            HillertSolver.MassBalanceReport mb = r.convergenceReport().massBalance();
            System.out.printf("    C: reason=%s  maxAbsMassResid=%.3e  (STEP-9 was ~1.7e-8)%n",
                    r.convergenceReport().reason(), mb.maxAbsResidual());
            checkMsg("C: Nb-Ti two-phase converges", r.overallConverged());
            checkMsg("C: final mass-balance residual is no worse than STEP-9 (<= 1e-7)",
                    mb.maxAbsResidual() <= 1e-7);
        }

        // ---- TEST D: Nb-Zr phase-removal case (same as Part 6 B / Part 7 D).
        //      The solver does NOT silently reset the target; after removal it
        //      recovers the target; converges with ~0 residual. ----
        {
            ModelSession sZ = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel mZ = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", sZ.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.98, 0.02 })));
            ph.add(new HillertSolver.Phase("beta", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.01, 0.99 })));
            double[] tgt = mbTarget(ph);
            double tgt0 = tgt[0], tgt1 = tgt[1];
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, 1073.0, 200, 30, 1.0e-8, null);
            HillertSolver.ConvergenceReport rep = r.convergenceReport();
            HillertSolver.PhaseResult alpha = r.phases().stream()
                    .filter(p -> p.label().equals("alpha")).findFirst().orElseThrow();
            System.out.printf("    D: reason=%s  alpha.x=%s  targetOverall=%s%n",
                    rep.reason(), fmt(alpha.composition()), fmt(rep.massBalance().targetOverall()));
            checkMsg("D: target NOT silently reset (targetOverall == input normalized)",
                    Math.abs(rep.massBalance().targetOverall()[0] - tgt0 / (tgt0 + tgt1)) < 1e-12);
            checkMsg("D: beta removed",
                    rep.phaseSetEvents().stream().anyMatch(e ->
                            e.type() == HillertSolver.PhaseSetEventType.PHASE_REMOVED));
            checkMsg("D: after removal the reduced set recovers the target and CONVERGES",
                    rep.reason() == HillertSolver.ConvergenceReason.CONVERGED
                            && rep.massBalance().maxAbsResidual() < 1e-6);
            checkMsg("D: alpha ends AT the overall composition (not its own free minimum)",
                    Math.abs(alpha.composition()[0] - tgt0) < 1e-4);
        }

        // ---- TEST E: single-phase target exactly equal to equilibrium
        //      composition -> convergence with essentially zero mass residual. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("bcc", s3, m3, 1.0, m3.randomStateFull(REF_X)));
            HillertSolver.Result r = HillertSolver.solve(ph, REF_X.clone(), T, 80, 20, 1.0e-9, null);
            System.out.printf("    E: reason=%s  maxAbsMassResid=%.3e%n",
                    r.convergenceReport().reason(), r.convergenceReport().massBalance().maxAbsResidual());
            checkMsg("E: np=1, target == composition converges",
                    r.convergenceReport().reason() == HillertSolver.ConvergenceReason.CONVERGED);
            checkMsg("E: mass residual is essentially zero",
                    r.convergenceReport().massBalance().maxAbsResidual() < 1e-12);
        }

        // ---- TEST F: non-normalized target [1,1] -> same normalized
        //      equilibrium composition, correct total phase amount, absolute
        //      target conserved. ----
        {
            List<HillertSolver.Phase> phN = new ArrayList<>();
            phN.add(new HillertSolver.Phase("a", s2, m2, 1.0, m2.randomStateFull(new double[] { 0.35, 0.65 })));
            phN.add(new HillertSolver.Phase("b", s2, m2, 1.0, m2.randomStateFull(new double[] { 0.65, 0.35 })));
            double[] tgtN = mbTarget(phN); // [1.0, 1.0]
            HillertSolver.Result rN = HillertSolver.solve(phN, tgtN, T, 50, 20, 1.0e-6, null);

            // reference: the same physical problem normalized to total 1
            List<HillertSolver.Phase> ph1 = new ArrayList<>();
            ph1.add(new HillertSolver.Phase("a", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.35, 0.65 })));
            ph1.add(new HillertSolver.Phase("b", s2, m2, 0.5, m2.randomStateFull(new double[] { 0.65, 0.35 })));
            HillertSolver.Result r1 = HillertSolver.solve(ph1, mbTarget(ph1), T, 50, 20, 1.0e-6, null);

            double totalN = rN.phases().stream().mapToDouble(HillertSolver.PhaseResult::amount).sum();
            double[] xN = rN.convergenceReport().massBalance().calculatedOverall();
            double[] x1 = r1.convergenceReport().massBalance().calculatedOverall();
            System.out.printf("    F: total phase amount=%.6f (target sum %.1f)  "
                            + "normalized comp N=%s vs unit=%s  absMassResid=%.3e%n",
                    totalN, tgtN[0] + tgtN[1], fmt(xN), fmt(x1),
                    rN.convergenceReport().massBalance().maxAbsResidual());
            checkMsg("F: total phase amount == sum of target (2.0)",
                    Math.abs(totalN - (tgtN[0] + tgtN[1])) < 1e-6);
            checkMsg("F: normalized equilibrium composition matches the unit-total run",
                    Math.abs(xN[0] - x1[0]) < 1e-6);
            checkMsg("F: absolute target conserved (residual in MOLES, not hidden by normalization)",
                    rN.convergenceReport().massBalance().maxAbsResidual() < 1e-6);
        }
    }

    // ====================================================================
    // PART 9 : phase addition (STEP 12)
    //
    // Two halves:
    //  (1) SYNTHETIC -- a self-contained regular-solution two-phase surface
    //      with a real gap, and a hand-rolled outer Hillert loop mirroring
    //      the production STEP-10 target-aware RHS. This validates the
    //      mathematical contract the production addition logic depends on:
    //      insert at epsilon, initialise at the relaxed argmin, let the
    //      target-aware Newton system redistribute -> exact common tangent.
    //      (No verified real-CVM miscibility-gap case is available -- see the
    //      class doc / STEP 11 PART 11 -- so the physics is validated here.)
    //  (2) REAL CVM -- the production HillertSolver on Nb-Ti / Nb-Zr, which
    //      have no CVM common-tangent gap, so the candidate scan must run,
    //      evaluate the inactive phase, and DECLINE. Verifies the control
    //      logic: candidate immutability, no spurious PHASE_ADDED, removal
    //      still works, no phase-set oscillation, bounded iterations.
    // ====================================================================

    // regular-solution synthetic surface (same construction as PART 2)
    private static final double SYN_R = 8.3144598;
    private static final double SYN_T = 1000.0;
    private static final double SYN_RT = SYN_R * SYN_T;
    private static final double SYN_AA = -30000.0, SYN_BA = -20000.0;
    private static final double SYN_AB = -24000.0, SYN_BB = -26000.0;
    private static final double SYN_OMA = 2.6 * SYN_RT, SYN_OMB = 2.6 * SYN_RT;

    private static void part9PhaseAddition() throws Exception {
        System.out.printf("%n--- PART 9 : phase addition (STEP 12) ---%n");

        double[] tie = synCommonTangent();
        double xAlpha = tie[0], xBeta = tie[1];
        double slope = synDGabs(xAlpha, SYN_AA, SYN_BA, SYN_OMA);
        double muRef1 = synGabs(xAlpha, SYN_AA, SYN_BA, SYN_OMA) + (1 - xAlpha) * slope;
        double muRef2 = synGabs(xAlpha, SYN_AA, SYN_BA, SYN_OMA) - xAlpha * slope;
        System.out.printf("    synthetic tie line: x_alpha=%.8f x_beta=%.8f  mu=[%.4f, %.4f]%n",
                xAlpha, xBeta, muRef1, muRef2);

        // ---- TEST A: clear phase addition -> analytic tie line ----
        // Start alpha alone on its CONVEX arc near the alpha-rich binodal
        // (r_i large -- alpha cannot represent the mid-gap target); beta
        // inactive. The synthetic loop scores beta, inserts it at epsilon =
        // 1e-6 * total from its relaxed argmin, and the STEP-10 RHS
        // redistributes. Must reach the analytic common tangent.
        {
            double xOverall = 0.5 * xAlpha + 0.5 * xBeta;
            double[] target = { xOverall, 1.0 - xOverall };
            double xAlphaSeed = xAlpha - 0.01;
            SynResult r = synOuterWithAddition(target, xAlphaSeed, 1e-6, 300, 1e-10);
            double nB = (xOverall - xAlpha) / (xBeta - xAlpha);
            double nA = 1.0 - nB;
            System.out.printf("    A: added=%s iters=%d reason=%s  x=[%.9f, %.9f]  N=[%.9f, %.9f]%n",
                    r.added, r.iters, r.reason, r.xA, r.xB, r.nA, r.nB);
            System.out.printf("       mu=[%.5f, %.5f]  G_total=%.6f  massResid=%.3e%n",
                    r.mu[0], r.mu[1], r.gTotal, r.massResid);
            checkMsg("A: a PHASE_ADDED event occurred", r.added);
            checkMsg("A: converged to the two-phase equilibrium", r.reason.equals("CONVERGED"));
            checkMsg(String.format("A: x_alpha matches analytic tie point (%.2e)", Math.abs(r.xA - xAlpha)),
                    Math.abs(r.xA - xAlpha) < 1e-6);
            checkMsg(String.format("A: x_beta matches analytic tie point (%.2e)", Math.abs(r.xB - xBeta)),
                    Math.abs(r.xB - xBeta) < 1e-6);
            checkMsg(String.format("A: N_alpha matches lever rule %.6f (%.2e)", nA, Math.abs(r.nA - nA)),
                    Math.abs(r.nA - nA) < 1e-5);
            checkMsg(String.format("A: N_beta matches lever rule %.6f (%.2e)", nB, Math.abs(r.nB - nB)),
                    Math.abs(r.nB - nB) < 1e-5);
            checkMsg(String.format("A: mu matches analytic common tangent (%.2e, %.2e)",
                            Math.abs(r.mu[0] - muRef1), Math.abs(r.mu[1] - muRef2)),
                    Math.abs(r.mu[0] - muRef1) < 0.5 && Math.abs(r.mu[1] - muRef2) < 0.5);
            checkMsg(String.format("A: mass residual near machine precision (%.2e)", r.massResid),
                    r.massResid < 1e-9);
        }

        // ---- TEST D: selected phase begins from the RELAXED state, not the
        //      stale seed. beta is seeded at x=0.60 (far from its argmin
        //      ~0.035); after insertion its composition must be the relaxed
        //      argmin, not 0.60. ----
        {
            double xOverall = 0.5 * xAlpha + 0.5 * xBeta;
            double[] target = { xOverall, 1.0 - xOverall };
            SynResult r = synOuterWithAddition(target, xAlpha - 0.01, 1e-6, 300, 1e-10);
            System.out.printf("    D: beta stale seed x=0.60 -> after insertion+solve x_beta=%.6f%n", r.xB);
            checkMsg("D: inserted phase did NOT keep its stale seed composition (0.60)",
                    Math.abs(r.xB - 0.60) > 0.4);
            checkMsg("D: inserted phase converged to the analytic tie point",
                    Math.abs(r.xB - xBeta) < 1e-6);
        }

        // ---- TEST E: epsilon sensitivity -- final equilibrium insensitive
        //      across several decades around 1e-6. ----
        {
            double xOverall = 0.5 * xAlpha + 0.5 * xBeta;
            double[] target = { xOverall, 1.0 - xOverall };
            double xAlphaSeed = xAlpha - 0.01;
            double[] fracs = { 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 };
            boolean allMatch = true;
            for (double f : fracs) {
                SynResult r = synOuterWithAddition(target, xAlphaSeed, f, 300, 1e-10);
                boolean ok = r.reason.equals("CONVERGED")
                        && Math.abs(r.xA - xAlpha) < 1e-6 && Math.abs(r.xB - xBeta) < 1e-6;
                System.out.printf("    E: epsFrac=%.0e -> reason=%s x=[%.9f,%.9f] match=%s%n",
                        f, r.reason, r.xA, r.xB, ok);
                allMatch &= ok;
            }
            checkMsg("E: final equilibrium is epsilon-insensitive across 1e-4..1e-8", allMatch);
        }

        // ================= REAL CVM control-logic tests =================
        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());

        ModelSession sTi = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel mTi = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", sTi.cecEntry, null);

        // ---- TEST B: no addition when the system is single-phase stable.
        //      Nb-Ti at 1000 K, x=0.5 -- alpha alone is the equilibrium; the
        //      inactive beta candidate must be evaluated and DECLINED. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", sTi, mTi, 1.0,
                    seededFull(mTi, new org.ce.model.equilibrium.CvmNewtonSolver(mTi), T,
                            new double[] { 0.5, 0.5 })));
            ph.add(new HillertSolver.Phase("beta", sTi, mTi, 0.0,
                    seededFull(mTi, new org.ce.model.equilibrium.CvmNewtonSolver(mTi), T,
                            new double[] { 0.15, 0.85 })));
            double[] tgt = { 0.5, 0.5 };
            List<String> log = new ArrayList<>();
            HillertSolver.Result r = HillertSolver.solve(ph, tgt, T, 80, 20, 1.0e-8, log::add);
            boolean anyAdded = r.convergenceReport().phaseSetEvents().stream()
                    .anyMatch(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_ADDED);
            boolean evaluated = log.stream().anyMatch(l -> l.contains("candidate 'beta' evaluated"));
            System.out.printf("    B: reason=%s  addedEvent=%s  candidateEvaluated=%s%n",
                    r.convergenceReport().reason(), anyAdded, evaluated);
            checkMsg("B: single-phase CONVERGED",
                    r.convergenceReport().reason() == HillertSolver.ConvergenceReason.CONVERGED);
            checkMsg("B: NO PHASE_ADDED event (dGf below threshold)", !anyAdded);
            checkMsg("B: the inactive candidate WAS evaluated (scan ran, then declined)", evaluated);
        }

        // ---- TEST C: candidate immutability -- the inactive Phase object is
        //      byte-for-byte unchanged when it is evaluated but not selected. ----
        {
            HillertSolver.Phase alpha = new HillertSolver.Phase("alpha", sTi, mTi, 1.0,
                    seededFull(mTi, new org.ce.model.equilibrium.CvmNewtonSolver(mTi), T,
                            new double[] { 0.5, 0.5 }));
            HillertSolver.Phase beta = new HillertSolver.Phase("beta", sTi, mTi, 0.0,
                    seededFull(mTi, new org.ce.model.equilibrium.CvmNewtonSolver(mTi), T,
                            new double[] { 0.15, 0.85 }));
            double[] uBefore = beta.uFull.clone();
            double amtBefore = beta.amount;
            boolean activeBefore = beta.active;
            List<HillertSolver.Phase> ph = new ArrayList<>(List.of(alpha, beta));
            HillertSolver.solve(ph, new double[] { 0.5, 0.5 }, T, 80, 20, 1.0e-8, null);
            System.out.printf("    C: beta uFull bitwise-unchanged=%s  amount %.3e->%.3e  active %s->%s%n",
                    bitwiseEqual(uBefore, beta.uFull), amtBefore, beta.amount, activeBefore, beta.active);
            checkMsg("C: candidate uFull is byte-for-byte unchanged by the scan",
                    bitwiseEqual(uBefore, beta.uFull));
            checkMsg("C: candidate amount unchanged (still 0)", beta.amount == amtBefore);
            checkMsg("C: candidate active flag unchanged (still false)", beta.active == activeBefore);
        }

        // ---- TEST F: phase REMOVAL still works after the addition logic is in
        //      place. Nb-Zr two-phase seed -> beta removed -> exactly one
        //      PHASE_REMOVED, no PHASE_ADDED, CONVERGED, mass conserved. ----
        {
            ModelSession sZ = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel mZ = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", sZ.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.98, 0.02 })));
            ph.add(new HillertSolver.Phase("beta", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.01, 0.99 })));
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), 1073.0, 200, 30, 1.0e-8, null);
            long removed = r.convergenceReport().phaseSetEvents().stream()
                    .filter(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_REMOVED).count();
            long added = r.convergenceReport().phaseSetEvents().stream()
                    .filter(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_ADDED).count();
            System.out.printf("    F: reason=%s  removed=%d  added=%d  massResid=%.3e%n",
                    r.convergenceReport().reason(), removed, added,
                    r.convergenceReport().massBalance().maxAbsResidual());
            checkMsg("F: exactly one PHASE_REMOVED", removed == 1);
            checkMsg("F: zero PHASE_ADDED (removed phase does not oscillate back)", added == 0);
            checkMsg("F: still CONVERGED after removal",
                    r.convergenceReport().reason() == HillertSolver.ConvergenceReason.CONVERGED);
            checkMsg("F: mass conserved", r.convergenceReport().massBalance().maxAbsResidual() < 1e-6);
        }

        // ---- TEST G: no phase-set infinite loop -- across the real-CVM runs
        //      the total phase-set event count is bounded and iterations stay
        //      under the cap. ----
        {
            ModelSession sZ = builder.build(new SystemId("Nb-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
            CVMGibbsModel mZ = CVMGibbsModel.of("Nb-Zr", "BCC_A2", "T", sZ.cecEntry, null);
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.98, 0.02 })));
            ph.add(new HillertSolver.Phase("beta", sZ, mZ, 0.5, mZ.randomStateFull(new double[] { 0.01, 0.99 })));
            int cap = 200;
            HillertSolver.Result r = HillertSolver.solve(ph, mbTarget(ph), 1073.0, cap, 30, 1.0e-8, null);
            int events = r.convergenceReport().phaseSetEvents().size();
            System.out.printf("    G: iters=%d (cap %d)  phaseSetEvents=%d  reason=%s%n",
                    r.outerIterations(), cap, events, r.convergenceReport().reason());
            checkMsg("G: phase-set event count is bounded (<= 2 * nPhases)", events <= 2 * ph.size());
            checkMsg("G: terminated before the iteration cap", r.outerIterations() < cap);
        }
    }

    // ---- synthetic surface helpers (x = mole fraction of component 0) ----
    private static double synGmix(double x, double om) {
        return SYN_RT * (x * Math.log(x) + (1 - x) * Math.log(1 - x)) + om * x * (1 - x);
    }
    private static double synDGmix(double x, double om) {
        return SYN_RT * (Math.log(x) - Math.log(1 - x)) + om * (1 - 2 * x);
    }
    private static double synD2Gmix(double x, double om) {
        return SYN_RT * (1.0 / x + 1.0 / (1 - x)) - 2 * om;
    }
    private static double synGabs(double x, double a, double b, double om) {
        return a * x + b * (1 - x) + synGmix(x, om);
    }
    private static double synDGabs(double x, double a, double b, double om) {
        return (a - b) + synDGmix(x, om);
    }

    /** argmin_x ( G_phase(x) - mu . [x, 1-x] ). */
    private static double synArgminGap(double[] mu, double a, double b, double om) {
        double x = 0.5;
        for (int i = 0; i < 200; i++) {
            double g = synDGabs(x, a, b, om) - (mu[0] - mu[1]);
            double h = synD2Gmix(x, om);
            double dx = -g / h;
            double s = 1.0;
            while (s > 1e-12) {
                double nx = x + s * dx;
                if (nx > 1e-12 && nx < 1 - 1e-12) { x = nx; break; }
                s *= 0.5;
            }
            if (Math.abs(dx) < 1e-14) break;
        }
        return x;
    }

    private static double[] synCommonTangent() {
        double xA = 0.9, xB = 0.05;
        for (int it = 0; it < 500; it++) {
            double f1 = synDGabs(xA, SYN_AA, SYN_BA, SYN_OMA) - synDGabs(xB, SYN_AB, SYN_BB, SYN_OMB);
            double f2 = synGabs(xB, SYN_AB, SYN_BB, SYN_OMB) - synGabs(xA, SYN_AA, SYN_BA, SYN_OMA)
                    - synDGabs(xA, SYN_AA, SYN_BA, SYN_OMA) * (xB - xA);
            double j11 = synD2Gmix(xA, SYN_OMA);
            double j12 = -synD2Gmix(xB, SYN_OMB);
            double j21 = -synD2Gmix(xA, SYN_OMA) * (xB - xA);
            double j22 = synDGabs(xB, SYN_AB, SYN_BB, SYN_OMB) - synDGabs(xA, SYN_AA, SYN_BA, SYN_OMA);
            double det = j11 * j22 - j12 * j21;
            if (Math.abs(det) < 1e-40) break;
            double dxA = -(j22 * f1 - j12 * f2) / det;
            double dxB = -(-j21 * f1 + j11 * f2) / det;
            double s = 1.0;
            while (s > 1e-10) {
                double nxA = xA + s * dxA, nxB = xB + s * dxB;
                if (nxA > 1e-10 && nxA < 1 - 1e-10 && nxB > 1e-10 && nxB < 1 - 1e-10) {
                    xA = nxA; xB = nxB; break;
                }
                s *= 0.5;
            }
            if (Math.abs(f1) < 1e-9 && Math.abs(f2) < 1e-7) break;
        }
        return new double[] { xA, xB };
    }

    private static final class SynResult {
        boolean added;
        int iters;
        String reason;
        double xA, xB, nA, nB, gTotal, massResid;
        double[] mu;
    }

    /**
     * Hand-rolled outer Hillert loop with STEP-10 target-aware RHS AND STEP-12
     * phase addition, on the synthetic surface. Slot 0 = alpha (SYN_AA/BA/OMA),
     * slot 1 = beta (SYN_AB/BB/OMB, inactive at start). Mirrors the production
     * algorithm: score the inactive candidate at the current mu, insert it at
     * epsilon from its relaxed argmin, then let the target-aware Newton system
     * redistribute. This is the STEP-11 PART 7/10 experiment, kept as a gate.
     */
    private static SynResult synOuterWithAddition(double[] target, double xAlphaSeed,
            double epsFrac, int maxIter, double tol) {
        double total = target[0] + target[1];
        double epsilon = epsFrac * total;

        double nA = total, xA = xAlphaSeed;
        double nB = 0.0, xB = 0.60;   // deliberately stale seed (TEST D)
        boolean betaActive = false;
        double[] mu = { 0, 0 };
        boolean added = false;
        String reason = "MAX_ITERATIONS";
        int it = 0;

        for (; it < maxIter; it++) {
            // --- BEFORE beta is active: hold alpha on its convex seed arc and
            //     immediately try to add beta. alpha must NOT be driven toward
            //     the (concave) mid-gap target while it is the sole phase --
            //     that is the unstable np=1 point STEP 11 PART 10 warns about.
            //     The production PhaseStep handles a real np=1 correctly; this
            //     synthetic just skips straight to the addition scan. ---
            if (!betaActive) {
                double gA = synGabs(xA, SYN_AA, SYN_BA, SYN_OMA);
                double dgA = synDGabs(xA, SYN_AA, SYN_BA, SYN_OMA);
                mu = new double[] { gA + (1 - xA) * dgA, gA - xA * dgA };
                double xbRelaxed = synArgminGap(mu, SYN_AB, SYN_BB, SYN_OMB);
                double gB = synGabs(xbRelaxed, SYN_AB, SYN_BB, SYN_OMB);
                double dgf = mu[0] * xbRelaxed + mu[1] * (1 - xbRelaxed) - gB;
                double addThreshold = Math.max(1.0, 1e-6 * Math.abs(gB));
                if (dgf > addThreshold) {
                    xB = xbRelaxed;              // RELAXED state, not the stale 0.60
                    nB = epsilon;
                    betaActive = true;
                    added = true;
                    continue;                   // widened system solves next iter
                }
                reason = "NO_FAVOURABLE_CANDIDATE";
                break;
            }

            // --- solve the current (two-phase) active set ---
            double hA = synD2Gmix(xA, SYN_OMA);
            double sensA0 = 1.0 / hA, sensA1 = -1.0 / hA;
            double dxA_0 = -synDGabs(xA, SYN_AA, SYN_BA, SYN_OMA) / hA;

            double[] muNew;
            double dNa, dNb, dxA, dxB;
            {
                double hB = synD2Gmix(xB, SYN_OMB);
                double sensB0 = 1.0 / hB, sensB1 = -1.0 / hB;
                double dxB_0 = -synDGabs(xB, SYN_AB, SYN_BB, SYN_OMB) / hB;

                double[][] A = new double[4][4];
                double[] b = new double[4];
                A[0][0] = xA; A[0][1] = 1 - xA; b[0] = synGabs(xA, SYN_AA, SYN_BA, SYN_OMA);
                A[1][0] = xB; A[1][1] = 1 - xB; b[1] = synGabs(xB, SYN_AB, SYN_BB, SYN_OMB);
                double rep0 = nA * xA + nB * xB;
                A[2][2] = xA; A[2][3] = xB;
                A[2][0] = nA * sensA0 + nB * sensB0;
                A[2][1] = nA * sensA1 + nB * sensB1;
                b[2] = (target[0] - rep0) - (nA * dxA_0 + nB * dxB_0);
                double rep1 = nA * (1 - xA) + nB * (1 - xB);
                A[3][2] = 1 - xA; A[3][3] = 1 - xB;
                A[3][0] = -(nA * sensA0 + nB * sensB0);
                A[3][1] = -(nA * sensA1 + nB * sensB1);
                b[3] = (target[1] - rep1) - (nA * (-dxA_0) + nB * (-dxB_0));
                double[] sol = synSolve4(A, b);
                if (sol == null) { reason = "SINGULAR"; break; }
                muNew = new double[] { sol[0], sol[1] };
                dNa = sol[2]; dNb = sol[3];
                dxA = dxA_0 + sensA0 * muNew[0] + sensA1 * muNew[1];
                dxB = dxB_0 + sensB0 * muNew[0] + sensB1 * muNew[1];
            }

            // line search (amount feasibility + x in (0,1))
            double lamAmt = Double.POSITIVE_INFINITY;
            if (dNa < 0) lamAmt = Math.min(lamAmt, -nA / dNa);
            if (dNb < 0) lamAmt = Math.min(lamAmt, -nB / dNb);
            double lambda = Math.min(1.0, 0.5 * lamAmt);
            boolean accepted = false;
            double stepNorm = Double.POSITIVE_INFINITY;
            for (int t = 0; t < 60; t++) {
                double nAt = nA + lambda * dNa, nBt = nB + lambda * dNb;
                double xAt = xA + lambda * dxA, xBt = xB + lambda * dxB;
                if (nAt >= 0 && nBt >= 0 && xAt > 1e-12 && xAt < 1 - 1e-12
                        && xBt > 1e-12 && xBt < 1 - 1e-12) {
                    stepNorm = Math.max(Math.abs(lambda * dxA), Math.abs(lambda * dxB));
                    stepNorm = Math.max(stepNorm,
                            Math.max(Math.abs(lambda * dNa), Math.abs(lambda * dNb)));
                    nA = nAt; nB = nBt; xA = xAt; xB = xBt; mu = muNew;
                    accepted = true;
                    break;
                }
                lambda *= 0.5;
            }
            if (!accepted) { reason = "LINE_SEARCH_FAILED"; break; }

            // convergence: joint step small AND mass residual small
            double massResid = Math.max(
                    Math.abs(target[0] - (nA * xA + nB * xB)),
                    Math.abs(target[1] - (nA * (1 - xA) + nB * (1 - xB))));
            if (stepNorm <= tol && massResid <= Math.max(1e-9, 10 * tol)) {
                reason = "CONVERGED"; it++; break;
            }
        }

        SynResult r = new SynResult();
        r.added = added;
        r.iters = it;
        r.reason = reason;
        r.xA = xA; r.xB = xB; r.nA = nA; r.nB = nB;
        r.mu = mu;
        r.gTotal = nA * synGabs(xA, SYN_AA, SYN_BA, SYN_OMA)
                + nB * synGabs(xB, SYN_AB, SYN_BB, SYN_OMB);
        r.massResid = Math.max(
                Math.abs(target[0] - (nA * xA + nB * xB)),
                Math.abs(target[1] - (nA * (1 - xA) + nB * (1 - xB))));
        return r;
    }

    private static double[] synSolve4(double[][] Ain, double[] bin) {
        int n = 4;
        double[][] A = new double[n][n];
        double[] b = bin.clone();
        for (int i = 0; i < n; i++) A[i] = Ain[i].clone();
        for (int col = 0; col < n; col++) {
            int piv = col;
            for (int r = col + 1; r < n; r++) if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) piv = r;
            if (Math.abs(A[piv][col]) < 1e-30) return null;
            double[] tmp = A[col]; A[col] = A[piv]; A[piv] = tmp;
            double tb = b[col]; b[col] = b[piv]; b[piv] = tb;
            for (int r = col + 1; r < n; r++) {
                double f = A[r][col] / A[col][col];
                for (int c = col; c < n; c++) A[r][c] -= f * A[col][c];
                b[r] -= f * b[col];
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double s = b[i];
            for (int c = i + 1; c < n; c++) s -= A[i][c] * x[c];
            x[i] = s / A[i][i];
        }
        return x;
    }

    // ====================================================================
    // PART 10 : release-blocker fixes (independent review, FINAL FIX STEP)
    //
    //   N1  -- a non-finite outer solve is caught and reported as
    //          ConvergenceReason.NUMERICAL_BREAKDOWN, not propagated and not
    //          presented as a plausible mu.
    //   R2  -- CancellationException from the phase-addition candidate scan
    //          propagates unchanged out of HillertSolver.solve.
    //   E2  -- relaxCandidate no longer has a blanket catch: an EXPECTED
    //          per-grid-point failure (CvmNewtonSolver returns non-converged)
    //          is skipped; there is no path that swallows an unexpected
    //          RuntimeException.
    // ====================================================================

    private static void part10ReleaseBlockerFixes() throws Exception {
        System.out.printf("%n--- PART 10 : release-blocker fixes (N1 non-finite, R2/E2 cancellation) ---%n");

        Workspace ws = new Workspace();
        CEWorkbenchContext ctx = new CEWorkbenchContext(ws);
        ModelSession.Builder builder = new ModelSession.Builder(ctx.getHamiltonianStore());
        ModelSession s = builder.build(new SystemId("Nb-Ti", "BCC_A2", "T"), EngineConfig.CVM, null);
        CVMGibbsModel m = CVMGibbsModel.of("Nb-Ti", "BCC_A2", "T", s.cecEntry, null);
        final double T = 1000.0;

        // ---- N1: a phase seeded at a PURE-ELEMENT composition x=[1,0]. Its
        //      randomStateFull and its absolute G are finite, but guFull() has a
        //      log(0) in the x-block -> NaN. That poisons PhaseStep ->
        //      EquilibriumMatrix.solve returns a NaN mu with a NaN residual,
        //      which LinearAlgebra.solveChecked's exact-singularity guard does
        //      NOT catch. The N1 guard must catch it: NUMERICAL_BREAKDOWN, no
        //      exception, mu left finite, no further PhaseStep. ----
        {
            double[] pureFull = m.randomStateFull(new double[] { 1.0, 0.0 });
            check("N1: seed randomStateFull([1,0]) is itself finite (the breakdown is downstream)",
                    java.util.Arrays.stream(pureFull).allMatch(Double::isFinite));

            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("pure", s, m, 1.0, pureFull));
            double[] tgt = { 1.0, 0.0 };   // matches the seed -> entry checks pass

            List<String> log = new ArrayList<>();
            HillertSolver.Result r = null;
            boolean threw = false;
            try {
                r = HillertSolver.solve(ph, tgt, T, 30, 20, 1.0e-9, log::add);
            } catch (RuntimeException ex) {
                threw = true;
                System.out.printf("    N1: solve THREW %s%n", ex);
            }
            check("N1: solve does NOT throw an uncaught exception on a non-finite outer system", !threw);
            check("N1: a Result is returned", r != null);
            if (r != null) {
                HillertSolver.ConvergenceReport rep = r.convergenceReport();
                System.out.printf("    N1: reason=%s  overallConverged=%s  mu=%s  iters=%d  lastStepAccepted=%s%n",
                        rep.reason(), r.overallConverged(), fmt(r.mu()),
                        rep.iterationsRun(), rep.lastStepAccepted());
                check("N1: reason == NUMERICAL_BREAKDOWN (not SINGULAR_OUTER_SYSTEM, not CONVERGED)",
                        rep.reason() == HillertSolver.ConvergenceReason.NUMERICAL_BREAKDOWN);
                check("N1: overallConverged == false", !r.overallConverged());
                check("N1: reported mu is finite (last finite value, not the NaN outer result)",
                        java.util.Arrays.stream(r.mu()).allMatch(Double::isFinite));
                check("N1: diagnostics are internally consistent (lastStepAccepted == false)",
                        !rep.lastStepAccepted());
                check("N1: it stopped at the first iteration -- no later PhaseStep attempted",
                        rep.iterationsRun() == 1);
                check("N1: mass-balance report is still populated (structurally valid Result)",
                        rep.massBalance() != null
                                && rep.massBalance().targetOverall() != null);
                boolean logged = log.stream().anyMatch(l -> l.contains("non-finite")
                        || l.contains("numerical breakdown"));
                check("N1: the breakdown was logged to the progress sink", logged);
            }
        }

        // ---- R2: CancellationException from the candidate scan propagates. ----
        // Nb-Ti at x=0.5 is single-phase stable, so the active-set solve
        // converges on iteration 1 and the phase-addition candidate scan runs
        // that same iteration. Interrupting the worker thread BEFORE solve()
        // means the first CvmNewtonSolver call inside relaxCandidate observes
        // the interrupt at its top-of-loop check and throws
        // CancellationException -- which must propagate unchanged, NOT be
        // turned into "candidate could not be evaluated" and NOT let solve()
        // return normally. Deterministic: no sleeps, no timing window.
        {
            final Throwable[] caught = new Throwable[1];
            final boolean[] normalReturn = { false };
            Thread worker = new Thread(() -> {
                List<HillertSolver.Phase> ph = new ArrayList<>();
                ph.add(new HillertSolver.Phase("alpha", s, m, 1.0,
                        seededFull(m, new org.ce.model.equilibrium.CvmNewtonSolver(m), T,
                                new double[] { 0.5, 0.5 })));
                ph.add(new HillertSolver.Phase("beta", s, m, 0.0,
                        seededFull(m, new org.ce.model.equilibrium.CvmNewtonSolver(m), T,
                                new double[] { 0.15, 0.85 })));
                Thread.currentThread().interrupt();   // set BEFORE entering solve
                try {
                    HillertSolver.solve(ph, new double[] { 0.5, 0.5 }, T, 80, 20, 1.0e-8, null);
                    normalReturn[0] = true;
                } catch (Throwable t) {
                    caught[0] = t;
                }
            }, "hillert-r2-worker");
            worker.start();
            worker.join(30_000);

            System.out.printf("    R2: worker alive=%s  normalReturn=%s  caught=%s%n",
                    worker.isAlive(), normalReturn[0],
                    caught[0] == null ? "(none)" : caught[0].getClass().getName());
            check("R2: solve did not return normally after interruption", !normalReturn[0]);
            check("R2: an exception propagated out of solve",
                    caught[0] != null);
            check("R2: the propagated exception is CancellationException (cooperative cancellation preserved)",
                    caught[0] instanceof java.util.concurrent.CancellationException);
            check("R2: the worker thread finished (no deadlock)", !worker.isAlive());
        }

        // ---- E2: an EXPECTED per-grid-point failure is skipped, not swallowed
        //      as a bug. Nb-Ti at x=0.5 is single-phase stable, so the candidate
        //      scan runs, evaluates the inactive beta over the whole grid via
        //      CvmNewtonSolver, and DECLINES (dGf below threshold). The run must
        //      converge normally with no PHASE_ADDED event and no exception --
        //      i.e. the grid scan tolerated whatever non-converged / invalid
        //      grid points it hit without aborting. ----
        {
            List<HillertSolver.Phase> ph = new ArrayList<>();
            ph.add(new HillertSolver.Phase("alpha", s, m, 1.0,
                    seededFull(m, new org.ce.model.equilibrium.CvmNewtonSolver(m), T,
                            new double[] { 0.5, 0.5 })));
            ph.add(new HillertSolver.Phase("beta", s, m, 0.0,
                    seededFull(m, new org.ce.model.equilibrium.CvmNewtonSolver(m), T,
                            new double[] { 0.15, 0.85 })));
            List<String> log = new ArrayList<>();
            HillertSolver.Result r = HillertSolver.solve(ph, new double[] { 0.5, 0.5 },
                    T, 80, 20, 1.0e-8, log::add);
            boolean added = r.convergenceReport().phaseSetEvents().stream()
                    .anyMatch(e -> e.type() == HillertSolver.PhaseSetEventType.PHASE_ADDED);
            boolean evaluated = log.stream().anyMatch(l -> l.contains("candidate 'beta' evaluated"));
            System.out.printf("    E2: reason=%s  candidateEvaluated=%s  added=%s%n",
                    r.convergenceReport().reason(), evaluated, added);
            check("E2: the candidate scan completed and the run converged",
                    r.convergenceReport().reason() == HillertSolver.ConvergenceReason.CONVERGED);
            check("E2: the inactive candidate was evaluated over the grid (expected failures tolerated)",
                    evaluated);
            check("E2: no PHASE_ADDED (declined cleanly, no exception swallowed as a decline)", !added);
        }
    }

    /**
     * The overall component inventory the seed phases represent:
     * {@code target_i = sum_p N_p x^p_i}. Computed at the call site so the
     * solver is given an explicit target, never left to infer one.
     */
    private static double[] mbTarget(List<HillertSolver.Phase> phases) {
        int k = phases.get(0).numComponents;
        double[] t = new double[k];
        for (HillertSolver.Phase p : phases) {
            double[] x = p.composition();
            for (int i = 0; i < k; i++) t[i] += p.amount * x[i];
        }
        return t;
    }

    private static boolean bitwiseEqual(double[] a, double[] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) {
            if (Double.doubleToRawLongBits(a[i]) != Double.doubleToRawLongBits(b[i])) return false;
        }
        return true;
    }

    /** Full joint vector [u;x] with u = the converged single-phase CVM solution at x. */
    private static double[] seededFull(CVMGibbsModel m,
            org.ce.model.equilibrium.CvmNewtonSolver nr, double T, double[] x) {
        org.ce.model.equilibrium.CvmNewtonSolver.Result eq = nr.solve(T, x, 1e-10, null, null);
        double[] full = new double[m.ncf() + m.numComponents()];
        System.arraycopy(eq.u(), 0, full, 0, m.ncf());
        System.arraycopy(x, 0, full, m.ncf(), m.numComponents());
        return full;
    }

    // ---- analytic helpers (x = mole fraction of comp 1) -----------------

    private static double Gmix(double x, double om, double RT) {
        return RT * (x * Math.log(x) + (1 - x) * Math.log(1 - x)) + om * x * (1 - x);
    }
    private static double dGmix(double x, double om, double RT) {
        return RT * (Math.log(x) - Math.log(1 - x)) + om * (1 - 2 * x);
    }
    private static double d2Gmix(double x, double om, double RT) {
        return RT * (1.0 / x + 1.0 / (1 - x)) - 2 * om;
    }
    private static double G0(double x, double a, double b) { return a * x + b * (1 - x); }
    private static double Gabs(double x, double a, double b, double om, double RT) {
        return G0(x, a, b) + Gmix(x, om, RT);
    }
    private static double dGabs(double x, double a, double b, double om, double RT) {
        return (a - b) + dGmix(x, om, RT);
    }

    private static double[] commonTangentAbs(double aA, double bA, double omA,
            double aB, double bB, double omB, double RT) {
        double xA = 0.9, xB = 0.05;
        for (int it = 0; it < 500; it++) {
            double f1 = dGabs(xA, aA, bA, omA, RT) - dGabs(xB, aB, bB, omB, RT);
            double f2 = Gabs(xB, aB, bB, omB, RT) - Gabs(xA, aA, bA, omA, RT)
                    - dGabs(xA, aA, bA, omA, RT) * (xB - xA);
            double j11 = d2Gmix(xA, omA, RT);
            double j12 = -d2Gmix(xB, omB, RT);
            double j21 = -d2Gmix(xA, omA, RT) * (xB - xA);
            double j22 = dGabs(xB, aB, bB, omB, RT) - dGabs(xA, aA, bA, omA, RT);
            double det = j11 * j22 - j12 * j21;
            if (Math.abs(det) < 1e-40) break;
            double dxA = -(j22 * f1 - j12 * f2) / det;
            double dxB = -(-j21 * f1 + j11 * f2) / det;
            double s = 1.0;
            while (s > 1e-10) {
                double nxA = xA + s * dxA, nxB = xB + s * dxB;
                if (nxA > 1e-10 && nxA < 1 - 1e-10 && nxB > 1e-10 && nxB < 1 - 1e-10) {
                    xA = nxA; xB = nxB; break;
                }
                s *= 0.5;
            }
            if (Math.abs(f1) < 1e-9 && Math.abs(f2) < 1e-7) break;
        }
        return new double[] { xA, xB };
    }

    private static double[] commonTangentMix(double omA, double omB, double RT) {
        double xA = 0.9, xB = 0.05;
        for (int it = 0; it < 500; it++) {
            double f1 = dGmix(xA, omA, RT) - dGmix(xB, omB, RT);
            double f2 = Gmix(xB, omB, RT) - Gmix(xA, omA, RT) - dGmix(xA, omA, RT) * (xB - xA);
            double j11 = d2Gmix(xA, omA, RT);
            double j12 = -d2Gmix(xB, omB, RT);
            double j21 = -d2Gmix(xA, omA, RT) * (xB - xA);
            double j22 = dGmix(xB, omB, RT) - dGmix(xA, omA, RT);
            double det = j11 * j22 - j12 * j21;
            if (Math.abs(det) < 1e-40) break;
            double dxA = -(j22 * f1 - j12 * f2) / det;
            double dxB = -(-j21 * f1 + j11 * f2) / det;
            double s = 1.0;
            while (s > 1e-10) {
                double nxA = xA + s * dxA, nxB = xB + s * dxB;
                if (nxA > 1e-10 && nxA < 1 - 1e-10 && nxB > 1e-10 && nxB < 1 - 1e-10) {
                    xA = nxA; xB = nxB; break;
                }
                s *= 0.5;
            }
            if (Math.abs(f1) < 1e-9 && Math.abs(f2) < 1e-7) break;
        }
        return new double[] { xA, xB };
    }

    // ---- utilities ------------------------------------------------------

    private static String elementName(int i) {
        return ELEMENTS.split("-")[i];
    }

    private static double maxAbsDiff(double[] a, double[] b) {
        double m = 0;
        for (int i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
        return m;
    }

    private static void relCheck(String what, double got, double ref, double tol) {
        double rel = Math.abs((got - ref) / ref);
        checkMsg(String.format("%s  (%.6f vs %.4f, rel %.2e <= %.0e)", what, got, ref, rel, tol),
                rel <= tol);
    }

    private static void check(String what, boolean ok) {
        checkMsg(what, ok);
    }

    private static void checkMsg(String what, boolean ok) {
        if (ok) {
            System.out.println("    [ok]   " + what);
        } else {
            failures++;
            System.out.println("    [FAIL] " + what);
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

    private HillertAbsoluteMuValidation() {
    }
}
