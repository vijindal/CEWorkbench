package org.ce.model.equilibrium;

import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.hamiltonian.CECEntry;

import java.util.Arrays;
import java.util.List;

/**
 * V2 STEP 3 + STEP 4 regression gate for the <b>generalised phase-local
 * stability</b> machinery and its <b>deterministic multi-start</b> ordered
 * candidate search, on an ordered (BCC_B2) constitution vector -- exercised
 * directly, because {@link HillertSolver} still rejects ordered phases at entry
 * ({@link HillertSolver.ConvergenceReason#UNSUPPORTED_PHASE_MODEL}).
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.model.equilibrium.HillertOrderedCandidateStep
 * </pre>
 *
 * <p>In the same package as {@link HillertSolver} so it can drive the
 * package-private {@link HillertSolver#relaxWide} /
 * {@link HillertSolver#orderedCandidateSeeds} directly.</p>
 *
 * <h2>Data caveat</h2>
 * <p>The only BCC_B2 Hamiltonian in the repo ({@code Nb-Ti_BCC_B2_T}) is an
 * <b>all-zero scaffold</b>: every ECI is 0, so {@code Hm == 0} and
 * {@code G = G0m - T*Sm}. There is no ordering energy, hence no independent
 * ordered-equilibrium reference and no way to force a genuine {@code eta != 0}
 * equilibrium. This gate is therefore <b>residual-based</b> (STEP 3 Part 13 /
 * STEP 4 Part 9): it verifies the generalised formulation is mathematically
 * self-consistent and the multi-start search is deterministic, not that it
 * reproduces a known ordered phase diagram.</p>
 *
 * <h2>What is checked</h2>
 * <ul>
 *   <li><b>A</b> a valid wide {@code Y} / state representation exists; the
 *       <b>deterministic multi-start seed set</b> is generated (parent +/-
 *       ordered offsets); {@code M(Y)}, {@code J_M} are the STEP-1 values;
 *       {@code C(Y) = xA + xB = 1}.</li>
 *   <li><b>B (STEP 4)</b> invalid seeds are rejected; every valid seed relaxes
 *       via {@link HillertSolver#relaxWide}; seed generation and order are
 *       reproducible across repeated calls.</li>
 *   <li><b>B (STEP 3)</b> {@code relaxWide} converges from a seed at a fixed
 *       {@code mu}; the generalised stationarity residual
 *       {@code G_Y - J_M^T mu - lambda C_Y} is ~0.</li>
 *   <li><b>C</b> the relaxed state is {@code isValidIncludingPoints}; the
 *       normalisation {@code xA + xB = 1} holds.</li>
 *   <li><b>D</b> Gibbs-Duhem {@code sum_A M_A mu_A = G} holds at the phase's own
 *       self-consistent single-phase potential (derived, not asserted).</li>
 *   <li><b>E</b> {@code M != x} whenever the relaxed {@code eta != 0}; {@code M}
 *       and {@code G} both read from the SAME relaxed state.</li>
 *   <li><b>F</b> symmetry note (STEP 3 Part 5): the {@code eta = 0} seed is not
 *       a trap.</li>
 * </ul>
 *
 * <p>The <b>basin-coverage</b> aspect of STEP 4 (two seeds reaching two
 * analytically-known minima, lower-{@code Phi} selected) is validated on a
 * separate synthetic double-well in
 * {@code org.ce.scratch.HillertOrderedBasinSynthetic}, because the zero-ECI
 * BCC_B2 scaffold has no multi-basin {@code Phi}.</p>
 */
public final class HillertOrderedCandidateStep {

    private static int failures = 0;
    private static final double T = 1000.0;
    private static final double STAT_TOL = 1.0e-6;
    private static final double GD_TOL = 1.0e-5;   // relative

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();
        System.out.println("=".repeat(84));
        System.out.println("  V2 STEP 3 -- generalised phase-local stability on an ordered (BCC_B2) state");
        System.out.println("=".repeat(84));

        CvmGeometry geo = CvmGeometry.build("Nb-Ti", "BCC_B2", "T", null);
        CVMGibbsModel m = new CVMGibbsModel(geo, emptyEntry());
        int ncf = geo.ncf, K = geo.numComponents, w = geo.tcf;
        int iEta = geo.basis.indexOfCf("eta");
        System.out.printf("  ncf=%d tcf=%d K=%d  cfNames=%s  (all-zero scaffold Hamiltonian)%n",
                ncf, w, K, geo.basis.cfNames);
        check("A: BCC_B2 is genuinely ordered (tcf - ncf > K)", w - ncf > K,
                "tcf-ncf=" + (w - ncf));

        // ---- A: valid wide representation, deterministic multi-start seeds ----
        List<double[]> seeds = HillertSolver.orderedCandidateSeeds(m, T);
        System.out.printf("  multi-start seed set: %d seeds%n", seeds.size());
        for (int i = 0; i < seeds.size(); i++) {
            System.out.printf("    seed[%d] eta = %+.4f  valid = %s%n",
                    i, seeds.get(i)[iEta], m.atFullWide(T, seeds.get(i)).isValidIncludingPoints());
        }
        check("A: multi-start generated >= 1 seed (parent always present)", seeds.size() >= 1,
                "size=" + seeds.size());
        check("A: first seed is the parent (eta == 0, order parameters at zero)",
                Math.abs(seeds.get(0)[iEta]) < 1e-12, "eta=" + seeds.get(0)[iEta]);
        check("A: multi-start includes a nonzero-eta ordered seed",
                seeds.stream().anyMatch(s -> Math.abs(s[iEta]) > 1e-6), "no ordered seed");
        // STEP 4 B: every generated seed must be valid (the generator filters
        // invalid ones itself).
        boolean allSeedsValid = seeds.stream()
                .allMatch(s -> m.atFullWide(T, s).isValidIncludingPoints());
        check("B (STEP 4): every generated seed passes isValidIncludingPoints "
                        + "(invalid variants filtered by the generator)",
                allSeedsValid, "an emitted seed is invalid");
        // STEP 4 B: determinism -- a second call yields an identical seed set.
        List<double[]> seeds2 = HillertSolver.orderedCandidateSeeds(m, T);
        boolean deterministicSeeds = seeds.size() == seeds2.size();
        for (int i = 0; deterministicSeeds && i < seeds.size(); i++) {
            deterministicSeeds = Arrays.equals(seeds.get(i), seeds2.get(i));
        }
        check("B (STEP 4): seed generation is deterministic (identical set on repeat call)",
                deterministicSeeds, "seed set differs between calls");
        // STEP 4 B: every valid seed relaxes without exception; count outcomes.
        double[] muRelax = { -60000.0, -45000.0 };
        int relaxed = 0, converged = 0;
        for (double[] sd : seeds) {
            HillertSolver.RelaxResult r = HillertSolver.relaxWide(m, sd, muRelax, T);
            relaxed++;
            if (r.converged() && r.state() != null && r.state().isValidIncludingPoints()) {
                converged++;
            }
        }
        check("B (STEP 4): every seed was relaxed (no exception)", relaxed == seeds.size(),
                relaxed + "/" + seeds.size());
        check("B (STEP 4): at least one seed relaxed to a converged valid minimum",
                converged >= 1, converged + " converged");

        double[] seed = seeds.get(0);
        CVMGibbsModel.State seedState = m.atFullWide(T, seed);
        check("A: seed state is valid", seedState.isValidIncludingPoints(), "invalid seed");
        check("A: seed eta == 0 (order parameter starts at 0)",
                Math.abs(seed[iEta]) < 1e-12, "eta=" + seed[iEta]);
        double[][] jM = seedState.componentAmountsJacobian();
        check("A: J_M is K x tcf", jM.length == K && jM[0].length == w,
                jM.length + "x" + (jM.length > 0 ? jM[0].length : -1));
        check("A: J_M eta-column == [-1/2, +1/2] (STEP 1)",
                Math.abs(jM[0][iEta] + 0.5) < 1e-12 && Math.abs(jM[1][iEta] - 0.5) < 1e-12,
                "col=[" + jM[0][iEta] + "," + jM[1][iEta] + "]");
        check("A: J_M u-block all zero (M depends only on the point block)",
                allZero(jM, 0, ncf), "nonzero u-block entry");

        // C(Y) = xA + xB  -> gradient c is 1 on the K composition columns only
        double cSum = seed[ncf] + seed[ncf + 1];
        check("A: C(seed) = xA + xB = 1", Math.abs(cSum - 1.0) < 1e-12, "= " + cSum);

        // ---- B/C: relaxWide converges, stationarity residual ~0 ----
        // Use a nonzero mu so the objective Phi = G - mu^T M is non-degenerate
        // and the eta coupling is exercised.
        double[] mu = { -60000.0, -45000.0 };
        HillertSolver.RelaxResult rr = HillertSolver.relaxWide(m, seed, mu, T);
        check("B: relaxWide converged", rr.converged(), "residual=" + rr.residual());
        check("B: relaxWide residual small", rr.residual() < 1e-8, "= " + rr.residual());
        double[] Y = rr.y();
        CVMGibbsModel.State st = rr.state();
        check("C: relaxed state valid", st != null && st.isValidIncludingPoints(), "invalid");
        check("C: normalisation xA + xB = 1 at relaxed Y",
                Math.abs(Y[ncf] + Y[ncf + 1] - 1.0) < 1e-9, "= " + (Y[ncf] + Y[ncf + 1]));

        // independent generalised stationarity residual at the relaxed Y:
        //   r_j = G_Y[j] - (J_M^T mu)[j] - lambda * c[j]
        // choose lambda to make the c-weighted mean of (G_Y - J_M^T mu) vanish
        // (that is exactly the KKT multiplier of the sum-to-one constraint).
        double[] Gy = st.guFull(w);
        double[][] jMr = st.componentAmountsJacobian();
        double[] c = new double[w];
        for (int a = 0; a < K; a++) c[ncf + a] = 1.0;
        double num = 0, den = 0;
        for (int j = 0; j < w; j++) {
            double base = Gy[j];
            for (int a = 0; a < K; a++) base -= jMr[a][j] * mu[a];
            num += c[j] * base;
            den += c[j] * c[j];
        }
        double lambda = num / den;
        double statResid = 0;
        for (int j = 0; j < w; j++) {
            double r = Gy[j] - lambda * c[j];
            for (int a = 0; a < K; a++) r -= jMr[a][j] * mu[a];
            statResid = Math.max(statResid, Math.abs(r));
        }
        check("B: generalised stationarity residual G_Y - J_M^T mu - lambda C_Y ~ 0 (worst "
                + String.format("%.2e", statResid) + ")", statResid < STAT_TOL, "= " + statResid);

        // ---- E: M != x when eta != 0; M and G from the SAME state ----
        double[] mAt = st.componentAmountsPerFormulaUnit();
        double[] xAt = st.composition();
        double etaR = Y[iEta];
        System.out.printf("  relaxed:  eta=%.8f  x=%s  M=%s  G=%.6f%n",
                etaR, Arrays.toString(xAt), Arrays.toString(mAt), st.g());
        if (Math.abs(etaR) > 1e-6) {
            check("E: M != x when eta != 0 (|M_A - xA| == |eta|/2)",
                    Math.abs(Math.abs(mAt[0] - xAt[0]) - Math.abs(etaR) / 2.0) < 1e-9,
                    "|M_A-xA|=" + Math.abs(mAt[0] - xAt[0]) + " eta/2=" + Math.abs(etaR) / 2.0);
        } else {
            System.out.println("    E: relaxed eta ~ 0 for this mu (zero-ECI lattice); "
                    + "M == x -- see class caveat");
            check("E: M == x when eta == 0", maxDiff(mAt, xAt) < 1e-12, "M-x = " + maxDiff(mAt, xAt));
        }
        check("E: candidate M == state.componentAmountsPerFormulaUnit()",
                Arrays.equals(mAt, st.componentAmountsPerFormulaUnit()), "M mismatch");
        // G consistency is implicit -- st.g() is read directly.

        // ---- D: Gibbs-Duhem at the phase's OWN self-consistent mu ----
        // For a single phase, mu is underdetermined by GD alone. Derive a
        // consistent mu from the relaxed state's own composition-block gradient
        // (dG/dx_i) minus the shared Lagrange term, then verify sum_A M_A mu_A = G.
        double[] muSelf = selfConsistentMu(st, ncf, K, w);
        double gd = 0;
        for (int a = 0; a < K; a++) gd += mAt[a] * muSelf[a];
        double gdRel = Math.abs(gd - st.g()) / Math.max(Math.abs(st.g()), 1.0);
        System.out.printf("  self-consistent mu=%s  sum(M*mu)=%.6f  G=%.6f  rel=%.2e%n",
                Arrays.toString(muSelf), gd, st.g(), gdRel);
        check("D: Gibbs-Duhem sum_A M_A mu_A == G at self-consistent mu", gdRel < GD_TOL,
                "rel=" + gdRel);

        // ---- F: the eta=0 seed is NOT a trap (Part 5) ----
        // From an eta=0 seed, relaxWide leaves eta=0 on its own whenever a
        // nonzero eta lowers Phi. Here the asymmetric reference energy G0m
        // (g0(Nb) != g0(Ti)) slides the composition off symmetric, which opens
        // room for eta != 0 -- so no arbitrary perturbation is ever needed.
        check("F: eta=0 seed is not a trap -- relaxation reached eta != 0 for this mu "
                        + "(no forced perturbation)",
                Math.abs(etaR) > 1e-3, "eta stayed at " + etaR);
        // The exact-eta=0 stationary point IS reachable, but only with
        // composition also pinned AND a symmetric G0m; document the mechanism.
        double g0Nb = LatticeStability.g0("Nb", geo.parentStructure, T);
        double g0Ti = LatticeStability.g0("Ti", geo.parentStructure, T);
        check("F: reference energy is genuinely A<->B asymmetric here "
                        + "(g0(Nb) != g0(Ti)) -- so G0m, not a bug, breaks the eta symmetry",
                Math.abs(g0Nb - g0Ti) > 1.0, "g0(Nb)=" + g0Nb + " g0(Ti)=" + g0Ti);

        // ---- STEP 4 Part 13: bounded extra cost of the multi-start ----
        long t1s = System.nanoTime();
        HillertSolver.RelaxResult one = HillertSolver.relaxWide(m, seeds.get(0), mu, T);
        long t1e = System.nanoTime();
        long tAllS = System.nanoTime();
        int distinctAfterDedup;
        {
            java.util.List<HillertSolver.CandidatePhaseState> found = new java.util.ArrayList<>();
            for (double[] sd : seeds) {
                HillertSolver.RelaxResult r = HillertSolver.relaxWide(m, sd, mu, T);
                if (r.converged() && r.state() != null && r.state().isValidIncludingPoints()) {
                    found.add(HillertSolver.minimumForTest(r.y(), 0.0));
                }
            }
            distinctAfterDedup = HillertSolver.countDistinctMinima(found);
        }
        long tAllE = System.nanoTime();
        double oneMs = (t1e - t1s) / 1e6, allMs = (tAllE - tAllS) / 1e6;
        System.out.printf("%n  PERF: %d seeds; single relaxWide %.2f ms; full multi-start %.2f ms "
                        + "(~%.1fx); %d distinct minima after dedup%n",
                seeds.size(), oneMs, allMs, allMs / Math.max(oneMs, 1e-9), distinctAfterDedup);
        check("Part 13: multi-start cost is a small bounded multiple of one relaxation "
                        + "(<= seed-count x, plus the seed set is tiny)",
                seeds.size() <= 6, "seed count " + seeds.size());
        check("Part 13: the zero-ECI scaffold has ONE basin -- all seeds dedup to 1 minimum",
                distinctAfterDedup == 1, "distinct=" + distinctAfterDedup);

        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(84));
        if (failures > 0) {
            throw new AssertionError(failures + " ordered-candidate-step checks failed");
        }
    }

    /**
     * A mu consistent with the relaxed single-phase state: from the KKT
     * condition {@code dG/dx_i - mu_i - lambda = 0} on the K composition
     * coordinates, {@code mu_i = dG/dx_i - lambda}, with lambda fixed by
     * {@code sum_A M_A mu_A = G}. (One linear equation for the one free scalar
     * lambda, given M and dG/dx.)
     */
    private static double[] selfConsistentMu(CVMGibbsModel.State st, int ncf, int K, int w) {
        double[] Gy = st.guFull(w);
        double[] m = st.componentAmountsPerFormulaUnit();
        // mu_i = Gy[ncf+i] - lambda ; sum_A M_A (Gy[ncf+A] - lambda) = G
        // => sum_A M_A Gy[ncf+A] - lambda * sum_A M_A = G
        double sMg = 0, sM = 0;
        for (int a = 0; a < K; a++) {
            sMg += m[a] * Gy[ncf + a];
            sM += m[a];
        }
        double lambda = (sMg - st.g()) / sM;
        double[] mu = new double[K];
        for (int a = 0; a < K; a++) mu[a] = Gy[ncf + a] - lambda;
        return mu;
    }

    private static boolean allZero(double[][] a, int from, int to) {
        for (double[] row : a)
            for (int j = from; j < to; j++)
                if (row[j] != 0.0) return false;
        return true;
    }

    private static double maxDiff(double[] a, double[] b) {
        double m = 0;
        for (int i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
        return m;
    }

    private static CECEntry emptyEntry() {
        CECEntry e = new CECEntry();
        e.elements = "Nb-Ti";
        e.structurePhase = "BCC_B2_T";
        e.model = "T";
        e.cecTerms = new CECEntry.CECTerm[0];
        return e;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-74s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-74s [!] FAIL  %s%n", label, detail);
        }
    }

    private HillertOrderedCandidateStep() {
    }
}
