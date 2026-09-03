package org.ce.scratch;

import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.hamiltonian.CECEntry;

import java.util.Arrays;

/**
 * V2 STEP 1 regression gate for the {@code M_A(Y)} thermodynamic contract:
 * {@link CVMGibbsModel.State#componentAmountsPerFormulaUnit()} and
 * {@link CVMGibbsModel.State#componentAmountsJacobian()}.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.HillertMAIdentity
 * </pre>
 *
 * <p>{@code composition()} is the phase's composition <em>coordinate</em>
 * {@code x_A}; {@code componentAmountsPerFormulaUnit()} is the conserved
 * quantity {@code M_A} = moles of component A per formula unit. For a
 * single-site one-atom-per-site disordered phase they are numerically
 * identical; for BCC_B2 (an ordered phase carrying the LRO parameter
 * {@code eta}) {@code M_A = xA - eta/2 != xA}.</p>
 *
 * <p><b>HillertSolver is deliberately not used</b> -- it still rejects BCC_B2 at
 * entry, and this gate must exercise the model contract directly.</p>
 *
 * <h2>Cases</h2>
 * <ol>
 *   <li><b>A</b> BCC_A2 (K=2,3,4): {@code M == x} exactly at random and skewed
 *       non-random valid states; Jacobian is exactly the selector matrix
 *       {@code e_{ncf+A}}.</li>
 *   <li><b>B</b> FCC_A1 (K=2,3): {@code M == x}; selector Jacobian.</li>
 *   <li><b>C</b> HCP_A3 (K=2,3): {@code M == x}; selector Jacobian.</li>
 *   <li><b>D</b> BCC_B2 (K=2): a valid ordered state with {@code eta != 0};
 *       {@code M_A = xA - eta/2}, {@code M_B = xB + eta/2}, independently
 *       computed; {@code M != x}.</li>
 *   <li><b>E</b> analytic Jacobian vs central finite differences of
 *       {@code componentAmountsPerFormulaUnit}, BCC_A2 and BCC_B2, at more than
 *       one state (the exact derivative is structurally constant); max abs
 *       discrepancy {@code < 1e-7}.</li>
 *   <li><b>F</b> immutability: mutating a returned {@code M} or Jacobian does
 *       not change a subsequently fetched copy or the state.</li>
 *   <li><b>Invariants</b> for every case: {@code M} finite, {@code M_A >= 0} for
 *       a valid state, {@code sum_A M_A == 1} (no silent renormalisation).</li>
 * </ol>
 */
public final class HillertMAIdentity {

    private static int failures = 0;
    private static final double T = 1000.0;
    private static final double EXACT = 0.0;         // M == x must be bit-exact
    private static final double FD_TOL = 1.0e-7;
    private static final double SUM_TOL = 1.0e-9;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(84));
        System.out.println("  V2 STEP 1 -- M_A(Y) thermodynamic contract");
        System.out.println("=".repeat(84));

        // M_A is a pure geometry quantity -- ECI-independent -- so every model
        // here is built straight from CvmGeometry + an empty CECEntry, with no
        // dependence on a stored Hamiltonian (FCC_A1 / HCP_A3 have none stored).

        // ---- A: BCC_A2 ----
        disorderedIdentity("Nb-Ti", "BCC_A2", "T", new double[][] {
                { 0.5, 0.5 }, { 0.2, 0.8 }, { 0.05, 0.95 } });
        disorderedIdentity("Nb-Ti-V", "BCC_A2", "T", new double[][] {
                { 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 0.2, 0.3, 0.5 }, { 0.1, 0.1, 0.8 } });
        disorderedIdentity("Nb-Ti-V-Zr", "BCC_A2", "T", new double[][] {
                { 0.25, 0.25, 0.25, 0.25 }, { 0.1, 0.2, 0.3, 0.4 } });

        // ---- B: FCC_A1 ----
        disorderedIdentity("Nb-Ti", "FCC_A1", "T", new double[][] {
                { 0.5, 0.5 }, { 0.3, 0.7 } });
        disorderedIdentity("Nb-Ti-V", "FCC_A1", "T", new double[][] {
                { 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 0.2, 0.3, 0.5 } });

        // ---- C: HCP_A3 ----
        disorderedIdentity("Nb-Ti", "HCP_A3", "T", new double[][] {
                { 0.5, 0.5 }, { 0.25, 0.75 } });
        disorderedIdentity("Nb-Ti-V", "HCP_A3", "T", new double[][] {
                { 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 0.2, 0.3, 0.5 } });

        // ---- D + E(B2) + F: BCC_B2 ----
        orderedBccB2();

        // ---- E(A2): finite-difference Jacobian for a disordered phase ----
        finiteDifferenceDisordered("Nb-Ti-V", "BCC_A2", "T");

        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(84));
        if (failures > 0) {
            throw new AssertionError(failures + " M_A contract checks failed");
        }
    }

    // ------------------------------------------------------------------
    // A / B / C : disordered phases -- M == x, selector Jacobian
    // ------------------------------------------------------------------

    private static void disorderedIdentity(
            String elements, String structure, String model, double[][] comps) {

        System.out.printf("%n--- %s / %s / %s ---%n", elements, structure, model);
        CvmGeometry geo = CvmGeometry.build(elements, structure, model, null);
        CVMGibbsModel m = new CVMGibbsModel(geo, emptyEntry(elements, structure + "_" + model));
        int K = geo.numComponents;
        int ncf = geo.ncf;
        int tcf = geo.tcf;

        // Jacobian is state-independent; check it once against the selector.
        double[][] jac = geo.componentAmountsJacobian();
        boolean selectorOk = jac.length == K;
        for (int a = 0; a < K && selectorOk; a++) {
            if (jac[a].length != tcf) { selectorOk = false; break; }
            for (int j = 0; j < tcf; j++) {
                double expect = (j == ncf + a) ? 1.0 : 0.0;
                if (jac[a][j] != expect) { selectorOk = false; break; }
            }
        }
        check(structure + " K=" + K + " Jacobian is exactly the selector e_{ncf+A}", selectorOk,
                "jac=" + Arrays.deepToString(jac));
        check(structure + " K=" + K + " tcf == ncf + K (disordered)", tcf == ncf + K,
                "tcf=" + tcf + " ncf=" + ncf + " K=" + K);

        for (double[] x : comps) {
            // random state at x, plus a skewed non-random state (perturb u).
            double[] uRand = m.randomStateU(x);
            checkOne(structure, "random", m, geo, T, x, uRand);

            double[] uSkew = uRand.clone();
            for (int i = 0; i < uSkew.length; i++) {
                uSkew[i] *= (i % 2 == 0) ? 0.95 : 1.05;
            }
            // only use the skewed state if it is still physically valid
            CVMGibbsModel.State st = m.at(T, x, uSkew);
            if (st.isValidIncludingPoints()) {
                checkOne(structure, "skewed", m, geo, T, x, uSkew);
            }
        }
    }

    private static void checkOne(String tag, String kind, CVMGibbsModel m, CvmGeometry geo,
            double t, double[] x, double[] u) {
        CVMGibbsModel.State st = m.at(t, x, u);
        double[] mA = st.componentAmountsPerFormulaUnit();
        double[] comp = st.composition();

        double worst = 0.0;
        for (int a = 0; a < x.length; a++) {
            worst = Math.max(worst, Math.abs(mA[a] - comp[a]));
        }
        check(tag + " (" + kind + " x=" + fmt(x) + "): M == composition() exactly",
                worst == EXACT, "worst |M-x| = " + worst);
        invariants(tag + " " + kind, mA, st.isValidIncludingPoints());
    }

    // ------------------------------------------------------------------
    // D : BCC_B2 ordered phase -- M_A = xA - eta/2
    // ------------------------------------------------------------------

    private static void orderedBccB2() {
        System.out.printf("%n--- Nb-Ti / BCC_B2 / T  (ordered: M_A != x_A) ---%n");
        CvmGeometry geo = CvmGeometry.build("Nb-Ti", "BCC_B2", "T", null);
        int K = geo.numComponents;
        int ncf = geo.ncf;
        int tcf = geo.tcf;
        System.out.printf("    K=%d ncf=%d tcf=%d  cfNames=%s%n", K, ncf, tcf, geo.basis.cfNames);
        check("BCC_B2 is genuinely ordered (tcf - ncf > K)", tcf - ncf > K,
                "tcf-ncf=" + (tcf - ncf) + " K=" + K);

        CECEntry empty = emptyEntry("Nb-Ti", "BCC_B2_T");
        CVMGibbsModel m = new CVMGibbsModel(geo, empty);

        // Build a valid ordered state directly: start from the random state at a
        // skewed composition (basis fills the full width incl. eta), then push
        // eta away from 0 by a hand construction that stays inside (0,1) on
        // every point occupation.
        double xA = 0.4, xB = 0.6;
        double[] uRand = m.randomStateFull(new double[] { xA, xB }); // width tcf, eta ~ 0
        int idxEta = geo.basis.indexOfCf("eta");
        int idxXA = geo.basis.indexOfCf("xA");
        int idxXB = geo.basis.indexOfCf("xB");
        check("BCC_B2 basis exposes xA, xB, eta columns",
                idxEta == ncf + 2 && idxXA == ncf && idxXB == ncf + 1,
                "idxXA=" + idxXA + " idxXB=" + idxXB + " idxEta=" + idxEta + " ncf=" + ncf);

        for (double eta : new double[] { 0.10, -0.15, 0.20 }) {
            double[] uFull = uRand.clone();
            uFull[idxXA] = xA;
            uFull[idxXB] = xB;
            uFull[idxEta] = eta;

            CVMGibbsModel.State st = m.atFullWide(T, uFull);

            double[] mA = st.componentAmountsPerFormulaUnit();
            double[] comp = st.composition(); // [xA, xB]

            double expA = xA - eta / 2.0;
            double expB = xB + eta / 2.0;
            double dA = Math.abs(mA[0] - expA);
            double dB = Math.abs(mA[1] - expB);
            check(String.format("BCC_B2 eta=%+.2f : M_A = xA - eta/2  (got %.10f, expect %.10f)",
                    eta, mA[0], expA), dA < 1e-12, "|delta| = " + dA);
            check(String.format("BCC_B2 eta=%+.2f : M_B = xB + eta/2  (got %.10f, expect %.10f)",
                    eta, mA[1], expB), dB < 1e-12, "|delta| = " + dB);

            double devA = Math.abs(mA[0] - comp[0]);
            check(String.format("BCC_B2 eta=%+.2f : M != composition() (|M_A - xA| = %.4f = |eta|/2)",
                    eta, devA), Math.abs(devA - Math.abs(eta) / 2.0) < 1e-12,
                    "devA=" + devA);

            invariants(String.format("BCC_B2 eta=%+.2f", eta), mA, st.isValidIncludingPoints());
            check(String.format("BCC_B2 eta=%+.2f : (xA-eta/2)+(xB+eta/2) == 1", eta),
                    Math.abs((expA + expB) - 1.0) < 1e-12, "sum=" + (expA + expB));
        }

        // ---- E (BCC_B2): analytic Jacobian vs finite differences ----
        double[] uBase = uRand.clone();
        uBase[idxXA] = xA;
        uBase[idxXB] = xB;
        uBase[idxEta] = 0.12;
        finiteDifferenceCheck("BCC_B2", m, geo, uBase);
        // second state -- the exact derivative is constant, so it must still match
        double[] uBase2 = uRand.clone();
        uBase2[idxXA] = 0.55;
        uBase2[idxXB] = 0.45;
        uBase2[idxEta] = -0.05;
        finiteDifferenceCheck("BCC_B2 (2nd state)", m, geo, uBase2);

        // ---- F: immutability ----
        immutability(m.atFullWide(T, uBase), geo);
    }

    // ------------------------------------------------------------------
    // E : finite-difference Jacobian verification
    // ------------------------------------------------------------------

    private static void finiteDifferenceDisordered(
            String elements, String structure, String model) {
        System.out.printf("%n--- Finite-difference Jacobian: %s / %s / %s ---%n",
                elements, structure, model);
        CvmGeometry geo = CvmGeometry.build(elements, structure, model, null);
        CVMGibbsModel m = new CVMGibbsModel(geo, emptyEntry(elements, structure + "_" + model));
        double[] x1 = { 0.2, 0.3, 0.5 };
        double[] uFull1 = m.randomStateFull(x1);
        finiteDifferenceCheck(structure + " state1", m, geo, uFull1);
        double[] x2 = { 0.45, 0.25, 0.30 };
        double[] uFull2 = m.randomStateFull(x2);
        finiteDifferenceCheck(structure + " state2", m, geo, uFull2);
    }

    /**
     * Central finite difference of {@code componentAmountsPerFormulaUnit} in
     * each {@code uFull} component vs the analytic Jacobian. Only requires that
     * {@code atFullWide} accepts a width-{@code tcf} vector (true for every
     * registered phase).
     */
    private static void finiteDifferenceCheck(
            String tag, CVMGibbsModel m, CvmGeometry geo, double[] uFullBase) {
        int tcf = geo.tcf;
        int K = geo.numComponents;
        double[][] jac = geo.componentAmountsJacobian();
        double h = 1.0e-6;

        double worst = 0.0;
        int worstA = -1, worstJ = -1;
        for (int j = 0; j < tcf; j++) {
            double[] up = uFullBase.clone();
            double[] dn = uFullBase.clone();
            up[j] += h;
            dn[j] -= h;
            double[] mUp = geo.componentAmounts(up);
            double[] mDn = geo.componentAmounts(dn);
            for (int a = 0; a < K; a++) {
                double fd = (mUp[a] - mDn[a]) / (2.0 * h);
                double d = Math.abs(fd - jac[a][j]);
                if (d > worst) { worst = d; worstA = a; worstJ = j; }
            }
        }
        check(tag + " : analytic Jacobian matches central FD (worst |delta| = "
                        + String.format("%.2e", worst) + " at [A=" + worstA + ", j=" + worstJ + "])",
                worst < FD_TOL, "worst=" + worst);
    }

    // ------------------------------------------------------------------
    // F : immutability
    // ------------------------------------------------------------------

    private static void immutability(CVMGibbsModel.State st, CvmGeometry geo) {
        System.out.printf("%n--- Immutability ---%n");

        double[] m1 = st.componentAmountsPerFormulaUnit();
        double[] pristineM = m1.clone();
        Arrays.fill(m1, Double.NaN);
        double[] m2 = st.componentAmountsPerFormulaUnit();
        check("mutating returned M does not affect a re-fetched M",
                Arrays.equals(m2, pristineM), "m2=" + Arrays.toString(m2));
        check("componentAmountsPerFormulaUnit() returns a distinct array",
                m1 != m2, "same ref");

        double[][] j1 = st.componentAmountsJacobian();
        double[][] pristineJ = new double[j1.length][];
        for (int a = 0; a < j1.length; a++) pristineJ[a] = j1[a].clone();
        for (double[] row : j1) Arrays.fill(row, Double.NaN);
        double[][] j2 = st.componentAmountsJacobian();
        boolean jOk = j2.length == pristineJ.length;
        for (int a = 0; a < j2.length && jOk; a++) jOk = Arrays.equals(j2[a], pristineJ[a]);
        check("mutating returned Jacobian does not affect a re-fetched Jacobian", jOk,
                "j2=" + Arrays.deepToString(j2));
        check("componentAmountsJacobian() returns a distinct matrix",
                j1 != j2 && (j1.length == 0 || j1[0] != j2[0]), "same ref");

        // geometry's own stored map unchanged
        double[][] geoMap = geo.componentAmountsJacobian();
        boolean geoOk = true;
        for (double[] row : geoMap) for (double v : row) if (Double.isNaN(v)) geoOk = false;
        check("CvmGeometry.componentAmountsMap not corrupted by caller mutation", geoOk,
                "NaN leaked into geometry map");
    }

    // ------------------------------------------------------------------
    // Invariants
    // ------------------------------------------------------------------

    private static void invariants(String tag, double[] mA, boolean validState) {
        boolean finite = true;
        for (double v : mA) finite &= Double.isFinite(v);
        check(tag + " : M components finite", finite, Arrays.toString(mA));

        if (validState) {
            boolean nonNeg = true;
            for (double v : mA) nonNeg &= v >= 0.0;
            check(tag + " : M components >= 0 for a valid state", nonNeg, Arrays.toString(mA));
        }

        double sum = 0.0;
        for (double v : mA) sum += v;
        check(tag + " : sum_A M_A == 1 (no renormalisation)",
                Math.abs(sum - 1.0) < SUM_TOL, "sum=" + sum);
    }

    // ------------------------------------------------------------------

    private static CECEntry emptyEntry(String elements, String structurePhase) {
        CECEntry e = new CECEntry();
        e.elements = elements;
        e.structurePhase = structurePhase;
        e.model = "T";
        e.cecTerms = new CECEntry.CECTerm[0];
        return e;
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-78s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-78s [!] FAIL  %s%n", label, detail);
        }
    }

    private static String fmt(double[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%.3f", v[i]));
            if (i < v.length - 1) sb.append(",");
        }
        return sb.append(']').toString();
    }

    private HillertMAIdentity() {
    }
}
