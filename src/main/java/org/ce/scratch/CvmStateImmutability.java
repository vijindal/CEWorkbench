package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CvmGeometry;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

/**
 * Regression gate for {@link CVMGibbsModel.State}'s defensive-copy contract.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.CvmStateImmutability
 * </pre>
 *
 * <p>{@code State} advertises immutability. {@link CVMGibbsModel.State#clusterVariables()}
 * used to return the internal {@code cv} 3-D array by reference, so a caller
 * could rewrite a state's cluster probabilities in place -- and the only
 * production caller ({@code CvmNewtonSolver}'s step-3 degeneracy check) happens
 * only to read it, so nothing caught the hole. This checks that a mutation of
 * any nesting level of the returned array leaves a subsequently fetched copy at
 * the original thermodynamic values, and re-confirms the already-defensive
 * {@code composition()}, {@code u()} and {@code eci()}.</p>
 */
public final class CvmStateImmutability {

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  CVMGibbsModel.State defensive-copy contract");
        System.out.println("=".repeat(78));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Ti-V", "BCC_A2", "T"), EngineConfig.CVM, null);

        CvmGeometry geo = CvmGeometry.build("Nb-Ti-V", "BCC_A2", "T", null);
        CVMGibbsModel model = new CVMGibbsModel(geo, session.cecEntry);

        double[] x = { 0.2, 0.3, 0.5 };
        double[] u = model.randomStateU(x);
        CVMGibbsModel.State state = model.at(1000.0, x, u);

        checkClusterVariables(state);
        checkFlatAccessor("composition()", state.composition(), state::composition);
        checkFlatAccessor("u()", state.u(), state::u);
        checkFlatAccessor("eci()", state.eci(), state::eci);

        // A property computed from cv must be unaffected by a caller's mutation
        // of a clusterVariables() result -- proves the deep copy really is
        // detached from the field the entropy sum reads.
        double smBefore = state.sm();
        double gmBefore = state.gm();
        double[][][] scribble = state.clusterVariables();
        for (double[][] tt : scribble) {
            if (tt == null) continue;
            for (double[] jj : tt) {
                if (jj == null) continue;
                java.util.Arrays.fill(jj, Double.NaN);
            }
        }
        check("sm() unchanged after scribbling on a clusterVariables() copy",
                state.sm() == smBefore, smBefore + " vs " + state.sm());
        check("gm() unchanged after scribbling on a clusterVariables() copy",
                state.gm() == gmBefore, gmBefore + " vs " + state.gm());

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " immutability checks failed");
        }
    }

    private static void checkClusterVariables(CVMGibbsModel.State state) {
        System.out.println("\n--- clusterVariables() ---");

        double[][][] cv1 = state.clusterVariables();

        // Snapshot the pristine values by value.
        double[][][] pristine = new double[cv1.length][][];
        for (int t = 0; t < cv1.length; t++) {
            if (cv1[t] == null) continue;
            pristine[t] = new double[cv1[t].length][];
            for (int j = 0; j < cv1[t].length; j++) {
                if (cv1[t][j] != null) pristine[t][j] = cv1[t][j].clone();
            }
        }

        // Distinct identity: a second fetch must not be the same array object
        // at any level.
        double[][][] cvA = state.clusterVariables();
        double[][][] cvB = state.clusterVariables();
        boolean sharedOuter = (cvA == cvB);
        boolean sharedInner = false;
        boolean sharedLeaf = false;
        for (int t = 0; t < cvA.length; t++) {
            if (cvA[t] == null) continue;
            if (cvA[t] == cvB[t]) sharedInner = true;
            for (int j = 0; j < cvA[t].length; j++) {
                if (cvA[t][j] != null && cvA[t][j] == cvB[t][j]) sharedLeaf = true;
            }
        }
        check("two fetches share no outer array", !sharedOuter, "same ref");
        check("two fetches share no mid-level array", !sharedInner, "same ref");
        check("two fetches share no leaf array", !sharedLeaf, "same ref");

        // Mutate every level of cv1.
        int mutated = 0;
        for (int t = 0; t < cv1.length; t++) {
            if (cv1[t] == null) continue;
            for (int j = 0; j < cv1[t].length; j++) {
                if (cv1[t][j] == null) continue;
                for (int v = 0; v < cv1[t][j].length; v++) {
                    cv1[t][j][v] = -999.0 - t - j - v;
                    mutated++;
                }
                // Also replace a whole row reference.
                cv1[t][j] = new double[] { 123.0 };
            }
            // And a whole block reference.
            cv1[t] = new double[][] { { 42.0 } };
        }
        check("test actually mutated some entries", mutated > 0, "mutated=" + mutated);

        // Fresh fetch must equal the pristine snapshot exactly.
        double[][][] cv2 = state.clusterVariables();
        double worst = 0;
        String at = "-";
        for (int t = 0; t < pristine.length; t++) {
            if (pristine[t] == null) continue;
            if (cv2[t] == null || cv2[t].length != pristine[t].length) {
                check("cv2 shape at t=" + t, false, "row-count mismatch");
                continue;
            }
            for (int j = 0; j < pristine[t].length; j++) {
                if (pristine[t][j] == null) continue;
                if (cv2[t][j] == null || cv2[t][j].length != pristine[t][j].length) {
                    check("cv2 shape at t=" + t + " j=" + j, false, "len mismatch");
                    continue;
                }
                for (int v = 0; v < pristine[t][j].length; v++) {
                    double d = Math.abs(cv2[t][j][v] - pristine[t][j][v]);
                    if (d > worst) { worst = d; at = "t=" + t + " j=" + j + " v=" + v; }
                }
            }
        }
        check("cv2 identical to pristine values (worst |delta| " + worst + " at " + at + ")",
                worst == 0.0, "worst=" + worst);

        // Independence within a copy: mutating one entry of a fresh copy must
        // not touch a sibling entry, a sibling row, or a sibling block of that
        // same copy.
        double[][][] cv3 = state.clusterVariables();
        int[] multiValLeaf = findMultiValueLeaf(cv3);
        if (multiValLeaf != null) {
            int t = multiValLeaf[0], j = multiValLeaf[1];
            double sibling1 = cv3[t][j][1];
            cv3[t][j][0] = 7777.0;
            check("mutating one entry does not affect a sibling entry in the same row",
                    cv3[t][j][1] == sibling1, sibling1 + " vs " + cv3[t][j][1]);
        } else {
            System.out.println("    (no leaf row with two entries; sibling-entry check skipped)");
        }

        int[] twoLeaves = findTwoLeaves(cv3);
        if (twoLeaves != null) {
            int t = twoLeaves[0], j1 = twoLeaves[1], j2 = twoLeaves[2];
            double sibling0 = cv3[t][j2][0];
            cv3[t][j1] = new double[] { -1.0 };
            check("replacing one row does not affect a sibling row",
                    cv3[t][j2][0] == sibling0, sibling0 + " vs " + cv3[t][j2][0]);
        } else {
            System.out.println("    (no cluster type with two rows; sibling-row check skipped)");
        }

        int[] twoBlocks = findTwoBlocks(cv3);
        if (twoBlocks != null) {
            int t1 = twoBlocks[0], t2 = twoBlocks[1];
            double other0 = cv3[t2][0][0];
            cv3[t1] = new double[][] { { -2.0 } };
            check("replacing one block does not affect a sibling block",
                    cv3[t2][0][0] == other0, other0 + " vs " + cv3[t2][0][0]);
        }
    }

    /** Finds a (t, j) whose leaf row has at least two entries. */
    private static int[] findMultiValueLeaf(double[][][] cv) {
        for (int t = 0; t < cv.length; t++) {
            if (cv[t] == null) continue;
            for (int j = 0; j < cv[t].length; j++) {
                if (cv[t][j] != null && cv[t][j].length >= 2) return new int[] { t, j };
            }
        }
        return null;
    }

    /** Finds a (t, j1, j2) with two non-empty leaf rows under one cluster type. */
    private static int[] findTwoLeaves(double[][][] cv) {
        for (int t = 0; t < cv.length; t++) {
            if (cv[t] == null) continue;
            int firstJ = -1;
            for (int j = 0; j < cv[t].length; j++) {
                if (cv[t][j] != null && cv[t][j].length > 0) {
                    if (firstJ < 0) firstJ = j;
                    else return new int[] { t, firstJ, j };
                }
            }
        }
        return null;
    }

    /** Finds two distinct cluster types, each with a non-empty first row. */
    private static int[] findTwoBlocks(double[][][] cv) {
        int first = -1;
        for (int t = 0; t < cv.length; t++) {
            if (cv[t] != null && cv[t].length > 0 && cv[t][0] != null && cv[t][0].length > 0) {
                if (first < 0) first = t;
                else return new int[] { first, t };
            }
        }
        return null;
    }

    private interface Fetch { double[] get(); }

    private static void checkFlatAccessor(String label, double[] first, Fetch refetch) {
        System.out.println("\n--- " + label + " ---");
        double[] pristine = first.clone();
        double[] a = refetch.get();
        check(label + " returns a distinct array", a != first && a != refetch.get(), "same ref");
        for (int i = 0; i < first.length; i++) first[i] = -12345.0 - i;
        double[] b = refetch.get();
        double worst = 0;
        for (int i = 0; i < pristine.length; i++) {
            worst = Math.max(worst, Math.abs(b[i] - pristine[i]));
        }
        check(label + " unaffected by caller mutation (worst |delta| " + worst + ")",
                worst == 0.0, "worst=" + worst);
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-64s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-64s [!] FAIL  %s%n", label, detail);
        }
    }

    private CvmStateImmutability() {
    }
}
