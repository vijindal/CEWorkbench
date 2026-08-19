package org.ce.scratch;

import org.ce.model.cvm.CvmGeometry;

import java.util.Arrays;

/**
 * Standalone smoke test for {@link CvmGeometry}: builds the geometry for
 * several system identities and runs {@link CvmGeometry#validate()} on each.
 *
 * <p>Deliberately touches no Hamiltonian, no temperature, and no correlation
 * functions -- if this passes, the lattice combinatorics feeding every
 * free-energy expression are structurally sound, independently of whether the
 * energy expressions themselves are.</p>
 */
public final class CvmGeometryProbe {

    public static void main(String[] args) {
        String[][] cases = {
                {"Nb-Ti", "BCC_A2", "T"},
                {"Nb-Ti-V", "BCC_A2", "T"},
                {"Nb-Ti-V-Zr", "BCC_A2", "T"},
        };

        int pass = 0, fail = 0;
        for (String[] c : cases) {
            System.out.println("\n==================================================");
            System.out.printf("BUILD  %s / %s / %s%n", c[0], c[1], c[2]);
            System.out.println("==================================================");
            try {
                CvmGeometry geo = CvmGeometry.build(c[0], c[1], c[2], null);
                System.out.println("  built:  " + geo);
                geo.validate();
                System.out.println("  validate(): PASS");
                dump(geo);
                if (!randomStateCheck(geo)) {
                    fail++;
                    continue;
                }
                pass++;
            } catch (Throwable t) {
                System.out.println("  FAILED: " + t.getClass().getSimpleName() + ": " + t.getMessage());
                for (StackTraceElement e : t.getStackTrace()) {
                    if (e.getClassName().startsWith("org.ce")) {
                        System.out.println("      at " + e);
                    }
                }
                fail++;
            }
        }

        System.out.println("\n==================================================");
        System.out.printf("RESULT: %s   (%d passed, %d failed)%n",
                fail == 0 ? "PASS" : "FAIL", pass, fail);
        System.out.println("==================================================");
    }

    /**
     * Random-state check, written entirely against {@link CvmGeometry}'s
     * public cluster algebra -- no field poking, no pipeline re-entry.
     *
     * <p>Lives here rather than on {@code CvmGeometry} because it needs a
     * composition to probe at: that is a test fixture, not cluster algebra,
     * and the geometry class deliberately exposes no method taking a
     * thermodynamic input other than as a plain vector argument.</p>
     *
     * <p>At equiatomic composition every configuration of an n-site cluster is
     * equally likely, so each cluster probability must be exactly
     * {@code (1/K)^n} -- an analytic ground truth independent of any
     * Hamiltonian. Normalisation {@code sum_v wcv*cv = 1} must hold at any
     * composition.</p>
     */
    private static boolean randomStateCheck(CvmGeometry geo) {
        int K = geo.numComponents;
        double[] x = new double[K];
        Arrays.fill(x, 1.0 / K);

        // The random-state CVCF vector. The pipeline caches it only when built
        // via runFullWorkflow; under plain run() it is null, so derive it from
        // the orthogonal random CFs through the basis transform. Note a
        // zero-filled u is NOT the random state -- CVCF CFs at random equal the
        // cluster probabilities themselves, not zero.
        double[] uFull = geo.pipelineResult.getEquiatomicCVCF();
        if (uFull == null) {
            uFull = geo.basis.computeRandomCvcfCFs(x, geo.pipelineResult);
        }
        double[][][] cv = geo.evaluateCVsFull(uFull);

        System.out.println("  --- random state (equiatomic) ---");
        int bad = 0;
        for (int t = 0; t < cv.length; t++) {
            int n = geo.pipelineResult.getDisClusData().getClusCoordList().get(t).getAllSites().size();
            double expected = Math.pow(1.0 / K, n);
            for (int j = 0; j < cv[t].length; j++) {
                int[] w = geo.wcv.get(t).get(j);
                double sum = 0, maxErr = 0;
                for (int v = 0; v < cv[t][j].length; v++) {
                    sum += w[v] * cv[t][j][v];
                    maxErr = Math.max(maxErr, Math.abs(cv[t][j][v] - expected));
                }
                boolean ok = maxErr <= 1e-9 && Math.abs(sum - 1.0) <= 1e-9;
                if (!ok) bad++;
                System.out.printf("    t=%d j=%d (n=%d): expect %.8f, maxErr %.2e, "
                                + "sum(wcv*cv)=%.12f  %s%n",
                        t, j, n, expected, maxErr, sum, ok ? "OK" : "[!] FAIL");
            }
        }
        System.out.println("    random-state check: " + (bad == 0 ? "PASS" : "FAIL (" + bad + " blocks)"));
        return bad == 0;
    }

    private static void dump(CvmGeometry geo) {
        System.out.println("  --- entropy prefactors ---");
        System.out.println("    tcdis  = " + geo.tcdis);
        System.out.println("    kb     = " + Arrays.toString(geo.kb));
        System.out.println("    mhdis  = " + Arrays.toString(geo.mhdis));
        System.out.println("    lc     = " + Arrays.toString(geo.lc));
        for (int t = 0; t < geo.tcdis; t++) {
            System.out.printf("    mh[%d]  = %s%n", t, Arrays.toString(geo.mh[t]));
            System.out.printf("    lcv[%d] = %s%n", t, Arrays.toString(geo.lcv[t]));
        }
        System.out.println("  --- basis ---");
        System.out.println("    ncf    = " + geo.ncf + "   tcf = " + geo.tcf
                + "   (tcf - ncf = " + (geo.tcf - geo.ncf) + ", K = " + geo.numComponents + ")");
        System.out.println("  --- cmat shape ---");
        for (int t = 0; t < geo.tcdis; t++) {
            for (int j = 0; j < geo.lc[t]; j++) {
                double[][] b = geo.cmat.get(t).get(j);
                System.out.printf("    cmat[%d][%d]: %d rows x %d cols, sum(wcv)=%d%n",
                        t, j, b.length, b[0].length,
                        Arrays.stream(geo.wcv.get(t).get(j)).sum());
            }
        }
    }
}
