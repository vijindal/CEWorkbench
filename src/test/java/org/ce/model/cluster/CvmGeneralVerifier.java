package org.ce.model.cluster;

import org.ce.model.cluster.cvcf.CvCfBasis;
import java.util.*;
import java.util.function.Consumer;

public class CvmGeneralVerifier {

    private static final double TOL = 1e-9;
    private static final double R = 8.314;

    private final Consumer<String> log;

    public CvmGeneralVerifier(Consumer<String> log) {
        this.log = log;
    }

    // ============================================================
    // MASTER ENTRY
    // ============================================================
    public void verifyAll(AllClusterData data, String name,
                          double T, double[][] compositions) {

        log.accept("\n========================================");
        log.accept("SYSTEM: " + name);
        log.accept("========================================");

        checkStructure(data);

        for (double[] x : compositions) {
            log.accept("\n--- Composition: " + Arrays.toString(x) + " ---");

            verifyOrthogonal(data, x);
            verifyTransformation(data, x);
            verifyCVCF(data, x);
            verifyEdgeCases(x);
        }

        log.accept("\nSYSTEM PASSED: " + name);
    }

    // ============================================================
    // STEP 0: STRUCTURE
    // ============================================================
    private void checkStructure(AllClusterData data) {
        if (!data.hasOrthogonalFoundation()) {
            throw new RuntimeException("Missing orthogonal foundation");
        }
        log.accept("✔ Structure OK");
    }

    // ============================================================
    // STEP 1: ORTHOGONAL BASIS
    // ============================================================
    private void verifyOrthogonal(AllClusterData data, double[] x) {

        var orth = data.getOrthogonalCMatrixResult();

        double[] uRand = ClusterVariableEvaluator.computeRandomCFs(
                x, x.length, orth.getCfBasisIndices(),
                data.getDisorderedCFResult().getNcf(),
                data.getDisorderedCFResult().getTcf());

        double[] uFull = ClusterVariableEvaluator.buildFullCFVector(
                uRand, x, x.length, orth.getCfBasisIndices(),
                data.getDisorderedCFResult().getNcf(),
                data.getDisorderedCFResult().getTcf());

        double[][][] cv = ClusterVariableEvaluator.evaluate(
                uFull, orth.getCmat(), orth.getLcv(),
                data.getDisorderedClusterResult().getTcdis(),
                data.getDisorderedClusterResult().getLc());

        checkProbability(cv, orth.getWcv());
        checkEntropy(cv, orth.getWcv(), data, x);
    }

    // ============================================================
    // STEP 2: TRANSFORMATION (T)
    // ============================================================
    private void verifyTransformation(AllClusterData data, double[] x) {

        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get(
                "BCC_A2", "T", x.length);

        double[] v = basis.computeRandomState(x,
                data.getOrthogonalCMatrixResult().getCfBasisIndices());

        double[][] T = basis.T;

        double[] uBack = new double[T.length];

        for (int i = 0; i < T.length; i++) {
            for (int j = 0; j < v.length; j++) {
                uBack[i] += T[i][j] * v[j];
            }
        }

        if (Math.abs(uBack[uBack.length - 1] - 1.0) > 1e-6) {
            throw new RuntimeException("Transformation inconsistency");
        }

        log.accept("✔ Transformation OK");
    }

    // ============================================================
    // STEP 3: CVCF BASIS
    // ============================================================
    private void verifyCVCF(AllClusterData data, double[] x) {

        var cmat = data.getCMatrixResult();

        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get(
                "BCC_A2", "T", x.length);

        double[] v = basis.computeRandomState(x,
                data.getOrthogonalCMatrixResult().getCfBasisIndices());

        double[][][] cv = ClusterVariableEvaluator.evaluate(
                v, cmat.getCmat(), cmat.getLcv(),
                data.getDisorderedClusterResult().getTcdis(),
                data.getDisorderedClusterResult().getLc());

        checkProbability(cv, cmat.getWcv());
        checkEntropy(cv, cmat.getWcv(), data, x);
    }

    // ============================================================
    // STEP 4: EDGE CASES
    // ============================================================
    private void verifyEdgeCases(double[] x) {

        double sum = Arrays.stream(x).sum();
        if (Math.abs(sum - 1.0) > 1e-9) {
            throw new RuntimeException("Composition not normalized");
        }

        for (double xi : x) {
            if (xi < -1e-12) {
                throw new RuntimeException("Negative composition");
            }
        }
    }

    // ============================================================
    // UTILITIES
    // ============================================================

    private void checkProbability(double[][][] cv, List<List<int[]>> wcv) {
        for (int t = 0; t < cv.length; t++) {
            for (int j = 0; j < cv[t].length; j++) {

                double sum = 0;
                for (int v = 0; v < cv[t][j].length; v++) {
                    double rho = cv[t][j][v];

                    if (rho < -1e-12) {
                        throw new RuntimeException("Negative probability");
                    }

                    sum += rho * wcv.get(t).get(j)[v];
                }

                if (Math.abs(sum - 1.0) > TOL) {
                    throw new RuntimeException("Probability mismatch");
                }
            }
        }

        log.accept("✔ Probability OK");
    }

    private void checkEntropy(double[][][] cv, List<List<int[]>> wcv,
                              AllClusterData data, double[] x) {

        double s = computeEntropy(cv, wcv, data);
        double sIdeal = computeIdealEntropy(x);

        double delta = Math.abs(s - sIdeal);

        log.accept(String.format("Entropy Δ = %.3e", delta));

        if (delta > 1e-7) {
            throw new RuntimeException("Entropy mismatch");
        }
    }

    private double computeEntropy(double[][][] cv, List<List<int[]>> wcv,
                                 AllClusterData data) {

        double s = 0;

        var kb = data.getDisorderedClusterResult().getKbCoefficients();
        var m  = data.getDisorderedClusterResult().getDisClusterData().getMultiplicities();
        int[] lc = data.getDisorderedClusterResult().getLc();

        for (int t = 0; t < cv.length; t++) {
            double sType = 0;

            for (int j = 0; j < lc[t]; j++) {
                for (int v = 0; v < cv[t][j].length; v++) {

                    double rho = cv[t][j][v];

                    if (rho > 1e-15) {
                        sType += wcv.get(t).get(j)[v] * rho * Math.log(rho);
                    }
                }
            }

            s += kb[t] * m.get(t) * sType;
        }

        return -R * s;
    }

    private double computeIdealEntropy(double[] x) {
        double s = 0;
        for (double xi : x) {
            if (xi > 1e-15) {
                s += xi * Math.log(xi);
            }
        }
        return -R * s;
    }
}
