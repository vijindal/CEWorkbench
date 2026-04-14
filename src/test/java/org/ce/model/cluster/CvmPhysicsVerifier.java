package org.ce.model.cluster;

import org.ce.model.cvm.CvCfBasis;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cluster.ClusterVariableEvaluator;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

/**
 * Utility for verifying CVM physics consistency.
 * Implements checks for random CFs, CV consistency, and entropy cross-basis validation.
 */
public class CvmPhysicsVerifier {

    private static final double TOL = 1e-12;
    private static final double R_GAS = 8.3144598;

    private final Consumer<String> logger;

    /**
     * Creates a verifier using the default System.out logger.
     */
    public CvmPhysicsVerifier() {
        this(System.out::println);
    }

    /**
     * Creates a verifier with a custom logger.
     */
    public CvmPhysicsVerifier(Consumer<String> logger) {
        this.logger = logger;
    }

    /**
     * Executes the full 4-step verification plan for default binary compositions.
     */
    public void verifyAll(AllClusterData data) {
        double[][] binaryCompositions = {
            {0.5, 0.5}, // Equiatomic Binary
            {0.3, 0.7}, // Off-stoichiometric Binary
            {0.1, 0.9}  // Dilute limit
        };
        verifyAll(data, "Default Binary", 0.0, binaryCompositions);
    }

    /**
     * Executes the full 4-step verification plan for custom system and compositions.
     * 
     * @param data           identification data
     * @param systemName    display name (e.g. Nb-Ti-V)
     * @param temperature   temp in K (for labeling)
     * @param compositions  list of mole fraction arrays to test
     */
    public void verifyAll(AllClusterData data, String systemName, double temperature, double[][] compositions) {
        log("================================================================================");
        log("                     CVM PHYSICS VERIFICATION REPORT");
        log(" SYSTEM: " + systemName);
        if (temperature > 1e-3) {
            log(" TEMPERATURE: " + temperature + " K");
        }
        log("================================================================================");
        
        boolean allPassed = true;
        for (double[] x : compositions) {
            try {
                verifyStep1(data, x);
                verifyStep4(data, x);
            } catch (Exception e) {
                log("  [ERROR] Steps failed for x=" + Arrays.toString(x) + ": " + e.getMessage());
                allPassed = false;
            }
        }
        
        log("================================================================================");
        if (allPassed) {
            log("                     VERIFICATION [" + systemName + "] - PASSED");
        } else {
            log("                     VERIFICATION [" + systemName + "] - FAILED");
        }
        log("================================================================================");
        
        if (!allPassed) {
             throw new RuntimeException("Physics verification failed for " + systemName + ". See report for details.");
        }
    }

    /**
     * Step 1: Random CF Generation, Orthogonal C-Matrix, & Entropy
     */
    public void verifyStep1(AllClusterData data, double[] x) {
        log("\n[STEP 1] Orthogonal Basis Verification at x=" + Arrays.toString(x));
        
        CMatrix.Result orth = data.getOrthogonalCMatrixResult();
        if (orth == null) {
            log("  FAILED: Orthogonal foundation missing from ClusterData");
            return;
        }

        int numComponents = x.length;
        int ncf = data.getDisorderedCFResult().getNcf();
        int tcf = data.getDisorderedCFResult().getTcf();
        
        // 1. Random CFs
        double[] uRand = ClusterVariableEvaluator.computeRandomCFs(
                x, numComponents, orth.getCfBasisIndices(), ncf, tcf);
        
        double[] uFull = ClusterVariableEvaluator.buildFullCFVector(
                uRand, x, numComponents, orth.getCfBasisIndices(), ncf, tcf);

        // 2. Evaluate CVs
        double[][][] cv = ClusterVariableEvaluator.evaluate(
                uFull, orth.getCmat(), orth.getLcv(), 
                data.getDisorderedClusterResult().getTcdis(), data.getDisorderedClusterResult().getLc());

        // 3. Check Consistency (Probabilities)
        checkCvProbabilityConsistency(cv, orth.getWcv(), x, "Orthogonal");

        // 4. Entropy
        double sOrth = computeEntropy(cv, orth.getWcv(), data.getDisorderedClusterResult().getKbCoefficients(), 
                                     data.getDisorderedClusterResult().getDisClusterData().getMultiplicities(),
                                     data.getDisorderedClusterResult().getLc());
        
        double sIdeal = computeIdealEntropy(x);
        double delta = Math.abs(sOrth - sIdeal);
        
        log(String.format("  - Entropy: S(orth)=%.12f, S(ideal)=%.12f, Δ=%.2e", sOrth, sIdeal, delta));
        if (delta > 1e-9) { 
            throw new RuntimeException("Orthogonal entropy mismatch! Δ=" + delta);
        }
        log("  - STEP 1 PASSED.");
    }

    /**
     * Step 4: CVCF Entropy & Cross-Basis Consistency
     */
    public void verifyStep4(AllClusterData data, double[] x) {
        log("\n[STEP 4] CVCF Basis Verification at x=" + Arrays.toString(x));
        
        CMatrix.Result cvcfCmat = data.getCMatrixResult();
        CMatrix.Result orthCmatResult = data.getOrthogonalCMatrixResult();
        int numComponents = x.length;

        // Fetch the CVCF basis from the registry for computing the proper random state
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", "T", numComponents);
        
        // 1. Compute random state in CVCF basis correctly
        // Instead of a zero vector, we must transform the orthogonal random state according to T.
        double[] vRandCVCF = basis.computeRandomState(x, orthCmatResult.getCfBasisIndices());
        
        log("  - [DIAGNOSTIC] CVCF random vector construction:");
        log("    vRandCVCF=" + Arrays.toString(vRandCVCF));

        // [Verification] Transform v back to u orthogonal to verify physical consistency
        double[][] T = basis.T; // mapping CVCF (v) to Orthogonal (u): u_i = sum_j T[i][j] * v_j
        double[] uBack = new double[T.length];
        for (int i = 0; i < T.length; i++) {
            for (int j = 0; j < vRandCVCF.length; j++) {
                uBack[i] += T[i][j] * vRandCVCF[j];
            }
        }
        double uConst = uBack[uBack.length - 1]; // orthogonal basis has constant at the end
        log(String.format("  - [VERIFICATION] u_const (after v->u transform) = %.12f", uConst));
        if (Math.abs(uConst - 1.0) > 1e-9) {
             log("    WARNING: Non-physical constant term detected in CVCF-basis random state!");
        }

        // 2. Evaluate CVs
        double[][][] cv = ClusterVariableEvaluator.evaluate(
                vRandCVCF, cvcfCmat.getCmat(), cvcfCmat.getLcv(), 
                data.getDisorderedClusterResult().getTcdis(), data.getDisorderedClusterResult().getLc());

        // [DIAGNOSTIC] Print probabilities for maximal cluster (Type 0)
        printDiagnosticCVs(cv, cvcfCmat.getWcv(), 0, "Maximal Cluster (Type 0)");

        // 3. Check Consistency (Probabilities)
        checkCvProbabilityConsistency(cv, cvcfCmat.getWcv(), x, "CVCF");

        // 4. Entropy
        double sCVCF = computeEntropy(cv, cvcfCmat.getWcv(), data.getDisorderedClusterResult().getKbCoefficients(), 
                                      data.getDisorderedClusterResult().getDisClusterData().getMultiplicities(),
                                      data.getDisorderedClusterResult().getLc());
        
        double sIdeal = computeIdealEntropy(x);
        double delta = Math.abs(sCVCF - sIdeal);
        
        log(String.format("  - Entropy: S(cvcf)=%.12f, S(ideal)=%.12f, Δ=%.2e", sCVCF, sIdeal, delta));
        
        if (delta > 1e-7) { // Using a slightly larger tolerance for CVCF numerical accumulation
            throw new RuntimeException("CVCF entropy mismatch in Step 4! Δ=" + delta);
        }
        log("  - STEP 4 PASSED.");
    }

    private void printDiagnosticCVs(double[][][] cv, List<List<int[]>> wcv, int type, String label) {
        if (type >= cv.length) return;
        log("  - [DIAGNOSTIC] Probabilities for " + label + ":");
        log(String.format("    %-10s %-15s %-15s %-15s", "Class v", "w_v", "rho_v (config)", "P_v (total)"));
        
        double totalP = 0;
        for (int j = 0; j < cv[type].length; j++) {
            double[] cvGroup = cv[type][j];
            int[] wcvGroup = wcv.get(type).get(j);
            for (int v = 0; v < cvGroup.length; v++) {
                double rho = cvGroup[v];
                double pTotal = rho * wcvGroup[v];
                totalP += pTotal;
                log(String.format("    %-10d %-15d %-15.8f %-15.8f", v, wcvGroup[v], rho, pTotal));
            }
        }
        log(String.format("    - Total Probability Sum: %.12f", totalP));
    }

    private void checkCvProbabilityConsistency(double[][][] cv, List<List<int[]>> wcv, double[] x, String label) {
        for (int t = 0; t < cv.length; t++) {
            for (int j = 0; j < cv[t].length; j++) {
                double[] cvGroup = cv[t][j];
                int[] wcvGroup = wcv.get(t).get(j);
                double sum = 0;
                for (int v = 0; v < cvGroup.length; v++) {
                    sum += cvGroup[v] * wcvGroup[v];
                }
                if (Math.abs(sum - 1.0) > 1e-7) {
                    log(String.format("  [CONSISTENCY FAIL] %s at (t=%d, j=%d): sum=%.12f", label, t, j, sum));
                    throw new RuntimeException(label + " probability sum mismatch: sum=" + sum);
                }
            }
        }
        log("  - " + label + " CV consistency OK.");
    }

    private double computeEntropy(double[][][] cv, List<List<int[]>> wcv, double[] kb, List<Double> m, int[] lc) {
        double s = 0;
        int tcdis = cv.length;
        for (int t = 0; t < tcdis; t++) {
            if (kb[t] == 0) continue;
            double sType = 0;
            for (int j = 0; j < lc[t]; j++) {
                double[] cvGroup = cv[t][j];
                int[] wcvGroup = wcv.get(t).get(j);
                for (int v = 0; v < cvGroup.length; v++) {
                    double rho = cvGroup[v];
                    if (rho > 1e-15) {
                        sType += wcvGroup[v] * rho * Math.log(rho);
                    }
                }
            }
            s += kb[t] * m.get(t) * sType;
        }
        return -R_GAS * s;
    }

    private double computeIdealEntropy(double[] x) {
        double s = 0;
        for (double xi : x) {
            if (xi > 1e-15) {
                s += xi * Math.log(xi);
            }
        }
        return -R_GAS * s;
    }

    private void log(String msg) {
        logger.accept(msg);
    }
}
