package org.ce.domain.cluster.cvcf;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.CMatrix;
import org.ce.model.cluster.ClusterVariableEvaluator;
import org.ce.model.cluster.cvcf.CvCfBasis;
import org.ce.model.cluster.cvcf.CvCfBasisTransformer;
import org.ce.model.cluster.cvcf.BccA2TModelCvCfTransformations;
import org.ce.calculation.workflow.ClusterIdentificationRequest;
import org.ce.calculation.workflow.ClusterIdentificationWorkflow;
import org.junit.jupiter.api.Test;

/**
 * Diagnostic test to verify the numerical consistency of Cluster Variables (CVs)
 * between the legacy orthogonal basis and the new CVCF basis.
 *
 * This test uses a specific ternary Nb-Ti-V state vector provided by the user.
 */
public class CvCfTransformationTest {

    private static final double TOLERANCE = 1e-10;

    @Test
    public void testTernaryConsistency() {
        // 1. Setup the ternary BCC_A2 basis (Nb-Ti-V)
        CvCfBasis basis = BccA2TModelCvCfTransformations.ternaryBasis();
        int ncf = basis.numNonPointCfs; // 18
        int tcf = basis.totalCfs();     // 21 (18 + 3 point variables)

        // 2. Define the user-provided CVCF non-point variables (v)
        // Order matches TERNARY_CF_NAMES: v4s, v3s, v22s, v21s
        double[] v = new double[ncf];
        v[0]  = 0.0116752;   // v4AB
        v[1]  = 0.00953038;  // v4AC
        v[2]  = 0.0192938;   // v4BC
        v[3]  = 0.00674946;  // v4ABC1
        v[4]  = 0.0222047;   // v4ABC2
        v[5]  = 0.0215488;   // v4ABC3
        v[6]  = 0.012542;    // v3AB
        v[7]  = 0.0149128;   // v3AC
        v[8]  = 0.00225646;  // v3BC
        v[9]  = 0.029726;    // v3ABC1
        v[10] = 0.0504037;   // v3ABC2
        v[11] = 0.048904;    // v3ABC3
        v[12] = 0.118221;    // v22AB (2nd neighbor)
        v[13] = 0.118642;    // v22AC (2nd neighbor)
        v[14] = 0.0930865;   // v22BC (2nd neighbor)
        v[15] = 0.111417;    // v21AB (1st neighbor)
        v[16] = 0.107397;    // v21AC (1st neighbor)
        v[17] = 0.133421;    // v21BC (1st neighbor)

        // Point variables: xA, xB, xC (User-provided)
        double[] x = { 0.33, 0.33, 0.34 };

        // Full CVCF vector: [v... | xA, xB, xC]
        double[] vFull = ClusterVariableEvaluator.buildFullCVCFVector(v, x, ncf, tcf);

        // 3. Get proper Orthogonal C-Matrix by running Stage 3 of the Type-1a workflow
        // We use the actual resource files for BCC_A2 T-model
        ClusterIdentificationRequest request = ClusterIdentificationRequest.builder()
                .disorderedClusterFile("clus/BCC_A2-T.txt")
                .orderedClusterFile("clus/BCC_A2-T.txt")
                .disorderedSymmetryGroup("BCC_A2-SG")
                .orderedSymmetryGroup("BCC_A2-SG")
                .numComponents(3)
                .structurePhase("BCC_A2")
                .model("T")
                .build();

        AllClusterData fullData = ClusterIdentificationWorkflow.identify(request, System.out::println);
        CMatrix.Result orthoResult = fullData.getCMatrixResult();
        
        // Extract block for printing
        double[][] orthoCmatBlock = orthoResult.getCmat().get(0).get(0);
        int[][] lcv = orthoResult.getLcv();
        int[] lc = new int[lcv.length];
        for (int i = 0; i < lcv.length; i++) lc[i] = orthoResult.getCmat().get(i).size();
        int nv = orthoCmatBlock.length;


        // 4. Transform to CVCF basis
        CMatrix.Result cvcfResult = CvCfBasisTransformer.transform(orthoResult, basis);
        double[][] cvcfCmatBlock = cvcfResult.getCmat().get(0).get(0);


        // 5. Calculate CVs in CVCF basis: rho_cvcf = C_cvcf * vFull + k
        double[][][] rhosCvcf = ClusterVariableEvaluator.evaluate(vFull, cvcfResult.getCmat(), lcv, 1, lc);
        double[] resultsCvcf = rhosCvcf[0][0];

        // 6. Calculate CVs in Orthogonal basis
        // First convert CVCF vector vFull to Orthogonal vector uFull using T matrix: u = T * v
        double[] uFull = new double[tcf];
        double[][] T = basis.T;
        for (int i = 0; i < tcf; i++) {
            double sum = 0.0;
            for (int j = 0; j < tcf; j++) {
                sum += T[i][j] * vFull[j];
            }
            uFull[i] = sum;
        }

        double[][][] rhosOrth = ClusterVariableEvaluator.evaluate(uFull, orthoResult.getCmat(), lcv, 1, lc);
        double[] resultsOrth = rhosOrth[0][0];

        // 7. Systematic Printing
        System.out.println("\n" + "=".repeat(80));
        System.out.println("CVCF BASIS TRANSFORMATION CONSISTENCY DIAGNOSTIC (Nb-Ti-V)");
        System.out.println("=".repeat(80));

        System.out.println("\n[STAGE 1] CVCF STATE VECTOR (v)");
        System.out.println("-".repeat(80));
        for (int i = 0; i < tcf; i++) {
            String name = (i < basis.cfNames.size()) ? basis.cfNames.get(i) : "const";
            System.out.printf("  v[%02d] (%-8s) : %16.10f\n", i, name, vFull[i]);
        }

        System.out.println("\n[STAGE 2] TRANSFORMED ORTHOGONAL STATE VECTOR (u = T \u00b7 v)");
        System.out.println("-".repeat(80));
        for (int i = 0; i < tcf; i++) {
            System.out.printf("  u[%02d]           : %16.10f\n", i + 1, uFull[i]);
        }

        System.out.println("\n[STAGE 3] C-MATRIX COMPARISON (Coefficients)");
        System.out.println("-".repeat(80));
        printMatrix("ORTHOGONAL C-MATRIX (Block 0,0)", orthoCmatBlock);
        printMatrix("CVCF C-MATRIX (Transformed)", cvcfCmatBlock);

        System.out.println("\n[STAGE 4] CLUSTER VARIABLE (\u03c1) CALCULATION (Dot Product)");
        System.out.println("-".repeat(80));
        System.out.println("  Basis Identity: \u03c1_cvcf = C_cvcf \u00b7 v_full  ==  \u03c1_orth = C_orth \u00b7 u_full");
        System.out.println("");
        System.out.printf("  %-4s  %-20s  %-20s  %-14s\n", "Idx", "CVCF Result", "Orthogonal Result", "Difference");
        for (int i = 0; i < nv; i++) {
            double diff = Math.abs(resultsCvcf[i] - resultsOrth[i]);
            System.out.printf("  [%02d]  %20.12f  %20.12f  %14.4e\n", 
                i, resultsCvcf[i], resultsOrth[i], diff);
            assertEquals(resultsOrth[i], resultsCvcf[i], TOLERANCE, "CV consistency failed at index " + i);
        }

        System.out.println("\n[STATUS] Basis Transformation numerically verified (Tolerance: " + TOLERANCE + ")");
        System.out.println("=".repeat(80));
    }

    private void printMatrix(String title, double[][] matrix) {
        System.out.println("  " + title + ":");
        int rows = matrix.length;
        int cols = matrix[0].length;
        for (int i = 0; i < rows; i++) {
            System.out.print("    Row " + i + " : [ ");
            for (int k = 0; k < cols; k++) {
                System.out.printf("%10.4f ", matrix[i][k]);
            }
            System.out.println("]");
        }
        System.out.println("");
    }
}
