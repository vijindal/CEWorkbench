package org.ce.domain.engine.cvm;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test entropy calculation with direct CVCF input values.
 *
 * Tests the CVM entropy formula with specific, simple CVCF values
 * to validate the entropy computation against Mathematica reference.
 *
 * Test case:
 *   Binary BCC_A2 system at T=1000K, equimolar (x_A=0.5, x_B=0.5)
 *   CVCF input:
 *     v4AB  = 0.0625
 *     v3AB  = 0.0
 *     v22AB = 0.25
 *     v21AB = 0.25
 */
@DisplayName("CVM Entropy — Direct CVCF Input Test")
class CVMEntropyDirectTest {

    @Test
    @DisplayName("Entropy with CVCF values: v4AB=0.0625, v3AB=0.0, v22AB=0.25, v21AB=0.25")
    void entropyWithDirectCVCFValues() {
        // === Binary BCC_A2 Parameters ===
        int ncf = 4;   // v4AB, v3AB, v22AB, v21AB
        int tcf = 6;   // ncf + K (K=2 for binary)
        int tcdis = 1; // 1 HSP type (point cluster)

        // Direct CVCF values as specified
        double[] v = {0.0625, 0.0, 0.25, 0.25};  // [v4AB, v3AB, v22AB, v21AB]
        double[] x = {0.5, 0.5};  // equimolar
        double temperature = 1000.0;

        // Enthalpy (not relevant for pure entropy test, set ECIs to zero)
        double[] eci = {0.0, 0.0, 0.0, 0.0};

        // === CVM parameters ===
        List<Double> mhdis = List.of(1.0);  // 1 HSP type with multiplicity 1
        double[] kb = {-1.0};  // Kikuchi-Barker coeff for point cluster (-R is applied separately)
        double[][] mh = {{1.0}};  // normalized multiplicity
        int[] lc = {1};  // 1 ordered cluster (the point cluster)

        // === C-matrix for point cluster ===
        // Point cluster has 1 correlation variable (the composition itself, xB)
        // CV[0] = sum of weights * xB = 1 * xB = 0.5
        int[][] lcv = {{1}};  // 1 CV for point cluster
        List<List<int[]>> wcv = List.of(List.of(new int[]{1}));  // weight = 1

        // C-matrix: cmat[t][j][v][k] where t=0, j=0, v=0 (point cluster has 1 CV)
        // [cmat[0..3] | cmat[4] cmat[5] | const]
        // For point cluster: CV = xB = 0.5
        // cmat[i] * v[i] for i=0..3 contribute to CV, cmat[4] * xA + cmat[5] * xB + const = xB
        // So: cmat[0..3] = [0, 0, 0, 0], cmat[4] = 0, cmat[5] = 1, const = 0
        double[][] cmPointGroup = {{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};  // 1 CV x 6 CF columns
        List<double[][]> groupList = List.of(cmPointGroup);
        List<List<double[][]>> cmat = List.of(groupList);

        // === Evaluate entropy ===
        CVMFreeEnergy.EvalResult eval = CVMFreeEnergy.evaluate(
                v, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        // === Expected entropy calculation ===
        // S = -R * kb[0] * mhdis[0] * mh[0][0] * w[0] * (xB * ln(xB) + xA * ln(xA))
        // With xA=0.5, xB=0.5:
        //   xB * ln(xB) = 0.5 * ln(0.5) = 0.5 * (-0.693147) = -0.346574
        //   xA * ln(xA) = 0.5 * ln(0.5) = -0.346574
        //   sum = -0.693147
        //
        // S = -R * (-1) * 1.0 * 1.0 * 1 * (-0.693147)
        //   = R * 0.693147
        //   = 8.3144598 * 0.693147
        //   ≈ 5.763
        double xA = x[0];
        double xB = x[1];
        double R = CVMFreeEnergy.R_GAS;
        double expectedS = -R * (-1.0) * 1.0 * 1.0 * 1.0 * (xA * Math.log(xA) + xB * Math.log(xB));

        System.out.println("=== Entropy Test: Direct CVCF Input ===");
        System.out.println("CVCF values: v4AB=" + v[0] + ", v3AB=" + v[1] + ", v22AB=" + v[2] + ", v21AB=" + v[3]);
        System.out.println("Composition: x_A=" + xA + ", x_B=" + xB);
        System.out.println("Temperature: " + temperature + " K");
        System.out.println("R_GAS: " + R);
        System.out.println("Expected entropy S: " + expectedS);
        System.out.println("Computed entropy S: " + eval.S);
        System.out.println("Difference: " + Math.abs(eval.S - expectedS));

        // Check entropy value
        assertEquals(expectedS, eval.S, 1.0e-6,
                "Entropy does not match expected value from direct CVCF inputs");

        // Check that entropy derivatives (gradient and Hessian) are computed
        // For point cluster only, entropy doesn't depend on non-point CVCFs,
        // so Scu[l] should be zero for all non-point l.
        for (int l = 0; l < ncf; l++) {
            double expectedScu = 0.0;  // Point cluster entropy doesn't depend on non-point CVCFs
            assertEquals(expectedScu, eval.Gcu[l], 1.0e-12,
                    "dS/dv[" + l + "] should be 0 (point cluster only)");
        }

        // Hessian d²S/dv² should also be zero (entropy from point cluster is independent of v)
        for (int l1 = 0; l1 < ncf; l1++) {
            for (int l2 = 0; l2 < ncf; l2++) {
                double expectedScuu = 0.0;
                assertEquals(expectedScuu, eval.Gcuu[l1][l2], 1.0e-12,
                        "d²S/dv²[" + l1 + "][" + l2 + "] should be 0 (point cluster only)");
            }
        }
    }

    @Test
    @DisplayName("Entropy finite-difference validation")
    void entropyFiniteDifferenceValidation() {
        int ncf = 4;
        int tcf = 6;
        int tcdis = 1;

        double[] v = {0.0625, 0.0, 0.25, 0.25};
        double[] x = {0.5, 0.5};
        double temperature = 1000.0;
        double[] eci = {0.0, 0.0, 0.0, 0.0};

        List<Double> mhdis = List.of(1.0);
        double[] kb = {-1.0};
        double[][] mh = {{1.0}};
        int[] lc = {1};
        int[][] lcv = {{1}};
        List<List<int[]>> wcv = List.of(List.of(new int[]{1}));
        double[][] cmPointGroup = {{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
        List<List<double[][]>> cmat = List.of(List.of(cmPointGroup));

        // Evaluate at v
        CVMFreeEnergy.EvalResult eval = CVMFreeEnergy.evaluate(
                v, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        // Finite difference check on entropy (should be zero derivative for point cluster only)
        double eps = 1.0e-8;
        for (int i = 0; i < ncf; i++) {
            double[] vp = v.clone();
            double[] vm = v.clone();
            vp[i] += eps;
            vm[i] -= eps;

            double sp = CVMFreeEnergy.evaluate(
                    vp, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf).S;
            double sm = CVMFreeEnergy.evaluate(
                    vm, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf).S;

            double fdDeriv = (sp - sm) / (2.0 * eps);
            System.out.println("dS/dv[" + i + "] analytic=" + eval.Gcu[i] + ", finite-diff=" + fdDeriv);

            assertEquals(eval.Gcu[i], fdDeriv, 1.0e-6,
                    "Entropy derivative dS/dv[" + i + "] finite-difference validation failed");
        }
    }

    @Test
    @DisplayName("Entropy value consistency across different CVCF values")
    void entropyConsistencyAcrossValues() {
        // Point cluster entropy should be independent of the CVCF values (v)
        // Only composition matters
        int ncf = 4;
        int tcf = 6;
        int tcdis = 1;
        double[] x = {0.5, 0.5};
        double temperature = 1000.0;
        double[] eci = {0.0, 0.0, 0.0, 0.0};
        List<Double> mhdis = List.of(1.0);
        double[] kb = {-1.0};
        double[][] mh = {{1.0}};
        int[] lc = {1};
        int[][] lcv = {{1}};
        List<List<int[]>> wcv = List.of(List.of(new int[]{1}));
        double[][] cmPointGroup = {{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
        List<List<double[][]>> cmat = List.of(List.of(cmPointGroup));

        // Compute entropy with first set of CVCF values
        double[] v1 = {0.0625, 0.0, 0.25, 0.25};
        CVMFreeEnergy.EvalResult eval1 = CVMFreeEnergy.evaluate(
                v1, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        // Compute entropy with different CVCF values
        double[] v2 = {-0.1, 0.15, -0.05, 0.33};
        CVMFreeEnergy.EvalResult eval2 = CVMFreeEnergy.evaluate(
                v2, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        // Both should give the same entropy (composition is the same)
        assertEquals(eval1.S, eval2.S, 1.0e-12,
                "Point cluster entropy should be independent of CVCF v values");

        System.out.println("Entropy consistency test: S1=" + eval1.S + ", S2=" + eval2.S);
    }
}
