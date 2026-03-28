package org.ce.domain.engine.cvm;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Direct entropy debug test for CVCF basis.
 *
 * Tests entropy calculation with specific CVCF values:
 *   v4AB  = 0.0625
 *   v3AB  = 0.0
 *   v22AB = 0.25
 *   v21AB = 0.25
 */
@DisplayName("CVM Entropy Debug — Direct CVCF Values")
class CVMEntropyDebugTest {

    @Test
    @DisplayName("Entropy with CVCF: v4AB=0.0625, v3AB=0.0, v22AB=0.25, v21AB=0.25")
    void debugEntropyDirectValues() {
        // === Setup ===
        int ncf = 4;   // v4AB, v3AB, v22AB, v21AB
        int tcf = 6;   // ncf + K (K=2 for binary)
        int tcdis = 1; // 1 HSP type (point cluster)

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

        // Point cluster C-matrix: CV = xB
        List<List<double[][]>> cmat = new ArrayList<>();
        List<double[][]> group0 = new ArrayList<>();
        group0.add(new double[][]{{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}});
        cmat.add(group0);

        // === Evaluate ===
        CVMFreeEnergy.EvalResult result = CVMFreeEnergy.evaluate(
            v, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        // === Expected entropy ===
        // For point cluster: S = -R * kb[0] * mhdis[0] * mh[0][0] * w[0] * (xA*ln(xA) + xB*ln(xB))
        // With xA=0.5, xB=0.5:
        //   xA*ln(xA) = 0.5*ln(0.5) = -0.346574
        //   xB*ln(xB) = 0.5*ln(0.5) = -0.346574
        //   sum = -0.693147
        // S = -8.3144598 * (-1) * 1 * 1 * 1 * (-0.693147)
        //   = 8.3144598 * 0.693147 ≈ 5.763

        double xA = x[0];
        double xB = x[1];
        double R = CVMFreeEnergy.R_GAS;
        double expectedS = -R * (-1.0) * 1.0 * 1.0 * 1.0 * (xA * Math.log(xA) + xB * Math.log(xB));

        // === Print results ===
        System.out.println("\n=== ENTROPY DEBUG TEST ===");
        System.out.println("Input CVCF values:");
        System.out.println("  v4AB  = " + v[0]);
        System.out.println("  v3AB  = " + v[1]);
        System.out.println("  v22AB = " + v[2]);
        System.out.println("  v21AB = " + v[3]);
        System.out.println("Composition: x_A=" + xA + ", x_B=" + xB);
        System.out.println("Temperature: " + temperature + " K");
        System.out.println();
        System.out.println("Result from CVMFreeEnergy.evaluate():");
        System.out.println("  G   = " + result.G);
        System.out.println("  H   = " + result.H);
        System.out.println("  S   = " + result.S);
        System.out.println();
        System.out.println("Expected entropy (from composition only):");
        System.out.println("  S_expected = " + expectedS);
        System.out.println();
        System.out.println("Gradient (dG/dv):");
        for (int i = 0; i < ncf; i++) {
            System.out.println("  Gcu[" + i + "] = " + result.Gcu[i]);
        }
        System.out.println();
        System.out.println("Hessian (d²G/dv²):");
        for (int i = 0; i < ncf; i++) {
            for (int j = 0; j < ncf; j++) {
                if (Math.abs(result.Gcuu[i][j]) > 1e-15) {
                    System.out.println("  Gcuu[" + i + "][" + j + "] = " + result.Gcuu[i][j]);
                }
            }
        }
        System.out.println();
        System.out.println("Comparison:");
        System.out.println("  S_computed  = " + result.S);
        System.out.println("  S_expected  = " + expectedS);
        System.out.println("  Difference  = " + Math.abs(result.S - expectedS));
        System.out.println("  Match? " + (Math.abs(result.S - expectedS) < 1e-6 ? "YES ✓" : "NO ✗"));
    }

    @Test
    @DisplayName("Entropy independence from non-point CVCF values")
    void debugEntropyIndependence() {
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
        List<List<double[][]>> cmat = new ArrayList<>();
        List<double[][]> group = new ArrayList<>();
        group.add(new double[][]{{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}});
        cmat.add(group);

        // Test with different CVCF values
        double[] v1 = {0.0625, 0.0, 0.25, 0.25};
        double[] v2 = {-0.1, 0.15, -0.05, 0.33};
        double[] v3 = {0.123, -0.456, 0.789, -0.012};

        CVMFreeEnergy.EvalResult eval1 = CVMFreeEnergy.evaluate(
            v1, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);
        CVMFreeEnergy.EvalResult eval2 = CVMFreeEnergy.evaluate(
            v2, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);
        CVMFreeEnergy.EvalResult eval3 = CVMFreeEnergy.evaluate(
            v3, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        System.out.println("\n=== ENTROPY INDEPENDENCE TEST ===");
        System.out.println("Test if entropy is independent of non-point CVCF values");
        System.out.println("(Should all be the same, depending only on composition)");
        System.out.println();
        System.out.println("v1 = " + java.util.Arrays.toString(v1) + " → S = " + eval1.S);
        System.out.println("v2 = " + java.util.Arrays.toString(v2) + " → S = " + eval2.S);
        System.out.println("v3 = " + java.util.Arrays.toString(v3) + " → S = " + eval3.S);
        System.out.println();
        System.out.println("S1 == S2? " + (Math.abs(eval1.S - eval2.S) < 1e-12 ? "YES ✓" : "NO ✗ Diff: " + Math.abs(eval1.S - eval2.S)));
        System.out.println("S1 == S3? " + (Math.abs(eval1.S - eval3.S) < 1e-12 ? "YES ✓" : "NO ✗ Diff: " + Math.abs(eval1.S - eval3.S)));
    }
}
