package org.ce.domain.engine.cvm;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class CVMFreeEnergyEnthalpyDerivativesTest {

    @Test
    @DisplayName("Enthalpy derivatives match H = sum(ECI_l * v_l)")
    void enthalpyDerivativesMatchExpression() {
        int ncf = 4;
        int tcf = 6;   // ncf + K for binary CVCF (K=2)
        int tcdis = 1;

        double[] v = {0.11, -0.07, 0.03, 0.19};
        double[] x = {0.5, 0.5};
        double[] eci = {3120.0, 6240.0, -260.0, -390.0};

        double temperature = 1000.0;
        List<Double> mhdis = List.of(1.0);
        double[] kb = {0.0}; // Disable entropy contribution for this derivative check
        double[][] mh = {{1.0}};
        int[] lc = {1};
        int[][] lcv = {{1}};
        List<List<int[]>> wcv = List.of(List.of(new int[]{1}));
        List<List<double[][]>> cmat = List.of(List.of((double[][]) new double[][]{
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
        }));

        CVMFreeEnergy.EvalResult eval = CVMFreeEnergy.evaluate(
                v, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        double expectedH = 0.0;
        for (int i = 0; i < ncf; i++) {
            expectedH += eci[i] * v[i];
        }
        assertEquals(expectedH, eval.H, 1.0e-12, "H must be linear in non-point CVCFs");

        for (int i = 0; i < ncf; i++) {
            assertEquals(eci[i], eval.Gcu[i], 1.0e-12, "dH/dv must equal ECI when S is disabled");
        }

        for (int i = 0; i < ncf; i++) {
            for (int j = 0; j < ncf; j++) {
                assertEquals(0.0, eval.Gcuu[i][j], 1.0e-12, "d2H/dv2 must be zero");
            }
        }

        double eps = 1.0e-8;
        for (int i = 0; i < ncf; i++) {
            double[] vp = v.clone();
            double[] vm = v.clone();
            vp[i] += eps;
            vm[i] -= eps;

            double hp = CVMFreeEnergy.evaluate(
                    vp, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf).H;
            double hm = CVMFreeEnergy.evaluate(
                    vm, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf).H;

            double fd = (hp - hm) / (2.0 * eps);
            assertEquals(eci[i], fd, 1.0e-6, "Finite-difference dH/dv mismatch at index " + i);
        }
    }
}
