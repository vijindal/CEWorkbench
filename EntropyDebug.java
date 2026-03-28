import org.ce.domain.engine.cvm.CVMFreeEnergy;
import java.util.ArrayList;
import java.util.List;

public class EntropyDebug {
    public static void main(String[] args) {
        System.out.println("=== ENTROPY DEBUG WITH DIRECT CVCF VALUES ===\n");

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
        double xA = x[0];
        double xB = x[1];
        double R = CVMFreeEnergy.R_GAS;
        double expectedS = -R * (-1.0) * 1.0 * 1.0 * 1.0 * (xA * Math.log(xA) + xB * Math.log(xB));

        // === Print results ===
        System.out.println("Input CVCF values:");
        System.out.println("  v4AB  = " + v[0]);
        System.out.println("  v3AB  = " + v[1]);
        System.out.println("  v22AB = " + v[2]);
        System.out.println("  v21AB = " + v[3]);
        System.out.println("Composition: x_A=" + xA + ", x_B=" + xB);
        System.out.println("Temperature: " + temperature + " K");
        System.out.println();

        System.out.println("Result from CVMFreeEnergy.evaluate():");
        System.out.println("  G   = " + String.format("%.12e", result.G));
        System.out.println("  H   = " + String.format("%.12e", result.H));
        System.out.println("  S   = " + String.format("%.12e", result.S));
        System.out.println();

        System.out.println("Expected entropy (from composition only):");
        System.out.println("  S_expected = " + String.format("%.12e", expectedS));
        System.out.println();

        System.out.println("Gradient (dG/dv):");
        for (int i = 0; i < ncf; i++) {
            System.out.println("  Gcu[" + i + "] = " + String.format("%.12e", result.Gcu[i]));
        }
        System.out.println();

        System.out.println("Hessian (d²G/dv²) [non-zero values only]:");
        boolean hasNonZero = false;
        for (int i = 0; i < ncf; i++) {
            for (int j = 0; j < ncf; j++) {
                if (Math.abs(result.Gcuu[i][j]) > 1e-15) {
                    System.out.println("  Gcuu[" + i + "][" + j + "] = " + String.format("%.12e", result.Gcuu[i][j]));
                    hasNonZero = true;
                }
            }
        }
        if (!hasNonZero) {
            System.out.println("  (all zero - as expected for point cluster only)");
        }
        System.out.println();

        System.out.println("Comparison:");
        System.out.println("  S_computed  = " + String.format("%.12e", result.S));
        System.out.println("  S_expected  = " + String.format("%.12e", expectedS));
        System.out.println("  Difference  = " + String.format("%.12e", Math.abs(result.S - expectedS)));
        System.out.println("  Match?      " + (Math.abs(result.S - expectedS) < 1e-6 ? "YES ✓" : "NO ✗"));
    }
}
