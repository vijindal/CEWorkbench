package org.ce.model.hamiltonian;

import java.util.*;
import org.ce.model.cluster.LinearAlgebra;

/**
 * Numerical CEC–CVM Transformer: basis transformation from lower-order (K_low) to higher-order (K_high) systems.
 *
 * Implements least-squares regression on random cluster configurations to properly transform
 * Cluster Expansion Coefficients (ECIs) from one basis to another.
 *
 * This is essential for cross-order usage of ECIs, e.g., transforming binary (K=2) ECIs
 * to ternary (K=3) basis for a ternary system.
 */
public class NumericalCECTransformer {

    // ============================================================
    // 🔷 Species Mapping
    // ============================================================

    /**
     * Maps species indices from a higher-order system to a lower-order subsystem.
     *
     * Example: For a ternary system (K=3: {A, B, C}) with binary subsystem Nb-V,
     * the mapping specifies which ternary species correspond to which binary species.
     */
    public static class SpeciesMapping {
        private final Map<Integer, Integer> map = new HashMap<>();

        /**
         * Add a mapping: higher-order index → lower-order index.
         *
         * @param highIndex species index in the higher-order system
         * @param lowIndex  species index in the lower-order subsystem
         */
        public void addMapping(int highIndex, int lowIndex) {
            map.put(highIndex, lowIndex);
        }

        public int getLowIndex(int highIndex) {
            return map.getOrDefault(highIndex, -1);
        }
    }

    // ============================================================
    // 🔷 Inden Basis
    // ============================================================

    /**
     * Builds the Inden basis for indene clusters with K species.
     *
     * For K even: {-K/2, -K/2+1, ..., -1, 1, 2, ..., K/2}
     * For K odd:  {-(K-1)/2, -(K-3)/2, ..., 0, ..., (K-3)/2, (K-1)/2}
     */
    public static double[] buildBasis(int K) {
        double[] basis = new double[K];

        if (K % 2 == 0) {
            int half = K / 2;
            for (int i = 0; i < half; i++) basis[i] = -half + i;
            for (int i = 0; i < half; i++) basis[half + i] = 1 + i;
        } else {
            int start = -((K - 1) / 2);
            for (int i = 0; i < K; i++) basis[i] = start + i;
        }

        return basis;
    }

    // ============================================================
    // 🔷 Index Encoding/Decoding
    // ============================================================

    /**
     * Decodes a single index into multi-index array in base-K representation.
     *
     * @param index  the encoded index
     * @param base   the base (K)
     * @param length the length of the output array
     * @param out    output array to fill (must be length long)
     */
    public static void decodeIndex(int index, int base, int length, int[] out) {
        for (int i = length - 1; i >= 0; i--) {
            out[i] = index % base;
            index /= base;
        }
    }

    /**
     * Encodes a multi-index array into a single index in base-K representation.
     *
     * @param arr  the multi-index array
     * @param base the base (K)
     * @return the encoded index
     */
    public static int encodeIndex(int[] arr, int base) {
        int idx = 0;
        for (int value : arr) {
            idx = idx * base + value;
        }
        return idx;
    }

    // ============================================================
    // 🔷 Configuration Generator
    // ============================================================

    /**
     * Generates random cluster configurations for regression.
     *
     * @param K        number of species
     * @param n        cluster size (e.g., 2 for pairs, 3 for triangles)
     * @param maxCount maximum number of configurations to generate
     * @return list of random configurations
     */
    public static List<int[]> generateConfigs(int K, int n, int maxCount) {
        List<int[]> configs = new ArrayList<>();
        Random rand = new Random(1);  // Deterministic seed for reproducibility

        for (int i = 0; i < maxCount; i++) {
            int[] c = new int[n];
            for (int j = 0; j < n; j++) {
                c[j] = rand.nextInt(K);
            }
            configs.add(c);
        }
        return configs;
    }

    // ============================================================
    // 🔷 Linear System Builder
    // ============================================================

    /**
     * Builds the linear system A·x = b for least-squares transformation.
     *
     * For each configuration:
     *   - Project it from high-order to low-order via SpeciesMapping
     *   - Use low-order ECI values as RHS (b)
     *   - Build basis rows for high-order correlation functions (A)
     *
     * @param configs   list of random configurations
     * @param J_low     ECI values in low-order system
     * @param K_low     number of species in low-order system
     * @param K_high    number of species in high-order system
     * @param n         cluster size
     * @param basisHigh high-order basis
     * @param mapping   species mapping from high to low
     * @param A         output matrix (rows = configs.size(), cols = K_high^n)
     * @param b         output RHS vector
     */
    public static void buildSystem(
            List<int[]> configs,
            double[] J_low,
            int K_low,
            int K_high,
            int n,
            double[] basisHigh,
            SpeciesMapping mapping,
            double[][] A,
            double[] b
    ) {

        int sizeHigh = (int) Math.pow(K_high, n);

        for (int i = 0; i < configs.size(); i++) {

            int[] sHigh = configs.get(i);
            int[] sLow = new int[n];

            // Project high-order configuration to low-order
            boolean valid = true;
            for (int k = 0; k < n; k++) {
                int lowIdx = mapping.getLowIndex(sHigh[k]);
                if (lowIdx == -1) {
                    valid = false;
                    break;
                }
                sLow[k] = lowIdx;
            }

            if (!valid) continue;

            // RHS: ECI value at the low-order configuration
            int lowIndex = encodeIndex(sLow, K_low);
            b[i] = J_low[lowIndex];

            // LHS: basis rows for all high-order configurations
            int[] m = new int[n];

            for (int j = 0; j < sizeHigh; j++) {

                decodeIndex(j, K_high, n, m);

                // Basis product: φ(m) = ∏_k basis_high[m_k]
                double phi = 1.0;
                for (int k = 0; k < n; k++) {
                    phi *= basisHigh[m[k]];
                }

                A[i][j] = phi;
            }
        }
    }

    // ============================================================
    // 🔷 Linear Solver (Normal Equations)
    // ============================================================

    /**
     * Solves the least-squares problem A·x = b using normal equations: (A^T·A)·x = A^T·b
     *
     * @param A matrix
     * @param b RHS vector
     * @return solution vector x
     */
    public static double[] solve(double[][] A, double[] b) {

        int rows = A.length;
        int cols = A[0].length;

        // Compute normal equations: A^T·A and A^T·b
        double[][] ATA = new double[cols][cols];
        double[] ATb = new double[cols];

        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < cols; j++) {
                for (int k = 0; k < rows; k++) {
                    ATA[i][j] += A[k][i] * A[k][j];
                }
            }
            for (int k = 0; k < rows; k++) {
                ATb[i] += A[k][i] * b[k];
            }
        }

        return LinearAlgebra.solve(ATA, ATb);
    }

    // ============================================================
    // 🔷 Main Transformation
    // ============================================================

    /**
     * Transforms ECIs from a lower-order basis to a higher-order basis.
     *
     * @param J_low    ECI values in lower-order system
     * @param K_low    number of species in lower-order system
     * @param K_high   number of species in higher-order system
     * @param n        cluster size (number of sites)
     * @param mapping  species mapping from high to low
     * @return ECI values in higher-order system
     */
    public static double[] transform(
            double[] J_low,
            int K_low,
            int K_high,
            int n,
            SpeciesMapping mapping
    ) {

        double[] basisHigh = buildBasis(K_high);

        // Number of random configurations for regression
        int numConfigs = 200;
        List<int[]> configs = generateConfigs(K_high, n, numConfigs);

        int sizeHigh = (int) Math.pow(K_high, n);

        double[][] A = new double[numConfigs][sizeHigh];
        double[] b = new double[numConfigs];

        buildSystem(
                configs,
                J_low,
                K_low,
                K_high,
                n,
                basisHigh,
                mapping,
                A,
                b
        );

        return solve(A, b);
    }

    // ============================================================
    // 🔷 Utility: Print Significant ECIs
    // ============================================================

    /**
     * Prints ECIs above a threshold (for debugging).
     *
     * @param J   ECI array
     * @param K   number of species
     * @param n   cluster size
     * @param tol threshold
     */
    public static void printSignificant(double[] J, int K, int n, double tol) {

        int[] idx = new int[n];

        for (int i = 0; i < J.length; i++) {
            if (Math.abs(J[i]) > tol) {
                decodeIndex(i, K, n, idx);
                System.out.println(Arrays.toString(idx) + " -> " + J[i]);
            }
        }
    }
}
