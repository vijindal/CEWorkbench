package org.ce.model.cluster;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;

import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

/**
 * JUnit 5 test suite for CMatrixGenerator.
 *
 * Tests the numerical C-matrix generation pipeline:
 * 1. Basis construction (K=2,3,4)
 * 2. R-matrix generation and validation
 * 3. Configuration enumeration
 * 4. Configuration expansion to coefficient vectors
 * 5. Clustering of identical rows and weight tallying
 */
@DisplayName("CMatrixGenerator")
class CMatrixGeneratorTest {

    private static final double TOLERANCE = 1e-10;

    // =========================================================================
    // Test Case 1: buildBasis(K) for K=2,3,4
    // =========================================================================

    @Test
    @DisplayName("buildBasis(2) returns {-1, 1}")
    void testBuildBasisK2() {
        double[] basis = CMatrixGenerator.buildBasis(2);
        assertEquals(2, basis.length);
        assertEquals(-1.0, basis[0], TOLERANCE);
        assertEquals(1.0, basis[1], TOLERANCE);
    }

    @Test
    @DisplayName("buildBasis(3) returns {-1, 0, 1}")
    void testBuildBasisK3() {
        double[] basis = CMatrixGenerator.buildBasis(3);
        assertEquals(3, basis.length);
        assertEquals(-1.0, basis[0], TOLERANCE);
        assertEquals(0.0, basis[1], TOLERANCE);
        assertEquals(1.0, basis[2], TOLERANCE);
    }

    @Test
    @DisplayName("buildBasis(4) returns {-2, -1, 1, 2}")
    void testBuildBasisK4() {
        double[] basis = CMatrixGenerator.buildBasis(4);
        assertEquals(4, basis.length);
        assertEquals(-2.0, basis[0], TOLERANCE);
        assertEquals(-1.0, basis[1], TOLERANCE);
        assertEquals(1.0, basis[2], TOLERANCE);
        assertEquals(2.0, basis[3], TOLERANCE);
    }

    // =========================================================================
    // Test Case 2: buildRMatrix(K=2) — verify exact values
    // =========================================================================

    @Test
    @DisplayName("buildRMatrix(2) produces [[0.5, -0.5], [0.5, 0.5]]")
    void testBuildRMatrixK2() {
        double[][] R = CMatrixGenerator.buildRMatrix(2);
        assertEquals(2, R.length);
        assertEquals(2, R[0].length);

        // R[0] = [0.5, -0.5]  (alpha=0, element A)
        assertEquals(0.5, R[0][0], TOLERANCE);
        assertEquals(-0.5, R[0][1], TOLERANCE);

        // R[1] = [0.5, 0.5]   (alpha=1, element B)
        assertEquals(0.5, R[1][0], TOLERANCE);
        assertEquals(0.5, R[1][1], TOLERANCE);
    }

    // =========================================================================
    // Test Case 3: buildRMatrix(K=3) — verify R*V = I
    // =========================================================================

    @Test
    @DisplayName("buildRMatrix(3) satisfies R * V = I")
    void testBuildRMatrixK3VandermondeInverse() {
        int K = 3;
        double[] basis = CMatrixGenerator.buildBasis(K);
        double[][] R = CMatrixGenerator.buildRMatrix(K);

        // Build Vandermonde matrix V
        double[][] V = new double[K][K];
        for (int row = 0; row < K; row++) {
            for (int col = 0; col < K; col++) {
                V[row][col] = (row == 0) ? 1.0 : Math.pow(basis[col], row);
            }
        }

        // Compute R * V
        double[][] product = multiplyMatrices(R, V, K);

        // Verify product = I
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                double expected = (i == j) ? 1.0 : 0.0;
                assertEquals(expected, product[i][j], TOLERANCE,
                    "R*V[" + i + "][" + j + "] should be " + expected);
            }
        }
    }

    // =========================================================================
    // Test Case 4: allConfigurations(N, K) — count and uniqueness
    // =========================================================================

    @Test
    @DisplayName("allConfigurations(1, 2) produces K^N = 2 configs")
    void testAllConfigurationsN1K2() {
        List<int[]> configs = CMatrixGenerator.allConfigurations(1, 2);
        assertEquals(2, configs.size());

        // Expect [0] and [1]
        assertArrayEquals(new int[]{0}, configs.get(0));
        assertArrayEquals(new int[]{1}, configs.get(1));
    }

    @Test
    @DisplayName("allConfigurations(2, 2) produces K^N = 4 configs")
    void testAllConfigurationsN2K2() {
        List<int[]> configs = CMatrixGenerator.allConfigurations(2, 2);
        assertEquals(4, configs.size());

        // Expect all 2-element binary sequences
        Set<String> seen = new HashSet<>();
        for (int[] cfg : configs) {
            assertEquals(2, cfg.length);
            assertTrue(cfg[0] >= 0 && cfg[0] < 2);
            assertTrue(cfg[1] >= 0 && cfg[1] < 2);
            seen.add(Arrays.toString(cfg));
        }
        assertEquals(4, seen.size(), "All 4 configs should be unique");
    }

    @Test
    @DisplayName("allConfigurations(3, 3) produces K^N = 27 configs")
    void testAllConfigurationsN3K3() {
        List<int[]> configs = CMatrixGenerator.allConfigurations(3, 3);
        assertEquals(27, configs.size());

        // Verify each config has the right length and valid element indices
        Set<String> seen = new HashSet<>();
        for (int[] cfg : configs) {
            assertEquals(3, cfg.length);
            for (int elem : cfg) {
                assertTrue(elem >= 0 && elem < 3);
            }
            seen.add(Arrays.toString(cfg));
        }
        assertEquals(27, seen.size(), "All 27 configs should be unique");
    }

    // =========================================================================
    // Test Case 5: expandConfiguration() for K=2, N=1 single-site cluster
    // =========================================================================

    @Test
    @DisplayName("expandConfiguration(K=2, N=1): config [0] yields [-0.5, 0.5] (CF0, const)")
    void testExpandConfigurationK2N1ElementA() {
        int K = 2;
        double[][] R = CMatrixGenerator.buildRMatrix(K);

        // Single-site cluster: siteIndices = [0]
        int[] siteIndices = {0};
        int[] config = {0};  // element A at site 0

        // CF table with 1 CF: basis function 1 at site 0
        Map<CMatrixGenerator.SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
        cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1})), 0);

        int totalCFs = 1;

        double[] row = CMatrixGenerator.expandConfiguration(siteIndices, config, R, cfTable, totalCFs);

        // Expected: [-0.5, 0.5]  (CF0 coeff, const)
        assertEquals(2, row.length);  // 1 CF + 1 const
        assertEquals(-0.5, row[0], TOLERANCE, "CF0 coefficient");
        assertEquals(0.5, row[1], TOLERANCE, "constant term");
    }

    @Test
    @DisplayName("expandConfiguration(K=2, N=1): config [1] yields [+0.5, 0.5] (CF0, const)")
    void testExpandConfigurationK2N1ElementB() {
        int K = 2;
        double[][] R = CMatrixGenerator.buildRMatrix(K);

        int[] siteIndices = {0};
        int[] config = {1};  // element B at site 0

        Map<CMatrixGenerator.SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
        cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1})), 0);

        int totalCFs = 1;

        double[] row = CMatrixGenerator.expandConfiguration(siteIndices, config, R, cfTable, totalCFs);

        // Expected: [+0.5, 0.5]  (CF0 coeff, const)
        assertEquals(2, row.length);
        assertEquals(0.5, row[0], TOLERANCE, "CF0 coefficient");
        assertEquals(0.5, row[1], TOLERANCE, "constant term");
    }

    // =========================================================================
    // Test Case 6: buildForCluster(K=2, binary pair)
    // =========================================================================

    @Test
    @DisplayName("buildForCluster(K=2, 2-site): lcv=4, sum(wcv)=4, const=0.25 each row")
    void testBuildForClusterK2BinaryPair() {
        int K = 2;
        double[][] R = CMatrixGenerator.buildRMatrix(K);

        // 2-site asymmetric pair: each site maps to different CFs
        int[] siteIndices = {0, 1};

        // CF table: 3 CFs (point 0, point 1, pair)
        Map<CMatrixGenerator.SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
        cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1})), 0);
        cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{1, 1})), 1);
        cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1}, new int[]{1, 1})), 2);

        int totalCFs = 3;

        CMatrixGenerator.ClusterCMatrix block =
            CMatrixGenerator.buildForCluster(siteIndices, K, R, cfTable, totalCFs);

        // Verify lcv
        assertEquals(4, block.lcv, "lcv should be 4 for 2-site asymmetric pair");

        // Verify sum of wcv equals K^N = 4
        int wSum = 0;
        for (int w : block.wcv) wSum += w;
        assertEquals(4, wSum, "sum(wcv) should equal K^N = 4");

        // Verify each row's constant term equals 0.25
        int constCol = totalCFs;  // column 3
        for (int v = 0; v < block.lcv; v++) {
            assertEquals(0.25, block.cmat[v][constCol], TOLERANCE,
                "Row " + v + " constant term should be 0.25");
        }
    }

    // =========================================================================
    // Test Case 7: buildForCluster(K=3, ternary single-site)
    // =========================================================================

    @Test
    @DisplayName("buildForCluster(K=3, 1-site): lcv=3, sum(wcv)=3, partition of unity")
    void testBuildForClusterK3TernarySingleSite() {
        int K = 3;
        double[][] R = CMatrixGenerator.buildRMatrix(K);

        // Single-site cluster
        int[] siteIndices = {0};

        // CF table: 2 CFs for a single-site ternary cluster
        // For K=3, we need basis functions 1 and 2 to cover the expansion
        Map<CMatrixGenerator.SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
        cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1})), 0);
        cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 2})), 1);

        int totalCFs = 2;

        CMatrixGenerator.ClusterCMatrix block =
            CMatrixGenerator.buildForCluster(siteIndices, K, R, cfTable, totalCFs);

        // Verify lcv = 3 (one row per element)
        assertEquals(3, block.lcv, "lcv should equal K = 3");

        // Verify sum of wcv equals K^N = 3
        int wSum = 0;
        for (int w : block.wcv) wSum += w;
        assertEquals(3, wSum, "sum(wcv) should equal K^N = 3");

        // Verify partition of unity: sum of rows (weighted by wcv) = [0, 0, ..., 1.0]
        // The weighted constant terms should sum to 1.0 (partition of unity)
        int constCol = totalCFs;  // column 2 (last column)
        double weightedConstSum = 0.0;
        for (int v = 0; v < block.lcv; v++) {
            weightedConstSum += block.wcv[v] * block.cmat[v][constCol];
        }
        assertEquals(1.0, weightedConstSum, TOLERANCE,
            "Weighted sum of constant terms should equal 1.0");
    }

    // =========================================================================
    // Helper: matrix multiplication
    // =========================================================================

    /**
     * Multiply two n×n matrices: C = A * B.
     */
    private static double[][] multiplyMatrices(double[][] A, double[][] B, int n) {
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += A[i][k] * B[k][j];
                }
                C[i][j] = sum;
            }
        }
        return C;
    }
}
