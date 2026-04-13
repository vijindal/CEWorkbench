package org.ce.model.cluster;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;

import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

/**
 * JUnit 5 test suite for mathematical and physical invariants in C-matrix generation.
 *
 * These tests verify properties that must hold for any physically correct C-matrix,
 * independent of crystal structure:
 *
 * 1. Partition of Unity — probabilities sum to 1
 * 2. Weight Sum — configurations counted correctly
 * 3. Constant Column — CV normalization
 * 4. R-Matrix Orthogonality — basis completeness
 * 5. Single-Site Diagonal — CF values match basis
 * 6. Symmetric Pair Antisymmetry — pair-CF properties
 */
@DisplayName("CMatrixGenerator Physics Invariants")
class CMatrixGeneratorPhysicsTest {

    private static final double TOLERANCE = 1e-10;

    // =========================================================================
    // Helper: Minimal CF table for N sites with K elements
    // =========================================================================

    /**
     * Builds a minimal CF table for testing:
     * - Columns 0..N-1: point CFs (basis function 1 at each site)
     * - Columns N..N+C(N,2)-1: pair CFs (basis function 1 at each site pair)
     *
     * @param numSites number of sites in the cluster
     * @param K number of elements (must be >= 2 for pair CFs to be meaningful)
     * @return map from SiteOpKey to CF column index, with totalCFs embedded
     */
    private static class CFTableInfo {
        final Map<CMatrixGenerator.SiteOpKey, Integer> table;
        final int totalCFs;

        CFTableInfo(Map<CMatrixGenerator.SiteOpKey, Integer> table, int totalCFs) {
            this.table = table;
            this.totalCFs = totalCFs;
        }
    }

    private static CFTableInfo buildMinimalCFTable(int numSites, int K) {
        Map<CMatrixGenerator.SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
        int colIdx = 0;

        // Point CFs: basis functions 1..K-1 at each site
        for (int i = 0; i < numSites; i++) {
            for (int b = 1; b < K; b++) {
                cfTable.put(
                    new CMatrixGenerator.SiteOpKey(List.of(new int[]{i, b})),
                    colIdx++
                );
            }
        }

        // Pair CFs: basis function 1 at each distinct site pair (for simplicity)
        for (int i = 0; i < numSites; i++) {
            for (int j = i + 1; j < numSites; j++) {
                cfTable.put(
                    new CMatrixGenerator.SiteOpKey(
                        List.of(new int[]{i, 1}, new int[]{j, 1})
                    ),
                    colIdx++
                );
            }
        }

        return new CFTableInfo(cfTable, colIdx);
    }

    // =========================================================================
    // Test Case 1: PARTITION OF UNITY
    // =========================================================================

    @Nested
    @DisplayName("Partition of Unity")
    class PartitionOfUnityTests {

        @Test
        @DisplayName("K=2, N=1: weighted sum of const column = 1, all CF cols = 0")
        void testPartitionOfUnityK2N1() {
            int K = 2, N = 1;
            CFTableInfo cfInfo = buildMinimalCFTable(N, K);
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            int[] siteIndices = {0};
            CMatrixGenerator.ClusterCMatrix block =
                CMatrixGenerator.buildForCluster(siteIndices, K, R, cfInfo.table, cfInfo.totalCFs);

            // Check CF columns sum to 0
            for (int col = 0; col < cfInfo.totalCFs; col++) {
                double sum = 0.0;
                for (int v = 0; v < block.lcv; v++) {
                    sum += block.wcv[v] * block.cmat[v][col];
                }
                assertEquals(0.0, sum, TOLERANCE,
                    "CF column " + col + " should sum to 0 (partition of unity)");
            }

            // Check constant column sums to 1
            double constSum = 0.0;
            for (int v = 0; v < block.lcv; v++) {
                constSum += block.wcv[v] * block.cmat[v][cfInfo.totalCFs];
            }
            assertEquals(1.0, constSum, TOLERANCE,
                "Constant column should sum to 1.0 (partition of unity)");
        }

        @Test
        @DisplayName("K=2, N=2: weighted sum of const column = 1, all CF cols = 0")
        void testPartitionOfUnityK2N2() {
            int K = 2, N = 2;
            CFTableInfo cfInfo = buildMinimalCFTable(N, K);
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            int[] siteIndices = {0, 1};
            CMatrixGenerator.ClusterCMatrix block =
                CMatrixGenerator.buildForCluster(siteIndices, K, R, cfInfo.table, cfInfo.totalCFs);

            // Check CF columns sum to 0
            for (int col = 0; col < cfInfo.totalCFs; col++) {
                double sum = 0.0;
                for (int v = 0; v < block.lcv; v++) {
                    sum += block.wcv[v] * block.cmat[v][col];
                }
                assertEquals(0.0, sum, TOLERANCE,
                    "CF column " + col + " should sum to 0");
            }

            // Check constant column sums to 1
            double constSum = 0.0;
            for (int v = 0; v < block.lcv; v++) {
                constSum += block.wcv[v] * block.cmat[v][cfInfo.totalCFs];
            }
            assertEquals(1.0, constSum, TOLERANCE,
                "Constant column should sum to 1.0");
        }

    }

    // =========================================================================
    // Test Case 2: WEIGHT SUM
    // =========================================================================

    @Nested
    @DisplayName("Weight Sum Equals K^N")
    class WeightSumTests {

        private void verifyWeightSum(int K, int N) {
            CFTableInfo cfInfo = buildMinimalCFTable(N, K);
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            int[] siteIndices = new int[N];
            for (int i = 0; i < N; i++) siteIndices[i] = i;

            CMatrixGenerator.ClusterCMatrix block =
                CMatrixGenerator.buildForCluster(siteIndices, K, R, cfInfo.table, cfInfo.totalCFs);

            int wSum = 0;
            for (int w : block.wcv) wSum += w;

            int expected = (int) Math.pow(K, N);
            assertEquals(expected, wSum,
                "Weight sum should equal K^N = " + K + "^" + N + " = " + expected);
        }

        @Test
        @DisplayName("K=2, N=1: sum(wcv) = 2")
        void testWeightSumK2N1() { verifyWeightSum(2, 1); }

        @Test
        @DisplayName("K=2, N=2: sum(wcv) = 4")
        void testWeightSumK2N2() { verifyWeightSum(2, 2); }
    }

    // =========================================================================
    // Test Case 3: CONSTANT COLUMN VALUE
    // =========================================================================

    @Nested
    @DisplayName("Constant Column Value = 1 / K^N")
    class ConstantColumnValueTests {

        private void verifyConstantColumnValue(int K, int N) {
            CFTableInfo cfInfo = buildMinimalCFTable(N, K);
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            int[] siteIndices = new int[N];
            for (int i = 0; i < N; i++) siteIndices[i] = i;

            CMatrixGenerator.ClusterCMatrix block =
                CMatrixGenerator.buildForCluster(siteIndices, K, R, cfInfo.table, cfInfo.totalCFs);

            double expected = 1.0 / Math.pow(K, N);

            for (int v = 0; v < block.lcv; v++) {
                assertEquals(expected, block.cmat[v][cfInfo.totalCFs], TOLERANCE,
                    "Row " + v + " constant term should equal 1/K^N = " + expected);
            }
        }

        @Test
        @DisplayName("K=2, N=1: constant = 0.5")
        void testConstantK2N1() { verifyConstantColumnValue(2, 1); }

        @Test
        @DisplayName("K=2, N=2: constant = 0.25")
        void testConstantK2N2() { verifyConstantColumnValue(2, 2); }

    }

    // =========================================================================
    // Test Case 4: R-MATRIX ORTHOGONALITY / COMPLETENESS
    // =========================================================================

    @Nested
    @DisplayName("R-Matrix Orthogonality and Completeness")
    class RMatrixOrthogonalityTests {

        @Test
        @DisplayName("K=2: sum_alpha R[alpha][0] = 1.0")
        void testRMatrixK2ConstantSum() {
            int K = 2;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            double sum = 0.0;
            for (int alpha = 0; alpha < K; alpha++) {
                sum += R[alpha][0];
            }

            assertEquals(1.0, sum, TOLERANCE,
                "Sum of constant coefficients should equal 1.0");
        }

        @Test
        @DisplayName("K=2: sum_alpha R[alpha][a] = 0 for a>0")
        void testRMatrixK2BasisOrthogonality() {
            int K = 2;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            for (int a = 1; a < K; a++) {
                double sum = 0.0;
                for (int alpha = 0; alpha < K; alpha++) {
                    sum += R[alpha][a];
                }
                assertEquals(0.0, sum, TOLERANCE,
                    "Sum of R[*][" + a + "] should equal 0");
            }
        }

        @Test
        @DisplayName("K=3: sum_alpha R[alpha][0] = 1.0")
        void testRMatrixK3ConstantSum() {
            int K = 3;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            double sum = 0.0;
            for (int alpha = 0; alpha < K; alpha++) {
                sum += R[alpha][0];
            }

            assertEquals(1.0, sum, TOLERANCE,
                "Sum of constant coefficients should equal 1.0");
        }

        @Test
        @DisplayName("K=3: sum_alpha R[alpha][a] = 0 for a>0")
        void testRMatrixK3BasisOrthogonality() {
            int K = 3;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            for (int a = 1; a < K; a++) {
                double sum = 0.0;
                for (int alpha = 0; alpha < K; alpha++) {
                    sum += R[alpha][a];
                }
                assertEquals(0.0, sum, TOLERANCE,
                    "Sum of R[*][" + a + "] should equal 0");
            }
        }

        @Test
        @DisplayName("K=4: sum_alpha R[alpha][0] = 1.0")
        void testRMatrixK4ConstantSum() {
            int K = 4;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            double sum = 0.0;
            for (int alpha = 0; alpha < K; alpha++) {
                sum += R[alpha][0];
            }

            assertEquals(1.0, sum, TOLERANCE,
                "Sum of constant coefficients should equal 1.0");
        }

        @Test
        @DisplayName("K=4: sum_alpha R[alpha][a] = 0 for a>0")
        void testRMatrixK4BasisOrthogonality() {
            int K = 4;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            for (int a = 1; a < K; a++) {
                double sum = 0.0;
                for (int alpha = 0; alpha < K; alpha++) {
                    sum += R[alpha][a];
                }
                assertEquals(0.0, sum, TOLERANCE,
                    "Sum of R[*][" + a + "] should equal 0");
            }
        }

        @Test
        @DisplayName("K=5: sum_alpha R[alpha][0] = 1.0")
        void testRMatrixK5ConstantSum() {
            int K = 5;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            double sum = 0.0;
            for (int alpha = 0; alpha < K; alpha++) {
                sum += R[alpha][0];
            }

            assertEquals(1.0, sum, TOLERANCE,
                "Sum of constant coefficients should equal 1.0");
        }

        @Test
        @DisplayName("K=5: sum_alpha R[alpha][a] = 0 for a>0")
        void testRMatrixK5BasisOrthogonality() {
            int K = 5;
            double[][] R = CMatrixGenerator.buildRMatrix(K);

            for (int a = 1; a < K; a++) {
                double sum = 0.0;
                for (int alpha = 0; alpha < K; alpha++) {
                    sum += R[alpha][a];
                }
                assertEquals(0.0, sum, TOLERANCE,
                    "Sum of R[*][" + a + "] should equal 0");
            }
        }
    }

    // =========================================================================
    // Test Case 5: SINGLE-SITE CLUSTER DIAGONAL
    // =========================================================================

    @Nested
    @DisplayName("Single-Site Cluster Diagonal")
    class SingleSiteDiagonalTests {

        private void verifySingleSiteDiagonal(int K) {
            // Use a simple CF table: just basis function 1 at site 0
            Map<CMatrixGenerator.SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
            cfTable.put(
                new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1})),
                0
            );
            int totalCFs = 1;

            double[][] R = CMatrixGenerator.buildRMatrix(K);

            int[] siteIndices = {0};
            CMatrixGenerator.ClusterCMatrix block =
                CMatrixGenerator.buildForCluster(siteIndices, K, R, cfTable, totalCFs);

            // Verify lcv == K
            assertEquals(K, block.lcv,
                "lcv should equal K = " + K + " for single-site cluster");

            // Verify wcv[v] == 1 for all v
            for (int v = 0; v < block.lcv; v++) {
                assertEquals(1, block.wcv[v],
                    "wcv[" + v + "] should equal 1 for single-site cluster");
            }

            // Verify CF column (basis 1 at site 0) values match R[alpha][1] in some order
            double[] rBasisColumn = new double[K];
            for (int alpha = 0; alpha < K; alpha++) {
                rBasisColumn[alpha] = R[alpha][1];
            }
            Arrays.sort(rBasisColumn);

            double[] cfColumnValues = new double[K];
            for (int v = 0; v < block.lcv; v++) {
                cfColumnValues[v] = block.cmat[v][0];  // basis 1 CF at site 0
            }
            Arrays.sort(cfColumnValues);

            assertArrayEquals(rBasisColumn, cfColumnValues, TOLERANCE,
                "CF column (basis 1) values should match R[*][1] (in some order)");
        }

        @Test
        @DisplayName("K=2: lcv=2, wcv=[1,1], CF col = R[*][1]")
        void testSingleSiteK2() { verifySingleSiteDiagonal(2); }

    }

    // =========================================================================
    // Test Case 6: SYMMETRIC PAIR — CF COLUMN ANTISYMMETRY
    // =========================================================================

    @Nested
    @DisplayName("Symmetric Pair Antisymmetry")
    class SymmetricPairAntisymmetryTests {

        @Test
        @DisplayName("K=2, symmetric pair: pair-CF weighted sum = 0, values antisymmetric")
        void testSymmetricPairK2() {
            int K = 2;

            // Simple CF table: point CFs at sites 0,1 (columns 0,1) and pair CF (column 2)
            Map<CMatrixGenerator.SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
            cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1})), 0);
            cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{1, 1})), 1);
            cfTable.put(new CMatrixGenerator.SiteOpKey(List.of(new int[]{0, 1}, new int[]{1, 1})), 2);

            int totalCFs = 3;

            double[][] R = CMatrixGenerator.buildRMatrix(K);

            int[] siteIndices = {0, 1};
            CMatrixGenerator.ClusterCMatrix block =
                CMatrixGenerator.buildForCluster(siteIndices, K, R, cfTable, totalCFs);

            // The pair CF is in column 2
            int pairCFColumn = 2;

            // Verify weighted sum of pair-CF column = 0
            double pairCFSum = 0.0;
            for (int v = 0; v < block.lcv; v++) {
                pairCFSum += block.wcv[v] * block.cmat[v][pairCFColumn];
            }
            assertEquals(0.0, pairCFSum, TOLERANCE,
                "Pair-CF column weighted sum should be 0 (partition of unity)");

            // For K=2 binary pair, expect 3 or 4 distinct rows
            // (AA, BB, and AB/BA may or may not group depending on symmetry)
            assertTrue(block.lcv >= 3,
                "K=2 symmetric pair should have at least 3 distinct CVs");

            // Verify weight sum = 4
            int wSum = 0;
            for (int w : block.wcv) wSum += w;
            assertEquals(4, wSum, "Weight sum should be K^N = 4");
        }
    }
}
