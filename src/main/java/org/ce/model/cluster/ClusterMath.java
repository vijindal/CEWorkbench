package org.ce.model.cluster;

import static org.ce.model.cluster.ClusterPrimitives.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
import org.ce.model.cluster.SpaceGroup.SymmetryOperation;

/**
 * Unified container for cluster-related mathematical operations and rules.
 *
 * <p>Consolidates four core CVM computation stages into a single pipeline:</p>
 * <ol>
 *   <li><b>Basis & R-Matrix:</b> Symmetric integer sequences and Inden (1992) R-matrices.</li>
 *   <li><b>P-Operator Rules:</b> Expansion coefficients for site-occupation operators.</li>
 *   <li><b>Nij Table:</b> Geometric containment counts (cluster hierarchy).</li>
 *   <li><b>Kikuchi-Baker:</b> Entropy coefficients to prevent double-counting.</li>
 * </ol>
 */
public final class ClusterMath {

    private static final double EPS = 1e-12;
    private static final double DELTA = 1e-6;

    private ClusterMath() {}

    // =========================================================================
    // 1. Basis & R-Matrix (formerly RMatrixCalculator)
    // =========================================================================

    /**
     * Builds the Inden (1992) R-matrix for the given number of elements.
     * The matrix is the inverse of a Vandermonde-like matrix from a symmetric basis.
     */
    public static double[][] buildRMatrix(int numElements) {
        if (numElements < 2) {
            throw new IllegalArgumentException("numElements must be >= 2");
        }

        double[] basis = buildBasis(numElements);
        double[][] matM = new double[numElements][numElements];

        for (int i = 0; i < numElements; i++) {
            for (int j = 0; j < numElements; j++) {
                matM[i][j] = (i == 0) ? 1.0 : Math.pow(basis[j], i);
            }
        }

        return LinearAlgebra.invert(matM);
    }

    /**
     * Returns the symmetric integer basis sequence for K components.
     * Even K: {-K/2, ..., -1, 1, ..., K/2}. Odd K: {-(K-1)/2, ..., 0, ..., (K-1)/2}.
     */
    public static double[] buildBasis(int numElements) {
        double[] basis = new double[numElements];
        if (numElements % 2 == 0) {
            int half = numElements / 2;
            for (int i = 0; i < half; i++) basis[i] = -half + i;
            for (int i = 0; i < half; i++) basis[half + i] = 1 + i;
        } else {
            int start = -((numElements - 1) / 2);
            for (int i = 0; i < numElements; i++) basis[i] = start + i;
        }
        return basis;
    }

    // =========================================================================
    // 2. P-Operator Rules (formerly PRules)
    // =========================================================================

    /**
     * Immutable container for p-operator expansion coefficients.
     */
    public static final class PRules {
        private final int numSites;
        private final int numElements;
        private final double[][] rMat;

        private PRules(int numSites, int numElements, double[][] rMat) {
            this.numSites = numSites;
            this.numElements = numElements;
            this.rMat = copyMatrix(rMat);
        }

        public static PRules build(int numSites, int numElements) {
            if (numSites <= 0) throw new IllegalArgumentException("numSites must be > 0");
            if (numElements < 2) throw new IllegalArgumentException("numElements must be >= 2");
            return new PRules(numSites, numElements, buildRMatrix(numElements));
        }

        public int getNumSites() { return numSites; }
        public int getNumElements() { return numElements; }

        /** Returns coefficients for p(site, element) in [1, s1, s2, ...] basis. */
        public double[] coefficientsFor(int siteIndex, int elementIndex) {
            if (siteIndex < 0 || siteIndex >= numSites) throw new IllegalArgumentException("siteIndex out of range");
            if (elementIndex < 0 || elementIndex >= numElements) throw new IllegalArgumentException("elementIndex out of range");
            return Arrays.copyOf(rMat[elementIndex], numElements);
        }

        public double[][] getRMatrix() { return copyMatrix(rMat); }

        private static double[][] copyMatrix(double[][] src) {
            double[][] copy = new double[src.length][];
            for (int i = 0; i < src.length; i++) copy[i] = Arrays.copyOf(src[i], src[i].length);
            return copy;
        }
    }

    // =========================================================================
    // 3. Nij Table (formerly NijTableCalculator)
    // =========================================================================

    /**
     * Computes the Nij containment table: nij[i][j] is the number of times cluster type j
     * appears as a geometrically distinct sub-cluster inside cluster type i.
     */
    public static int[][] computeNijTable(List<Cluster> disClusCoordList, List<List<Cluster>> disClusOrbitList) {
        int numClus = disClusCoordList.size();
        int[][] nijTable = new int[numClus][numClus];

        for (int i = 0; i < numClus; i++) {
            List<Cluster> subClusters = ClusterCFIdentificationPipeline.genSubClusCoord(disClusCoordList.get(i));
            for (int k = subClusters.size() - 1; k >= 0; k--) {
                Cluster subClus = subClusters.get(k);
                int subSize = subClus.getAllSites().size();
                if (subSize == 0) continue;

                for (int j = i; j < numClus; j++) {
                    if (subSize != disClusCoordList.get(j).getAllSites().size()) continue;
                    if (isContained(disClusOrbitList.get(j), subClus)) {
                        nijTable[i][j]++;
                    }
                }
            }
        }
        return nijTable;
    }

    public static void printNijTable(int[][] nijTable) {
        int n = nijTable.length;
        System.out.println("[ClusterMath] nijTable (" + n + "Ã—" + n + "):");
        System.out.print("       ");
        for (int j = 0; j < n; j++) System.out.printf(" j=%-3d", j);
        System.out.println();
        for (int i = 0; i < n; i++) {
            System.out.printf("  i=%-3d", i);
            for (int j = 0; j < n; j++) System.out.printf(" %-5d", nijTable[i][j]);
            System.out.println();
        }
    }

    // =========================================================================
    // 4. Kikuchi-Baker (formerly KikuchiBakerCalculator)
    // =========================================================================

    /**
     * Computes KB entropy coefficients via the inclusion-exclusion recurrence:
     * kb[j] = 1 - Î£_{i < j} (m[i] * nij[i][j] * kb[i]) / m[j]
     */
    public static double[] computeKikuchiBaker(double[] multiplicities, int[][] nijTable) {
        int tcdis = multiplicities.length;
        if (nijTable.length != tcdis || nijTable[0].length != tcdis) {
            throw new IllegalArgumentException("Dimension mismatch between multiplicities and Nij table.");
        }

        double[] kb = new double[tcdis];
        for (int j = 0; j < tcdis; j++) {
            double sumTerm = 0.0;
            for (int i = 0; i < j; i++) {
                sumTerm += multiplicities[i] * nijTable[i][j] * kb[i];
            }
            kb[j] = (multiplicities[j] - sumTerm) / multiplicities[j];
        }
        return kb;
    }

    public static void printKikuchiBaker(double[] multiplicities, double[] kb) {
        System.out.println("[ClusterMath] Kikuchi-Baker Coefficients:");
        System.out.printf("  %-6s %-12s %-14s%n", "Type", "mhdis", "KB coeff");
        for (int t = 0; t < kb.length; t++) {
            System.out.printf("  t=%-4d %-12.4f %-14.8f%n", t, multiplicities[t], kb[t]);
        }
    }

    // =========================================================================
    // 5. Geometric Utilities (formerly ClusterUtils)
    // =========================================================================

    public static boolean isTranslated(Cluster c1, Cluster c2) {
        if (c1.getSublattices().size() != c2.getSublattices().size())
            return false;

        Set<Position> diffSet = new HashSet<>();
        for (int i = 0; i < c1.getSublattices().size(); i++) {
            Sublattice sub1 = c1.getSublattices().get(i);
            Sublattice sub2 = c2.getSublattices().get(i);
            List<Site> s1 = sub1.getSites();
            List<Site> s2 = sub2.getSites();
            if (s1.size() != s2.size()) return false;
            for (int j = 0; j < s1.size(); j++) {
                Site site1 = s1.get(j);
                Site site2 = s2.get(j);
                if (!site1.getSymbol().equals(site2.getSymbol())) return false;
                diffSet.add(site2.getPosition().subtract(site1.getPosition()));
            }
        }
        if (diffSet.size() > 1) return false;
        if (diffSet.isEmpty())  return true;
        Position d = diffSet.iterator().next();
        return isIntegerShift(d.getX()) && isIntegerShift(d.getY()) && isIntegerShift(d.getZ());
    }

    private static boolean isIntegerShift(double value) {
        return Math.abs(value - Math.round(value)) < DELTA;
    }

    public static boolean isContained(List<Cluster> orbit, Cluster cluster) {
        for (Cluster existing : orbit) {
            if (isTranslated(existing, cluster)) return true;
        }
        return false;
    }

    public static List<Cluster> generateOrbit(Cluster cluster, List<SymmetryOperation> spaceGroup) {
        List<Cluster> orbit = new ArrayList<>();
        for (SymmetryOperation op : spaceGroup) {
            Cluster transformed = op.applyToCluster(cluster);
            if (!isContained(orbit, transformed))
                orbit.add(transformed);
        }
        return orbit;
    }

    // =========================================================================
    // 6. List/Summary operations
    // =========================================================================

    public static int countNonEmpty(ClusCoordListData data) {
        if (data == null || data.getClusCoordList().isEmpty()) return 0;
        List<Cluster> list = data.getClusCoordList();
        int count = 0;
        for (Cluster c : list) {
            if (c.getAllSites().size() > 0) count++;
        }
        return count;
    }

    public static void printOrbitDebug(Cluster cluster, List<SymmetryOperation> spaceGroup) {
        List<Cluster> orbit = generateOrbit(cluster, spaceGroup);
        System.out.println("[ClusterMath.generateOrbit]");
        System.out.println("  seed cluster       : " + cluster);
        System.out.println("  space group ops    : " + spaceGroup.size());
        System.out.println("  generated orbit size : " + orbit.size());
        for (int i = 0; i < orbit.size(); i++) {
            System.out.println("  orbit[" + i + "] : " + orbit.get(i));
        }
    }
}
