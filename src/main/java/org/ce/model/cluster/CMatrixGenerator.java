package org.ce.model.cluster;

import java.util.*;
import java.util.function.Consumer;
import static org.ce.model.cluster.ClusterPrimitives.*;
import static org.ce.model.cluster.ClusterKeys.*;

/**
 * C-Matrix Generator — pure numerical implementation, general K elements.
 *
 * No symbolic math. Every operation is a concrete numerical array or map lookup.
 *
 * INPUTS  (all computed upstream before calling this class):
 *   K                — number of chemical elements
 *   clusters         — ordered list of clusters, each described as an array
 *                      of site indices (into the global siteList)
 *   cfTable          — the CF lookup table: maps a (siteIndex[], basisIndex[])
 *                      pattern to a CF column number (built from CF identification)
 *   totalCFs         — total number of CF columns in the C-matrix
 *
 * OUTPUTS for each cluster:
 *   cmat[row][col]   — coefficient of CF column 'col' in CV row 'row'
 *                      last column (index totalCFs) is the constant term
 *   lcv              — number of distinct CV rows (= cmat.length)
 *   wcv[row]         — number of configurations that produced CV row 'row'
 */
public final class CMatrixGenerator {

    private CMatrixGenerator() {}

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }

    // =========================================================================
    // STEP 1 — Build the R-matrix (basis-change table)
    // =========================================================================
    //
    // For K elements, define a symmetric integer sequence called the "basis":
    //
    //   K even: basis = { -K/2, -(K/2)+1, ..., -1,   1, 2, ..., K/2 }
    //   K odd:  basis = { -(K-1)/2, ..., -1, 0, 1, ..., (K-1)/2 }
    //
    // Example K=2: basis = {-1,  1}
    // Example K=3: basis = {-1,  0,  1}
    // Example K=4: basis = {-2, -1,  1,  2}
    //
    // Build a K×K Vandermonde-like matrix V where:
    //   V[row][col] = basis[col]^row     (row=0 gives all 1s)
    //
    // The R-matrix is the numerical inverse of V:
    //   R = V^{-1}
    //
    // Meaning of R:
    //   R[alpha][a] is the coefficient of basis-function a in the expansion
    //   of "atom alpha occupies this site".
    //
    //   p(site, alpha) = R[alpha][0]*1 + R[alpha][1]*s1(site)
    //                                  + R[alpha][2]*s2(site) + ...
    //
    //   where s0 = 1 (constant), s1, s2, ... are the site-operator basis functions.

    /**
     * Returns basis[] of length K.
     */
    public static double[] buildBasis(int K) {
        double[] b = new double[K];
        if (K % 2 == 0) {
            int half = K / 2;
            for (int i = 0; i < half; i++) b[i]        = -half + i;   // -K/2 .. -1
            for (int i = 0; i < half; i++) b[half + i] =  1    + i;   //  1   .. K/2
        } else {
            int start = -((K - 1) / 2);
            for (int i = 0; i < K; i++) b[i] = start + i;             // symmetric around 0
        }
        return b;
    }

    /**
     * Returns the K×K R-matrix.
     * R[alpha][a] = coefficient of basis-function a for element alpha.
     * a=0 is the constant term; a=1,2,...,K-1 are the site-operator powers.
     */
    public static double[][] buildRMatrix(int K) {
        return buildRMatrix(K, null);
    }

    /**
     * Internal overload with diagnostic output to sink.
     */
    public static double[][] buildRMatrix(int K, Consumer<String> sink) {
        emit(sink, "\n[STAGE 1] Building R-Matrix (K=" + K + ")");
        double[] basis = buildBasis(K);
        emit(sink, "  Basis values: " + Arrays.toString(basis));

        // Build Vandermonde matrix V[row][col] = basis[col]^row
        double[][] V = new double[K][K];
        for (int row = 0; row < K; row++) {
            for (int col = 0; col < K; col++) {
                V[row][col] = (row == 0) ? 1.0 : Math.pow(basis[col], row);
            }
        }

        emit(sink, "  Vandermonde Matrix V:");
        for (int i = 0; i < K; i++) {
            StringBuilder sb = new StringBuilder("    [");
            for (int j = 0; j < K; j++) {
                sb.append(String.format(" %8.4f", V[i][j]));
            }
            sb.append(" ]");
            emit(sink, sb.toString());
        }

        // R = V^{-1}  using Gaussian elimination with partial pivoting
        double[][] R = invertMatrix(V, K);

        emit(sink, "  Resulting R-Matrix (V^-1):");
        for (int i = 0; i < K; i++) {
            StringBuilder sb = new StringBuilder("    [");
            for (int j = 0; j < K; j++) {
                sb.append(String.format(" %8.4f", R[i][j]));
            }
            sb.append(" ]");
            emit(sink, sb.toString());
        }

        return R;
    }

    // =========================================================================
    // STEP 2 — Build the CF lookup table
    // =========================================================================
    //
    // The CF lookup table maps a sorted list of (siteIndex, basisFunctionIndex)
    // pairs to a column number in the C-matrix.
    //
    // This table is built from the CF identification stage. Each CF is defined
    // by a representative cluster of sites decorated with basis-function indices
    // (e.g. site 2 carries s1, site 5 carries s2 → CF column 7).
    //
    // Key rule: two site-operator products map to the SAME CF column if and
    // only if they are related by the symmetry of the highest-symmetric phase
    // (same geometry = same inter-site distances in the same order, AND same
    // basis-function decoration pattern).
    //
    // In Java: represent each product as int[][] where each row is {siteIndex, basisIndex},
    // sort rows canonically, compute a geometry signature, look up the column.
    //
    // This table is an INPUT to this class. See buildCFTable() below for how
    // to construct it from the CF identification result.

    /**
     * Canonical key for a site-operator product.
     * Sorted by siteIndex then basisIndex so order of multiplication doesn't matter.
     */
    public static final class SiteOpKey {
        final int[] siteIndices;    // sorted
        final int[] basisIndices;   // aligned with siteIndices

        SiteOpKey(List<int[]> ops) {
            // ops: each element is {siteIndex, basisIndex}
            List<int[]> sorted = new ArrayList<>(ops);
            sorted.sort((a, b) -> a[0] != b[0] ? Integer.compare(a[0], b[0])
                                                : Integer.compare(a[1], b[1]));
            siteIndices  = new int[sorted.size()];
            basisIndices = new int[sorted.size()];
            for (int i = 0; i < sorted.size(); i++) {
                siteIndices[i]  = sorted.get(i)[0];
                basisIndices[i] = sorted.get(i)[1];
            }
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof SiteOpKey)) return false;
            SiteOpKey k = (SiteOpKey) o;
            return Arrays.equals(siteIndices, k.siteIndices)
                && Arrays.equals(basisIndices, k.basisIndices);
        }

        @Override
        public int hashCode() {
            return 31 * Arrays.hashCode(siteIndices) + Arrays.hashCode(basisIndices);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder("[");
            for (int i = 0; i < siteIndices.length; i++) {
                sb.append("(site=").append(siteIndices[i])
                  .append(",b=").append(basisIndices[i]).append(")");
                if (i < siteIndices.length - 1) sb.append(", ");
            }
            return sb.append("]").toString();
        }
    }

    // =========================================================================
    // STEP 3 — Expand one configuration into a coefficient vector
    // =========================================================================
    //
    // Given:
    //   siteIndices[s] = which global site is position s in this cluster
    //   config[s]      = which element (0..K-1) sits at position s
    //   R              = the K×K R-matrix from Step 1
    //   cfTable        = map from SiteOpKey → CF column number
    //   totalCFs       = total number of CF columns
    //
    // We compute a double[] row of length (totalCFs + 1) where:
    //   row[col]       = coefficient of CF 'col' in this CV
    //   row[totalCFs]  = constant term
    //
    // Algorithm:
    //   We maintain a running polynomial as Map<SiteOpKey, double>.
    //   Initially: { empty_key → 1.0 }   (the number 1)
    //
    //   For each site s:
    //     The occupation operator for (siteIndices[s], config[s]) expands as:
    //       p = R[config[s]][0]*1 + R[config[s]][1]*s1(siteIndices[s])
    //                             + R[config[s]][2]*s2(siteIndices[s]) + ...
    //
    //     Multiply every existing term in the polynomial by this expansion:
    //       for each existing (key, coeff):
    //         for each basis power a = 0..K-1:
    //           new_ops = key.ops + (a>0 ? {siteIndices[s], a} : nothing)
    //           new_key = SiteOpKey(new_ops)
    //           nextPoly[new_key] += coeff * R[config[s]][a]
    //
    //   After all sites, map each monomial to a CF column:
    //     empty monomial → constant term
    //     non-empty      → cfTable.get(key) gives the CF column

    /**
     * Computes the coefficient row for one configuration of one cluster.
     *
     * @param siteIndices  global site index of each cluster site (length = cluster size N)
     * @param config       element index (0..K-1) at each cluster site (length N)
     * @param R            K×K R-matrix from buildRMatrix(K)
     * @param cfTable      map from SiteOpKey to CF column index (0-based)
     * @param totalCFs     total number of CF columns (constant goes in col totalCFs)
     * @return             double[] of length (totalCFs + 1)
     */
    public static double[] expandConfiguration(
            int[]                    siteIndices,
            int[]                    config,
            double[][]               R,
            Map<SiteOpKey, Integer>  cfTable,
            int                      totalCFs) {

        int N = siteIndices.length;
        int K = R.length;

        // Polynomial: monomial → accumulated coefficient
        // A monomial is represented as an ordered list of {siteIndex, basisIndex} pairs.
        // We use List<int[]> as the mutable form and SiteOpKey as the lookup key.
        Map<SiteOpKey, Double> poly = new LinkedHashMap<>();
        poly.put(new SiteOpKey(Collections.emptyList()), 1.0);  // start: constant 1

        for (int s = 0; s < N; s++) {
            int globalSite  = siteIndices[s];
            int elementIdx  = config[s];
            double[] coeffs = R[elementIdx];   // length K: coeffs[0]=constant, coeffs[a]=coeff of s_a

            Map<SiteOpKey, Double> next = new LinkedHashMap<>();

            for (Map.Entry<SiteOpKey, Double> entry : poly.entrySet()) {
                SiteOpKey existingKey  = entry.getKey();
                double    existingCoeff = entry.getValue();

                for (int a = 0; a < K; a++) {
                    double c = coeffs[a];
                    if (Math.abs(c) < 1e-14) continue;

                    // Build new ops list = existing ops + (siteIndex, a) if a > 0
                    List<int[]> newOps = new ArrayList<>();
                    for (int i = 0; i < existingKey.siteIndices.length; i++) {
                        newOps.add(new int[]{existingKey.siteIndices[i],
                                             existingKey.basisIndices[i]});
                    }
                    if (a > 0) {
                        newOps.add(new int[]{globalSite, a});
                    }

                    SiteOpKey newKey = new SiteOpKey(newOps);
                    next.merge(newKey, existingCoeff * c, Double::sum);
                }
            }

            poly = next;
        }

        // Assemble result row
        double[] row = new double[totalCFs + 1];

        for (Map.Entry<SiteOpKey, Double> entry : poly.entrySet()) {
            SiteOpKey key   = entry.getKey();
            double    coeff = entry.getValue();

            if (Math.abs(coeff) < 1e-14) continue;

            if (key.siteIndices.length == 0) {
                // Empty monomial = constant term
                row[totalCFs] += coeff;
            } else {
                Integer col = cfTable.get(key);
                if (col == null) {
                    throw new IllegalStateException(
                        "No CF column found for site-op product: " + key
                        + "\nCheck that cfTable contains all symmetry-equivalent products.");
                }
                row[col] += coeff;
            }
        }

        return row;
    }

    // =========================================================================
    // STEP 4 — Generate all K^N configurations
    // =========================================================================
    //
    // For a cluster of N sites and K elements, enumerate all K^N assignments.
    // config[i] = element index at site i, values 0..K-1.
    // This is just counting in base K.

    /**
     * Returns all K^N element assignments for an N-site cluster.
     * Each config is an int[] of length N with values in 0..K-1.
     */
    public static List<int[]> allConfigurations(int N, int K) {
        int total = 1;
        for (int i = 0; i < N; i++) total *= K;   // K^N

        List<int[]> configs = new ArrayList<>(total);
        for (int i = 0; i < total; i++) {
            int[] cfg = new int[N];
            int   x   = i;
            for (int s = N - 1; s >= 0; s--) {
                cfg[s] = x % K;
                x     /= K;
            }
            configs.add(cfg);
        }
        return configs;
    }

    // =========================================================================
    // STEP 5 — Group identical rows (Tally), build cmat, lcv, wcv
    // =========================================================================
    //
    // After expanding all K^N configurations, many produce identical coefficient
    // rows. We group them:
    //   - ncv (= lcv) = number of distinct rows
    //   - wcv[v]      = count of configurations that produced row v
    //   - cmat[v]     = the coefficient array for row v
    //
    // Two rows are "identical" if they agree within numerical tolerance.
    // We use rounded long values (multiply by 1e10, round to long) as the hash key.

    /**
     * Canonical key for a coefficient row (tolerance-aware).
     */
    private static final class RowKey {
        private final long[] rounded;
        private static final double SCALE = 1e10;

        RowKey(double[] row) {
            rounded = new long[row.length];
            for (int i = 0; i < row.length; i++) {
                rounded[i] = Math.round(row[i] * SCALE);
            }
        }

        @Override public boolean equals(Object o) {
            return o instanceof RowKey && Arrays.equals(rounded, ((RowKey) o).rounded);
        }

        @Override public int hashCode() { return Arrays.hashCode(rounded); }
    }

    /**
     * Result for one cluster: the C-matrix block, CV count, and CV weights.
     */
    public static final class ClusterCMatrix {
        public final double[][] cmat;  // [ncv][totalCFs+1]
        public final int        lcv;   // number of distinct CVs
        public final int[]      wcv;   // weight of each CV row

        ClusterCMatrix(double[][] cmat, int lcv, int[] wcv) {
            this.cmat = cmat;
            this.lcv  = lcv;
            this.wcv  = wcv;
        }
    }

    /**
     * Builds the C-matrix block for one cluster.
     * Public entry point that takes R-matrix explicitly (for backward compatibility).
     *
     * @param siteIndices  global site indices for the N cluster sites
     * @param K            number of chemical elements
     * @param R            K×K R-matrix
     * @param cfTable      SiteOpKey → CF column index
     * @param totalCFs     total number of CF columns
     * @return             ClusterCMatrix containing cmat, lcv, wcv
     */
    public static ClusterCMatrix buildForCluster(
            int[]                   siteIndices,
            int                     K,
            double[][]              R,
            Map<SiteOpKey, Integer> cfTable,
            int                     totalCFs) {
        return buildForClusterInternal(siteIndices, K, R, cfTable, totalCFs, null);
    }

    /**
     * Enhanced version of buildForCluster with diagnostic output.
     */
    public static ClusterCMatrix buildForCluster(
            int[]                   siteIndices,
            int                     K,
            double[][]              R,
            Map<SiteOpKey, Integer> cfTable,
            int                     totalCFs,
            Consumer<String>        sink) {
        return buildForClusterInternal(siteIndices, K, R, cfTable, totalCFs, sink);
    }

    /**
     * Builds the C-matrix block for one cluster.
     * Overload that derives R-matrix from K automatically.
     *
     * @param siteIndices  global site indices for the N cluster sites
     * @param K            number of chemical elements
     * @param cfTable      SiteOpKey → CF column index
     * @param totalCFs     total number of CF columns
     * @return             ClusterCMatrix containing cmat, lcv, wcv
     */
    public static ClusterCMatrix buildForClusterWithK(
            int[]                   siteIndices,
            int                     K,
            Map<SiteOpKey, Integer> cfTable,
            int                     totalCFs) {
        return buildForClusterWithK(siteIndices, K, cfTable, totalCFs, null);
    }

    /**
     * Enhanced version of buildForClusterWithK with diagnostic output.
     */
    public static ClusterCMatrix buildForClusterWithK(
            int[]                   siteIndices,
            int                     K,
            Map<SiteOpKey, Integer> cfTable,
            int                     totalCFs,
            Consumer<String>        sink) {
        double[][] R = buildRMatrix(K, sink);
        return buildForClusterInternal(siteIndices, K, R, cfTable, totalCFs, sink);
    }

    /**
     * Builds the C-matrix block for one cluster.
     * Overload that extracts siteIndices from a Cluster and siteList.
     * Uses inline position-lookup logic (does not call CMatrix.flattenSiteIndices).
     *
     * @param cluster      the cluster to process
     * @param K            number of chemical elements
     * @param siteList     list of site positions for index lookup
     * @param cfTable      SiteOpKey → CF column index
     * @param totalCFs     total number of CF columns
     * @return             ClusterCMatrix containing cmat, lcv, wcv
     */
    public static ClusterCMatrix buildForCluster(
            Cluster                       cluster,
            int                           K,
            List<Position>                siteList,
            Map<SiteOpKey, Integer>       cfTable,
            int                           totalCFs) {

        // Inline position-lookup logic: extract siteIndices from cluster
        List<Integer> siteIndicesList = new ArrayList<>();
        for (Sublattice sub : cluster.getSublattices()) {
            for (Site site : sub.getSites()) {
                Position pos = site.getPosition();
                int idx = -1;
                for (int i = 0; i < siteList.size(); i++) {
                    if (siteList.get(i).equals(pos)) {
                        idx = i;
                        break;
                    }
                }
                if (idx < 0) {
                    throw new IllegalStateException("Site position not found in siteList: " + pos);
                }
                siteIndicesList.add(idx);
            }
        }

        // Convert List<Integer> to int[]
        int[] siteIndices = new int[siteIndicesList.size()];
        for (int i = 0; i < siteIndicesList.size(); i++) {
            siteIndices[i] = siteIndicesList.get(i);
        }

        return buildForClusterWithK(siteIndices, K, cfTable, totalCFs, null);
    }

    /**
     * Enhanced version of buildForCluster with diagnostic output.
     */
    public static ClusterCMatrix buildForCluster(
            Cluster                       cluster,
            int                           K,
            List<Position>                siteList,
            Map<SiteOpKey, Integer>       cfTable,
            int                           totalCFs,
            Consumer<String>              sink) {

        // Inline position-lookup logic: extract siteIndices from cluster
        List<Integer> siteIndicesList = new ArrayList<>();
        for (Sublattice sub : cluster.getSublattices()) {
            for (Site site : sub.getSites()) {
                Position pos = site.getPosition();
                int idx = -1;
                for (int i = 0; i < siteList.size(); i++) {
                    if (siteList.get(i).equals(pos)) {
                        idx = i;
                        break;
                    }
                }
                if (idx < 0) {
                    throw new IllegalStateException("Site position not found in siteList: " + pos);
                }
                siteIndicesList.add(idx);
            }
        }

        int[] siteIndices = new int[siteIndicesList.size()];
        for (int i = 0; i < siteIndicesList.size(); i++) {
            siteIndices[i] = siteIndicesList.get(i);
        }

        return buildForClusterWithK(siteIndices, K, cfTable, totalCFs, sink);
    }

    /**
     * Internal implementation of C-matrix generation.
     * Kept private so existing tests still compile.
     *
     * @param siteIndices  global site indices for the N cluster sites
     * @param K            number of chemical elements
     * @param R            K×K R-matrix
     * @param cfTable      SiteOpKey → CF column index
     * @param totalCFs     total number of CF columns
     * @return             ClusterCMatrix containing cmat, lcv, wcv
     */
    private static ClusterCMatrix buildForClusterInternal(
            int[]                   siteIndices,
            int                     K,
            double[][]              R,
            Map<SiteOpKey, Integer> cfTable,
            int                     totalCFs,
            Consumer<String>        sink) {

        int N = siteIndices.length;
        int KN = (int) Math.pow(K, N);
        emit(sink, String.format("  Processing cluster with %d sites (K=%d, configs=%d)...", N, K, KN));

        List<int[]> configs = allConfigurations(N, K);

        // Tally: row key → (row, count)
        Map<RowKey, double[]> seenRows   = new LinkedHashMap<>();
        Map<RowKey, Integer>  seenCounts = new LinkedHashMap<>();

        for (int[] config : configs) {
            double[] row = expandConfiguration(siteIndices, config, R, cfTable, totalCFs);
            RowKey   key = new RowKey(row);
            seenRows.putIfAbsent(key, row);
            seenCounts.merge(key, 1, Integer::sum);
        }

        int ncv  = seenRows.size();
        emit(sink, String.format("    -> Tally complete. Identified %d distinct CV rows.", ncv));
        double[][] cmat = new double[ncv][];
        int[]      wcv  = new int[ncv];

        int idx = 0;
        for (Map.Entry<RowKey, double[]> e : seenRows.entrySet()) {
            cmat[idx] = e.getValue();
            wcv[idx]  = seenCounts.get(e.getKey());
            idx++;
        }

        return new ClusterCMatrix(cmat, ncv, wcv);
    }

    // =========================================================================
    // STEP 6 — Top-level: process all cluster types and groups
    // =========================================================================
    //
    // clusterGroups[t][j] = int[] of site indices for cluster (type t, group j)
    // Returns parallel arrays of ClusterCMatrix blocks.

    /**
     * Builds C-matrix blocks for every cluster type and group.
     *
     * @param clusterGroups  [t][j] → int[] of global site indices
     * @param K              number of chemical elements
     * @param cfTable        SiteOpKey → CF column index (from CF identification)
     * @param totalCFs       total number of CF columns
     * @return               [t][j] → ClusterCMatrix
     */
    public static ClusterCMatrix[][] buildAll(
            int[][][]               clusterGroups,
            int                     K,
            Map<SiteOpKey, Integer> cfTable,
            int                     totalCFs) {
        return buildAll(clusterGroups, K, cfTable, totalCFs, null);
    }

    /**
     * Enhanced version of buildAll with diagnostic output.
     */
    public static ClusterCMatrix[][] buildAll(
            int[][][]               clusterGroups,
            int                     K,
            Map<SiteOpKey, Integer> cfTable,
            int                     totalCFs,
            Consumer<String>        sink) {

        double[][] R = buildRMatrix(K, sink);

        int T = clusterGroups.length;
        ClusterCMatrix[][] result = new ClusterCMatrix[T][];

        emit(sink, String.format("\n[STAGE 3] Processing %d Cluster Types...", T));
        for (int t = 0; t < T; t++) {
            int J = clusterGroups[t].length;
            result[t] = new ClusterCMatrix[J];
            emit(sink, String.format("  Type %d has %d groups...", t, J));
            for (int j = 0; j < J; j++) {
                emit(sink, String.format("    Group %d:", j));
                result[t][j] = buildForClusterInternal(clusterGroups[t][j], K, R, cfTable, totalCFs, sink);
            }
        }

        return result;
    }

    // =========================================================================
    // STEP 7 — Build the CF lookup table from the CF identification result
    // =========================================================================
    //
    // The CF identification stage produces, for each CF (t, j, k), a list of
    // representative site-op pairs: [{siteIndex_1, basisIndex_1}, ...].
    //
    // We also need geometry-equivalent entries: two different site assignments
    // that have the same inter-site distance pattern AND same basis indices map
    // to the same CF column.
    //
    // This method builds the table in two passes:
    //   Pass 1 — register the representative directly from CF identification.
    //   Pass 2 — for every combination of sites and basis assignments up to
    //            maxBody sites, compute the geometry signature and add any
    //            match found in Pass 1.
    //
    // Geometry signature of a product {(site_i, b_i)}:
    //   For each site i:  pair = (b_i, sorted list of rounded distances to all other sites)
    //   Sort all pairs canonically.
    //   The signature is the string representation of the sorted list.
    //
    // IMPORTANT: distances must be Cartesian (apply the lattice matrix to
    // fractional coordinates first), otherwise non-cubic structures fail.

    /**
     * Builds the CF lookup table.
     *
     * @param cfRepresentatives  list of CF definitions; each entry is a
     *                           List<int[]> where int[] = {siteIndex, basisIndex},
     *                           and the list position is the CF column number
     * @param sitePositions      [numSites][3] Cartesian coordinates of each site
     * @return                   SiteOpKey → CF column index
     */
    public static Map<SiteOpKey, Integer> buildCFTable(
            List<List<int[]>> cfRepresentatives,
            double[][]        sitePositions) {
        return buildCFTable(cfRepresentatives, sitePositions, null);
    }

    /**
     * Enhanced version of buildCFTable with diagnostic output.
     */
    public static Map<SiteOpKey, Integer> buildCFTable(
            List<List<int[]>> cfRepresentatives,
            double[][]        sitePositions,
            Consumer<String>  sink) {

        emit(sink, "\n[STAGE 2] Building CF Lookup Table (Geometry Augmentation)");
        Map<SiteOpKey, Integer> table   = new LinkedHashMap<>();
        Map<String, Integer>    geoMap  = new LinkedHashMap<>();   // signature → column

        int maxBody  = 0;
        int maxBasis = 0;

        // Pass 1 — register representatives
        emit(sink, "  Pass 1: Registering " + cfRepresentatives.size() + " representatives...");
        for (int col = 0; col < cfRepresentatives.size(); col++) {
            List<int[]> ops = cfRepresentatives.get(col);
            SiteOpKey   key = new SiteOpKey(ops);
            table.put(key, col);

            String sig = geometrySignature(ops, sitePositions);
            geoMap.putIfAbsent(sig, col);

            maxBody  = Math.max(maxBody,  ops.size());
            for (int[] op : ops) maxBasis = Math.max(maxBasis, op[1]);
        }
        emit(sink, "    -> Registered " + table.size() + " unique SiteOpKeys.");

        // Pass 2 — geometry augmentation
        emit(sink, "  Pass 2: Geometry augmentation (maxBody=" + maxBody + ", maxBasis=" + maxBasis + ")...");
        int countBefore = table.size();
        int nSites = sitePositions.length;
        for (int n = 1; n <= Math.min(maxBody, nSites); n++) {
            for (int[] combo : combinations(nSites, n)) {
                for (int[] bases : basisTuples(n, maxBasis)) {
                    List<int[]> ops = new ArrayList<>();
                    for (int i = 0; i < n; i++) ops.add(new int[]{combo[i], bases[i]});

                    SiteOpKey key = new SiteOpKey(ops);
                    if (table.containsKey(key)) continue;

                    String sig = geometrySignature(ops, sitePositions);
                    Integer col = geoMap.get(sig);
                    if (col != null) table.put(key, col);
                }
            }
        }
        emit(sink, "    -> Added " + (table.size() - countBefore) + " equivalent keys via geometry.");
        emit(sink, "    -> Final CF table size: " + table.size() + " entries.");

        return table;
    }

    /**
     * Geometry signature of a site-op product.
     * Two products have the same signature iff they belong to the same CF orbit.
     *
     * Concretely: for each site in the product, record its basis index and its
     * sorted list of Cartesian distances to all other sites in the product.
     * Sort these per-site records lexicographically → canonical signature string.
     */
    private static String geometrySignature(List<int[]> ops, double[][] cartPos) {
        int n = ops.size();
        List<String> perSite = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            int basisIdx = ops.get(i)[1];
            int siteI    = ops.get(i)[0];
            List<Long> dists = new ArrayList<>(n - 1);
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                int siteJ = ops.get(j)[0];
                dists.add(Math.round(cartesianDistance(cartPos[siteI], cartPos[siteJ]) * 1_000_000));
            }
            Collections.sort(dists);
            perSite.add(basisIdx + ":" + dists.toString());
        }
        Collections.sort(perSite);
        return perSite.toString();
    }

    private static double cartesianDistance(double[] a, double[] b) {
        double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
        return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }

    // =========================================================================
    // Matrix inversion — Gauss-Jordan with partial pivoting
    // =========================================================================

    private static double[][] invertMatrix(double[][] A, int n) {
        double[][] M = new double[n][2 * n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n + i] = 1.0;   // augment with identity
        }

        for (int col = 0; col < n; col++) {
            // Partial pivot
            int pivot = col;
            for (int row = col + 1; row < n; row++)
                if (Math.abs(M[row][col]) > Math.abs(M[pivot][col])) pivot = row;
            double[] tmp = M[col]; M[col] = M[pivot]; M[pivot] = tmp;

            double diag = M[col][col];
            if (Math.abs(diag) < 1e-14) throw new ArithmeticException("Singular matrix at col " + col);
            for (int j = 0; j < 2 * n; j++) M[col][j] /= diag;

            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double factor = M[row][col];
                for (int j = 0; j < 2 * n; j++) M[row][j] -= factor * M[col][j];
            }
        }

        double[][] inv = new double[n][n];
        for (int i = 0; i < n; i++) System.arraycopy(M[i], n, inv[i], 0, n);
        return inv;
    }

    // =========================================================================
    // Combinatorics helpers
    // =========================================================================

    /** All k-element combinations from {0..n-1}. */
    private static List<int[]> combinations(int n, int k) {
        List<int[]> result = new ArrayList<>();
        generateCombinations(new int[k], 0, 0, n, k, result);
        return result;
    }

    private static void generateCombinations(int[] combo, int depth, int start, int n, int k, List<int[]> out) {
        if (depth == k) { out.add(combo.clone()); return; }
        for (int i = start; i < n; i++) {
            combo[depth] = i;
            generateCombinations(combo, depth + 1, i + 1, n, k, out);
        }
    }

    /** All k-tuples with values in {1..maxBasis}. */
    private static List<int[]> basisTuples(int k, int maxBasis) {
        List<int[]> result = new ArrayList<>();
        generateTuples(new int[k], 0, k, maxBasis, result);
        return result;
    }

    private static void generateTuples(int[] tuple, int depth, int k, int maxBasis, List<int[]> out) {
        if (depth == k) { out.add(tuple.clone()); return; }
        for (int b = 1; b <= maxBasis; b++) {
            tuple[depth] = b;
            generateTuples(tuple, depth + 1, k, maxBasis, out);
        }
    }

    // =========================================================================
    // DIAGNOSTIC — print a ClusterCMatrix block
    // =========================================================================

    public static void printBlock(ClusterCMatrix block, int totalCFs) {
        System.out.printf("lcv = %d%n", block.lcv);
        for (int v = 0; v < block.lcv; v++) {
            System.out.printf("  row %2d (w=%2d): ", v, block.wcv[v]);
            for (int c = 0; c <= totalCFs; c++) {
                if (Math.abs(block.cmat[v][c]) > 1e-12)
                    System.out.printf("col[%d]=% .6f  ", c, block.cmat[v][c]);
            }
            System.out.println();
        }
    }

    // =========================================================================
    // STEP 3 — runAndPrint: Full pipeline entry point
    // =========================================================================

    /**
     * Runs the C-matrix generation pipeline end-to-end and prints results.
     * Uses only existing public methods from cluster builders; reads but does not modify.
     *
     * @param clusterResult    cluster identification result
     * @param cfResult         CF identification result
     * @param maxClusters      maximal clusters (for building siteList)
     * @param K                number of chemical elements
     * @param sink             output sink for all printed lines
     */
    public static void runAndPrint(
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult      cfResult,
            List<Cluster>               maxClusters,
            int                         K,
            java.util.function.Consumer<String> sink) {

        // Step a: Build siteList
        List<Position> siteList = ClusterBuilders.buildSiteList(maxClusters);
        sink.accept("[CMatrixGenerator.runAndPrint]");
        sink.accept("  siteList size: " + siteList.size());

        // Step b: Build CF site-operator list
        List<List<List<List<SiteOp>>>> cfSiteOpList = ClusterBuilders.buildCFSiteOpList(
            cfResult.getGroupedCFData(), siteList);
        sink.accept("  cfSiteOpList size: " + cfSiteOpList.size());

        // Step c: Build cfTable, then apply geometry augmentation
        List<List<int[]>> cfRepresentatives = new ArrayList<>();

        // Convert cfSiteOpList to cfRepresentatives (list of int[][] per CF)
        for (int t = 0; t < cfSiteOpList.size(); t++) {
            List<List<List<int[]>>> typeGroups = new ArrayList<>();
            for (int j = 0; j < cfSiteOpList.get(t).size(); j++) {
                List<List<int[]>> groupList = new ArrayList<>();
                for (int k = 0; k < cfSiteOpList.get(t).get(j).size(); k++) {
                    List<SiteOp> siteOpList = cfSiteOpList.get(t).get(j).get(k);
                    List<int[]> siteOpInts = new ArrayList<>();
                    for (SiteOp op : siteOpList) {
                        siteOpInts.add(new int[]{op.getSiteIndex(), op.getBasisIndex()});
                    }
                    groupList.add(siteOpInts);
                }
                typeGroups.add(groupList);
            }

            // Flatten: t -> j -> k -> siteOp list
            for (List<List<int[]>> typeGroup : typeGroups) {
                for (List<int[]> cfOps : typeGroup) {
                    cfRepresentatives.add(cfOps);
                }
            }
        }

        // Build cfTable: first register representatives
        Map<SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
        double[][] siteCartesian = new double[siteList.size()][3];
        for (int i = 0; i < siteList.size(); i++) {
            Position p = siteList.get(i);
            siteCartesian[i][0] = p.getX();
            siteCartesian[i][1] = p.getY();
            siteCartesian[i][2] = p.getZ();
            // Note: assuming cubic structure; for non-cubic, apply lattice matrix first
        }

        // Use buildCFTable for geometry augmentation
        cfTable = buildCFTable(cfRepresentatives, siteCartesian, sink);

        // Step d: Get totalCFs
        int totalCFs = cfResult.getTcf();
        sink.accept("  totalCFs: " + totalCFs);

        // Step e: Get ordinal clusters
        List<List<Cluster>> ordClusters = clusterResult.getOrdClusterData().getCoordList();

        // Step f: Process each cluster and print
        for (int t = 0; t < ordClusters.size(); t++) {
            List<Cluster> clusterGroup = ordClusters.get(t);
            for (int j = 0; j < clusterGroup.size(); j++) {
                Cluster cluster = clusterGroup.get(j);
                ClusterCMatrix block = buildForCluster(cluster, K, siteList, cfTable, totalCFs, sink);

                int N = cluster.getAllSites().size();
                int KN = (int) Math.pow(K, N);

                sink.accept("");
                sink.accept(String.format("=== Cluster t=%d j=%d (N=%d sites, K=%d) ===", t, j, N, K));
                sink.accept("lcv = " + block.lcv);
                StringBuilder wcvLine = new StringBuilder("wcv = [");
                for (int v = 0; v < block.lcv; v++) {
                    if (v > 0) wcvLine.append(", ");
                    wcvLine.append(block.wcv[v]);
                }
                wcvLine.append("]");
                sink.accept(wcvLine.toString());

                sink.accept("cmat:");
                for (int v = 0; v < block.lcv; v++) {
                    StringBuilder rowLine = new StringBuilder(String.format("  row %d (w=%d): ", v, block.wcv[v]));
                    for (int c = 0; c <= totalCFs; c++) {
                        double val = block.cmat[v][c];
                        if (Math.abs(val) > 1e-12) {
                            String label = (c < totalCFs) ? "CF" + c : "CONST";
                            rowLine.append(String.format("%s=% .6f  ", label, val));
                        }
                    }
                    sink.accept(rowLine.toString());
                }

                int wSum = 0;
                for (int w : block.wcv) wSum += w;
                sink.accept(String.format("weight_sum = %d  (expected K^N = %d)", wSum, KN));

                // Check partition of unity
                boolean partitionOK = true;
                for (int c = 0; c < totalCFs; c++) {
                    double sum = 0.0;
                    for (int v = 0; v < block.lcv; v++) {
                        sum += block.wcv[v] * block.cmat[v][c];
                    }
                    if (Math.abs(sum) > 1e-9) {
                        partitionOK = false;
                        break;
                    }
                }
                double constSum = 0.0;
                for (int v = 0; v < block.lcv; v++) {
                    constSum += block.wcv[v] * block.cmat[v][totalCFs];
                }
                if (Math.abs(constSum - 1.0) > 1e-9) {
                    partitionOK = false;
                }

                sink.accept("partition_of_unity: " + (partitionOK ? "OK" : "FAIL"));
            }
        }
    }

    // =========================================================================
    // MINIMAL SELF-TEST — binary BCC pair cluster, K=2
    // =========================================================================
    //
    // For K=2, basis={-1,1}, R = 0.5*[[-1,1],[1,1]] (transposed convention).
    //
    // A pair cluster has 2 sites. There are 4 configurations:
    //   (A,A), (A,B), (B,A), (B,B)
    //
    // Expected CVs (for binary): 3 distinct (AA+BB group, AB group, BA group
    // — but AB=BA by symmetry → 2 after geometry, or 3 before).
    // With just numerical ops, you get 4 rows before tally, then grouping reduces them.

    public static void main(String[] args) {
        int K = 2;

        // R-matrix
        double[][] R = buildRMatrix(K);
        System.out.println("R-matrix for K=" + K + ":");
        for (double[] row : R) System.out.println("  " + Arrays.toString(row));

        // Basis
        System.out.println("Basis: " + Arrays.toString(buildBasis(K)));

        // Minimal CF table for a 2-site cluster:
        // CF 0 = point CF at site 0 with basis function 1  → {(0,1)}
        // CF 1 = point CF at site 1 with basis function 1  → {(1,1)}
        // CF 2 = pair CF s1(0)*s1(1)                       → {(0,1),(1,1)}
        Map<SiteOpKey, Integer> cfTable = new LinkedHashMap<>();
        cfTable.put(new SiteOpKey(List.of(new int[]{0, 1})),                    0);  // point site 0
        cfTable.put(new SiteOpKey(List.of(new int[]{1, 1})),                    1);  // point site 1
        cfTable.put(new SiteOpKey(List.of(new int[]{0, 1}, new int[]{1, 1})),   2);  // pair

        int totalCFs = 3;  // columns 0,1,2 + constant at 3

        // Pair cluster: sites 0 and 1
        int[] siteIndices = {0, 1};

        ClusterCMatrix block = buildForClusterWithK(siteIndices, K, cfTable, totalCFs, System.out::println);

        System.out.println("\nC-matrix block for binary pair cluster:");
        System.out.printf("  Columns: CF0(pt0)  CF1(pt1)  CF2(pair)  CONST%n");
        printBlock(block, totalCFs);

        // Verify: sum of wcv must equal K^N = 4
        int wSum = 0;
        for (int w : block.wcv) wSum += w;
        System.out.printf("%nWeight sum = %d (expected %d)%n", wSum, (int)Math.pow(K, siteIndices.length));
    }
}
