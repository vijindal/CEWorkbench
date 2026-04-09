package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterResults.*;

import static org.ce.domain.cluster.SpaceGroup.SymmetryOperation;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Generates the complete list of distinct cluster coordinate types for a
 * given set of maximal clusters and a crystallographic space group.
 *
 * <p>This class is the Java equivalent of the Mathematica function
 * {@code genClusCoordList[maxClusCoord, spaceGroup]}.  Starting from a list
 * of <em>maximal</em> (largest) clusters, it:</p>
 * <ol>
 *   <li>Enumerates all geometric sub-clusters of each maximal cluster
 *       (or, in the decorated variant, all decorated sub-clusters).</li>
 *   <li>Filters out clusters that are symmetry-equivalent to one already seen
 *       (using {@link org.ce.domain.identification.symmetry.OrbitUtils#isContained}).</li>
 *   <li>Generates the full symmetry orbit for each new cluster type.</li>
 *   <li>Normalises multiplicities by the total point-cluster orbit count.</li>
 *   <li>Sorts the final list by descending cluster size.</li>
 * </ol>
 *
 * <h2>Usage</h2>
 * <pre>{@code
 * List<Cluster> maxClusters = InputLoader.parseClusterFile("cluster/A2-T.txt");
 * List<SymmetryOperation> symOps = InputLoader.parseSymmetryFile("A2-SG");
 *
 * // Geometric (undecorated) version
 * ClusCoordListResult result = ClusCoordListGenerator.generate(maxClusters, symOps);
 *
 * // Decorated (multicomponent CF) version
 * List<String> symbols = List.of("s1", "s2");
 * ClusCoordListResult cfResult = ClusCoordListGenerator.generate(maxClusters, symOps, symbols);
 * }</pre>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     ClusCoordListResult
 * @see     SubClusterGenerators
 */
public class ClusCoordListGenerator {

    /** Private constructor â€” all methods are static utilities. */
    private ClusCoordListGenerator() {}

    // -------------------------------------------------------------------------
    // Public API
    // -------------------------------------------------------------------------

    /**
     * Generates the cluster-coordinate list for a <em>geometric</em>
     * (undecorated) structure.
     *
     * <p>Equivalent to the Mathematica call:
     * {@code genClusCoordList[maxClusCoord, spaceGroup]}</p>
     *
     * @param maxClusCoord list of maximal cluster representatives;
     *                     must not be {@code null} or empty
     * @param spaceGroup   symmetry operations of the space group;
     *                     must not be {@code null}
     * @return fully populated {@link ClusCoordListResult}
     */
    public static ClusCoordListResult generate(
            List<Cluster> maxClusCoord,
            List<SymmetryOperation> spaceGroup) {
        return generateInternal(maxClusCoord, spaceGroup, null);
    }

    /**
     * Generates the cluster-coordinate list for a <em>decorated</em>
     * (multicomponent correlation function) structure.
     *
     * <p>Equivalent to the Mathematica call:
     * {@code genClusCoordList[maxClusCoord, spaceGroup, basisSymbolList]}</p>
     *
     * @param maxClusCoord    list of maximal cluster representatives;
     *                        must not be {@code null} or empty
     * @param spaceGroup      symmetry operations of the space group;
     *                        must not be {@code null}
     * @param basisSymbolList species symbols used for decoration
     *                        (e.g. {@code ["s1", "s2"]} for binary);
     *                        must not be {@code null}
     * @return fully populated {@link ClusCoordListResult}
     */
    public static ClusCoordListResult generate(
            List<Cluster> maxClusCoord,
            List<SymmetryOperation> spaceGroup,
            List<String> basisSymbolList) {
        return generateInternal(maxClusCoord, spaceGroup, basisSymbolList);
    }

    // -------------------------------------------------------------------------
    // Internal implementation
    // -------------------------------------------------------------------------

    /**
     * Core implementation shared by both public overloads.
     *
     * @param maxClusCoord    maximal clusters
     * @param spaceGroup      space-group operations
     * @param basisSymbolList {@code null} for geometric mode, non-null for decorated mode
     * @return {@link ClusCoordListResult}
     */
    private static ClusCoordListResult generateInternal(
            List<Cluster>           maxClusCoord,
            List<SymmetryOperation> spaceGroup,
            List<String>            basisSymbolList) {

        List<Cluster>          clusCoordList   = new ArrayList<>();
        List<List<Cluster>>    subClusOrbitList = new ArrayList<>();
        List<Integer>          subClusMList    = new ArrayList<>();
        List<List<Integer>>    rc              = new ArrayList<>();

        int numMaxClus = maxClusCoord.size();

        // ------------------------------------------------------------------
        // STEP 1: Iterate maximal clusters and collect unique sub-cluster types
        // ------------------------------------------------------------------
        for (int k = 0; k < numMaxClus; k++) {

            List<Cluster> subClusCoord;

            if (basisSymbolList == null) {
                // Geometric (undecorated) sub-clusters
                subClusCoord = SubClusterGenerators.generateSubClusters(
                        maxClusCoord.get(k));
            } else {
                // Decorated sub-clusters for CFs
                subClusCoord = SubClusterGenerators.generateDecorated(
                        maxClusCoord.get(k), basisSymbolList);
            }

            // Sort descending by size (Mathematica convention)
            subClusCoord.sort(
                Comparator.comparingInt((Cluster c) -> c.getAllSites().size())
                          .reversed());

            // Mathematica loop: For[i=numSubClus, i>=1, i--]
            for (int i = subClusCoord.size() - 1; i >= 0; i--) {

                Cluster candidate = subClusCoord.get(i);
                boolean isNew = true;

                for (List<Cluster> orbit : subClusOrbitList) {
                    if (ClusterUtils.isContained(orbit, candidate)) {
                        isNew = false;
                        break;
                    }
                }

                if (isNew) {
                    clusCoordList.add(candidate);

                    List<Cluster> orbit =
                            ClusterUtils.generateOrbit(candidate, spaceGroup);

                    subClusOrbitList.add(orbit);
                    subClusMList.add(orbit.size());

                    // Per-sublattice site counts
                    List<Integer> rcEntry = new ArrayList<>();
                    for (Sublattice sub : candidate.getSublattices()) {
                        rcEntry.add(sub.getSites().size());
                    }
                    rc.add(rcEntry);
                }
            }
        }

        int tc = clusCoordList.size();

        // ------------------------------------------------------------------
        // STEP 2: Normalise multiplicities by total point-cluster orbit count
        // ------------------------------------------------------------------
        int    numPointSubClusFound  = 0;
        double pointM                = 0;
        List<String> pointClusKeys   = new ArrayList<>();

        for (int i = 0; i < tc; i++) {
            List<Site> flat = clusCoordList.get(i).getAllSites();
            if (flat.size() == 1) {
                String key = flat.get(0).getPosition().toString();
                if (!pointClusKeys.contains(key)) {
                    pointClusKeys.add(key);
                    pointM += subClusOrbitList.get(i).size();
                }
                numPointSubClusFound++;
            }
        }

        List<Double> normalizedM = new ArrayList<>();
        for (Integer m : subClusMList) {
            normalizedM.add(m / pointM);
        }

        // ------------------------------------------------------------------
        // STEP 3: Final sort â€” descending by total site count
        // ------------------------------------------------------------------
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < tc; i++) indices.add(i);

        indices.sort((i1, i2) ->
            compareClusterTypesCanonical(clusCoordList.get(i1), clusCoordList.get(i2)));

        List<Cluster>          finalClusList = new ArrayList<>();
        List<Double>           finalM        = new ArrayList<>();
        List<List<Cluster>>    finalOrbit    = new ArrayList<>();
        List<List<Integer>>    finalRc       = new ArrayList<>();

        for (int idx : indices) {
            finalClusList.add(clusCoordList.get(idx));
            finalM.add(normalizedM.get(idx));
            finalOrbit.add(subClusOrbitList.get(idx));
            finalRc.add(rc.get(idx));
        }

        return new ClusCoordListResult(
                finalClusList, finalM, finalOrbit, finalRc,
                tc, numPointSubClusFound);
    }

    /**
     * Canonical cluster-type comparator used for final ordering.
     *
     * <p>Primary key: descending number of sites (existing behavior).</p>
     * <p>Secondary key: descending pair-distance signature for equal-size clusters
     * to keep shell ordering stable (e.g. BCC 2nd-NN pair before 1st-NN pair).</p>
     * <p>Tertiary key: lexicographic site coordinates for full determinism.</p>
     */
    private static int compareClusterTypesCanonical(Cluster a, Cluster b) {
        int nCmp = Integer.compare(b.getAllSites().size(), a.getAllSites().size());
        if (nCmp != 0) {
            return nCmp;
        }

        long[] sa = distanceSignatureDesc(a);
        long[] sb = distanceSignatureDesc(b);
        int m = Math.min(sa.length, sb.length);
        for (int i = 0; i < m; i++) {
            if (sa[i] != sb[i]) {
                return Long.compare(sb[i], sa[i]); // descending
            }
        }
        if (sa.length != sb.length) {
            return Integer.compare(sb.length, sa.length);
        }

        List<Site> fa = a.sorted().getAllSites();
        List<Site> fb = b.sorted().getAllSites();
        int s = Math.min(fa.size(), fb.size());
        for (int i = 0; i < s; i++) {
            int cmp = Cluster.compareSites(fa.get(i), fb.get(i));
            if (cmp != 0) {
                return cmp;
            }
        }
        return Integer.compare(fa.size(), fb.size());
    }

    /**
     * Returns sorted (descending) squared pair distances as an integer signature.
     */
    private static long[] distanceSignatureDesc(Cluster c) {
        List<Site> sites = c.getAllSites();
        if (sites.size() < 2) {
            return new long[0];
        }
        List<Long> d2 = new ArrayList<>();
        for (int i = 0; i < sites.size(); i++) {
            for (int j = i + 1; j < sites.size(); j++) {
                double dx = sites.get(i).getPosition().getX() - sites.get(j).getPosition().getX();
                double dy = sites.get(i).getPosition().getY() - sites.get(j).getPosition().getY();
                double dz = sites.get(i).getPosition().getZ() - sites.get(j).getPosition().getZ();
                long q = Math.round((dx * dx + dy * dy + dz * dz) * 1_000_000_000L);
                d2.add(q);
            }
        }
        d2.sort((x, y) -> Long.compare(y, x)); // descending
        long[] out = new long[d2.size()];
        for (int i = 0; i < d2.size(); i++) {
            out[i] = d2.get(i);
        }
        return out;
    }

    // -------------------------------------------------------------------------
    // Debug static helper
    // -------------------------------------------------------------------------

    /**
     * Prints a detailed trace of the inputs passed to the generator and the
     * resulting {@link ClusCoordListResult} to standard output.
     *
     * <p>Call this after {@link #generate} to understand the full pipeline:</p>
     * <pre>{@code
     * ClusCoordListResult result = ClusCoordListGenerator.generate(maxClus, symOps);
     * ClusCoordListGenerator.printDebug(maxClus, symOps, null, result);
     * }</pre>
     *
     * @param maxClusCoord    input maximal clusters
     * @param spaceGroup      input symmetry operations
     * @param basisSymbolList input basis symbols ({@code null} for geometric mode)
     * @param result          the computed result
     */
    public static void printDebug(
            List<Cluster>           maxClusCoord,
            List<SymmetryOperation> spaceGroup,
            List<String>            basisSymbolList,
            ClusCoordListResult     result) {

        System.out.println("==============================");
        System.out.println("[ClusCoordListGenerator] DEBUG");
        System.out.println("==============================");
        System.out.println("  INPUT:");
        System.out.println("    maximal cluster count : " + maxClusCoord.size());
        System.out.println("    space group ops       : " + spaceGroup.size());
        System.out.println("    mode                  : " +
                (basisSymbolList == null ? "geometric (undecorated)"
                                        : "decorated " + basisSymbolList));
        System.out.println("  OUTPUT:");
        result.printDebug();
        System.out.println("==============================");
    }
}



