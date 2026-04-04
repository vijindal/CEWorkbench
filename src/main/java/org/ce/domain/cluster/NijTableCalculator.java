package org.ce.domain.cluster;

import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.SubClusterGenerator;
import org.ce.domain.cluster.OrbitUtils;

import java.util.List;

/**
 * Computes the Nij containment table for a list of cluster types.
 *
 * <h2>Definition</h2>
 * <p>{@code nijTable[i][j]} is the number of times HSP cluster type {@code j}
 * appears as a geometrically distinct sub-cluster inside HSP cluster type
 * {@code i}, where "distinct" is measured by containment in the orbit of
 * cluster type {@code j}.</p>
 *
 * <p>Only entries with {@code j >= i} (i.e. cluster {@code j} is the same
 * size as or smaller than cluster {@code i}) are non-zero; larger clusters
 * cannot be sub-clusters of smaller ones.</p>
 *
 * <h2>Mathematica equivalent</h2>
 * <pre>
 * nijTable = getNijTable[disClusCoordList, mhdis, disClusOrbitList]
 *
 * getNijTable[clusCoordList_, clusMList_, clusOrbitList_] := Module[...
 *   For[i=1, i<=numClus, i++,
 *     subClusCoord = genSubClusCoord[clusCoordList[[i]]];
 *     For[k=numSubClus, k>=1, k--,
 *       For[j=1, j<=numClus, j++,
 *         If[j>=i,
 *           If[Map[Length,...]==Map[Length,...],  (* size check *)
 *             If[isContained[clusOrbitList[[j]], subClusCoord[[k]]],
 *               nijTable[[i]][[j]] += 1;
 *             ]
 *           ]
 *         ]
 *       ]
 *     ]
 *   ]
 * ]
 * </pre>
 *
 * <h2>Role in CVM</h2>
 * <p>The Nij table is needed to compute the Kikuchi-Baker entropy coefficients
 * via the inclusion-exclusion recurrence in {@link KikuchiBakerCalculator}.
 * It captures how many times each smaller cluster is "covered" by each
 * larger cluster in the CVM cluster hierarchy.</p>
 *
 * @see KikuchiBakerCalculator
 * @see ClusterIdentifier
 */
public class NijTableCalculator {

    private NijTableCalculator() {}

    /**
     * Computes the full Nij containment table.
     *
     * <p>The table is of size {@code tcdis Ã— tcdis} where {@code tcdis} is the
     * number of HSP cluster types (excluding the empty cluster).
     * Entry {@code [i][j]} counts how many sub-clusters of cluster {@code i}
     * are orbit-equivalent to cluster type {@code j}.</p>
     *
     * @param disClusCoordList canonical representative for each HSP cluster type;
     *                         list is indexed {@code 0..tcdis-1} in descending
     *                         size order (largest cluster first)
     * @param disClusOrbitList full symmetry orbit for each HSP cluster type;
     *                         must be parallel to {@code disClusCoordList}
     * @return square integer matrix of size {@code tcdis Ã— tcdis};
     *         entry {@code [i][j]} is the Nij containment count
     */
    public static int[][] compute(
            List<Cluster>          disClusCoordList,
            List<List<Cluster>>    disClusOrbitList) {

        int numClus = disClusCoordList.size();
        int[][] nijTable = new int[numClus][numClus];

        for (int i = 0; i < numClus; i++) {

            // Generate all geometric sub-clusters of cluster i
            // Equivalent to: subClusCoord = genSubClusCoord[clusCoordList[[i]]]
            List<Cluster> subClusters =
                    SubClusterGenerator.generateSubClusters(
                            disClusCoordList.get(i));

            int numSubClus = subClusters.size();

            // Mathematica: For[k=numSubClus, k>=1, k--]  (descending, skips empty)
            for (int k = numSubClus - 1; k >= 0; k--) {

                Cluster subClus = subClusters.get(k);
                int subSize = subClus.getAllSites().size();

                // Skip empty cluster (size == 0)
                if (subSize == 0) continue;

                // Mathematica: For[j=1, j<=numClus, j++] If[j>=i, ...]
                // In Java: j >= i means j is the same index or further
                // (clusters are sorted descending by size, so j >= i means
                //  cluster[j] is the same size or smaller than cluster[i])
                for (int j = i; j < numClus; j++) {

                    Cluster clusterJ = disClusCoordList.get(j);
                    int sizeJ = clusterJ.getAllSites().size();

                    // Size check: sub-cluster must match cluster j's site count
                    // Mathematica: Map[Length, subClusCoord[[k]]] == Map[Length, clusCoordList[[j]]]
                    // Here we compare total site count (for single-sublattice disordered clusters)
                    if (subSize != sizeJ) continue;

                    // Orbit containment check
                    // Mathematica: isContained[clusOrbitList[[j]], subClusCoord[[k]]]
                    List<Cluster> orbitJ = disClusOrbitList.get(j);

                    if (OrbitUtils.isContained(orbitJ, subClus)) {
                        nijTable[i][j]++;
                    }
                }
            }
        }

        return nijTable;
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints the Nij table to standard output in a readable grid format.
     *
     * @param nijTable the table to print
     */
    public static void printDebug(int[][] nijTable) {
        int n = nijTable.length;
        System.out.println("[NijTableCalculator] nijTable (" + n + "Ã—" + n + "):");
        System.out.print("       ");
        for (int j = 0; j < n; j++) System.out.printf(" j=%-3d", j);
        System.out.println();
        for (int i = 0; i < n; i++) {
            System.out.printf("  i=%-3d", i);
            for (int j = 0; j < n; j++) System.out.printf(" %-5d", nijTable[i][j]);
            System.out.println();
        }
    }
}



