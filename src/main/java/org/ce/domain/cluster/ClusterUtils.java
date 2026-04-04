package org.ce.domain.cluster;

import java.util.List;

/**
 * Shared utilities for handling CVM cluster configurations and summaries.
 */
public final class ClusterUtils {

    private ClusterUtils() {}

    /**
     * Counts the number of non-empty clusters in a result set.
     * CVM results in CEWorkbench are typically ordered from largest cluster 
     * down to point clusters, with the "empty" cluster at the very end.
     *
     * @param data the cluster configuration data
     * @return the number of non-point, non-empty clusters
     */
    public static int countNonEmpty(ClusCoordListResult data) {
        if (data == null || data.getClusCoordList().isEmpty()) {
            return 0;
        }
        List<Cluster> list = data.getClusCoordList();
        int count = 0;
        for (Cluster c : list) {
            if (c.getAllSites().size() > 0) {
                count++;
            }
        }
        return count;
    }
    /**
     * Returns a new result set containing only the first {@code n} entries.
     * This is typically used to strip the "empty" cluster from the end of a list.
     *
     * @param data the full result set
     * @param n the number of entries to keep
     * @return a trimmed copy of the result set
     */
    public static ClusCoordListResult trimToNonEmpty(ClusCoordListResult data, int n) {
        if (data == null) return null;
        return new ClusCoordListResult(
                data.getClusCoordList().subList(0, n),
                data.getMultiplicities().subList(0, n),
                data.getOrbitList().subList(0, n),
                data.getRcList().subList(0, n),
                n,
                data.getNumPointSubClusFound());
    }
}
