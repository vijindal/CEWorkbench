package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Builds a unique list of site coordinates from maximal clusters.
 */
public final class SiteListBuilder {

    private SiteListBuilder() {}

    /**
     * Returns the unique site coordinates in the order first encountered.
     *
     * @param maxClusters maximal clusters used in the CVM approximation
     * @return list of unique coordinates
     */
    public static List<Position> buildSiteList(List<Cluster> maxClusters) {
        if (maxClusters == null) {
            throw new IllegalArgumentException("maxClusters must not be null");
        }

        List<Position> siteList = new ArrayList<>();
        for (Cluster cluster : maxClusters) {
            for (Sublattice sub : cluster.getSublattices()) {
                for (Site site : sub.getSites()) {
                    Position pos = site.getPosition();
                    if (!contains(siteList, pos)) {
                        siteList.add(pos);
                    }
                }
            }
        }

        return siteList;
    }

    private static boolean contains(List<Position> list, Position pos) {
        for (Position existing : list) {
            if (existing.equals(pos)) {
                return true;
            }
        }
        return false;
    }
}


