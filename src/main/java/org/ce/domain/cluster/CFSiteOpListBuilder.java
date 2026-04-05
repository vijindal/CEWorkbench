package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterResults.*;

import static org.ce.domain.cluster.ClusterKeys.*;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Builds per-CF site-operator lists from grouped CF clusters.
 */
public final class CFSiteOpListBuilder {

    private CFSiteOpListBuilder() {}

    /**
     * Returns a nested list matching grouped CF structure:
     * [t][j][k] -> list of SiteOp for CF k in group j of type t.
     */
    public static List<List<List<List<SiteOp>>>> build(
            GroupedCFResult groupedCFData,
            List<Position> siteList) {

        if (groupedCFData == null) {
            throw new IllegalArgumentException("groupedCFData must not be null");
        }
        if (siteList == null) {
            throw new IllegalArgumentException("siteList must not be null");
        }

        List<List<List<Cluster>>> coordData = groupedCFData.getCoordData();
        List<List<List<List<SiteOp>>>> result = new ArrayList<>();

        for (List<List<Cluster>> typeGroups : coordData) {
            List<List<List<SiteOp>>> typeOut = new ArrayList<>();
            for (List<Cluster> group : typeGroups) {
                List<List<SiteOp>> groupOut = new ArrayList<>();
                for (Cluster cfCluster : group) {
                    groupOut.add(buildForCluster(cfCluster, siteList));
                }
                typeOut.add(groupOut);
            }
            result.add(typeOut);
        }

        return result;
    }

    private static List<SiteOp> buildForCluster(Cluster cluster, List<Position> siteList) {
        List<SiteOp> ops = new ArrayList<>();
        for (Sublattice sub : cluster.getSublattices()) {
            for (Site site : sub.getSites()) {
                int siteIndex = indexOf(siteList, site.getPosition());
                if (siteIndex < 0) {
                    throw new IllegalStateException("Site position not found in site list: "
                            + site.getPosition());
                }
                int basisIndex = parseBasisIndex(site.getSymbol());
                ops.add(new SiteOp(siteIndex, basisIndex));
            }
        }
        return ops;
    }

    private static int indexOf(List<Position> siteList, Position pos) {
        for (int i = 0; i < siteList.size(); i++) {
            if (siteList.get(i).equals(pos)) {
                return i;
            }
        }
        return -1;
    }

    private static int parseBasisIndex(String symbol) {
        if (symbol == null) {
            throw new IllegalArgumentException("Site symbol is null; CF sites must be decorated");
        }
        if (!symbol.startsWith("s")) {
            throw new IllegalArgumentException("Unexpected site symbol: " + symbol);
        }
        try {
            return Integer.parseInt(symbol.substring(1));
        } catch (NumberFormatException ex) {
            throw new IllegalArgumentException("Invalid basis symbol: " + symbol, ex);
        }
    }
}



