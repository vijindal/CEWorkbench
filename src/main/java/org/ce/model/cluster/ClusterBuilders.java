package org.ce.model.cluster;

import static org.ce.model.cluster.ClusterPrimitives.*;
import static org.ce.model.cluster.ClusterResults.*;
import static org.ce.model.cluster.ClusterKeys.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Unified utility class for building structural lists and grouping data.
 */
public final class ClusterBuilders {

    private ClusterBuilders() {}

    // =========================================================================
    // SiteListBuilder
    // =========================================================================

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

    // =========================================================================
    // CFSiteOpListBuilder
    // =========================================================================

    /**
     * Builds per-CF site-operator lists from grouped CF clusters.
     * Returns a nested list matching grouped CF structure:
     * [t][j][k] -> list of SiteOp for CF k in group j of type t.
     */
    public static List<List<List<List<SiteOp>>>> buildCFSiteOpList(
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

    // =========================================================================
    // CFGroupGenerator
    // =========================================================================

    public static GroupedCFResult groupCFData(
            ClusCoordListResult disClusData,
            ClusCoordListResult disCFData,
            ClassifiedClusterResult ordCFData,
            List<String> basisSymbolList) {

        List<Cluster> disClusCoordData = disClusData.getClusCoordList();
        List<Cluster> disCFCoordData = disCFData.getClusCoordList();
        List<List<Cluster>> disClusOrbitList = disClusData.getOrbitList();

        List<List<Cluster>> ordCFCoordData = ordCFData.getCoordList();
        List<List<Double>> ordCFMData = ordCFData.getMultiplicityList();
        List<List<List<Cluster>>> ordCFOrbitData = ordCFData.getOrbitList();
        List<List<List<Integer>>> ordCFRData = ordCFData.getRcList();

        // ---------------------------------------------------
        // STEP 1: Remove decorations from disordered CF clusters
        // Equivalent to transDisCFCoordData in Mathematica
        // ---------------------------------------------------
        List<Cluster> transDisCFCoordData = new ArrayList<>();
        for (Cluster cf : disCFCoordData) {
            List<Sublattice> newSubs = new ArrayList<>();
            for (Sublattice sub : cf.getSublattices()) {
                List<Site> newSites = new ArrayList<>();
                for (Site s : sub.getSites()) {
                    newSites.add(new Site(s.getPosition(), basisSymbolList.get(0)));
                }
                newSubs.add(new Sublattice(newSites));
            }
            transDisCFCoordData.add(new Cluster(newSubs));
        }

        // ---------------------------------------------------
        // STEP 2: Group ordered CFs into disordered clusters
        // ---------------------------------------------------
        List<List<List<Cluster>>> groupedCFCoordData = new ArrayList<>();
        List<List<List<Double>>> groupedCFMData = new ArrayList<>();
        List<List<List<List<Cluster>>>> groupedCFOrbitData = new ArrayList<>();
        List<List<List<List<Integer>>>> groupedCFRData = new ArrayList<>();

        int tcdis = disClusCoordData.size();
        for (int i = 0; i < tcdis; i++) {
            groupedCFCoordData.add(new ArrayList<>());
            groupedCFMData.add(new ArrayList<>());
            groupedCFOrbitData.add(new ArrayList<>());
            groupedCFRData.add(new ArrayList<>());

            List<Cluster> disOrbit = disClusOrbitList.get(i);
            for (int j = 0; j < transDisCFCoordData.size(); j++) {
                Cluster undecoratedCF = transDisCFCoordData.get(j);
                if (ClusterUtils.isContained(disOrbit, undecoratedCF)) {
                    groupedCFCoordData.get(i).add(ordCFCoordData.get(j));
                    groupedCFMData.get(i).add(ordCFMData.get(j));
                    groupedCFOrbitData.get(i).add(ordCFOrbitData.get(j));
                    groupedCFRData.get(i).add(ordCFRData.get(j));
                }
            }
        }

        return new GroupedCFResult(
                groupedCFCoordData,
                groupedCFMData,
                groupedCFOrbitData,
                groupedCFRData
        );
    }
}
