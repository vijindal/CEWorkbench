package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Unified utility for generating sub-clusters (both undecorated and decorated).
 */
public class SubClusterGenerators {

    private SubClusterGenerators() {}

    // =========================================
    // 1. Undecorated sub-clusters
    // Exact Mathematica translation
    // Works for both disordered (A2) and ordered (B2)
    // =========================================
    public static List<Cluster> generateSubClusters(Cluster clusCoord) {

        List<Cluster> result = new ArrayList<>();

        // 1. Number of sublattices
        List<Sublattice> originalSubs = clusCoord.getSublattices();
        int numSubLattice = originalSubs.size();

        // 2. Flatten cluster across sublattices
        List<Site> disClusCoord = clusCoord.getAllSites();

        // 3. Sort flattened list (Mathematica sortClusCoord)
        disClusCoord = sortFlatSites(disClusCoord);

        // 4. Generate all subsets (INCLUDING EMPTY)
        List<List<Site>> subsets = generateAllSubsets(disClusCoord);

        // 5. Reconstruct sublattice grouping
        for (List<Site> subset : subsets) {

            // Initialize empty sublattices
            List<Sublattice> subClusterSublattices = new ArrayList<>();

            for (int i = 0; i < numSubLattice; i++) {
                subClusterSublattices.add(new Sublattice(new ArrayList<>()));
            }

            // Distribute subset sites back to original sublattices
            for (Site site : subset) {
                for (int k = 0; k < numSubLattice; k++) {
                    List<Site> originalSubSites = originalSubs.get(k).getSites();
                    if (originalSubSites.contains(site)) {
                        subClusterSublattices.get(k).getSites().add(site);
                    }
                }
            }

            result.add(new Cluster(subClusterSublattices));
        }

        return result;
    }

    private static List<Site> sortFlatSites(List<Site> sites) {
        List<Site> sorted = new ArrayList<>(sites);
        for (int i = 1; i < sorted.size(); i++) {
            Site x = sorted.get(i);
            int j = i - 1;
            while (j >= 0 && Cluster.compareSites(sorted.get(j), x) > 0) {
                sorted.set(j + 1, sorted.get(j));
                j--;
            }
            sorted.set(j + 1, x);
        }
        return sorted;
    }

    private static List<List<Site>> generateAllSubsets(List<Site> sites) {
        List<List<Site>> subsets = new ArrayList<>();
        int n = sites.size();
        int total = 1 << n;

        for (int mask = 0; mask < total; mask++) {
            List<Site> subset = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    subset.add(sites.get(i));
                }
            }
            subsets.add(subset);
        }
        return subsets;
    }

    // =========================================================
    // 2. Decorated sub-clusters
    // genSubClusCoord[clusCoord, basisSymbolList]
    // Exact Mathematica translation
    // =========================================================
    public static List<Cluster> generateDecorated(Cluster clusCoord, List<String> basisSymbolList) {

        int numSubLattice = clusCoord.getSublattices().size();

        // 1. Flatten geometry
        List<Site> flatSites = new ArrayList<>(clusCoord.getAllSites());
        flatSites.sort(Cluster::compareSites);

        // 2. Build disClus structure
        // For each site:
        //   option 0 -> empty
        //   option j -> (coord, basisSymbol[j])
        List<List<Site>> disClus = new ArrayList<>();

        for (Site site : flatSites) {
            List<Site> siteOptions = new ArrayList<>();
            siteOptions.add(null); // empty option
            for (String symbol : basisSymbolList) {
                siteOptions.add(new Site(site.getPosition(), symbol));
            }
            disClus.add(siteOptions);
        }

        // 3. Cartesian product (Tuples equivalent)
        List<List<Site>> tuples = cartesianProduct(disClus);

        // 4. Rebuild sublattice grouping
        List<Cluster> result = new ArrayList<>();

        for (List<Site> tuple : tuples) {
            List<Sublattice> newSubs = new ArrayList<>();
            for (int k = 0; k < numSubLattice; k++) {
                newSubs.add(new Sublattice(new ArrayList<>()));
            }

            for (Site site : tuple) {
                if (site == null) continue;

                for (int k = 0; k < numSubLattice; k++) {
                    List<Site> original = clusCoord.getSublattices().get(k).getSites();
                    if (containsPosition(original, site.getPosition())) {
                        newSubs.get(k).getSites().add(site);
                    }
                }
            }

            result.add(new Cluster(newSubs));
        }

        return result;
    }

    private static List<List<Site>> cartesianProduct(List<List<Site>> lists) {
        List<List<Site>> result = new ArrayList<>();
        cartesianRecursive(lists, 0, new ArrayList<>(), result);
        return result;
    }

    private static void cartesianRecursive(List<List<Site>> lists, int depth, List<Site> current, List<List<Site>> result) {
        if (depth == lists.size()) {
            result.add(new ArrayList<>(current));
            return;
        }

        for (Site s : lists.get(depth)) {
            current.add(s);
            cartesianRecursive(lists, depth + 1, current, result);
            current.remove(current.size() - 1);
        }
    }

    private static boolean containsPosition(List<Site> list, Position pos) {
        for (Site s : list) {
            if (s.getPosition().equalsWithTolerance(pos, 1e-8)) return true;
        }
        return false;
    }
}
