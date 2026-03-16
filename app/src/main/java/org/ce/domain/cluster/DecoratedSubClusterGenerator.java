package org.ce.domain.cluster;

import java.util.ArrayList;
import java.util.List;

public class DecoratedSubClusterGenerator {

    // =========================================================
    // genSubClusCoord[clusCoord, basisSymbolList]
    // Exact Mathematica translation
    // =========================================================
    public static List<Cluster> generate(
            Cluster clusCoord,
            List<String> basisSymbolList) {

        int numSubLattice =
                clusCoord.getSublattices().size();

        // 1ï¸âƒ£ Flatten geometry
        List<Site> flatSites =
                new ArrayList<>(clusCoord.getAllSites());

        flatSites.sort(Cluster::compareSites);

        // 2ï¸âƒ£ Build disClus structure
        // For each site:
        //   option 0 â†’ empty
        //   option j â†’ (coord, basisSymbol[j])
        List<List<Site>> disClus =
                new ArrayList<>();

        for (Site site : flatSites) {

            List<Site> siteOptions =
                    new ArrayList<>();

            // empty option
            siteOptions.add(null);

            // decoration options
            for (String symbol : basisSymbolList) {
                siteOptions.add(
                        new Site(site.getPosition(), symbol));
            }

            disClus.add(siteOptions);
        }

        // 3ï¸âƒ£ Cartesian product (Tuples equivalent)
        List<List<Site>> tuples =
                cartesianProduct(disClus);

        // 4ï¸âƒ£ Rebuild sublattice grouping
        List<Cluster> result =
                new ArrayList<>();

        for (List<Site> tuple : tuples) {

            List<Sublattice> newSubs =
                    new ArrayList<>();

            for (int k = 0; k < numSubLattice; k++) {
                newSubs.add(
                        new Sublattice(new ArrayList<>()));
            }

            for (Site site : tuple) {

                if (site == null)
                    continue;

                for (int k = 0;
                     k < numSubLattice;
                     k++) {

                    List<Site> original =
                            clusCoord.getSublattices()
                                     .get(k)
                                     .getSites();

                    if (containsPosition(original,
                            site.getPosition())) {

                        newSubs.get(k)
                               .getSites()
                               .add(site);
                    }
                }
            }

            result.add(new Cluster(newSubs));
        }

        return result;
    }

    // =========================================================
    // Cartesian product (Tuples equivalent)
    // =========================================================
    private static List<List<Site>> cartesianProduct(
            List<List<Site>> lists) {

        List<List<Site>> result =
                new ArrayList<>();

        cartesianRecursive(lists, 0,
                new ArrayList<>(), result);

        return result;
    }

    private static void cartesianRecursive(
            List<List<Site>> lists,
            int depth,
            List<Site> current,
            List<List<Site>> result) {

        if (depth == lists.size()) {
            result.add(new ArrayList<>(current));
            return;
        }

        for (Site s : lists.get(depth)) {
            current.add(s);
            cartesianRecursive(
                    lists, depth + 1,
                    current, result);
            current.remove(current.size() - 1);
        }
    }

    // =========================================================
    // Compare by position only (ignore symbol)
    // Mathematica MemberQ behavior
    // =========================================================
    private static boolean containsPosition(
            List<Site> list,
            Position pos) {

        for (Site s : list) {
            if (s.getPosition()
                 .equalsWithTolerance(pos, 1e-8))
                return true;
        }

        return false;
    }
}



