package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;

public class SubClusterGenerator {

    // =========================================
    // genSubClusCoord
    // Exact Mathematica translation
    // Works for both disordered (A2) and ordered (B2)
    // =========================================
    public static List<Cluster> generateSubClusters(Cluster clusCoord) {

        List<Cluster> result = new ArrayList<>();

        // 1ï¸âƒ£ Number of sublattices
        List<Sublattice> originalSubs =
                clusCoord.getSublattices();

        int numSubLattice = originalSubs.size();

        // 2ï¸âƒ£ Flatten cluster across sublattices
        List<Site> disClusCoord =
                clusCoord.getAllSites();

        // 3ï¸âƒ£ Sort flattened list (Mathematica sortClusCoord)
        disClusCoord = sortFlatSites(disClusCoord);

        // 4ï¸âƒ£ Generate all subsets (INCLUDING EMPTY)
        List<List<Site>> subsets =
                generateAllSubsets(disClusCoord);

        // 5ï¸âƒ£ Reconstruct sublattice grouping
        for (List<Site> subset : subsets) {

            // Initialize empty sublattices
            List<Sublattice> subClusterSublattices =
                    new ArrayList<>();

            for (int i = 0; i < numSubLattice; i++) {
                subClusterSublattices.add(
                        new Sublattice(new ArrayList<>()));
            }

            // Distribute subset sites back to original sublattices
            for (Site site : subset) {

                for (int k = 0; k < numSubLattice; k++) {

                    List<Site> originalSubSites =
                            originalSubs.get(k).getSites();

                    if (originalSubSites.contains(site)) {

                        subClusterSublattices.get(k)
                                .getSites()
                                .add(site);
                    }
                }
            }

            result.add(new Cluster(subClusterSublattices));
        }

        return result;
    }

    // =========================================
    // Equivalent to sortClusCoord on flat list
    // Uses insertion sort to match Mathematica
    // =========================================
    private static List<Site> sortFlatSites(List<Site> sites) {

        List<Site> sorted =
                new ArrayList<>(sites);

        for (int i = 1; i < sorted.size(); i++) {

            Site x = sorted.get(i);
            int j = i - 1;

            while (j >= 0 &&
                    Cluster.compareSites(sorted.get(j), x) > 0) {

                sorted.set(j + 1, sorted.get(j));
                j--;
            }

            sorted.set(j + 1, x);
        }

        return sorted;
    }

    // =========================================
    // Equivalent to Mathematica Subsets
    // Includes empty subset
    // =========================================
    private static List<List<Site>> generateAllSubsets(List<Site> sites) {

        List<List<Site>> subsets =
                new ArrayList<>();

        int n = sites.size();
        int total = 1 << n;

        for (int mask = 0; mask < total; mask++) {

            List<Site> subset =
                    new ArrayList<>();

            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    subset.add(sites.get(i));
                }
            }

            subsets.add(subset);
        }

        return subsets;
    }
}



