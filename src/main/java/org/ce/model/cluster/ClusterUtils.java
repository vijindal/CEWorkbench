package org.ce.domain.cluster;

import static org.ce.domain.cluster.SpaceGroup.SymmetryOperation;
import static org.ce.domain.cluster.ClusterPrimitives.*;
import static org.ce.domain.cluster.ClusterResults.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Shared utilities for handling CVM cluster configurations, orbits, and summaries.
 */
public final class ClusterUtils {

    /** Tolerance for checking whether a coordinate difference is an integer shift. */
    private static final double DELTA = 1e-6;

    private ClusterUtils() {}

    // =========================================================================
    // List/Summary operations
    // =========================================================================

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

    // =========================================================================
    // Orbit operations (formerly OrbitUtils)
    // =========================================================================

    public static boolean isTranslated(Cluster c1, Cluster c2) {
        if (c1.getSublattices().size() != c2.getSublattices().size())
            return false;

        Set<Position> diffSet = new HashSet<>();

        for (int i = 0; i < c1.getSublattices().size(); i++) {
            Sublattice sub1 = c1.getSublattices().get(i);
            Sublattice sub2 = c2.getSublattices().get(i);

            List<Site> s1 = sub1.getSites();
            List<Site> s2 = sub2.getSites();

            if (s1.size() != s2.size()) return false;

            for (int j = 0; j < s1.size(); j++) {
                Site site1 = s1.get(j);
                Site site2 = s2.get(j);

                // Species must match
                if (!site1.getSymbol().equals(site2.getSymbol()))
                    return false;

                diffSet.add(site2.getPosition().subtract(site1.getPosition()));
            }
        }

        if (diffSet.size() > 1) return false;
        if (diffSet.isEmpty())  return true;  // both empty

        Position d = diffSet.iterator().next();
        return isIntegerShift(d.getX())
            && isIntegerShift(d.getY())
            && isIntegerShift(d.getZ());
    }

    public static boolean isContained(List<Cluster> orbit, Cluster cluster) {
        for (Cluster existing : orbit) {
            if (isTranslated(existing, cluster)) return true;
        }
        return false;
    }

    public static List<Cluster> generateOrbit(
            Cluster                 cluster,
            List<SymmetryOperation> spaceGroup) {

        List<Cluster> orbit = new ArrayList<>();

        for (SymmetryOperation op : spaceGroup) {
            Cluster transformed = op.applyToCluster(cluster);
            if (!isContained(orbit, transformed))
                orbit.add(transformed);
        }

        return orbit;
    }

    private static boolean isIntegerShift(double value) {
        return Math.abs(value - Math.round(value)) < DELTA;
    }

    // -------------------------------------------------------------------------
    // Debug helpers
    // -------------------------------------------------------------------------

    public static void printContainmentDebug(List<Cluster> orbit, Cluster cluster) {
        System.out.println("[ClusterUtils.isContained]");
        System.out.println("  orbit size : " + orbit.size());
        System.out.println("  candidate  : " + cluster);
        System.out.println("  result     : " + isContained(orbit, cluster));
    }

    public static void printOrbitDebug(
            Cluster                 cluster,
            List<SymmetryOperation> spaceGroup) {

        List<Cluster> orbit = generateOrbit(cluster, spaceGroup);

        System.out.println("[ClusterUtils.generateOrbit]");
        System.out.println("  seed cluster       : " + cluster);
        System.out.println("  space group ops    : " + spaceGroup.size());
        System.out.println("  generated orbit size : " + orbit.size());
        for (int i = 0; i < orbit.size(); i++) {
            System.out.println("  orbit[" + i + "] : " + orbit.get(i));
        }
    }
}
