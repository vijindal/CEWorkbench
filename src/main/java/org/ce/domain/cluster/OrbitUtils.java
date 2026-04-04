package org.ce.domain.cluster;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Utility class providing core orbit-related operations on {@link Cluster}
 * objects under crystallographic symmetry.
 *
 * <p>The three operations provided here â€” translation-equivalence check,
 * orbit-containment test, and orbit generation â€” are the building blocks
 * of the CVM cluster enumeration algorithm.</p>
 *
 * <h2>Translation equivalence</h2>
 * <p>Two clusters are <em>translation-equivalent</em> if one can be obtained
 * from the other by adding the same integer lattice vector to every site
 * position.  This is the key equivalence relation used in the cluster
 * coordinate list algorithm ({@link org.ce.domain.identification.engine.ClusCoordListGenerator}).</p>
 *
 * <h2>Orbit containment</h2>
 * <p>A cluster {@code c} is <em>contained</em> in an orbit if any member of
 * the orbit is translation-equivalent to {@code c}.</p>
 *
 * <h2>Orbit generation</h2>
 * <p>The orbit of a cluster is the set of all distinct images produced by
 * applying each space-group operation and retaining those not already in
 * the orbit (translation-equivalence is used for deduplication).</p>
 *
 * @author  CVM Project
 * @version 1.0
 */
public class OrbitUtils {

    /** Tolerance for checking whether a coordinate difference is an integer shift. */
    private static final double DELTA = 1e-6;

    /** Private constructor â€” all methods are static utilities. */
    private OrbitUtils() {}

    // -------------------------------------------------------------------------
    // Core operations
    // -------------------------------------------------------------------------

    /**
     * Tests whether two clusters are related by a pure lattice translation.
     *
     * <p>The check is <em>decoration-aware</em>: sites with different species
     * symbols are never considered translation-equivalent, even if their
     * positions differ by an integer vector.</p>
     *
     * <p>Algorithm:</p>
     * <ol>
     *   <li>Assert equal sublattice counts and equal per-sublattice site counts.</li>
     *   <li>Assert that all corresponding site symbols match.</li>
     *   <li>Collect all pairwise position differences into a set; they must all
     *       be the same vector (i.e. the set has size â‰¤ 1).</li>
     *   <li>Assert that the unique difference vector has integer components.</li>
     * </ol>
     *
     * @param c1 first cluster; must not be {@code null}
     * @param c2 second cluster; must not be {@code null}
     * @return {@code true} iff {@code c1} and {@code c2} are related by a
     *         lattice translation and have identical species decoration
     */
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

    /**
     * Tests whether a given cluster is contained (translation-equivalent)
     * within a pre-computed orbit.
     *
     * @param orbit   the set of clusters forming the orbit; must not be {@code null}
     * @param cluster the cluster to test; must not be {@code null}
     * @return {@code true} iff any member of {@code orbit} is translation-equivalent
     *         to {@code cluster}
     */
    public static boolean isContained(List<Cluster> orbit, Cluster cluster) {
        for (Cluster existing : orbit) {
            if (isTranslated(existing, cluster)) return true;
        }
        return false;
    }

    /**
     * Generates the complete symmetry orbit of a cluster under the given
     * space-group operations.
     *
     * <p>Each operation is applied to {@code cluster}; the result is added to
     * the orbit only if it is not translation-equivalent to any cluster already
     * in the orbit.</p>
     *
     * @param cluster    the seed cluster; must not be {@code null}
     * @param spaceGroup the space-group operations; must not be {@code null}
     * @return the full orbit as an ordered list of distinct clusters
     */
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

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Returns {@code true} if {@code value} is within {@value #DELTA} of
     * the nearest integer.
     *
     * @param value the value to test
     * @return {@code true} iff {@code value} â‰ˆ integer
     */
    private static boolean isIntegerShift(double value) {
        return Math.abs(value - Math.round(value)) < DELTA;
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug report for a containment check to standard output.
     *
     * <p>Useful for tracing why a cluster is or is not added during orbit
     * generation.</p>
     *
     * <p>Output format:</p>
     * <pre>
     * [OrbitUtils.isContained]
     *   orbit size : 6
     *   candidate  : [(0.5,0.0,0.5),s1 ...]
     *   result     : true
     * </pre>
     *
     * @param orbit    the orbit to search in
     * @param cluster  the candidate cluster
     */
    public static void printContainmentDebug(List<Cluster> orbit, Cluster cluster) {
        System.out.println("[OrbitUtils.isContained]");
        System.out.println("  orbit size : " + orbit.size());
        System.out.println("  candidate  : " + cluster);
        System.out.println("  result     : " + isContained(orbit, cluster));
    }

    /**
     * Prints a structured debug report for an orbit generation call.
     *
     * <p>Output format:</p>
     * <pre>
     * [OrbitUtils.generateOrbit]
     *   seed cluster : [(0.0,0.0,0.0),s1]
     *   space group ops : 48
     *   generated orbit size : 8
     *   orbit[0] : ...
     *   orbit[1] : ...
     *   ...
     * </pre>
     *
     * @param cluster    the seed cluster
     * @param spaceGroup the space-group operations used
     */
    public static void printOrbitDebug(
            Cluster                 cluster,
            List<SymmetryOperation> spaceGroup) {

        List<Cluster> orbit = generateOrbit(cluster, spaceGroup);

        System.out.println("[OrbitUtils.generateOrbit]");
        System.out.println("  seed cluster       : " + cluster);
        System.out.println("  space group ops    : " + spaceGroup.size());
        System.out.println("  generated orbit size : " + orbit.size());
        for (int i = 0; i < orbit.size(); i++) {
            System.out.println("  orbit[" + i + "] : " + orbit.get(i));
        }
    }
}



