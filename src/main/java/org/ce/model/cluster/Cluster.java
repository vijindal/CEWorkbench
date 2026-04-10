package org.ce.model.cluster;

import static org.ce.model.cluster.ClusterPrimitives.*;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents a cluster of lattice sites, grouped by sublattice.
 *
 * <p>A {@code Cluster} is the central data structure of the CVM framework.
 * It contains one or more {@link Sublattice} objects, each holding a list
 * of {@link Site} objects.  This two-level hierarchy is required to correctly
 * handle ordered structures (e.g. B2) where different Wyckoff positions carry
 * independent occupation variables.</p>
 *
 * <h2>Sublattice convention</h2>
 * <ul>
 *   <li>A <em>disordered</em> structure (A2, A1) has exactly one sublattice.</li>
 *   <li>An <em>ordered</em> structure (B2, L1â‚‚) has two or more sublattices.</li>
 * </ul>
 *
 * <h2>Sorting and canonicalization</h2>
 * <p>Methods such as {@link #sorted()} and {@link #compareSites(Site, Site)}
 * produce a canonical ordering that is used when comparing clusters for
 * orbit equivalence.  The ordering is: ascending x, then y, then z.</p>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     Sublattice
 * @see     Site
 */
public class Cluster {

    /** Tolerance for floating-point coordinate comparison in site ordering. */
    private static final double TOL = 1e-10;

    /** Ordered list of sublattices that make up this cluster. */
    private final List<Sublattice> sublattices;

    /**
     * Constructs a cluster from an ordered list of sublattices.
     *
     * @param sublattices list of sublattices; must not be {@code null}
     */
    @JsonCreator
    public Cluster(@JsonProperty("sublattices") List<Sublattice> sublattices) {
        this.sublattices = sublattices;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the ordered list of sublattices in this cluster.
     *
     * @return list of {@link Sublattice} objects; never {@code null}
     */
    public List<Sublattice> getSublattices() { return sublattices; }

    /**
     * Flattens all sublattices into a single ordered list of sites.
     *
     * <p>This is used by geometry-only algorithms that do not distinguish
     * between sublattices (e.g. orbit generation for disordered structures).</p>
     *
     * @return a new mutable list containing every {@link Site} in this cluster,
     *         in sublattice order
     */
    public List<Site> getAllSites() {
        List<Site> all = new ArrayList<>();
        for (Sublattice sub : sublattices) {
            all.addAll(sub.getSites());
        }
        return all;
    }

    // -------------------------------------------------------------------------
    // Sorting / canonicalization
    // -------------------------------------------------------------------------

    /**
     * Returns a new {@code Cluster} in which every sublattice is individually
     * sorted in ascending coordinate order.
     *
     * <p>Equivalent to the Mathematica {@code sortClusCoord} function applied
     * at cluster level.  Used to produce canonical representatives for orbit
     * containment checks.</p>
     *
     * @return a new, sorted {@code Cluster}; the original is unchanged
     */
    public Cluster sorted() {
        List<Sublattice> sortedSubs = new ArrayList<>();
        for (Sublattice sub : sublattices) {
            sortedSubs.add(sub.sorted());
        }
        return new Cluster(sortedSubs);
    }

    /**
     * Comparator for {@link Site} objects by position (x â†’ y â†’ z, ascending).
     *
     * <p>Two sites whose coordinates differ by less than {@value #TOL} in a
     * given dimension are considered equal in that dimension and the next
     * dimension is checked.</p>
     *
     * @param a first site; must not be {@code null}
     * @param b second site; must not be {@code null}
     * @return negative, zero, or positive integer consistent with
     *         {@link java.util.Comparator} contract
     */
    public static int compareSites(Site a, Site b) {
        double dx = a.getPosition().getX() - b.getPosition().getX();
        if (Math.abs(dx) > TOL) return dx < 0 ? -1 : 1;

        double dy = a.getPosition().getY() - b.getPosition().getY();
        if (Math.abs(dy) > TOL) return dy < 0 ? -1 : 1;

        double dz = a.getPosition().getZ() - b.getPosition().getZ();
        if (Math.abs(dz) > TOL) return dz < 0 ? -1 : 1;

        return 0;
    }

    /**
     * Returns a compact string listing all sublattices.
     *
     * @return e.g. {@code "[[(0.0,0.0,0.0),s1], [(0.5,0.5,0.5),s1]]"}
     */
    @Override
    public String toString() { return sublattices.toString(); }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this cluster to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [Cluster]
     *   sublattice count : 1
     *   total sites      : 4
     *   Sublattice[0] (3 sites)
     *     [0] (0.000000, 0.000000, 0.000000), symbol=s1
     *     [1] (0.500000, -0.500000, 0.500000), symbol=s1
     *     ...
     * </pre>
     */
    public void printDebug() {
        List<Site> allSites = getAllSites();
        System.out.println("[Cluster]");
        System.out.println("  sublattice count : " + sublattices.size());
        System.out.println("  total sites      : " + allSites.size());
        for (int s = 0; s < sublattices.size(); s++) {
            List<Site> sites = sublattices.get(s).getSites();
            System.out.printf("  Sublattice[%d] (%d sites)%n", s, sites.size());
            for (int i = 0; i < sites.size(); i++) {
                Site site = sites.get(i);
                System.out.printf("    [%d] %s, symbol=%s%n",
                        i, site.getPosition(), site.getSymbol());
            }
        }
    }
}


