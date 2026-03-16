package org.ce.domain.cluster;

import java.util.ArrayList;
import java.util.List;

/**
 * An ordered collection of {@link Site} objects that belong to the same
 * crystallographic sublattice within a {@link Cluster}.
 *
 * <p>In multi-sublattice structures such as B2 (CsCl-type), each Wyckoff
 * position defines a separate sublattice.  A {@code Sublattice} groups all
 * cluster sites that share the same Wyckoff position, which is important for
 * applying decorated correlation functions correctly.</p>
 *
 * <p>Instances are mutable in the sense that the underlying {@code List<Site>}
 * can be modified by callers who obtain it via {@link #getSites()}.  Use
 * {@link #sorted()} when a canonical, stable ordering is required.</p>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     Cluster
 * @see     Site
 */
public class Sublattice {

    /** Ordered list of sites belonging to this sublattice. */
    private final List<Site> sites;

    /**
     * Constructs a sublattice backed by the given site list.
     *
     * @param sites mutable list of sites; ownership is transferred to this object
     */
    public Sublattice(List<Site> sites) {
        this.sites = sites;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the live, mutable list of sites in this sublattice.
     *
     * @return list of {@link Site} objects; never {@code null}
     */
    public List<Site> getSites() { return sites; }

    // -------------------------------------------------------------------------
    // Transformations
    // -------------------------------------------------------------------------

    /**
     * Returns a new {@code Sublattice} whose sites are sorted in ascending
     * coordinate order (x first, then y, then z) using
     * {@link Cluster#compareSites(Site, Site)}.
     *
     * <p>This corresponds to the Mathematica {@code sortClusCoord} function
     * applied at the sublattice level, and is used to produce canonical
     * cluster representatives for orbit generation.</p>
     *
     * @return a new, sorted {@code Sublattice}; the original is unchanged
     */
    public Sublattice sorted() {
        List<Site> sortedSites = new ArrayList<>(sites);
        // Insertion sort â€” matches Mathematica ordering behaviour
        for (int i = 1; i < sortedSites.size(); i++) {
            Site x = sortedSites.get(i);
            int j = i - 1;
            while (j >= 0 && Cluster.compareSites(sortedSites.get(j), x) > 0) {
                sortedSites.set(j + 1, sortedSites.get(j));
                j--;
            }
            sortedSites.set(j + 1, x);
        }
        return new Sublattice(sortedSites);
    }

    /**
     * Returns a compact string representation listing all sites.
     *
     * @return e.g. {@code "[(0.0,0.0,0.0),s1, (0.5,0.5,0.5),s1]"}
     */
    @Override
    public String toString() { return sites.toString(); }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this sublattice to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [Sublattice]
     *   site count : 2
     *   [0] (0.000000, 0.000000, 0.000000), symbol=s1
     *   [1] (0.500000, 0.500000, 0.500000), symbol=s1
     * </pre>
     */
    public void printDebug() {
        System.out.println("[Sublattice]");
        System.out.println("  site count : " + sites.size());
        for (int i = 0; i < sites.size(); i++) {
            Site s = sites.get(i);
            System.out.printf("  [%d] %s, symbol=%s%n",
                    i, s.getPosition(), s.getSymbol());
        }
    }
}


