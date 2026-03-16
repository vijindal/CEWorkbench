package org.ce.domain.cluster;

import java.util.Objects;

/**
 * Represents a single lattice site within a cluster.
 *
 * <p>A {@code Site} is the atomic unit of the cluster model. It pairs a
 * fractional-coordinate position with a chemical-species symbol (e.g.
 * {@code "s1"}, {@code "s2"}).  Sites are immutable; all transformation
 * methods return new instances.</p>
 *
 * <p>Equality is position- <em>and</em> symbol-aware: two sites that occupy
 * the same position but carry different species labels are considered
 * distinct.</p>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     Sublattice
 * @see     Cluster
 */
public class Site {

    /** Fractional coordinate of this site in the unit cell. */
    private final Position position;

    /**
     * Chemical-species symbol assigned to this site (e.g. {@code "s1"}).
     * {@code null} represents a geometrically empty (undecorated) site.
     */
    private final String symbol;

    /**
     * Constructs a site at the given position with the given species symbol.
     *
     * @param position fractional-coordinate position; must not be {@code null}
     * @param symbol   species label; may be {@code null} for undecorated sites
     */
    public Site(Position position, String symbol) {
        this.position = position;
        this.symbol   = symbol;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the fractional-coordinate position of this site.
     *
     * @return position vector; never {@code null}
     */
    public Position getPosition() { return position; }

    /**
     * Returns the chemical-species symbol of this site.
     *
     * @return species label, or {@code null} for an undecorated site
     */
    public String getSymbol() { return symbol; }

    // -------------------------------------------------------------------------
    // Transformations
    // -------------------------------------------------------------------------

    /**
     * Returns a new site translated by the given displacement vector.
     * The symbol is preserved.
     *
     * @param t translation vector (fractional coordinates); must not be {@code null}
     * @return a new {@code Site} at {@code position + t}
     */
    public Site translate(Position t) {
        return new Site(position.add(t), symbol);
    }

    // -------------------------------------------------------------------------
    // Equality / hashing
    // -------------------------------------------------------------------------

    /**
     * Two sites are equal if and only if they share the same position
     * (within {@code Position}'s built-in tolerance) <em>and</em> the
     * same species symbol.
     *
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Site)) return false;
        Site other = (Site) obj;
        return position.equals(other.position)
            && Objects.equals(symbol, other.symbol);
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return Objects.hash(position, symbol);
    }

    /**
     * Returns a compact string representation: {@code "position,symbol"}.
     *
     * @return e.g. {@code "(0.500000, 0.000000, 0.500000),s1"}
     */
    @Override
    public String toString() {
        return position + "," + symbol;
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this site to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [Site]
     *   position : (0.500000, 0.000000, 0.500000)
     *   symbol   : s1
     * </pre>
     */
    public void printDebug() {
        System.out.println("[Site]");
        System.out.println("  position : " + position);
        System.out.println("  symbol   : " + symbol);
    }
}


