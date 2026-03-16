package org.ce.domain.cluster;

import java.util.List;

/**
 * Represents a symmetry-equivalence class (orbit) of clusters.
 *
 * <p>Under the space-group symmetry of the crystal, many distinct clusters
 * in the supercell are related by symmetry and therefore share the same
 * CVM interaction coefficient.  A {@code ClusterType} bundles:</p>
 * <ul>
 *   <li>A canonical <em>representative</em> cluster (the first member of
 *       the orbit chosen by the generator).</li>
 *   <li>The full <em>orbit</em> â€” all symmetry-equivalent clusters.</li>
 *   <li>The <em>multiplicity</em> â€” the number of orbit members (equal to
 *       {@code orbit.size()}).</li>
 * </ul>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     Cluster
 */
public class ClusterType {

    /** Canonical representative of this symmetry class. */
    private final Cluster representative;

    /** All clusters related to the representative by space-group symmetry. */
    private final List<Cluster> orbit;

    /**
     * Number of symmetry-equivalent clusters in this class.
     * Equals {@code orbit.size()}.
     */
    private final int multiplicity;

    /**
     * Constructs a {@code ClusterType} with a given representative, orbit,
     * and multiplicity.
     *
     * @param representative canonical cluster for this type; must not be {@code null}
     * @param orbit          full list of symmetry-equivalent clusters; must not be {@code null}
     * @param multiplicity   size of the orbit (must equal {@code orbit.size()})
     */
    public ClusterType(Cluster representative,
                       List<Cluster> orbit,
                       int multiplicity) {
        this.representative = representative;
        this.orbit          = orbit;
        this.multiplicity   = multiplicity;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the canonical representative cluster for this type.
     *
     * @return representative {@link Cluster}; never {@code null}
     */
    public Cluster getRepresentative() { return representative; }

    /**
     * Returns all symmetry-equivalent clusters forming the orbit of this type.
     *
     * @return unmodifiable view of the orbit; never {@code null}
     */
    public List<Cluster> getOrbit() { return orbit; }

    /**
     * Returns the number of clusters in this symmetry orbit (the multiplicity).
     *
     * @return orbit size â‰¥ 1
     */
    public int getMultiplicity() { return multiplicity; }

    /**
     * Returns a compact string summarising this cluster type.
     *
     * @return e.g. {@code "ClusterType{mult=8, rep=[(0.0,0.0,0.0),s1]}"}
     */
    @Override
    public String toString() {
        return "ClusterType{mult=" + multiplicity
             + ", rep=" + representative + "}";
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this cluster type to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [ClusterType]
     *   multiplicity  : 8
     *   orbit size    : 8
     *   representative: [Cluster with 2 sites]
     *     Sublattice[0] (2 sites)
     *       [0] (0.000000, 0.000000, 0.000000), symbol=s1
     *       ...
     *   orbit[0] : ...
     *   orbit[1] : ...
     * </pre>
     */
    public void printDebug() {
        System.out.println("[ClusterType]");
        System.out.println("  multiplicity : " + multiplicity);
        System.out.println("  orbit size   : " + orbit.size());
        System.out.println("  representative:");
        for (Sublattice sub : representative.getSublattices()) {
            for (Site s : sub.getSites()) {
                System.out.printf("    %s, symbol=%s%n",
                        s.getPosition(), s.getSymbol());
            }
        }
        System.out.println("  orbit members:");
        for (int i = 0; i < orbit.size(); i++) {
            System.out.printf("    orbit[%d] : %s%n", i, orbit.get(i));
        }
    }
}


