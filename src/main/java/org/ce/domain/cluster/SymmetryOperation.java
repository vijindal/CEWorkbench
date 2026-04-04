package org.ce.domain.cluster;

import java.util.ArrayList;
import java.util.List;

/**
 * An element of the crystallographic space group: a combined rotation and
 * translation acting on fractional coordinates.
 *
 * <p>A symmetry operation transforms a point {@code r} via the affine map</p>
 * <pre>
 *   r' = R Â· r + t
 * </pre>
 * <p>where {@code R} is a 3Ã—3 rotation/reflection matrix and {@code t} is a
 * 3-component translation vector, both expressed in fractional coordinates.</p>
 *
 * <p>The class provides convenience methods to apply the operation to a single
 * {@link Site} or to an entire {@link Cluster}.</p>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     SpaceGroup
 */
public class SymmetryOperation {

    /** 3Ã—3 rotation/reflection matrix in fractional coordinates. */
    private final double[][] rotation;

    /** 3-component translation vector in fractional coordinates. */
    private final double[] translation;

    /**
     * Constructs a symmetry operation from an explicit rotation matrix and
     * translation vector.
     *
     * @param rotation    3Ã—3 matrix (row-major); must not be {@code null}
     * @param translation length-3 translation vector; must not be {@code null}
     */
    public SymmetryOperation(double[][] rotation, double[] translation) {
        this.rotation    = rotation;
        this.translation = translation;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the 3Ã—3 rotation matrix of this operation.
     *
     * @return rotation matrix (row-major); never {@code null}
     */
    public double[][] getRotation() { return rotation; }

    /**
     * Returns the translation vector of this operation.
     *
     * @return length-3 array; never {@code null}
     */
    public double[] getTranslation() { return translation; }

    // -------------------------------------------------------------------------
    // Application
    // -------------------------------------------------------------------------

    /**
     * Applies this symmetry operation to a single site.
     *
     * <p>Computes {@code r' = R Â· r + t} and returns a new site at {@code r'}
     * with the same species symbol as the input.  The symbol is preserved
     * because symmetry operations do not change chemical identity.</p>
     *
     * @param site the site to transform; must not be {@code null}
     * @return a new {@link Site} at the transformed position
     */
    public Site applyToSite(Site site) {
        Position r = site.getPosition();

        double x = rotation[0][0]*r.getX() + rotation[0][1]*r.getY()
                 + rotation[0][2]*r.getZ() + translation[0];

        double y = rotation[1][0]*r.getX() + rotation[1][1]*r.getY()
                 + rotation[1][2]*r.getZ() + translation[1];

        double z = rotation[2][0]*r.getX() + rotation[2][1]*r.getY()
                 + rotation[2][2]*r.getZ() + translation[2];

        return new Site(new Position(x, y, z), site.getSymbol());
    }

    /**
     * Applies this symmetry operation to every site in a cluster.
     *
     * <p>Each sublattice is transformed independently and then
     * {@link Sublattice#sorted()} is called so that the result is
     * in canonical coordinate order.</p>
     *
     * @param cluster the cluster to transform; must not be {@code null}
     * @return a new {@link Cluster} with all sites transformed
     */
    public Cluster applyToCluster(Cluster cluster) {
        List<Sublattice> newSublattices = new ArrayList<>();
        for (Sublattice sub : cluster.getSublattices()) {
            List<Site> newSites = new ArrayList<>();
            for (Site s : sub.getSites()) {
                newSites.add(applyToSite(s));
            }
            newSublattices.add(new Sublattice(newSites).sorted());
        }
        return new Cluster(newSublattices);
    }

    /**
     * Returns a human-readable representation of this symmetry operation.
     *
     * @return string of the form {@code "R=[[...]] t=[...]"}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("SymmetryOperation{R=[");
        for (int r = 0; r < 3; r++) {
            sb.append("[");
            for (int c = 0; c < 3; c++) {
                sb.append(String.format("%.4f", rotation[r][c]));
                if (c < 2) sb.append(", ");
            }
            sb.append("]");
            if (r < 2) sb.append(", ");
        }
        sb.append(String.format("], t=[%.4f, %.4f, %.4f]}",
                translation[0], translation[1], translation[2]));
        return sb.toString();
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this symmetry operation to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [SymmetryOperation]
     *   rotation matrix (R):
     *     [ 1.0000  0.0000  0.0000 ]
     *     [ 0.0000  1.0000  0.0000 ]
     *     [ 0.0000  0.0000  1.0000 ]
     *   translation (t):
     *     [ 0.0000  0.0000  0.0000 ]
     * </pre>
     */
    public void printDebug() {
        System.out.println("[SymmetryOperation]");
        System.out.println("  rotation matrix (R):");
        for (int r = 0; r < 3; r++) {
            System.out.printf("    [ %7.4f  %7.4f  %7.4f ]%n",
                    rotation[r][0], rotation[r][1], rotation[r][2]);
        }
        System.out.println("  translation (t):");
        System.out.printf("    [ %7.4f  %7.4f  %7.4f ]%n",
                translation[0], translation[1], translation[2]);
    }
}


