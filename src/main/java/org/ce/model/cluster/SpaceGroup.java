package org.ce.model.cluster;

import static org.ce.model.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Encapsulates the full crystallographic space group of a structure,
 * together with the rotation/translation pair that maps the ordered
 * (B2) supercell coordinates into the disordered (A2) reference frame.
 */
public class SpaceGroup {

    /** Human-readable identifier, e.g. {@code "A2-SG"}. */
    private final String name;

    /** All symmetry operations (rotations + translations) of this space group. */
    private final List<SymmetryOperation> operations;

    private final double[][] rotateMat;
    private final double[] translateMat;

    public SpaceGroup(String name,
                      List<SymmetryOperation> operations,
                      double[][] rotateMat,
                      double[] translateMat) {
        this.name         = name;
        this.operations   = operations;
        this.rotateMat    = rotateMat;
        this.translateMat = translateMat;
    }

    public String getName() { return name; }
    public List<SymmetryOperation> getOperations() { return operations; }
    public double[][] getRotateMat() { return rotateMat; }
    public double[] getTranslateMat() { return translateMat; }
    public int order() { return operations.size(); }

    @Override
    public String toString() {
        return "SpaceGroup{name=" + name + ", order=" + operations.size() + "}";
    }

    public void printDebug() {
        System.out.println("[SpaceGroup]");
        System.out.println("  name        : " + name);
        System.out.println("  order       : " + operations.size());
        System.out.println("  rotateMat   :");
        for (int r = 0; r < 3; r++) {
            System.out.printf("    [ %7.4f  %7.4f  %7.4f ]%n",
                    rotateMat[r][0], rotateMat[r][1], rotateMat[r][2]);
        }
        System.out.printf("  translateMat: [ %7.4f  %7.4f  %7.4f ]%n",
                translateMat[0], translateMat[1], translateMat[2]);
        for (int i = 0; i < Math.min(3, operations.size()); i++) {
            System.out.println("    [" + i + "] " + operations.get(i));
        }
    }

    /**
     * An element of the crystallographic space group: a combined rotation and
     * translation acting on fractional coordinates.
     */
    public static class SymmetryOperation {

        private final double[][] rotation;
        private final double[] translation;

        public SymmetryOperation(double[][] rotation, double[] translation) {
            this.rotation    = rotation;
            this.translation = translation;
        }

        public double[][] getRotation() { return rotation; }
        public double[] getTranslation() { return translation; }

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
}
