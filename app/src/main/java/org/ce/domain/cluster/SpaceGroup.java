package org.ce.domain.cluster;

import java.util.List;

/**
 * Encapsulates the full crystallographic space group of a structure,
 * together with the rotation/translation pair that maps the ordered
 * (B2) supercell coordinates into the disordered (A2) reference frame.
 *
 * <p>A {@code SpaceGroup} holds:</p>
 * <ul>
 *   <li>A list of {@link SymmetryOperation} objects (rotation + translation
 *       pairs in fractional coordinates) that generate the full site symmetry.</li>
 *   <li>An optional affine transformation ({@code rotateMat}, {@code translateMat})
 *       used to map ordered cluster coordinates into the disordered reference
 *       lattice â€” required for the B2â†’A2 orbit equivalence checks.</li>
 * </ul>
 *
 * <p>Instances are created by {@link org.ce.storage.input.SpaceGroupParser} which reads
 * the {@code symmetry/*.txt} resource files in Mathematica brace-notation
 * format.</p>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     SymmetryOperation
 * @see     org.ce.storage.input.SpaceGroupParser
 */
public class SpaceGroup {

    /** Human-readable identifier, e.g. {@code "A2-SG"}. */
    private final String name;

    /** All symmetry operations (rotations + translations) of this space group. */
    private final List<SymmetryOperation> operations;

    /**
     * 3Ã—3 rotation matrix (row-major) mapping ordered â†’ disordered coordinates.
     * Read from the corresponding {@code *_mat.txt} resource file.
     */
    private final double[][] rotateMat;

    /**
     * 3-component translation vector mapping ordered â†’ disordered coordinates.
     * Read from the corresponding {@code *_mat.txt} resource file.
     */
    private final double[] translateMat;

    /**
     * Constructs a {@code SpaceGroup} with all required data.
     *
     * @param name         human-readable label; must not be {@code null}
     * @param operations   list of symmetry operations; must not be {@code null}
     * @param rotateMat    3Ã—3 ordered-to-disordered rotation matrix; must not be {@code null}
     * @param translateMat length-3 ordered-to-disordered translation vector; must not be {@code null}
     */
    public SpaceGroup(String name,
                      List<SymmetryOperation> operations,
                      double[][] rotateMat,
                      double[] translateMat) {
        this.name         = name;
        this.operations   = operations;
        this.rotateMat    = rotateMat;
        this.translateMat = translateMat;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the human-readable name of this space group.
     *
     * @return non-null name string
     */
    public String getName() { return name; }

    /**
     * Returns all symmetry operations of this space group.
     *
     * @return list of {@link SymmetryOperation}; never {@code null}
     */
    public List<SymmetryOperation> getOperations() { return operations; }

    /**
     * Returns the 3Ã—3 rotation matrix for the ordered-to-disordered
     * coordinate transformation.
     *
     * @return 3Ã—3 double array (row-major); never {@code null}
     */
    public double[][] getRotateMat() { return rotateMat; }

    /**
     * Returns the translation vector for the ordered-to-disordered
     * coordinate transformation.
     *
     * @return length-3 double array; never {@code null}
     */
    public double[] getTranslateMat() { return translateMat; }

    /**
     * Returns the number of symmetry operations in this space group.
     *
     * @return order of the space group â‰¥ 1
     */
    public int order() { return operations.size(); }

    /**
     * Returns a compact string describing this space group.
     *
     * @return e.g. {@code "SpaceGroup{name=A2-SG, order=48}"}
     */
    @Override
    public String toString() {
        return "SpaceGroup{name=" + name + ", order=" + operations.size() + "}";
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this space group to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [SpaceGroup]
     *   name        : A2-SG
     *   order       : 48
     *   rotateMat   :
     *     [ 1.0  0.0  0.0 ]
     *     [ 0.0  1.0  0.0 ]
     *     [ 0.0  0.0  1.0 ]
     *   translateMat: [ 0.0  0.0  0.0 ]
     *   operations  : (first 3 shown)
     *     [0] SymmetryOperation{...}
     *     ...
     * </pre>
     */
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
        System.out.println("  operations  : (showing first " +
                Math.min(3, operations.size()) + " of " + operations.size() + ")");
        for (int i = 0; i < Math.min(3, operations.size()); i++) {
            System.out.println("    [" + i + "] " + operations.get(i));
        }
    }
}


