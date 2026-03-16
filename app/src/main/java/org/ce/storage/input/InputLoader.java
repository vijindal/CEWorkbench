package org.ce.storage.input;

import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.SpaceGroup;
import org.ce.domain.cluster.SymmetryOperation;

import java.util.List;

/**
 * Facade for loading CVM input data from classpath resources.
 *
 * <p>{@code InputLoader} wraps the lower-level parsers ({@link ClusterParser}
 * and {@link SpaceGroupParser}) with a simplified API and uniform exception
 * handling (checked exceptions are wrapped as {@link RuntimeException}).
 * All resources are expected to be present in the application JAR or on the
 * classpath under the {@code cluster/} and {@code symmetry/} directories.</p>
 *
 * <h2>Resource naming conventions</h2>
 * <ul>
 *   <li>Cluster files: {@code cluster/<name>.txt}  (e.g. {@code cluster/A2-T.txt})</li>
 *   <li>Space-group files: {@code symmetry/<name>.txt} and
 *       {@code symmetry/<name>_mat.txt}  (e.g. {@code symmetry/A2-SG.txt})</li>
 * </ul>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     ClusterParser
 * @see     SpaceGroupParser
 */
public class InputLoader {

    /** Private constructor â€” all methods are static utilities. */
    private InputLoader() {}

    // -------------------------------------------------------------------------
    // Cluster loading
    // -------------------------------------------------------------------------

    /**
     * Parses a cluster file from the classpath and returns the maximal clusters.
     *
     * <p>The file must be in the Mathematica nested-brace format used throughout
     * the CVM project (e.g. the output of {@code exportClusCoord} in Mathematica).
     * The path must include the folder prefix, e.g. {@code "cluster/A2-T.txt"}.</p>
     *
     * @param path classpath-relative path including folder and extension;
     *             must not be {@code null}
     * @return list of parsed {@link Cluster} objects; never {@code null}
     * @throws RuntimeException if the resource is not found or the file is malformed
     */
    public static List<Cluster> parseClusterFile(String path) {
        try {
            return ClusterParser.parseFromResources(path);
        } catch (Exception e) {
            throw new RuntimeException("Failed to load cluster file: " + path, e);
        }
    }

    // -------------------------------------------------------------------------
    // Space-group loading
    // -------------------------------------------------------------------------

    /**
     * Parses a space-group resource pair and returns the symmetry operations.
     *
     * <p>Two resource files are read: {@code symmetry/<baseName>.txt}
     * (the operations) and {@code symmetry/<baseName>_mat.txt}
     * (the rotation/translation matrix for coordinate transformations).
     * Supply only the base name without path prefix or extension,
     * e.g. {@code "A2-SG"}.</p>
     *
     * @param baseName base name of the space-group resource (no path, no extension);
     *                 must not be {@code null}
     * @return list of {@link SymmetryOperation} objects; never {@code null}
     * @throws RuntimeException if either resource is not found or is malformed
     */
    public static List<SymmetryOperation> parseSymmetryFile(String baseName) {
        try {
            SpaceGroup group = SpaceGroupParser.parseFromResources(baseName);
            return group.getOperations();
        } catch (Exception e) {
            throw new RuntimeException(
                    "Failed to load space-group file: " + baseName, e);
        }
    }

    /**
     * Parses a space-group resource pair and returns the full {@link SpaceGroup}
     * object, which also includes the ordered-to-disordered transformation matrix.
     *
     * @param baseName base name of the space-group resource; must not be {@code null}
     * @return fully populated {@link SpaceGroup}; never {@code null}
     * @throws RuntimeException if either resource is not found or is malformed
     */
    public static SpaceGroup parseSpaceGroup(String baseName) {
        try {
            return SpaceGroupParser.parseFromResources(baseName);
        } catch (Exception e) {
            throw new RuntimeException(
                    "Failed to load space-group file: " + baseName, e);
        }
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Loads and prints a structured debug summary of a cluster file to standard
     * output.  Useful for verifying that a resource file is parsed correctly
     * before passing it to the main algorithm.
     *
     * <p>Output format:</p>
     * <pre>
     * [InputLoader] cluster file: cluster/A2-T.txt
     *   maximal clusters loaded : 1
     *   cluster[0]:
     *     Sublattice[0] (4 sites)
     *       [0] (0.000000, ...) ...
     * </pre>
     *
     * @param path classpath-relative path to the cluster resource
     */
    public static void printClusterFileDebug(String path) {
        System.out.println("[InputLoader] cluster file: " + path);
        List<Cluster> clusters = parseClusterFile(path);
        System.out.println("  maximal clusters loaded : " + clusters.size());
        for (int i = 0; i < clusters.size(); i++) {
            System.out.println("  cluster[" + i + "]:");
            clusters.get(i).printDebug();
        }
    }

    /**
     * Loads and prints a structured debug summary of a space-group resource to
     * standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [InputLoader] space group: A2-SG
     *   [SpaceGroup]
     *     name  : A2-SG
     *     order : 48
     *     ...
     * </pre>
     *
     * @param baseName base name of the space-group resource (no path, no extension)
     */
    public static void printSpaceGroupDebug(String baseName) {
        System.out.println("[InputLoader] space group: " + baseName);
        SpaceGroup sg = parseSpaceGroup(baseName);
        sg.printDebug();
    }
}


