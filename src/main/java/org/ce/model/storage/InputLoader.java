package org.ce.model.storage;

import static org.ce.model.cluster.SpaceGroup.SymmetryOperation;

import org.ce.model.cluster.Cluster;
import static org.ce.model.cluster.ClusterPrimitives.*;
import org.ce.model.cluster.SpaceGroup;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Loads CVM input data (cluster geometry and symmetry) from the filesystem
 * workspace ({@code Workspace.inputsDir()}), so that adding a new structure
 * only requires dropping a file into {@code inputs/clus/} or {@code inputs/sym/}
 * — no rebuild/repackage needed.
 *
 * <h2>File naming conventions</h2>
 * <ul>
 *   <li>Cluster files:    {@code inputs/clus/<name>.txt}  (e.g. {@code clus/BCC_A2-T.txt})</li>
 *   <li>Space-group files: {@code inputs/sym/<baseName>.txt}</li>
 * </ul>
 */
public class InputLoader {

    private InputLoader() {}

    // =========================================================================
    // Public API — resolves against the default filesystem Workspace
    // =========================================================================

    /**
     * Parses a cluster file from the workspace filesystem.
     *
     * @param path inputs-dir-relative path (e.g. {@code "clus/BCC_A2-T.txt"})
     * @return list of parsed {@link Cluster} objects
     */
    public static List<Cluster> parseClusterFile(String path) {
        return parseClusterFileFromPath(new Workspace().inputsDir(), path);
    }

    /**
     * Parses a space-group file pair from the workspace filesystem and
     * returns the full {@link SpaceGroup}.
     *
     * @param baseName base name without path or extension (e.g. {@code "BCC_A2-SG"})
     * @return fully populated {@link SpaceGroup}
     */
    public static SpaceGroup parseSpaceGroup(String baseName) {
        return parseSpaceGroupFromPath(new Workspace().inputsDir(), baseName);
    }

    /**
     * Parses a space-group file from the workspace filesystem and returns
     * only the symmetry operations.
     *
     * @param baseName base name without path or extension
     * @return list of {@link SymmetryOperation} objects
     */
    public static List<SymmetryOperation> parseSymmetryFile(String baseName) {
        return parseSpaceGroup(baseName).getOperations();
    }

    // =========================================================================
    // Filesystem-based API — explicit inputsDir (for a non-default Workspace)
    // =========================================================================

    /**
     * Parses a cluster file from the given inputs directory.
     *
     * @param inputsDir   base inputs directory path
     * @param clusterFile relative path within inputs dir (e.g. "clus/BCC_A2-T.txt")
     * @return list of parsed {@link Cluster} objects
     */
    public static List<Cluster> parseClusterFileFromPath(Path inputsDir, String clusterFile) {
        Path filePath = inputsDir.resolve(clusterFile);
        try {
            return ClusterParser.parseFromPath(filePath);
        } catch (Exception e) {
            throw new RuntimeException("Failed to load cluster file: " + filePath, e);
        }
    }

    /**
     * Parses a space-group from the given inputs directory.
     *
     * @param inputsDir base inputs directory path
     * @param baseName  base name without path or extension (e.g. "BCC_A2-SG")
     * @return fully populated {@link SpaceGroup}
     */
    public static SpaceGroup parseSpaceGroupFromPath(Path inputsDir, String baseName) {
        try {
            return SpaceGroupParser.parseFromPath(inputsDir, baseName);
        } catch (Exception e) {
            throw new RuntimeException(
                    "Failed to load space-group file: " + inputsDir.resolve("sym").resolve(baseName + ".txt"), e);
        }
    }

    // -------------------------------------------------------------------------
    // Debug helpers
    // -------------------------------------------------------------------------

    public static void printClusterFileDebug(String path) {
        System.out.println("[InputLoader] cluster file: " + path);
        List<Cluster> clusters = parseClusterFile(path);
        System.out.println("  maximal clusters loaded : " + clusters.size());
        for (int i = 0; i < clusters.size(); i++) {
            System.out.println("  cluster[" + i + "]:");
            clusters.get(i).printDebug();
        }
    }

    public static void printSpaceGroupDebug(String baseName) {
        System.out.println("[InputLoader] space group: " + baseName);
        parseSpaceGroup(baseName).printDebug();
    }

    // =========================================================================
    // ClusterParser — parses Mathematica nested-brace cluster files
    // =========================================================================

    private static final class ClusterParser {

        private ClusterParser() {}

        static List<Cluster> parseFromPath(Path filePath) throws Exception {

            if (!Files.exists(filePath))
                throw new RuntimeException("File not found: " + filePath);

            String content = Files.readString(filePath);
            return parseClusterContent(content);
        }

        private static List<Cluster> parseClusterContent(String content) {
            content = stripLeadingComment(content);
            // Remove outermost braces
            content = content.trim();
            content = content.substring(1, content.length() - 1);

            List<Cluster> clusters = new ArrayList<>();
            int index = 0;
            while (index < content.length()) {
                if (content.charAt(index) == '{') {
                    int end = findMatchingBrace(content, index);
                    clusters.add(parseSingleCluster(content.substring(index + 1, end)));
                    index = end + 1;
                } else {
                    index++;
                }
            }
            return clusters;
        }

        private static Cluster parseSingleCluster(String block) {
            List<Sublattice> sublattices = new ArrayList<>();
            int index = 0;
            while (index < block.length()) {
                if (block.charAt(index) == '{') {
                    int end = findMatchingBrace(block, index);
                    sublattices.add(parseSublattice(block.substring(index + 1, end)));
                    index = end + 1;
                } else {
                    index++;
                }
            }
            return new Cluster(sublattices);
        }

        private static Sublattice parseSublattice(String block) {
            List<Site> sites = new ArrayList<>();
            int index = 0;
            while (index < block.length()) {
                if (block.charAt(index) == '{') {
                    int end = findMatchingBrace(block, index);
                    sites.add(parseSite(block.substring(index + 1, end)));
                    index = end + 1;
                } else {
                    index++;
                }
            }
            return new Sublattice(sites);
        }

        private static Site parseSite(String block) {
            String[] tokens = block.split(",");
            double x = parseMathematicaNumber(tokens[0]);
            double y = parseMathematicaNumber(tokens[1]);
            double z = parseMathematicaNumber(tokens[2]);
            return new Site(new Position(x, y, z), "s1");
        }

        /**
         * Parses a Mathematica numeric literal, including the fraction forms
         * used in some cluster-coordinate exports: {@code "1/3"}, {@code "-1/3"},
         * {@code "-(1/3)"}. Plain decimals ({@code "0.5"}, {@code "-1"}) parse
         * as before.
         */
        static double parseMathematicaNumber(String token) {
            String t = token.trim();
            boolean negated = false;
            if (t.startsWith("-(") && t.endsWith(")")) {
                negated = true;
                t = t.substring(2, t.length() - 1);
            }
            int slash = t.indexOf('/');
            double value;
            if (slash >= 0) {
                double num = Double.parseDouble(t.substring(0, slash).trim());
                double den = Double.parseDouble(t.substring(slash + 1).trim());
                value = num / den;
            } else {
                value = Double.parseDouble(t);
            }
            return negated ? -value : value;
        }

        /**
         * Strips a leading Mathematica {@code (* ... *)} comment block (used
         * by some plain-text notebook exports), including nested comments.
         * Leaves the content untouched if it doesn't start with a comment.
         */
        static String stripLeadingComment(String content) {
            String s = content.stripLeading();
            while (s.startsWith("(*")) {
                int depth = 0;
                int i = 0;
                int end = -1;
                while (i < s.length() - 1) {
                    if (s.charAt(i) == '(' && s.charAt(i + 1) == '*') {
                        depth++;
                        i += 2;
                    } else if (s.charAt(i) == '*' && s.charAt(i + 1) == ')') {
                        depth--;
                        i += 2;
                        if (depth == 0) {
                            end = i;
                            break;
                        }
                    } else {
                        i++;
                    }
                }
                if (end < 0)
                    throw new RuntimeException("Unterminated (* comment *) block");
                s = s.substring(end).stripLeading();
            }
            return s;
        }

        private static int findMatchingBrace(String s, int start) {
            int depth = 0;
            for (int i = start; i < s.length(); i++) {
                if (s.charAt(i) == '{') depth++;
                if (s.charAt(i) == '}') depth--;
                if (depth == 0) return i;
            }
            throw new RuntimeException("Unbalanced braces in cluster file.");
        }
    }

    // =========================================================================
    // SpaceGroupParser — parses symmetry operation and transformation matrix files
    // =========================================================================

    private static final class SpaceGroupParser {

        private SpaceGroupParser() {}

        static SpaceGroup parseFromPath(Path inputsDir, String baseName) throws Exception {

            Path symFile = inputsDir.resolve("sym").resolve(baseName + ".txt");

            if (!Files.exists(symFile))
                throw new RuntimeException("File not found: " + symFile);

            String content = Files.readString(symFile);

            return parseSpaceGroupContent(baseName, content);
        }

        private static SpaceGroup parseSpaceGroupContent(String baseName, String content) {

            // Remove all braces and parse numbers
            String cleanContent = content.replaceAll("\\{", "").replaceAll("}", "");
            String[] tokens = cleanContent.split(",");

            List<Double> numbers = new ArrayList<>();
            for (String t : tokens) {
                if (!t.trim().isEmpty()) numbers.add(ClusterParser.parseMathematicaNumber(t));
            }

            // --- Calculate number of symmetry operations ---
            // Format: (numOps * 12) numbers for ops + 12 numbers for matrix
            int matrixSize = 12;
            int totalNumbers = numbers.size();
            int numOps = (totalNumbers - matrixSize) / matrixSize;

            List<SymmetryOperation> ops = new ArrayList<>();

            for (int i = 0; i < numOps; i++) {
                double[][] rot = new double[3][3];
                double[] trans = new double[3];
                for (int r = 0; r < 3; r++) {
                    for (int c = 0; c < 3; c++) {
                        rot[r][c] = numbers.get(i * matrixSize + r * 4 + c);
                    }
                    trans[r] = numbers.get(i * matrixSize + r * 4 + 3);
                }
                ops.add(new SymmetryOperation(rot, trans));
            }

            // --- Extract rotation + translation matrix (last 12 numbers: 3x3 rotation + 3x1 translation) ---
            int matStartIndex = numOps * matrixSize;
            double[][] rotateMat = new double[3][3];
            double[] translateMat = new double[3];

            for (int i = 0; i < 9; i++) {
                rotateMat[i / 3][i % 3] = numbers.get(matStartIndex + i);
            }
            for (int i = 0; i < 3; i++) {
                translateMat[i] = numbers.get(matStartIndex + 9 + i);
            }

            return new SpaceGroup(baseName, ops, rotateMat, translateMat);
        }
    }
}
