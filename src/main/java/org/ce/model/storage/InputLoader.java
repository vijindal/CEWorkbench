package org.ce.storage;

import static org.ce.domain.cluster.SpaceGroup.SymmetryOperation;

import org.ce.domain.cluster.Cluster;
import static org.ce.domain.cluster.ClusterPrimitives.*;
import org.ce.domain.cluster.SpaceGroup;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Loads CVM input data (cluster geometry and symmetry) from classpath resources.
 *
 * <h2>Resource naming conventions</h2>
 * <ul>
 *   <li>Cluster files:    {@code clus/<name>.txt}  (e.g. {@code clus/BCC_A2-T.txt})</li>
 *   <li>Space-group files: {@code sym/<baseName>.txt} and {@code sym/<baseName>_mat.txt}</li>
 * </ul>
 */
public class InputLoader {

    private InputLoader() {}

    // =========================================================================
    // Public API
    // =========================================================================

    /**
     * Parses a cluster file from the classpath.
     *
     * @param path classpath-relative path (e.g. {@code "clus/BCC_A2-T.txt"})
     * @return list of parsed {@link Cluster} objects
     */
    public static List<Cluster> parseClusterFile(String path) {
        try {
            return ClusterParser.parseFromResources(path);
        } catch (Exception e) {
            throw new RuntimeException("Failed to load cluster file: " + path, e);
        }
    }

    /**
     * Parses a space-group resource pair and returns the full {@link SpaceGroup}.
     *
     * @param baseName base name without path or extension (e.g. {@code "A2-SG"})
     * @return fully populated {@link SpaceGroup}
     */
    public static SpaceGroup parseSpaceGroup(String baseName) {
        try {
            return SpaceGroupParser.parseFromResources(baseName);
        } catch (Exception e) {
            throw new RuntimeException("Failed to load space-group file: " + baseName, e);
        }
    }

    /**
     * Parses a space-group resource pair and returns only the symmetry operations.
     *
     * @param baseName base name without path or extension
     * @return list of {@link SymmetryOperation} objects
     */
    public static List<SymmetryOperation> parseSymmetryFile(String baseName) {
        return parseSpaceGroup(baseName).getOperations();
    }

    // =========================================================================
    // Filesystem-based API (for development/project-local workspace)
    // =========================================================================

    /**
     * Parses a cluster file from the filesystem (development workspace).
     *
     * @param inputsDir base inputs directory path
     * @param clusterFile relative path within inputs dir (e.g. "clus/A2-T.txt")
     * @return list of parsed {@link Cluster} objects
     */
    public static List<Cluster> parseClusterFileFromPath(Path inputsDir, String clusterFile) {
        try {
            Path filePath = inputsDir.resolve(clusterFile);
            return ClusterParser.parseFromPath(filePath);
        } catch (Exception e) {
            throw new RuntimeException("Failed to load cluster file: " + clusterFile, e);
        }
    }

    /**
     * Parses a space-group from the filesystem (development workspace).
     *
     * @param inputsDir base inputs directory path
     * @param baseName base name without path or extension (e.g. "A2-SG")
     * @return fully populated {@link SpaceGroup}
     */
    public static SpaceGroup parseSpaceGroupFromPath(Path inputsDir, String baseName) {
        try {
            return SpaceGroupParser.parseFromPath(inputsDir, baseName);
        } catch (Exception e) {
            throw new RuntimeException("Failed to load space-group file: " + baseName, e);
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

        static List<Cluster> parseFromResources(String resourcePath) throws Exception {

            InputStream is = ClusterParser.class
                    .getClassLoader()
                    .getResourceAsStream(resourcePath);

            if (is == null)
                throw new RuntimeException("File not found: " + resourcePath);

            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            StringBuilder sb = new StringBuilder();
            String line;
            while ((line = br.readLine()) != null) sb.append(line.trim());
            br.close();

            return parseClusterContent(sb.toString());
        }

        static List<Cluster> parseFromPath(Path filePath) throws Exception {

            if (!Files.exists(filePath))
                throw new RuntimeException("File not found: " + filePath);

            String content = Files.readString(filePath);
            return parseClusterContent(content);
        }

        private static List<Cluster> parseClusterContent(String content) {
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
            double x = Double.parseDouble(tokens[0]);
            double y = Double.parseDouble(tokens[1]);
            double z = Double.parseDouble(tokens[2]);
            return new Site(new Position(x, y, z), "s1");
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

        static SpaceGroup parseFromResources(String baseName) throws Exception {

            // Single file contains both symmetry operations and transformation matrix
            InputStream is = SpaceGroupParser.class
                    .getClassLoader()
                    .getResourceAsStream("sym/" + baseName + ".txt");

            if (is == null)
                throw new RuntimeException("File not found: " + baseName + ".txt");

            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            StringBuilder sb = new StringBuilder();
            String line;
            while ((line = br.readLine()) != null) sb.append(line.trim());
            br.close();

            return parseSpaceGroupContent(baseName, sb.toString());
        }

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
                if (!t.trim().isEmpty()) numbers.add(Double.parseDouble(t.trim()));
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
