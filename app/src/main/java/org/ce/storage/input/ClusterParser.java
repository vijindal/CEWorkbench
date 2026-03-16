package org.ce.storage.input;

import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.Sublattice;
import org.ce.domain.cluster.Site;
import org.ce.domain.cluster.Position;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class ClusterParser {

    // =========================================
    // Parses nested cluster structure
    // Preserves maximal cluster â†’ sublattice â†’ sites
    // =========================================
    public static List<Cluster> parseFromResources(String resourcePath)
            throws Exception {

        InputStream is = ClusterParser.class
                .getClassLoader()
                .getResourceAsStream(resourcePath);

        if (is == null)
            throw new RuntimeException("File not found: " + resourcePath);

        BufferedReader br =
                new BufferedReader(new InputStreamReader(is));

        StringBuilder sb = new StringBuilder();
        String line;

        while ((line = br.readLine()) != null) {
            sb.append(line.trim());
        }

        br.close();

        String content = sb.toString().trim();

        // Remove outermost braces
        content = content.substring(1, content.length() - 1);

        List<Cluster> clusters = new ArrayList<>();

        int index = 0;

        while (index < content.length()) {

            if (content.charAt(index) == '{') {

                int end = findMatchingBrace(content, index);
                String clusterBlock =
                        content.substring(index + 1, end);

                clusters.add(parseSingleCluster(clusterBlock));

                index = end + 1;
            } else {
                index++;
            }
        }

        return clusters;
    }

    // =========================================
    // Parses one maximal cluster
    // =========================================
    private static Cluster parseSingleCluster(String block) {

        List<Sublattice> sublattices = new ArrayList<>();

        int index = 0;

        while (index < block.length()) {

            if (block.charAt(index) == '{') {

                int end = findMatchingBrace(block, index);
                String subBlock =
                        block.substring(index + 1, end);

                sublattices.add(parseSublattice(subBlock));

                index = end + 1;
            } else {
                index++;
            }
        }

        return new Cluster(sublattices);
    }

    // =========================================
    // Parses one sublattice
    // =========================================
    private static Sublattice parseSublattice(String block) {

        List<Site> sites = new ArrayList<>();

        int index = 0;

        while (index < block.length()) {

            if (block.charAt(index) == '{') {

                int end = findMatchingBrace(block, index);
                String siteBlock =
                        block.substring(index + 1, end);

                sites.add(parseSite(siteBlock));

                index = end + 1;
            } else {
                index++;
            }
        }

        return new Sublattice(sites);
    }

    // =========================================
    // Parses one site coordinate
    // Default symbol = "s1"
    // =========================================
    private static Site parseSite(String block) {

        String[] tokens = block.split(",");

        double x = Double.parseDouble(tokens[0]);
        double y = Double.parseDouble(tokens[1]);
        double z = Double.parseDouble(tokens[2]);

        return new Site(new Position(x, y, z), "s1");
    }

    // =========================================
    // Utility to find matching brace
    // =========================================
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



