package org.ce.domain.cluster;

/**
 * Immutable bundle of cluster identification and correlation function identification results.
 */
public class ClusterData {

    private final ClusterIdentificationResult clusterResult;
    private final CFIdentificationResult cfResult;

    /**
     * Creates a ClusterData bundle from identification results.
     *
     * @param clusterResult result from Stage 1 (cluster identification)
     * @param cfResult result from Stage 2 (correlation function identification)
     */
    public ClusterData(
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult) {
        this.clusterResult = clusterResult;
        this.cfResult = cfResult;
    }

    public ClusterIdentificationResult getClusterResult() {
        return clusterResult;
    }

    public CFIdentificationResult getCfResult() {
        return cfResult;
    }

    /**
     * Prints a detailed summary of all cluster identification and CF results
     * to standard output.
     */
    public void printResults() {
        printClusterResults();
        System.out.println();
        printCFResults();
    }

    /**
     * Prints the cluster identification results (Stage 1).
     */
    public void printClusterResults() {
        if (clusterResult != null) {
            clusterResult.printDebug();
        } else {
            System.out.println("Cluster result is null");
        }
    }

    /**
     * Prints the correlation function identification results (Stage 2).
     */
    public void printCFResults() {
        if (cfResult == null) {
            System.out.println("CF result is null");
            return;
        }

        System.out.println("================================================================================");
        System.out.println("                   CORRELATION FUNCTION IDENTIFICATION RESULT");
        System.out.println("================================================================================");

        System.out.println("\n------ Stage 2a: HSP CFs ------");
        System.out.println("  tcfdis = " + cfResult.getTcfdis() + "  (HSP CF types, excl. empty)");

        System.out.println("\n------ Stage 2b: Phase CFs (classified) ------");
        System.out.println("  tcf   = " + cfResult.getTcf() + "  (total CF types)");
        System.out.println("  nxcf  = " + cfResult.getNxcf() + "  (point-cluster CFs)");
        System.out.println("  ncf   = " + cfResult.getNcf() + "  (non-point CFs)");

        System.out.println("\n------ lcf array (CFs per HSP cluster type) ------");
        if (cfResult.getLcf() != null) {
            int[][] lcf = cfResult.getLcf();
            for (int t = 0; t < lcf.length; t++) {
                System.out.printf("  t=%-4d lcf=", t);
                for (int j = 0; j < lcf[t].length; j++) {
                    System.out.printf("%d ", lcf[t][j]);
                }
                System.out.println();
            }
        }

        System.out.println("================================================================================");
    }

    /**
     * Returns a summary string of the cluster results.
     */
    public String getSummary() {
        if (clusterResult == null) {
            return "Cluster result: null";
        }
        return String.format(
                "Cluster Results: tcdis=%d, tc=%d, nxcf=%d",
                clusterResult.getTcdis(),
                clusterResult.getTc(),
                cfResult != null ? cfResult.getNxcf() : 0
        );
    }
}
