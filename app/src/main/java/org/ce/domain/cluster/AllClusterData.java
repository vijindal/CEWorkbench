package org.ce.domain.cluster;

/**
 * Complete bundle of cluster and correlation function identification results
 * for both disordered and ordered phases.
 *
 * <p>The cluster and CF identification workflows analyze both phases independently,
 * computing separate metrics and results for each. This class maintains the semantic
 * distinction between disordered and ordered phase results.</p>
 *
 * <p>Note: For practical reasons (when disordered and ordered phases are identical,
 * e.g., A2→A2), the disordered and ordered result objects may be the same instance.
 * However, they represent logically distinct phase calculations.</p>
 */
public class AllClusterData {

    private final ClusterIdentificationResult disorderedClusterResult;
    private final ClusterIdentificationResult orderedClusterResult;
    private final CFIdentificationResult disorderedCFResult;
    private final CFIdentificationResult orderedCFResult;
    private final CMatrixResult cMatrixResult;

    /**
     * Creates AllClusterData from identification results for both phases.
     *
     * <p>While a single {@code ClusterIdentificationResult} may internally contain
     * data for both disordered and ordered phases, this constructor accepts them
     * as distinct parameters to maintain architectural clarity about which phase
     * each result represents.</p>
     *
     * @param disorderedClusterResult cluster identification for disordered (HSP) phase
     * @param orderedClusterResult cluster identification for ordered phase
     * @param disorderedCFResult correlation function identification for disordered phase
     * @param orderedCFResult correlation function identification for ordered phase
     * @param cMatrixResult C-matrix identification result
     */
    public AllClusterData(
            ClusterIdentificationResult disorderedClusterResult,
            ClusterIdentificationResult orderedClusterResult,
            CFIdentificationResult disorderedCFResult,
            CFIdentificationResult orderedCFResult,
            CMatrixResult cMatrixResult) {

        this.disorderedClusterResult = disorderedClusterResult;
        this.orderedClusterResult = orderedClusterResult;
        this.disorderedCFResult = disorderedCFResult;
        this.orderedCFResult = orderedCFResult;
        this.cMatrixResult = cMatrixResult;
    }

    /**
     * Default no-arg constructor for deserialization (e.g., JSON loading).
     */
    public AllClusterData() {
        this.disorderedClusterResult = null;
        this.orderedClusterResult = null;
        this.disorderedCFResult = null;
        this.orderedCFResult = null;
        this.cMatrixResult = null;
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    public ClusterIdentificationResult getDisorderedClusterResult() {
        return disorderedClusterResult;
    }

    public ClusterIdentificationResult getOrderedClusterResult() {
        return orderedClusterResult;
    }

    public CFIdentificationResult getDisorderedCFResult() {
        return disorderedCFResult;
    }

    public CFIdentificationResult getOrderedCFResult() {
        return orderedCFResult;
    }

    public CMatrixResult getCMatrixResult() {
        return cMatrixResult;
    }

    // =========================================================================
    // Print Methods
    // =========================================================================

    /**
     * Prints a detailed summary of all identification results.
     */
    public void printResults() {
        System.out.println("\n================================================================================");
        System.out.println("                     CLUSTER IDENTIFICATION RESULTS");
        System.out.println("================================================================================");

        System.out.println("\n──── DISORDERED PHASE (HSP) ────");
        if (disorderedClusterResult != null) {
            disorderedClusterResult.printDebug();
        }

        System.out.println("\n──── ORDERED PHASE ────");
        if (orderedClusterResult != null) {
            orderedClusterResult.printDebug();
        }

        System.out.println("\n================================================================================");
        System.out.println("                  CORRELATION FUNCTION IDENTIFICATION RESULTS");
        System.out.println("================================================================================");

        System.out.println("\n──── DISORDERED PHASE CFs ────");
        printCFSummary(disorderedCFResult);

        System.out.println("\n──── ORDERED PHASE CFs ────");
        printCFSummary(orderedCFResult);
    }

    /**
     * Prints a summary of CF results.
     */
    private void printCFSummary(CFIdentificationResult cfResult) {
        if (cfResult == null) {
            System.out.println("CF result is null");
            return;
        }

        System.out.println("  tcfdis = " + cfResult.getTcfdis() + "  (CF types, excl. empty)");
        System.out.println("  tcf    = " + cfResult.getTcf() + "  (total CF types)");
        System.out.println("  nxcf   = " + cfResult.getNxcf() + "  (point-cluster CFs)");
        System.out.println("  ncf    = " + cfResult.getNcf() + "  (non-point CFs)");

        if (cfResult.getLcf() != null) {
            System.out.println("  lcf array:");
            int[][] lcf = cfResult.getLcf();
            for (int t = 0; t < lcf.length; t++) {
                System.out.printf("    t=%-4d lcf=[", t);
                for (int j = 0; j < lcf[t].length; j++) {
                    System.out.printf("%d", lcf[t][j]);
                    if (j < lcf[t].length - 1) System.out.print(", ");
                }
                System.out.println("]");
            }
        }
    }

    /**
     * Returns a summary string of all results.
     */
    public String getSummary() {
        return String.format(
                "AllClusterData: dis(tcdis=%d, tc=%d) ord(tcdis=%d, tc=%d)",
                disorderedClusterResult != null ? disorderedClusterResult.getTcdis() : 0,
                disorderedClusterResult != null ? disorderedClusterResult.getTc() : 0,
                orderedClusterResult != null ? orderedClusterResult.getTcdis() : 0,
                orderedClusterResult != null ? orderedClusterResult.getTc() : 0
        );
    }
}
