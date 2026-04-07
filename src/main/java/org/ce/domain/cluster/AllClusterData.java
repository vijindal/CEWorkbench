package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonAlias;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

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
@JsonIgnoreProperties(ignoreUnknown = true)
public class AllClusterData {

    @JsonProperty("disorderedClusterResult")
    private final ClusterIdentificationResult disorderedClusterResult;

    @JsonProperty("orderedClusterResult")
    private final ClusterIdentificationResult orderedClusterResult;

    @JsonProperty("disorderedCFResult")
    private final CFIdentificationResult disorderedCFResult;

    @JsonProperty("orderedCFResult")
    private final CFIdentificationResult orderedCFResult;

    @JsonProperty("orthogonalCMatrixResult")
    private final CMatrix.Result orthogonalCMatrixResult;

    @JsonProperty("cMatrixResult")
    private final CMatrix.Result cMatrixResult;

    // Basis-specific symbol lists
    @JsonProperty("uList")
    private final List<String> uList;  // Orthogonal CF symbols

    @JsonProperty("vList")
    private final List<String> vList;  // CVCF CF symbols

    @JsonProperty("eOList")
    private final List<String> eOList; // Orthogonal CEC symbols

    @JsonProperty("eList")
    private final List<String> eList;  // CVCF CEC symbols

    @JsonProperty("equiatomicCVCF")
    private final double[] equiatomicCVCF; // Random state in CVCF basis

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
     * @param equiatomicCVCF random state in CVCF basis at x=1/K
     */
    @JsonCreator
    public AllClusterData(
            @JsonProperty("disorderedClusterResult") ClusterIdentificationResult disorderedClusterResult,
            @JsonProperty("orderedClusterResult")    ClusterIdentificationResult orderedClusterResult,
            @JsonProperty("disorderedCFResult")      CFIdentificationResult      disorderedCFResult,
            @JsonProperty("orderedCFResult")         CFIdentificationResult      orderedCFResult,
            @JsonProperty("orthogonalCMatrixResult") CMatrix.Result              orthogonalCMatrixResult,
            @JsonProperty("cMatrixResult") @JsonAlias({"cMatrixResult", "cmatrixResult", "CMatrix.Result"}) CMatrix.Result cMatrixResult,
            @JsonProperty("uList") @JsonAlias("ulist") List<String> uList,
            @JsonProperty("vList") @JsonAlias("vlist") List<String> vList,
            @JsonProperty("eOList") @JsonAlias("eolist") List<String> eOList,
            @JsonProperty("eList") @JsonAlias("elist") List<String> eList,
            @JsonProperty("equiatomicCVCF") double[] equiatomicCVCF) {

        this.disorderedClusterResult = disorderedClusterResult;
        this.orderedClusterResult = orderedClusterResult;
        this.disorderedCFResult = disorderedCFResult;
        this.orderedCFResult = orderedCFResult;
        this.orthogonalCMatrixResult = orthogonalCMatrixResult;
        this.cMatrixResult = cMatrixResult;
        this.uList = uList;
        this.vList = vList;
        this.eOList = eOList;
        this.eList = eList;
        this.equiatomicCVCF = equiatomicCVCF;
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    @JsonProperty("disorderedClusterResult")
    public ClusterIdentificationResult getDisorderedClusterResult() {
        return disorderedClusterResult;
    }

    @JsonProperty("orderedClusterResult")
    public ClusterIdentificationResult getOrderedClusterResult() {
        return orderedClusterResult;
    }

    @JsonProperty("disorderedCFResult")
    public CFIdentificationResult getDisorderedCFResult() {
        return disorderedCFResult;
    }

    @JsonProperty("orderedCFResult")
    public CFIdentificationResult getOrderedCFResult() {
        return orderedCFResult;
    }

    /**
     * Returns true if this data bundle includes the Stage 3 orthogonal foundation.
     * Legacy cluster_data.json files may lack this information.
     */
    @JsonIgnore
    public boolean hasOrthogonalFoundation() {
        return orthogonalCMatrixResult != null;
    }

    @JsonProperty("orthogonalCMatrixResult")
    public CMatrix.Result getOrthogonalCMatrixResult() {
        return orthogonalCMatrixResult;
    }

    @JsonProperty("cMatrixResult")
    public CMatrix.Result getCMatrixResult() {
        return cMatrixResult;
    }

    @JsonProperty("uList")
    public List<String> getUList()  { return uList; }

    @JsonProperty("vList")
    public List<String> getVList()  { return vList; }

    @JsonProperty("eOList")
    public List<String> getEOList() { return eOList; }

    @JsonProperty("eList")
    public List<String> getEList()  { return eList; }

    @JsonProperty("equiatomicCVCF")
    public double[] getEquiatomicCVCF() { return equiatomicCVCF; }

    // =========================================================================
    // Print Methods
    // =========================================================================

    /**
     * Prints a detailed summary of all identification results.
     */
    public void printSummary(java.util.function.Consumer<String> sink) {
        sink.accept("================================================================================");
        sink.accept("                     ALL CLUSTER IDENTIFICATION DATA");
        sink.accept("================================================================================");

        if (disorderedClusterResult != null) {
            disorderedClusterResult.printSummary(sink);
        }

        if (disorderedCFResult != null) {
            disorderedCFResult.printSummary(sink);
        }

        if (orthogonalCMatrixResult != null) {
            orthogonalCMatrixResult.printSummary("Orthogonal Basis (Structural Foundation)", sink);
            printStage3Diagnostics(sink);
        }

        if (cMatrixResult != null) {
            cMatrixResult.printSummary("Final CVCF Result (Minimized Hamiltonian Basis)", sink);
        }

        printStage4Diagnostics(sink);

        sink.accept("\nSummary of Basis Specific Symbol Lists:");
        sink.accept("  - uList (OrthCFs): " + uList);
        sink.accept("  - vList (CVCFs):   " + vList);
        sink.accept("  - eOList (OrthCECs):" + eOList);
        sink.accept("  - eList (CVCFs):   " + eList);
        sink.accept("================================================================================");
    }

    public void printResults() {
        printSummary(System.out::println);
    }

    /**
     * Returns a summary string of all results.
     */
    @JsonIgnore
    public String getSummary() {
        return String.format(
                "AllClusterData: dis(tcdis=%d, tc=%d) ord(tcdis=%d, tc=%d)",
                disorderedClusterResult != null ? disorderedClusterResult.getTcdis() : 0,
                disorderedClusterResult != null ? disorderedClusterResult.getTc() : 0,
                orderedClusterResult != null ? orderedClusterResult.getTcdis() : 0,
                orderedClusterResult != null ? orderedClusterResult.getTc() : 0
        );
    }

    /**
     * Immutable bundle of cluster identification and correlation function identification results
     * for a single phase. Primarily used for legacy compatibility in MCS engine.
     */
    public static class ClusterData {
        private final ClusterIdentificationResult clusterResult;
        private final CFIdentificationResult cfResult;
        public ClusterData(ClusterIdentificationResult clusterResult, CFIdentificationResult cfResult) {
            this.clusterResult = clusterResult;
            this.cfResult = cfResult;
        }
        public ClusterIdentificationResult getClusterResult() { return clusterResult; }
        public CFIdentificationResult getCfResult() { return cfResult; }
    }

    private void printStage3Diagnostics(Consumer<String> sink) {
        if (sink == null || orthogonalCMatrixResult == null || orderedCFResult == null || disorderedClusterResult == null) {
            return;
        }

        sink.accept("\n  [STAGE 3 DIAGNOSTICS: Random state at equiatomic composition]");
        
        int numComponents = orderedCFResult.getNxcf() + 1;
        
        // 1. Equiatomic mole fractions
        double[] x = new double[numComponents];
        Arrays.fill(x, 1.0 / numComponents);
        sink.accept(String.format("    Composition: x = %s", Arrays.toString(x)));

        // 2. Random CFs (uRand)
        double[] uRand = ClusterVariableEvaluator.computeRandomCFs(
                x, numComponents, orthogonalCMatrixResult.getCfBasisIndices(), 
                orderedCFResult.getNcf(), orderedCFResult.getTcf());
        
        sink.accept("    Random CFs (u):");
        for (int i = 0; i < uRand.length; i++) {
            sink.accept(String.format("      u%-3d = %11.4e", i, uRand[i]));
        }

        // 3. Full CF vector (including point CFs)
        double[] uFull = ClusterVariableEvaluator.buildFullCFVector(
                uRand, x, numComponents, orthogonalCMatrixResult.getCfBasisIndices(), 
                orderedCFResult.getNcf(), orderedCFResult.getTcf());

        // 4. Evaluate CVs (C * uFull)
        double[][][] cv = ClusterVariableEvaluator.evaluate(
                uFull, orthogonalCMatrixResult.getCmat(), orthogonalCMatrixResult.getLcv(), 
                disorderedClusterResult.getTcdis(), disorderedClusterResult.getLc());

        sink.accept("    Resulting CVs (Cluster Variables):");
        for (int t = 0; t < cv.length; t++) {
            for (int j = 0; j < cv[t].length; j++) {
                sink.accept(String.format("      Block (t=%d, j=%d): %s", t, j, Arrays.toString(cv[t][j])));
            }
        }
    }

    private void printStage4Diagnostics(Consumer<String> sink) {
        if (sink == null || equiatomicCVCF == null || cMatrixResult == null || orderedCFResult == null || disorderedClusterResult == null) {
            return;
        }

        sink.accept("\n  [STAGE 4 DIAGNOSTICS: Random state at equiatomic composition (CVCF Basis)]");
        int numComponents = orderedCFResult.getNxcf() + 1;
        double[] x = new double[numComponents];
        Arrays.fill(x, 1.0 / numComponents);
        sink.accept(String.format("    Composition: x = %s", Arrays.toString(x)));

        // Random CFs in CVCF basis (v)
        sink.accept("    Random CFs (v):");
        for (int i = 0; i < equiatomicCVCF.length; i++) {
            sink.accept(String.format("      v%-3d = %11.4e", i, equiatomicCVCF[i]));
        }

        // Evaluate CVs (C * v)
        // Since equiatomicCVCF already includes point variables (x_i), 
        // we use it as the full vector.
        double[][][] cv = ClusterVariableEvaluator.evaluate(
                equiatomicCVCF, cMatrixResult.getCmat(), cMatrixResult.getLcv(), 
                disorderedClusterResult.getTcdis(), disorderedClusterResult.getLc());

        sink.accept("    Resulting CVs (Cluster Variables):");
        for (int t = 0; t < cv.length; t++) {
            for (int j = 0; j < cv[t].length; j++) {
                sink.accept(String.format("      Block (t=%d, j=%d): %s", t, j, Arrays.toString(cv[t][j])));
            }
        }
    }
}
