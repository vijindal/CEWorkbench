package org.ce.model.cluster;

import com.fasterxml.jackson.annotation.JsonAlias;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.storage.InputLoader;

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

    @JsonProperty("pipelineResult")
    private final ClusterCFIdentificationPipeline.PipelineResult pipelineResult;

    @JsonProperty("matrixData")
    private final CMatrixPipeline.CMatrixData matrixData;

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
     * @param pipelineResult pipeline result
     * @param matrixData matrix data
     * @param uList u list
     * @param vList v list
     * @param eOList eO list
     * @param eList e list
     * @param equiatomicCVCF random state in CVCF basis at x=1/K
     */
    @JsonCreator
    public AllClusterData(
            @JsonProperty("pipelineResult") ClusterCFIdentificationPipeline.PipelineResult pipelineResult,
            @JsonProperty("matrixData") CMatrixPipeline.CMatrixData matrixData,
            @JsonProperty("uList") @JsonAlias("ulist") List<String> uList,
            @JsonProperty("vList") @JsonAlias("vlist") List<String> vList,
            @JsonProperty("eOList") @JsonAlias("eolist") List<String> eOList,
            @JsonProperty("eList") @JsonAlias("elist") List<String> eList,
            @JsonProperty("equiatomicCVCF") double[] equiatomicCVCF) {

        this.pipelineResult = pipelineResult;
        this.matrixData = matrixData;
        this.uList = uList;
        this.vList = vList;
        this.eOList = eOList;
        this.eList = eList;
        this.equiatomicCVCF = equiatomicCVCF;
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    @JsonIgnore
    public ClusterCFIdentificationPipeline.PipelineResult getPipelineResult() {
        return pipelineResult;
    }

    @JsonIgnore
    public CMatrixPipeline.CMatrixData getMatrixData() {
        return matrixData;
    }

    @JsonProperty("orderedCFResult")
    public CFIdentificationResult getOrderedCFResult() {
        return pipelineResult != null ? pipelineResult.toCFIdentificationResult() : null;
    }

    @JsonProperty("disorderedClusterResult")
    public ClusterIdentificationResult getDisorderedClusterResult() {
        return pipelineResult != null ? pipelineResult.toClusterIdentificationResult() : null;
    }

    @JsonProperty("disorderedCFResult")
    public CFIdentificationResult getDisorderedCFResult() {
        return pipelineResult != null ? pipelineResult.toCFIdentificationResult() : null;
    }

    @JsonProperty("orderedClusterResult")
    public ClusterIdentificationResult getOrderedClusterResult() {
        return pipelineResult != null ? pipelineResult.toClusterIdentificationResult() : null;
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
    public void printSummary(Consumer<String> sink) {
        if (sink == null) return;
        sink.accept("================================================================================");
        sink.accept("                       CLUSTER IDENTIFICATION RESULT");
        sink.accept("================================================================================");

        if (pipelineResult != null) {
            sink.accept(String.format("\nIDENTIFICATION PIPELINE: ncf=%d, total-cfs=%d", 
                    pipelineResult.getNcf(), pipelineResult.getTcf()));
        }

        if (matrixData != null) {
            sink.accept(String.format("C-MATRIX PIPELINE: %d types, %d sites", 
                    matrixData.getCmat().size(), matrixData.getSiteList().size()));
        }

        printStage3Diagnostics(sink);
        printStage4Diagnostics(sink);
        
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
                "AllClusterData: tcdis=%d, tc=%d",
                pipelineResult != null ? pipelineResult.getTcdis() : 0,
                pipelineResult != null ? pipelineResult.getTc() : 0
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
        if (sink == null || matrixData == null || pipelineResult == null) return;
        
        sink.accept("\n  [STAGE 3 DIAGNOSTICS: Random state at equiatomic composition]");
        int K = pipelineResult.getNumComponents();
        double[] x = new double[K];
        java.util.Arrays.fill(x, 1.0 / K);
        CMatrixPipeline.verifyRandomCVs(x, pipelineResult, matrixData, sink);
    }

    private void printStage4Diagnostics(Consumer<String> sink) {
        if (sink == null || matrixData == null || pipelineResult == null || equiatomicCVCF == null) return;

        sink.accept("\n  [STAGE 4 DIAGNOSTICS: Random state equiatomic (CVCF Basis)]");
        // We use the already computed equiatomicCVCF provided to the constructor
        double[][][] cv = CMatrixPipeline.evaluateCVs(
                equiatomicCVCF, 
                matrixData.getCmat(), 
                matrixData.getLcv(), 
                pipelineResult.getTcdis(), 
                pipelineResult.getLc());

        sink.accept("    Resulting CVs (rho_v):");
        for (int t = 0; t < cv.length; t++) {
            for (int j = 0; j < cv[t].length; j++) {
                sink.accept(String.format("      Type %d, Group %d: %s", t, j, Arrays.toString(cv[t][j])));
            }
        }
    }

    // =========================================================================
    // Static Cluster Identification Factory (moved from ClusterIdentificationWorkflow)
    // =========================================================================

    /**
     * Convenience overload — identifies without progress reporting.
     */
    public static AllClusterData identify(ClusterIdentificationRequest config) {
        return identify(config, null);
    }

    /**
     * Orchestrates the full cluster and CF identification workflow.
     *
     * <p>Moved from {@code ClusterIdentificationWorkflow.identify()} to eliminate
     * upward dependency from model layer to calculation layer. All operations
     * are purely model-layer; this method is a static factory.</p>
     *
     * @param config identification request (file paths, matrices, numComponents)
     * @param progressSink optional sink for progress text; may be null
     * @return complete cluster identification result
     */
    public static AllClusterData identify(ClusterIdentificationRequest config, Consumer<String> progressSink) {
        java.util.logging.Logger LOG = java.util.logging.Logger.getLogger(AllClusterData.class.getName());
        LOG.info("AllClusterData.identify — START");

        int numComponents = config.getNumComponents();
        String structurePhase = config.getStructurePhase();
        String model = config.getModel();

        // 1. Load resources
        emit(progressSink, "\n[STAGE 0]: Loading Inputs...");
        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(config.getDisorderedClusterFile());
        disorderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup disorderedSpaceGroup = InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());
        
        List<Cluster> orderedClusters = InputLoader.parseClusterFile(config.getOrderedClusterFile());
        orderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup orderedSpaceGroup = InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());

        // 2. Stage 1 & 2: Cluster + CF Identification
        emit(progressSink, "\n[STAGE 1/2]: Running Identification Pipeline...");
        ClusterCFIdentificationPipeline.PipelineResult pipelineResult = ClusterCFIdentificationPipeline.run(
                disorderedClusters,
                disorderedSpaceGroup.getOperations(),
                orderedClusters,
                orderedSpaceGroup.getOperations(),
                config.getTransformationMatrix(),
                new double[] { config.getTranslationVector().getX(),
                        config.getTranslationVector().getY(),
                        config.getTranslationVector().getZ() },
                numComponents,
                progressSink);

        // 3. Stage 3: C-Matrix foundation
        emit(progressSink, "\n[STAGE 3]: Running C-Matrix Pipeline...");
        CMatrixPipeline.CMatrixData matrixData = CMatrixPipeline.run(
                pipelineResult.toClusterIdentificationResult(),
                pipelineResult.toCFIdentificationResult(),
                disorderedClusters,
                numComponents,
                progressSink);

        // 4. Stage 4: CVCF Transformation
        emit(progressSink, "\n[STAGE 4]: Basis Transformation...");
        CvCfBasis cvcfBasis = CvCfBasis.generate(structurePhase, 
                pipelineResult, 
                matrixData, 
                model, progressSink);

        // Final bundle
        double[] equiX = new double[numComponents];
        Arrays.fill(equiX, 1.0 / numComponents);
        double[] vRandEqui = cvcfBasis.computeRandomState(equiX, matrixData.getCfBasisIndices());

        AllClusterData finalResult = new AllClusterData(
                pipelineResult,
                matrixData,
                pipelineResult.toCFIdentificationResult().getUNames(),
                cvcfBasis.cfNames,
                pipelineResult.toCFIdentificationResult().getEONames(),
                cvcfBasis.eciNames,
                vRandEqui);

        LOG.info("AllClusterData.identify — EXIT");
        return finalResult;
    }

    private static void emit(Consumer<String> sink, String message) {
        if (sink != null) {
            sink.accept(message);
        }
    }
}
