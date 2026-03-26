package org.ce.workflow;

import org.ce.domain.cluster.*;
import org.ce.domain.cluster.cvcf.BccA2CvCfTransformations;
import org.ce.storage.InputLoader;

import java.util.List;
import java.util.logging.Logger;

/**
 * Workflow for cluster and correlation function identification.
 *
 * <p>This workflow orchestrates the complete identification pipeline for both
 * disordered (highest-symmetry) and ordered phases:</p>
 * <ol>
 *   <li>Load disordered cluster geometry and symmetry (InputLoader)</li>
 *   <li>Load ordered cluster geometry and symmetry (InputLoader)</li>
 *   <li>Identify clusters in both phases (ClusterIdentifier)</li>
 *   <li>Identify correlation functions in both phases (CFIdentifier)</li>
 *   <li>Return complete results as AllClusterData</li>
 * </ol>
 *
 * @author CVM Project
 * @version 2.0
 */
public class ClusterIdentificationWorkflow {

    private static final Logger LOG = Logger.getLogger(ClusterIdentificationWorkflow.class.getName());

    private ClusterIdentificationWorkflow() {}

    /**
     * Executes the complete cluster and correlation function identification workflow.
     *
     * @param config the identification request with all configuration parameters
     * @return bundle of identification results containing both disordered and ordered data
     * @throws RuntimeException if resources are not found or identification fails
     */
    public static AllClusterData identify(ClusterIdentificationRequest config) {
        LOG.info("ClusterIdentificationWorkflow.identify — START");

        // =====================================================================
        // 1. Load disordered phase (HSP) data
        // =====================================================================
        LOG.fine("Loading disordered phase data...");

        List<Cluster> disorderedClusters =
                InputLoader.parseClusterFile(config.getDisorderedClusterFile());

        SpaceGroup disorderedSpaceGroup =
                InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());

        List<SymmetryOperation> disorderedSymOps =
                disorderedSpaceGroup.getOperations();

        LOG.fine("Disordered clusters loaded: " + disorderedClusters.size());

        // =====================================================================
        // 2. Load ordered phase data
        // =====================================================================
        LOG.fine("Loading ordered phase data...");

        List<Cluster> orderedClusters =
                InputLoader.parseClusterFile(config.getOrderedClusterFile());

        SpaceGroup orderedSpaceGroup =
                InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());

        List<SymmetryOperation> orderedSymOps =
                orderedSpaceGroup.getOperations();

        LOG.fine("Ordered clusters loaded: " + orderedClusters.size());

        // =====================================================================
        // 3. Identify clusters (both disordered and ordered phases)
        // =====================================================================
        LOG.fine("Identifying clusters for both phases...");

        ClusterIdentificationResult clusterResult =
                ClusterIdentifier.identify(
                        disorderedClusters,
                        disorderedSymOps,
                        orderedClusters,
                        orderedSymOps,
                        config.getTransformationMatrix(),
                        new double[]{config.getTranslationVector().getX(),
                                     config.getTranslationVector().getY(),
                                     config.getTranslationVector().getZ()}
                );

        LOG.fine("Clusters identified: tcdis=" + clusterResult.getTcdis() + ", tc=" + clusterResult.getTc());

        // =====================================================================
        // 4. Identify CFs (both disordered and ordered phases)
        // =====================================================================
        LOG.fine("Identifying correlation functions for both phases...");

        CFIdentificationResult cfResult =
                CFIdentifier.identify(
                        clusterResult,
                        clusterResult.getDisClusterData().getClusCoordList(),
                        disorderedSymOps,
                        orderedClusters,
                        orderedSymOps,
                        orderedSpaceGroup.getRotateMat(),
                        orderedSpaceGroup.getTranslateMat(),
                        config.getNumComponents()
                );

        LOG.fine("CFs identified: tcfdis=" + cfResult.getTcfdis() + ", tcf=" + cfResult.getTcf());

        // =====================================================================
        // 5. Build C-Matrix (Stage 3)
        // =====================================================================
        LOG.fine("Building C-Matrix...");

        // Use orderedClusters as maxClusters (full geometry reference)
        CMatrixResult cMatrix = CMatrixBuilder.build(
                clusterResult,
                cfResult,
                orderedClusters,              // ✅ correct choice
                config.getNumComponents(),
                BccA2CvCfTransformations.binaryBasis()  // TODO: select basis by numComponents
        );

        LOG.fine("C-Matrix built successfully");

        // =====================================================================
        // 6. Package and return results
        // =====================================================================
        // Note: ClusterIdentificationResult and CFIdentificationResult each contain
        // analysis of both disordered and ordered phases. We pass them as both
        // the "disordered" and "ordered" parameters to maintain architectural clarity
        // about the 4 distinct inputs (dis cluster file, dis symmetry, ord cluster file,
        // ord symmetry), even though the results derive from a single analysis.
        AllClusterData result = new AllClusterData(
                clusterResult,    // disordered cluster result
                clusterResult,    // ordered cluster result
                cfResult,         // disordered CF result
                cfResult,         // ordered CF result
                cMatrix           // C-Matrix result
        );

        LOG.info("ClusterIdentificationWorkflow.identify — EXIT: " + result.getSummary());

        return result;
    }
}
