package org.ce.workflow;

import static org.ce.domain.cluster.SpaceGroup.SymmetryOperation;

import org.ce.domain.cluster.*;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.storage.InputLoader;
import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Workflow for cluster and correlation function identification.
 *
 * <p>
 * This workflow orchestrates the complete identification pipeline for both
 * disordered (highest-symmetry) and ordered phases:
 * </p>
 * <ol>
 * <li>Load disordered cluster geometry and symmetry (InputLoader)</li>
 * <li>Load ordered cluster geometry and symmetry (InputLoader)</li>
 * <li>Identify clusters in both phases (ClusterIdentifier)</li>
 * <li>Identify correlation functions in both phases (CFIdentifier)</li>
 * <li>Return complete results as AllClusterData</li>
 * </ol>
 *
 * @author CVM Project
 * @version 2.0
 */
public class ClusterIdentificationWorkflow {

        private static final Logger LOG = Logger.getLogger(ClusterIdentificationWorkflow.class.getName());

        private ClusterIdentificationWorkflow() {
        }

        /**
         * Executes the complete cluster and correlation function identification
         * workflow.
         *
         * @param config the identification request with all configuration parameters
         * @return bundle of identification results containing both disordered and
         *         ordered data
         * @throws RuntimeException if resources are not found or identification fails
         */
        public static AllClusterData identify(ClusterIdentificationRequest config) {
                return identify(config, null);
        }

        public static AllClusterData identify(ClusterIdentificationRequest config, Consumer<String> progressSink) {
                LOG.info("ClusterIdentificationWorkflow.identify — START");
                emit(progressSink, "TYPE-1a START");
                emit(progressSink, "  orderedClusterFile    : " + config.getOrderedClusterFile());
                emit(progressSink, "  orderedSymmetryGroup  : " + config.getOrderedSymmetryGroup());
                emit(progressSink, "  disorderedClusterFile : " + config.getDisorderedClusterFile());
                emit(progressSink, "  disorderedSymmetryGroup: " + config.getDisorderedSymmetryGroup());
                emit(progressSink, "  numComponents         : " + config.getNumComponents());

                // Validate CVCF basis is supported for this (structure, numComponents)
                String structurePhase = config.getStructurePhase();
                int numComponents = config.getNumComponents();
                if (!CvCfBasis.Registry.INSTANCE.isSupported(structurePhase, numComponents)) {
                        String msg = "CVCF basis not supported for " + structurePhase
                                        + " with " + numComponents + " components.\n"
                                        + CvCfBasis.Registry.INSTANCE.supportedSummary();
                        LOG.severe(msg);
                        emit(progressSink, "ERROR: " + msg);
                        throw new IllegalArgumentException(msg);
                }
                emit(progressSink, "  basis supported: " + structurePhase + ", K=" + numComponents);

                // =====================================================================
                // 1. Load disordered phase (HSP) data
                // =====================================================================
                LOG.fine("Loading disordered phase data...");
                emit(progressSink, "STAGE 1a: Load disordered phase geometry/symmetry");

                List<Cluster> disorderedClusters = InputLoader.parseClusterFile(config.getDisorderedClusterFile());

                SpaceGroup disorderedSpaceGroup = InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());

                List<SymmetryOperation> disorderedSymOps = disorderedSpaceGroup.getOperations();

                LOG.fine("Disordered clusters loaded: " + disorderedClusters.size());
                emit(progressSink, "  disordered max clusters: " + disorderedClusters.size());
                emit(progressSink, "  disordered sym ops     : " + disorderedSymOps.size());

                // =====================================================================
                // 2. Load ordered phase data
                // =====================================================================
                LOG.fine("Loading ordered phase data...");
                emit(progressSink, "STAGE 1b: Load ordered phase geometry/symmetry");

                List<Cluster> orderedClusters = InputLoader.parseClusterFile(config.getOrderedClusterFile());

                SpaceGroup orderedSpaceGroup = InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());

                List<SymmetryOperation> orderedSymOps = orderedSpaceGroup.getOperations();

                LOG.fine("Ordered clusters loaded: " + orderedClusters.size());
                emit(progressSink, "  ordered max clusters   : " + orderedClusters.size());
                emit(progressSink, "  ordered sym ops        : " + orderedSymOps.size());

                // =====================================================================
                // 3. Identify clusters (both disordered and ordered phases)
                // =====================================================================
                LOG.fine("Identifying clusters for both phases...");
                emit(progressSink, "STAGE 1c: Cluster identification");

                ClusterIdentificationResult clusterResult = ClusterIdentifier.identify(
                                disorderedClusters,
                                disorderedSymOps,
                                orderedClusters,
                                orderedSymOps,
                                config.getTransformationMatrix(),
                                new double[] { config.getTranslationVector().getX(),
                                                config.getTranslationVector().getY(),
                                                config.getTranslationVector().getZ() });

                LOG.fine("Clusters identified: tcdis=" + clusterResult.getTcdis() + ", tc=" + clusterResult.getTc());
                emit(progressSink, "  cluster result: tcdis=" + clusterResult.getTcdis()
                                + ", tc=" + clusterResult.getTc() + ", nxcdis=" + clusterResult.getNxcdis()
                                + ", nxc=" + clusterResult.getNxc());
                emit(progressSink, "  lc            : " + Arrays.toString(clusterResult.getLc()));
                emit(progressSink, "  kb            : " + Arrays.toString(clusterResult.getKbCoefficients()));

                // =====================================================================
                // 4. Identify CFs (both disordered and ordered phases)
                // =====================================================================
                LOG.fine("Identifying correlation functions for both phases...");
                emit(progressSink, "STAGE 2: CF identification");

                CFIdentificationResult cfResult = CFIdentifier.identify(
                                clusterResult,
                                clusterResult.getDisClusterData().getClusCoordList(),
                                disorderedSymOps,
                                orderedClusters,
                                orderedSymOps,
                                orderedSpaceGroup.getRotateMat(),
                                orderedSpaceGroup.getTranslateMat(),
                                config.getNumComponents());

                LOG.fine("CFs identified: tcfdis=" + cfResult.getTcfdis() + ", tcf=" + cfResult.getTcf());
                emit(progressSink, "  CF result: tcfdis=" + cfResult.getTcfdis()
                                + ", tcf=" + cfResult.getTcf()
                                + ", ncf=" + cfResult.getNcf()
                                + ", nxcf=" + cfResult.getNxcf());
                emit(progressSink, "  EXPECTED c-matrix columns = tcf + 1 = " + (cfResult.getTcf() + 1)
                                + " (last column is constant term)");
                emit(progressSink, "  lcf table:");
                dump2DInt("    lcf", cfResult.getLcf(), progressSink);

                // =====================================================================
                // 5. Build C-Matrix (Stage 3)
                // =====================================================================
                LOG.fine("Building C-Matrix...");
                emit(progressSink, "STAGE 3: Build C-matrix (orthogonal -> CVCF)");

                // Use orderedClusters as maxClusters (full geometry reference)
                CMatrix.Result cMatrix = CMatrix.build(
                                clusterResult,
                                cfResult,
                                orderedClusters, // ✅ correct choice
                                numComponents,
                                CvCfBasis.Registry.INSTANCE.get(config.getStructurePhase(), numComponents),
                                progressSink);

                LOG.fine("C-Matrix built successfully");
                emit(progressSink, "  c-matrix built");
                emit(progressSink, "  lcv table:");
                dump2DInt("    lcv", cMatrix.getLcv(), progressSink);
                emit(progressSink, "  (Above: both orthogonal and CVCF full c-matrix blocks)");

                // =====================================================================
                // 6. Generate basis-specific symbol lists
                // =====================================================================
                LOG.fine("Retrieving symbolic lists for CFs and CECs from Stage 2 result...");
                CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get(config.getStructurePhase(), numComponents);

                List<String> uList  = cfResult.getUNames();
                List<String> eOList = cfResult.getEONames();
                List<String> vList  = basis.cfNames;
                List<String> eList  = basis.eciNames;

                // =====================================================================
                // 7. Package and return results
                // =====================================================================
                // Note: ClusterIdentificationResult and CFIdentificationResult each contain
                // analysis of both disordered and ordered phases. We pass them as both
                // the "disordered" and "ordered" parameters to maintain architectural clarity
                // about the 4 distinct inputs (dis cluster file, dis symmetry, ord cluster
                // file,
                // ord symmetry), even though the results derive from a single analysis.
                AllClusterData finalResult = new AllClusterData(
                                clusterResult, // disordered cluster result
                                clusterResult, // ordered cluster result
                                cfResult, // disordered CF result
                                cfResult, // ordered CF result
                                cMatrix // C-Matrix result
                                , uList, vList, eOList, eList);

                LOG.info("ClusterIdentificationWorkflow.identify — EXIT: " + finalResult.getSummary());
                emit(progressSink, "TYPE-1a DONE: " + finalResult.getSummary());

                return finalResult;
        }

        private static void dump2DInt(String label, int[][] table, Consumer<String> sink) {
                if (table == null) {
                        emit(sink, label + " = null");
                        return;
                }
                for (int i = 0; i < table.length; i++) {
                        emit(sink, label + "[" + i + "] = " + Arrays.toString(table[i]));
                }
        }

        private static void emit(Consumer<String> sink, String line) {
                if (sink != null) {
                        sink.accept(line);
                }
        }
}
