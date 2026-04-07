package org.ce.workflow;

import java.util.Arrays;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.CMatrix;
import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.ClusterIdentifier;
import org.ce.domain.cluster.CFIdentifier;
import org.ce.domain.cluster.SpaceGroup;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.storage.InputLoader;

/**
 * Orchestrates the full cluster and CF identification workflow.
 */
public class ClusterIdentificationWorkflow {

    private static final Logger LOG = Logger.getLogger(ClusterIdentificationWorkflow.class.getName());

    private ClusterIdentificationWorkflow() {}

    public static AllClusterData identify(ClusterIdentificationRequest config) {
        return identify(config, null);
    }

    public static AllClusterData identify(ClusterIdentificationRequest config, Consumer<String> progressSink) {
        LOG.info("ClusterIdentificationWorkflow.identify — START");
        emit(progressSink, "TYPE-1a START");

        int numComponents = config.getNumComponents();
        String structurePhase = config.getStructurePhase();
        String model = config.getModel();

        // [STAGE 1a/b]: Load resources
        emit(progressSink, "\n[STAGE 1a]: Loading Disordered clusters...");
        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(config.getDisorderedClusterFile());
        SpaceGroup disorderedSpaceGroup = InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());
        List<SpaceGroup.SymmetryOperation> disorderedSymOps = disorderedSpaceGroup.getOperations();

        emit(progressSink, "[STAGE 1b]: Loading Ordered clusters...");
        List<Cluster> orderedClusters = InputLoader.parseClusterFile(config.getOrderedClusterFile());
        SpaceGroup orderedSpaceGroup = InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());
        List<SpaceGroup.SymmetryOperation> orderedSymOps = orderedSpaceGroup.getOperations();

        // [STAGE 1c]: Identify clusters
        emit(progressSink, "\n[STAGE 1c]: Identifying Clusters (Structural Symmetry)...");
        ClusterIdentificationResult clusterResult = ClusterIdentifier.identify(
                disorderedClusters,
                disorderedSymOps,
                orderedClusters,
                orderedSymOps,
                config.getTransformationMatrix(),
                new double[] { config.getTranslationVector().getX(),
                        config.getTranslationVector().getY(),
                        config.getTranslationVector().getZ() });

        // [STAGE 2]: Identify CFs
        emit(progressSink, "\n[STAGE 2]: Identifying Correlation Function (CF) Orbits...");
        CFIdentificationResult cfResult = CFIdentifier.identify(
                clusterResult,
                clusterResult.getDisClusterData().getClusCoordList(),
                disorderedSymOps,
                orderedClusters,
                orderedSymOps,
                orderedSpaceGroup.getRotateMat(),
                orderedSpaceGroup.getTranslateMat(),
                config.getNumComponents());

        // [STAGE 3]: Orthogonal C-Matrix
        emit(progressSink, "\n[STAGE 3]: Building Orthogonal C-Matrix foundation...");
        CMatrix.Result orthMatrix = CMatrix.buildOrthogonal(
                clusterResult,
                cfResult,
                orderedClusters,
                numComponents);

        // [STAGE 4]: CVCF Transformation
        emit(progressSink, "\n[STAGE 4]: Transforming to CVCF Basis (Thermodynamic Basis)...");
        CvCfBasis cvcfBasis;
        if ("T".equalsIgnoreCase(model)) {
            emit(progressSink, "  [NOTE] Generating transformation matrix dynamically from CV definitions...");
            cvcfBasis = CvCfBasis.dynamic(structurePhase, clusterResult, cfResult, orthMatrix, model);
        } else {
            cvcfBasis = CvCfBasis.Registry.INSTANCE.get(structurePhase, model, numComponents);
        }

        CMatrix.Result cMatrix = org.ce.domain.cluster.cvcf.CvCfBasisTransformer.transform(orthMatrix,
                cvcfBasis);

        List<String> uList = cfResult.getUNames();
        List<String> eOList = cfResult.getEONames();
        List<String> vList = cvcfBasis.cfNames;
        List<String> eList = cvcfBasis.eciNames;

        printTransformationMatrix(cvcfBasis.T, progressSink);

        double[] equiX = new double[numComponents];
        Arrays.fill(equiX, 1.0 / numComponents);
        double[] vRandEqui = cvcfBasis.computeRandomState(equiX, orthMatrix.getCfBasisIndices());

        AllClusterData finalResult = new AllClusterData(
                clusterResult,
                clusterResult,
                cfResult,
                cfResult,
                orthMatrix,
                cMatrix,
                uList,
                vList,
                eOList,
                eList,
                vRandEqui);

        LOG.info("ClusterIdentificationWorkflow.identify — EXIT");
        return finalResult;
    }

    private static void printTransformationMatrix(double[][] matrix, Consumer<String> sink) {
        if (sink == null || matrix == null) return;
        sink.accept("\n  TRANSFORMATION MATRIX (M):");
        for (double[] row : matrix) {
            sink.accept("    " + Arrays.toString(row));
        }
    }

    private static void emit(Consumer<String> sink, String message) {
        if (sink != null) {
            sink.accept(message);
        }
        LOG.fine(message);
    }
}
