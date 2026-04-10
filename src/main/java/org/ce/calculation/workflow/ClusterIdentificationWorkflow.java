package org.ce.calculation.workflow;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.CFIdentificationResult;
import org.ce.model.cluster.CMatrix;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterIdentificationResult;
import org.ce.model.cluster.ClusterIdentifier;
import org.ce.model.cluster.CFIdentifier;
import org.ce.model.cluster.SpaceGroup;
import org.ce.model.cluster.cvcf.CvCfBasis;
import org.ce.model.storage.InputLoader;

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

        // [STAGE 0a/b]: Load resources
        emit(progressSink, "\n[STAGE 0a]: Loading Disordered clusters...");
        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(config.getDisorderedClusterFile());
        emit(progressSink, "  Loaded: " + config.getDisorderedClusterFile() + " — OK (" + disorderedClusters.size() + " clusters)");
        SpaceGroup disorderedSpaceGroup = InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());
        List<SpaceGroup.SymmetryOperation> disorderedSymOps = disorderedSpaceGroup.getOperations();
        emit(progressSink, "  Loaded: " + config.getDisorderedSymmetryGroup() + " — OK (" + disorderedSymOps.size() + " sym ops)");

        emit(progressSink, "[STAGE 0b]: Loading Ordered clusters...");
        List<Cluster> orderedClusters = InputLoader.parseClusterFile(config.getOrderedClusterFile());
        emit(progressSink, "  Loaded: " + config.getOrderedClusterFile() + " — OK (" + orderedClusters.size() + " clusters)");
        SpaceGroup orderedSpaceGroup = InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());
        List<SpaceGroup.SymmetryOperation> orderedSymOps = orderedSpaceGroup.getOperations();
        emit(progressSink, "  Loaded: " + config.getOrderedSymmetryGroup() + " — OK (" + orderedSymOps.size() + " sym ops)");

        // [STAGE 1a/1b]: Identify clusters
        emit(progressSink, "\n[STAGE 1a/1b]: Identifying Clusters (Structural Symmetry)...");
        
        if (config.getTranslationVector() == null) {
            throw new IllegalStateException("translationVector is null in identify(config)");
        }

        ClusterIdentificationResult clusterResult = ClusterIdentifier.identify(
                disorderedClusters,
                disorderedSymOps,
                orderedClusters,
                orderedSymOps,
                config.getTransformationMatrix(),
                new double[] { config.getTranslationVector().getX(),
                        config.getTranslationVector().getY(),
                        config.getTranslationVector().getZ() });
        
        emit(progressSink, String.format("  [STAGE 2a OK] Clusters identified: tcdis=%d, nxcdis=%d, tc=%d, nxc=%d, kb[0]=%.4f",
                clusterResult.getTcdis(), clusterResult.getNxcdis(),
                clusterResult.getTc(), clusterResult.getNxc(),
                clusterResult.getKbCoefficients()[0]));

        // [STAGE 2a/2b]: Identify CFs
        emit(progressSink, "\n[STAGE 2a/2b]: Identifying Correlation Function (CF) Orbits...");
        CFIdentificationResult cfResult = CFIdentifier.identify(
                clusterResult,
                clusterResult.getDisClusterData().getClusCoordList(),
                disorderedSymOps,
                orderedClusters,
                orderedSymOps,
                orderedSpaceGroup.getRotateMat(),
                orderedSpaceGroup.getTranslateMat(),
                config.getNumComponents());
        
        emit(progressSink, String.format("  [STAGE 2b OK] CFs identified: tcfdis=%d, tcf=%d, nxcf=%d, ncf=%d, uNames=%d, eoNames=%d",
                cfResult.getTcfdis(), cfResult.getTcf(), cfResult.getNxcf(), cfResult.getNcf(),
                cfResult.getUNames().size(), cfResult.getEONames().size()));

        // [STAGE 3]: Orthogonal C-Matrix
        emit(progressSink, "\n[STAGE 3]: Building Orthogonal C-Matrix foundation...");
        CMatrix.Result orthMatrix = CMatrix.buildOrthogonal(
                clusterResult,
                cfResult,
                orderedClusters,
                numComponents);
        
        int[][] lcv = orthMatrix.getLcv();
        java.util.List<java.util.List<int[]>> wcv = orthMatrix.getWcv();
        int totalCvRows = 0;
        for (int[] lcvt : lcv) for (int v : lcvt) totalCvRows += v;
        emit(progressSink, String.format("  [STAGE 3 OK] C-Matrix built: %d cluster types, %d total CV rows", lcv.length, totalCvRows));

        // [STAGE 4]: CVCF Transformation
        emit(progressSink, "\n[STAGE 4]: Transforming to CVCF Basis (Thermodynamic Basis)...");
        CvCfBasis cvcfBasis = CvCfBasis.dynamic(structurePhase, clusterResult, cfResult, orthMatrix, model, progressSink);

        CMatrix.Result cMatrix = org.ce.model.cluster.cvcf.CvCfBasisTransformer.transform(orthMatrix,
                cvcfBasis);

        List<String> uList = cfResult.getUNames();
        List<String> eOList = cfResult.getEONames();
        List<String> vList = cvcfBasis.cfNames;
        List<String> eList = cvcfBasis.eciNames;

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

    private static void emit(Consumer<String> sink, String message) {
        if (sink != null) {
            sink.accept(message);
        }
        LOG.fine(message);
    }
}
