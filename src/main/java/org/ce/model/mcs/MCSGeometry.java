package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.ModelSession;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult;
import org.ce.model.storage.InputLoader;
import org.ce.model.cluster.SpaceGroup;
import org.ce.model.cvm.CvCfBasis;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Encapsulates the temperature-independent geometry and orbit structure of an MCS model.
 * 
 * <p>Building this object involves expensive lattice construction and cluster embedding
 * generation. By caching {@code MCSGeometry}, we can reuse the spatial structure across
 * multiple temperature points in a thermodynamic scan, only re-evaluating the ECIs.</p>
 */
public final class MCSGeometry {

    private static final Logger LOG = Logger.getLogger(MCSGeometry.class.getName());

    final ClusCoordListData       clusterData;
    final List<Vector3D>          positions;
    final Embeddings              emb;
    final List<List<Cluster>>     orbits;
    final int[]                   orbitSizes;
    final int[]                   multiSiteEmbedCounts;
    final CvCfBasis               basis;
    final List<List<Embeddings.Embedding>> cfEmbeddings;
    final double[][]              basisMatrix;
    final int                     L;
    final int                     numComp;

    private MCSGeometry(ModelSession session, int L, Consumer<String> progressSink) {
        this.L = L;
        this.numComp = session.numComponents();

        PipelineResult pr = session.clusterData;
        CvCfBasis basisRef = session.cvcfBasis;

        if (pr == null || basisRef == null) {
            emit(progressSink, "  [MCS Geometry] Session data incomplete, running local identification...");
            pr = runLocalIdentification(session, progressSink);
            basisRef = runLocalBasisResolution(session, pr, progressSink);
        }

        this.clusterData = pr.getDisorderedClusterResult().getDisClusterData();
        this.basis = basisRef;

        int tc = clusterData.getTc();

        // Build lattice geometry and embeddings
        this.positions = buildBCCPositions(L);
        this.emb = Embeddings.generate(positions, clusterData, L);

        this.orbits = clusterData.getOrbitList();
        this.orbitSizes = new int[tc];
        for (int t = 0; t < tc; t++) {
            orbitSizes[t] = orbits.get(t).size();
        }

        this.multiSiteEmbedCounts = emb.multiSiteEmbedCountsPerType(tc);

        CMatrixPipeline.CMatrixData matData = pr.getMatrixData();
        if (matData == null && basis != null && basis.cvcfCMatrixData != null) {
            matData = basis.cvcfCMatrixData;
        }

        if (matData == null && basis != null) {
            // Supplemental calculation of matrix data if it was missing 
            matData = CMatrixPipeline.run(
                    pr.toClusterIdentificationResult(),
                    pr.toCFIdentificationResult(),
                    pr.getDisClusData().getClusCoordList(),
                    numComp,
                    progressSink);
        }

        if (basis != null && matData != null) {
            this.cfEmbeddings = Embeddings.generateCfEmbeddings(
                    emb.getAllEmbeddings(), clusterData, matData.getCfBasisIndices(), pr.getLcf());
            this.basisMatrix = Embeddings.buildBasisValues(numComp);
        } else {
            this.cfEmbeddings = null;
            this.basisMatrix = null;
        }

        LOG.info(String.format("MCSGeometry built — L=%d, N=%d sites, numComp=%d", L, positions.size(), numComp));
    }

    /**
     * Builds a new geometry for the given session and supercell size.
     */
    public static MCSGeometry build(ModelSession session, int L, Consumer<String> progressSink) {
        return new MCSGeometry(session, L, progressSink);
    }

    private PipelineResult runLocalIdentification(ModelSession session, Consumer<String> progressSink) {
        String structure = session.systemId.structure();
        String model = session.systemId.model();
        int numComp = session.numComponents();

        String parentStructure = resolveParentStructure(structure);
        String disorderedFile = "clus/" + parentStructure + "-" + model + ".txt";
        String orderedFile = "clus/" + structure + "-" + model + ".txt";
        String disorderedSGName = parentStructure + "-SG";
        String orderedSGName = structure + "-SG";

        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(disorderedFile);
        SpaceGroup disorderedSG = InputLoader.parseSpaceGroup(disorderedSGName);
        List<Cluster> orderedClusters = InputLoader.parseClusterFile(orderedFile);
        SpaceGroup orderedSG = InputLoader.parseSpaceGroup(orderedSGName);

        return ClusterCFIdentificationPipeline.run(
                disorderedClusters, disorderedSG.getOperations(),
                orderedClusters, orderedSG.getOperations(),
                orderedSG.getRotateMat(), orderedSG.getTranslateMat(),
                numComp, progressSink);
    }

    private CvCfBasis runLocalBasisResolution(ModelSession session, PipelineResult pr, Consumer<String> progressSink) {
        CMatrixPipeline.CMatrixData cmatOrth = CMatrixPipeline.run(
                pr.toClusterIdentificationResult(),
                pr.toCFIdentificationResult(),
                pr.getDisClusData().getClusCoordList(),
                session.numComponents(),
                progressSink);

        return CvCfBasis.generate(session.systemId.structure(), pr, cmatOrth, session.systemId.model(), progressSink);
    }

    private String resolveParentStructure(String structure) {
        if (structure == null) return null;
        String base = structure.replace("_CVCF", "");
        if (base.equals("BCC_B2")) return "BCC_A2";
        if (base.equals("FCC_L12")) return "FCC_A1";
        return base;
    }

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }

    public int nSites() {
        return positions.size();
    }

    private static List<Vector3D> buildBCCPositions(int L) {
        if (L < 1) throw new IllegalArgumentException("L must be >= 1");
        List<Vector3D> pos = new ArrayList<>(2 * L * L * L);
        for (int ix = 0; ix < L; ix++)
            for (int iy = 0; iy < L; iy++)
                for (int iz = 0; iz < L; iz++) {
                    pos.add(new Vector3D(ix,       iy,       iz      ));
                    pos.add(new Vector3D(ix + 0.5, iy + 0.5, iz + 0.5));
                }
        return pos;
    }
}
