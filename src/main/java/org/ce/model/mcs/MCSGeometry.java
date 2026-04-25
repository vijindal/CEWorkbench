package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.debug.MCSDebug;
import org.ce.model.ModelSession;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterIdentificationRequest;
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

        // ── Stage 1a: Cluster & CF Identification ─────────────────────────────
        emit(progressSink, "  [Session] Stage 1a: Cluster identification...");

        String structure = session.systemId.structure();
        String model = session.systemId.model();

        // 0. Resolve Configuration
        ClusterIdentificationRequest config = ClusterIdentificationRequest.fromSystem(
                session.systemId.elements(), structure, model);

        // 1. Load resources
        emit(progressSink, "\n[STAGE 0]: Loading Inputs...");
        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(config.getDisorderedClusterFile());
        disorderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup disorderedSpaceGroup = InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());

        List<Cluster> orderedClusters = InputLoader.parseClusterFile(config.getOrderedClusterFile());
        orderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup orderedSpaceGroup = InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());

        // 2. Stage 1 & 2: Identification Pipeline
        emit(progressSink, "\n[STAGE 1/2]: Running Identification Pipeline...");
        PipelineResult pr = ClusterCFIdentificationPipeline.run(
                disorderedClusters,
                disorderedSpaceGroup.getOperations(),
                orderedClusters,
                orderedSpaceGroup.getOperations(),
                config.getTransformationMatrix(),
                new double[] { 
                    config.getTranslationVector().getX(),
                    config.getTranslationVector().getY(),
                    config.getTranslationVector().getZ() 
                },
                numComp,
                progressSink);

        // 3. Stage 3: C-Matrix foundation
        emit(progressSink, "\n[STAGE 3]: Running C-Matrix Pipeline...");
        CMatrixPipeline.CMatrixData matrixData = CMatrixPipeline.run(
                pr.toClusterIdentificationResult(),
                pr.toCFIdentificationResult(),
                disorderedClusters,
                numComp,
                progressSink);

        // ── Stage 1d: Basis Transformation ────────────────────────────────────
        emit(progressSink, "\n  [Session] Stage 1d: Basis transformation...");
        String parentStructure = ClusterCFIdentificationPipeline.resolveParentStructure(structure);
        CvCfBasis basisRef = CvCfBasis.generate(parentStructure, pr, matrixData, model, progressSink);

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

        if (basis != null && matData != null) {
            this.cfEmbeddings = Embeddings.generateCfEmbeddings(
                    emb.getAllEmbeddings(), clusterData, matData.getCfBasisIndices(), pr.getLcf());
            this.basisMatrix = Embeddings.buildBasisValues(numComp);
        } else {
            this.cfEmbeddings = null;
            this.basisMatrix = null;
        }

        LOG.info(String.format("MCSGeometry built — L=%d, N=%d sites, numComp=%d", L, positions.size(), numComp));

        // ── MCS-DBG: embedding + CF structure validation ──
        if (MCSDebug.ENABLED) {
            MCSDebug.separator("GEOMETRY BUILT — L=" + L + ", N=" + positions.size());
            MCSDebug.log("GEO", "Total embeddings: %d", emb.totalEmbeddingCount());
            MCSDebug.log("GEO", "Orbit types (tc): %d", tc);
            for (int t = 0; t < tc; t++) {
                MCSDebug.log("GEO", "  orbit[%d]: %d members, multiSiteEmbedCount=%d",
                        t, orbitSizes[t], multiSiteEmbedCounts[t]);
            }
            if (cfEmbeddings != null) {
                MCSDebug.log("GEO", "cfEmbeddings: %d CF columns", cfEmbeddings.size());
                for (int l = 0; l < cfEmbeddings.size(); l++) {
                    MCSDebug.log("GEO", "  cfEmb[%d]: %d embeddings", l, cfEmbeddings.get(l).size());
                }
            } else {
                MCSDebug.log("GEO", "cfEmbeddings: NULL (no CF measurement possible)");
            }
            if (basisMatrix != null) {
                MCSDebug.matrix("GEO", "basisMatrix (Embeddings.buildBasisValues)", basisMatrix);
            }
        }
    }

    /**
     * Builds a new geometry for the given session and supercell size.
     */
    public static MCSGeometry build(ModelSession session, int L, Consumer<String> progressSink) {
        return new MCSGeometry(session, L, progressSink);
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
