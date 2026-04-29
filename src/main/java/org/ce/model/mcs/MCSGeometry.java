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
 * Encapsulates the static, temperature-independent geometry and orbit structure of a 
 * Monte Carlo Simulation (MCS) model.
 * 
 * <p>
 * This class serves as the structural "blueprint" for the simulation. It defines the physical 
 * arrangement of sites, identifies all unique cluster orbits allowed by the space group 
 * symmetry, and pre-computes high-performance indices for energy evaluations.
 * </p>
 * 
 * <h3>Key Responsibilities:</h3>
 * <ul>
 *   <li><b>Lattice Construction:</b> Generates the 3D atomic positions for a supercell 
 *       of dimension {@code L} based on the crystal structure (e.g., BCC_A2).</li>
 *   <li><b>Orbit Resolution:</b> Identifies and deduplicates cluster instances into 
 *       symmetrically equivalent orbits.</li>
 *   <li><b>Embedding Generation:</b> Locates all instances (embeddings) of every cluster 
 *       type within the periodic supercell.</li>
 *   <li><b>Performance Indexing:</b> Builds the <b>CSR (Compressed Sparse Row)</b> 
 *       site-to-embedding index and <b>Flattened Embedding Data</b> structures. These indices 
 *       enable $O(1)$ local energy updates during the simulation's trial moves.</li>
 * </ul>
 * 
 * <h3>Computational Complexity:</h3>
 * <p>
 * Instantiating this class is an <b>expensive $O(N)$ operation</b> because it involves 
 * recursive cluster searching and symmetry checks across the entire supercell. However, 
 * because the spatial structure is independent of temperature, this object is designed 
 * to be <b>cached and reused</b> across multiple thermodynamic points in a temperature scan.
 * </p>
 * 
 * <h3>Thread Safety:</h3>
 * <p>
 * This class is <b>effectively immutable</b> after construction. All fields are final and 
 * represent fixed geometric properties. It can be safely shared across multiple 
 * {@link AlloyMC} engine instances.
 * </p>
 * 
 * @see org.ce.model.mcs.AlloyMC
 * @see org.ce.model.mcs.Embeddings
 * @see org.ce.model.mcs.LatticeConfig
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
    final Embeddings.CsrSiteToCfIndex     siteToCfIndex;
    final Embeddings.FlatEmbData          flatEmbData;
    final double[][]                      basisMatrix;
    final double[]                        flatBasisMatrix;
    final int                             L;
    final int                             numComp;
    final int                             ncf;

    private MCSGeometry(ModelSession session, int L, Consumer<String> progressSink) {
        this.L = L;
        this.numComp = session.numComponents();

        // ── Step 1a: Cluster & CF Identification ──────────────────────────────
        emit(progressSink, "  [Session] Step 1a: Cluster identification...");

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
        
        this.ncf = pr.getNcf();

        // 3. Stage 3: C-Matrix foundation
        emit(progressSink, "\n[STAGE 3]: Running C-Matrix Pipeline...");
        CMatrixPipeline.CMatrixData matrixData = CMatrixPipeline.run(
                pr.toClusterIdentificationResult(),
                pr.toCFIdentificationResult(),
                disorderedClusters,
                numComp,
                progressSink);

        // ── Step 1b: Basis Transformation ─────────────────────────────────────
        emit(progressSink, "\n  [Session] Step 1b: Basis transformation...");
        String parentStructure = ClusterCFIdentificationPipeline.resolveParentStructure(structure);
        CvCfBasis basisRef = CvCfBasis.generate(parentStructure, pr, matrixData, model, progressSink);

        this.clusterData = pr.getDisorderedClusterResult().getDisClusterData();
        this.basis = basisRef;

        int tc = clusterData.getTc();

        // ── Step 1c: Lattice & Embedding Generation ────────────────────────────
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
                    this.positions, clusterData, L, matData.getCfBasisIndices(), pr.getLcf());
            this.basisMatrix = Embeddings.buildBasisValues(numComp);
            
            // Build high-performance indices for ΔE
            this.siteToCfIndex = Embeddings.buildSiteToCfIndex(cfEmbeddings, nSites());
            this.flatEmbData   = Embeddings.FlatEmbData.build(cfEmbeddings);
            
            // Flatten basis matrix for high-speed indexing: [occ * numComp + alpha]
            this.flatBasisMatrix = new double[numComp * numComp];
            for (int s = 0; s < numComp; s++) {
                // Alpha 0 is the point function/constant (usually 1.0)
                flatBasisMatrix[s * numComp + 0] = 1.0; 
                for (int a = 1; a < numComp; a++) {
                    flatBasisMatrix[s * numComp + a] = basisMatrix[s][a - 1];
                }
            }
        } else {
            this.cfEmbeddings  = null;
            this.basisMatrix   = null;
            this.flatBasisMatrix = null;
            this.siteToCfIndex = null;
            this.flatEmbData   = null;
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

    public List<List<Embeddings.Embedding>> cfEmbeddings() {
        return cfEmbeddings;
    }

    public Embeddings.CsrSiteToCfIndex getSiteToCfIndex() { return siteToCfIndex; }
    public Embeddings.FlatEmbData getFlatEmbData() { return flatEmbData; }
    public double[][] getBasisMatrix() { return basisMatrix; }
    public double[] getFlatBasisMatrix() { return flatBasisMatrix; }
    public CvCfBasis getBasis() { return basis; }
    public int getNcf() { return ncf; }

    /**
     * Transforms an orthogonal CF state vector [uOrthNonPoint | uPoint | Empty]
     * to the CVCF basis [vCvcfNonPoint | moleFractions].
     * 
     * @param correlationFunctions The full orthogonal state vector.
     * @param composition          The mole fractions to append.
     * @return The full-length CVCF basis state vector.
     */
    public double[] getCvcfCorrelationFunctions(double[] correlationFunctions, double[] composition) {
        if (basis == null || basis.Tinv == null) {
            return java.util.Arrays.copyOf(correlationFunctions, ncf + numComp);
        }
        
        double[][] Tinv = basis.Tinv;
        double[] vFull = new double[ncf + numComp];
        
        // 1. Transform non-point terms [0..ncf-1]
        for (int i = 0; i < ncf; i++) {
            double sum = 0.0;
            for (int j = 0; j < Tinv[i].length && j < correlationFunctions.length; j++) {
                sum += Tinv[i][j] * correlationFunctions[j];
            }
            vFull[i] = sum;
        }
        
        // 2. Append mole fractions (compositions) [ncf..ncf+numComp-1]
        if (composition != null) {
            System.arraycopy(composition, 0, vFull, ncf, numComp);
        }
        
        return vFull;
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
