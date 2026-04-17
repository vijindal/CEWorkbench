package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;

import org.ce.model.cluster.Cluster;
import java.util.*;
import java.util.logging.Logger;

/** Generates all embeddings of abstract cluster types onto supercell lattice sites. */
public class EmbeddingGenerator {

    private static final Logger LOG = Logger.getLogger(EmbeddingGenerator.class.getName());

    private EmbeddingGenerator() {}

    public static EmbeddingData generateEmbeddings(
            List<Vector3D>      latticePositions,
            ClusCoordListData   clusterData,
            int                 L) {

        int N = latticePositions.size();

        Map<Vector3DKey, Integer> posToIndex = new HashMap<>();
        for (int i = 0; i < N; i++) {
            posToIndex.put(new Vector3DKey(reduceMod(latticePositions.get(i), L)), i);
        }

        List<ClusterTemplate> templates = buildTemplates(clusterData);

        List<EmbeddingData.Embedding>   allEmbeddings    = new ArrayList<>();
        @SuppressWarnings("unchecked")
        List<EmbeddingData.Embedding>[] siteToEmbeddings = new ArrayList[N];
        for (int i = 0; i < N; i++) siteToEmbeddings[i] = new ArrayList<>();

        for (int i = 0; i < N; i++) {
            Vector3D anchor = latticePositions.get(i);
            List<EmbeddingData.Embedding> raw = new ArrayList<>();

            for (ClusterTemplate template : templates) {
                Vector3D[] rel     = template.getRelativeVectors();
                int[]      indices = new int[rel.length];
                boolean    valid   = true;

                for (int k = 0; k < rel.length; k++) {
                    Vector3D target = reduceMod(anchor.add(rel[k]), L);
                    Integer  j      = posToIndex.get(new Vector3DKey(target));
                    if (j == null) { valid = false; break; }
                    indices[k] = j;
                }

                if (!valid) continue;

                int ttype = template.getClusterType();
                int omIdx = template.getOrbitMemberIndex();
                List<Site> sites     = clusterData.getOrbitList().get(ttype).get(omIdx).getAllSites();
                int        anchorIdx = template.getAnchorIndex();
                int[]      alphas    = new int[sites.size()];

                int slot = 0;
                alphas[slot++] = SiteOperatorBasis.alphaFromSymbol(sites.get(anchorIdx).getSymbol());
                for (int k = 0; k < sites.size(); k++) {
                    if (k != anchorIdx)
                        alphas[slot++] = SiteOperatorBasis.alphaFromSymbol(sites.get(k).getSymbol());
                }

                raw.add(new EmbeddingData.Embedding(ttype, omIdx, indices, alphas));
            }

            Set<String>     seen    = new LinkedHashSet<>();
            List<EmbeddingData.Embedding> deduped = new ArrayList<>();
            for (EmbeddingData.Embedding e : raw) {
                int[] sorted = e.getSiteIndices().clone();
                Arrays.sort(sorted);
                String key = e.getClusterType() + ":" + Arrays.toString(sorted);
                if (seen.add(key)) deduped.add(e);
            }

            siteToEmbeddings[i] = deduped;
            allEmbeddings.addAll(deduped);
        }

        LOG.fine("EmbeddingGenerator — " + allEmbeddings.size() + " total embeddings for N=" + N);
        return new EmbeddingData(allEmbeddings, siteToEmbeddings);
    }

    private static List<ClusterTemplate> buildTemplates(ClusCoordListData clusterData) {
        List<ClusterTemplate> templates = new ArrayList<>();
        List<List<Cluster>>   orbitList = clusterData.getOrbitList();

        for (int t = 0; t < orbitList.size(); t++) {
            List<Cluster> orbit = orbitList.get(t);
            int clusterSize = orbit.isEmpty() ? 0 : orbit.get(0).getAllSites().size();
            if (clusterSize < 2) continue;

            for (int o = 0; o < orbit.size(); o++) {
                List<Site> sites = orbit.get(o).getAllSites();
                if (sites.isEmpty()) {
                    templates.add(new ClusterTemplate(t, o, new Vector3D[0], 0));
                    continue;
                }
                int n = sites.size();
                for (int anchor = 0; anchor < n; anchor++) {
                    Vector3D anchorPos = toVector3D(sites.get(anchor).getPosition());
                    Vector3D[] rel     = new Vector3D[n];
                    rel[0] = new Vector3D(0, 0, 0);
                    int slot = 1;
                    for (int k = 0; k < n; k++) {
                        if (k == anchor) continue;
                        rel[slot++] = toVector3D(sites.get(k).getPosition()).subtract(anchorPos);
                    }
                    templates.add(new ClusterTemplate(t, o, rel, anchor));
                }
            }
        }
        return templates;
    }

    private static Vector3D toVector3D(Position p) {
        return new Vector3D(p.getX(), p.getY(), p.getZ());
    }

    static Vector3D reduceMod(Vector3D v, double L) {
        return new Vector3D(
                v.getX() - L * Math.floor(v.getX() / L),
                v.getY() - L * Math.floor(v.getY() / L),
                v.getZ() - L * Math.floor(v.getZ() / L));
    }

    // -------------------------------------------------------------------------
    // Private inner classes (only used by EmbeddingGenerator)
    // -------------------------------------------------------------------------

    /** HashMap key for Vector3D with tolerance-safe hashing. */
    private static final class Vector3DKey {

        private static final double ROUND = 1e6;
        private final long     hx, hy, hz;
        private final Vector3D v;

        Vector3DKey(Vector3D v) {
            this.v  = v;
            this.hx = Math.round(v.getX() * ROUND);
            this.hy = Math.round(v.getY() * ROUND);
            this.hz = Math.round(v.getZ() * ROUND);
        }

        @Override public int hashCode() {
            return Long.hashCode(hx) * 961 + Long.hashCode(hy) * 31 + Long.hashCode(hz);
        }

        @Override public boolean equals(Object obj) {
            if (!(obj instanceof Vector3DKey)) return false;
            return v.equals(((Vector3DKey) obj).v);
        }
    }

    /** Position-relative representation of one orbit member, used only during template building. */
    private static final class ClusterTemplate {

        private final int        clusterType;
        private final int        orbitMemberIndex;
        private final int        anchorIndex;
        private final Vector3D[] relativeVectors;

        ClusterTemplate(int clusterType, int orbitMemberIndex,
                Vector3D[] relativeVectors, int anchorIndex) {
            this.clusterType      = clusterType;
            this.orbitMemberIndex = orbitMemberIndex;
            this.relativeVectors  = relativeVectors;
            this.anchorIndex      = anchorIndex;
        }

        int        getClusterType()      { return clusterType; }
        int        getOrbitMemberIndex() { return orbitMemberIndex; }
        Vector3D[] getRelativeVectors()  { return relativeVectors; }
        int        getAnchorIndex()      { return anchorIndex; }
    }

    /**
     * Generates a list of embeddings for each orthogonal CF column index.
     *
     * <p>Uses {@code lcf[t][j]} to determine which orbit type {@code t} each CF column
     * {@code l} belongs to. This is necessary because multiple orbit types can share the
     * same cluster size (e.g., two distinct pair types in BCC_A2 T). Matching by size
     * alone would assign both pair CF columns to the same orbit type.</p>
     *
     * @param baseEmbeddings  all embeddings generated by {@link #generateEmbeddings}
     * @param clusterData     orbit list for the disordered structure
     * @param cfBasisIndices  alpha (σ-power) indices per CF column, as from {@code CMatrixPipeline.CMatrixData}
     * @param lcf             lcf[t][j] = number of decorations for orbit type t, group j
     */
    public static List<List<EmbeddingData.Embedding>> generateCfEmbeddings(
            List<EmbeddingData.Embedding> baseEmbeddings,
            ClusCoordListData   clusterData,
            int[][] cfBasisIndices,
            int[][] lcf) {

        if (cfBasisIndices == null || baseEmbeddings == null) return null;

        // Build orbit-type index for each CF column by iterating lcf in the same order
        // as buildCfColumnMap: t outer, j middle, k inner.
        int ncf = cfBasisIndices.length;
        int[] colToOrbitType = new int[ncf];
        int col = 0;
        for (int t = 0; lcf != null && t < lcf.length && col < ncf; t++) {
            for (int j = 0; j < lcf[t].length && col < ncf; j++) {
                for (int k = 0; k < lcf[t][j] && col < ncf; k++) {
                    colToOrbitType[col++] = t;
                }
            }
        }

        Map<Integer, List<EmbeddingData.Embedding>> typeMap = new HashMap<>();
        for (EmbeddingData.Embedding e : baseEmbeddings) {
            typeMap.computeIfAbsent(e.getClusterType(), k -> new ArrayList<>()).add(e);
        }

        List<List<EmbeddingData.Embedding>> cfEmbeddings = new ArrayList<>(ncf);
        for (int l = 0; l < ncf; l++) {
            int[] alphas = cfBasisIndices[l];
            List<EmbeddingData.Embedding> matched = new ArrayList<>();

            if (alphas == null || alphas.length < 2) {
                cfEmbeddings.add(matched);
                continue;
            }

            int t = colToOrbitType[l];
            List<EmbeddingData.Embedding> typeEmbs = typeMap.get(t);
            if (typeEmbs != null) {
                for (EmbeddingData.Embedding base : typeEmbs) {
                    matched.add(new EmbeddingData.Embedding(
                            t, base.getOrbitMemberIndex(),
                            base.getSiteIndices().clone(), alphas.clone()));
                }
            }
            cfEmbeddings.add(matched);
        }

        return cfEmbeddings;
    }
}
