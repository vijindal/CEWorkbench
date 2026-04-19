package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
import org.ce.model.cvm.CvCfBasis;

import java.util.*;
import java.util.logging.Logger;

/**
 * Unified Hamiltonian layer for MCS: holds all cluster embeddings for a supercell,
 * provides factory methods to generate them from cluster geometry, and supplies
 * static utilities for energy and observable computation.
 *
 * <p>Combines what were previously three separate files:
 * EmbeddingData (product), EmbeddingGenerator (producer), LocalEnergyCalc (computations).</p>
 */
public class Embeddings {

    private static final Logger LOG = Logger.getLogger(Embeddings.class.getName());

    private final List<Embedding>   allEmbeddings;
    private final List<Embedding>[] siteToEmbeddings;

    public Embeddings(List<Embedding> allEmbeddings, List<Embedding>[] siteToEmbeddings) {
        this.allEmbeddings    = allEmbeddings;
        this.siteToEmbeddings = siteToEmbeddings;
    }

    // ── Container query methods ───────────────────────────────────────────────

    public List<Embedding>   getAllEmbeddings()    { return allEmbeddings; }
    public List<Embedding>[] getSiteToEmbeddings() { return siteToEmbeddings; }
    public int totalEmbeddingCount()               { return allEmbeddings.size(); }
    public int siteCount()                         { return siteToEmbeddings.length; }

    public double[] computeHmixCoeff(double[] eci, int tc) {
        int N = siteCount();
        int[] counts = new int[tc];
        int[] sizes  = new int[tc];
        for (Embedding e : allEmbeddings) {
            int t = e.getClusterType();
            if (t < tc) { counts[t]++; sizes[t] = e.size(); }
        }
        double[] coeff = new double[tc];
        for (int t = 0; t < tc && t < eci.length; t++) {
            if (sizes[t] > 1) coeff[t] = eci[t] * counts[t] / ((double) sizes[t] * N);
        }
        return coeff;
    }

    public int[] multiSiteEmbedCountsPerType(int tc) {
        int[] counts = new int[tc];
        for (Embedding e : allEmbeddings) {
            int t = e.getClusterType();
            if (t < tc && e.size() > 1) counts[t]++;
        }
        return counts;
    }

    public List<Embedding> getEmbeddingsForTypeAtSite(int clusterType, int siteIndex) {
        List<Embedding> result = new ArrayList<>();
        for (Embedding e : siteToEmbeddings[siteIndex]) {
            if (e.getClusterType() == clusterType) result.add(e);
        }
        return result;
    }

    // ── Factory: generate embeddings from cluster geometry ────────────────────

    /** Generates all cluster embeddings for a supercell from cluster geometry definitions. */
    public static Embeddings generate(
            List<Vector3D>    latticePositions,
            ClusCoordListData clusterData,
            int               L) {

        int N = latticePositions.size();

        Map<Vector3DKey, Integer> posToIndex = new HashMap<>();
        for (int i = 0; i < N; i++)
            posToIndex.put(new Vector3DKey(reduceMod(latticePositions.get(i), L)), i);

        List<ClusterTemplate> templates = buildTemplates(clusterData);

        List<Embedding>   allEmbeddings    = new ArrayList<>();
        @SuppressWarnings("unchecked")
        List<Embedding>[] siteToEmbeddings = new ArrayList[N];
        for (int i = 0; i < N; i++) siteToEmbeddings[i] = new ArrayList<>();

        for (int i = 0; i < N; i++) {
            Vector3D anchor = latticePositions.get(i);
            List<Embedding> raw = new ArrayList<>();

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
                alphas[slot++] = alphaFromSymbol(sites.get(anchorIdx).getSymbol());
                for (int k = 0; k < sites.size(); k++) {
                    if (k != anchorIdx)
                        alphas[slot++] = alphaFromSymbol(sites.get(k).getSymbol());
                }

                raw.add(new Embedding(ttype, omIdx, indices, alphas));
            }

            Set<String>     seen    = new LinkedHashSet<>();
            List<Embedding> deduped = new ArrayList<>();
            for (Embedding e : raw) {
                int[] sorted = e.getSiteIndices().clone();
                Arrays.sort(sorted);
                String key = e.getClusterType() + ":" + Arrays.toString(sorted);
                if (seen.add(key)) deduped.add(e);
            }

            siteToEmbeddings[i] = deduped;
            allEmbeddings.addAll(deduped);
        }

        LOG.fine("Embeddings.generate — " + allEmbeddings.size() + " total embeddings for N=" + N);
        return new Embeddings(allEmbeddings, siteToEmbeddings);
    }

    /** Generates per-CF-column embedding lists for direct CVCF measurement. */
    public static List<List<Embedding>> generateCfEmbeddings(
            List<Embedding>   baseEmbeddings,
            ClusCoordListData clusterData,
            int[][]           cfBasisIndices,
            int[][]           lcf) {

        if (cfBasisIndices == null || baseEmbeddings == null) return null;

        int ncf = cfBasisIndices.length;
        int[] colToOrbitType = new int[ncf];
        int col = 0;
        for (int t = 0; lcf != null && t < lcf.length && col < ncf; t++)
            for (int j = 0; j < lcf[t].length && col < ncf; j++)
                for (int k = 0; k < lcf[t][j] && col < ncf; k++)
                    colToOrbitType[col++] = t;

        Map<Integer, List<Embedding>> typeMap = new HashMap<>();
        for (Embedding e : baseEmbeddings)
            typeMap.computeIfAbsent(e.getClusterType(), k -> new ArrayList<>()).add(e);

        List<List<Embedding>> cfEmbeddings = new ArrayList<>(ncf);
        for (int l = 0; l < ncf; l++) {
            int[] alphas = cfBasisIndices[l];
            List<Embedding> matched = new ArrayList<>();

            if (alphas == null || alphas.length < 2) {
                cfEmbeddings.add(matched);
                continue;
            }

            int t = colToOrbitType[l];
            List<Embedding> typeEmbs = typeMap.get(t);
            if (typeEmbs != null) {
                for (Embedding base : typeEmbs) {
                    matched.add(new Embedding(
                            t, base.getOrbitMemberIndex(),
                            base.getSiteIndices().clone(), alphas.clone()));
                }
            }
            cfEmbeddings.add(matched);
        }

        return cfEmbeddings;
    }

    // ── Energy and observable computation (Hamiltonian layer) ─────────────────

    public static double clusterProduct(Embedding e, LatticeConfig config) {
        double prod  = 1.0;
        int[]  idx    = e.getSiteIndices();
        int[]  alphas = e.getAlphaIndices();
        for (int k = 0; k < idx.length; k++)
            prod *= config.evaluateBasis(alphas[k], config.getOccupation(idx[k]));
        return prod;
    }

    public static double totalEnergy(LatticeConfig config,
                                     Embeddings emb,
                                     double[] eci,
                                     List<List<Cluster>> orbits) {
        double sum = 0.0;
        for (Embedding e : emb.getAllEmbeddings()) {
            int size = e.size();
            if (size > 0)
                sum += eci[e.getClusterType()] * clusterProduct(e, config) / size;
            else
                sum += eci[e.getClusterType()] * clusterProduct(e, config);
        }
        return sum;
    }

    public static double deltaESingleSite(int i,
                                          int newOcc,
                                          LatticeConfig config,
                                          Embeddings emb,
                                          double[] eci,
                                          List<List<Cluster>> orbits) {
        int oldOcc = config.getOccupation(i);
        if (oldOcc == newOcc) return 0.0;

        double dE = 0.0;
        for (Embedding e : emb.getSiteToEmbeddings()[i]) {
            int    t      = e.getClusterType();
            int[]  idx    = e.getSiteIndices();
            int[]  alphas = e.getAlphaIndices();
            double restProduct = 1.0;
            int    alphaI = -1;

            for (int k = 0; k < idx.length; k++) {
                if (idx[k] == i) alphaI = alphas[k];
                else restProduct *= config.evaluateBasis(alphas[k], config.getOccupation(idx[k]));
            }

            if (alphaI < 0) continue;
            double phiOld = config.evaluateBasis(alphaI, oldOcc);
            double phiNew = config.evaluateBasis(alphaI, newOcc);
            int size = e.size();
            double energyCont = eci[t] * (phiNew - phiOld) * restProduct;
            dE += (size > 0) ? (energyCont / size) : energyCont;
        }
        return dE;
    }

    public static double deltaEExchange(int i, int j,
                                        LatticeConfig config,
                                        Embeddings emb,
                                        double[] eci,
                                        List<List<Cluster>> orbits) {
        int occI = config.getOccupation(i);
        int occJ = config.getOccupation(j);
        if (occI == occJ) return 0.0;

        double dEi = deltaESingleSite(i, occJ, config, emb, eci, orbits);

        config.setOccupation(i, occJ);
        double dEj;
        try {
            dEj = deltaESingleSite(j, occI, config, emb, eci, orbits);
        } finally {
            config.setOccupation(i, occI);
        }
        return dEi + dEj;
    }

    /** Builds the basis value lookup table [occ][alpha-1] for direct CF measurement. */
    public static double[][] buildBasisValues(int numComp) {
        double[] sequence = org.ce.model.cluster.ClusterMath.buildBasis(numComp);
        double[][] values = new double[numComp][numComp - 1];
        for (int occ = 0; occ < numComp; occ++)
            for (int alpha = 1; alpha <= numComp - 1; alpha++)
                values[occ][alpha - 1] = Math.pow(sequence[occ], alpha);
        return values;
    }

    /** Measures CVCF cluster variables from configuration by averaging embedding products per CF column. */
    public static double[] measureCVsFromConfig(
            LatticeConfig config,
            List<List<Embedding>> cfEmbeddings,
            double[][] basisMatrix,
            int ncf) {
        double[] v = new double[ncf];
        for (int l = 0; l < ncf; l++) {
            List<Embedding> embs = cfEmbeddings.get(l);
            if (embs == null || embs.isEmpty()) { v[l] = 0.0; continue; }
            double sum = 0.0;
            for (Embedding e : embs) {
                int[] sites  = e.getSiteIndices();
                int[] alphas = e.getAlphaIndices();
                double prod  = 1.0;
                for (int k = 0; k < sites.length; k++)
                    prod *= basisMatrix[config.getOccupation(sites[k])][alphas[k] - 1];
                sum += prod;
            }
            v[l] = sum / embs.size();
        }
        return v;
    }

    /**
     * Transforms orthogonal CFs to CVCF basis: v_cvcf = Tinv * [uOrthNonPoint | uPoint].
     * Returns uOrthNonPoint unchanged when Tinv is null (caller handles logging).
     */
    public static double[] applyTinvTransform(
            double[] uOrthNonPoint, double[] composition, CvCfBasis basis) {
        double[][] Tinv = basis.Tinv;
        int ncf = basis.numNonPointCfs;
        if (Tinv == null) return uOrthNonPoint.clone();

        int K = basis.numComponents;
        double[] basisSeq = org.ce.model.cluster.ClusterMath.buildBasis(K);
        int nPoint = K - 1;
        double[] uPoint = new double[nPoint];
        for (int k = 0; k < nPoint; k++)
            for (int s = 0; s < K; s++)
                uPoint[k] += composition[s] * Math.pow(basisSeq[s], k + 1);

        int tRows = Tinv[0].length;
        double[] uFull = new double[tRows];
        int nonPt = Math.min(uOrthNonPoint.length, tRows);
        System.arraycopy(uOrthNonPoint, 0, uFull, 0, nonPt);
        for (int k = 0; k < nPoint && nonPt + k < tRows; k++)
            uFull[nonPt + k] = uPoint[k];

        double[] vCvcf = new double[ncf];
        for (int i = 0; i < ncf; i++) {
            double sum = 0.0;
            for (int j = 0; j < Tinv[i].length && j < tRows; j++)
                sum += Tinv[i][j] * uFull[j];
            vCvcf[i] = sum;
        }
        return vCvcf;
    }

    // ── Private generation helpers ────────────────────────────────────────────

    /** Parses alpha index from site symbol (e.g. "s1" → 1). */
    static int alphaFromSymbol(String symbol) {
        if (symbol == null || !symbol.startsWith("s"))
            throw new IllegalArgumentException("Site symbol must start with 's', got: " + symbol);
        try {
            return Integer.parseInt(symbol.substring(1));
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Cannot parse alpha from symbol: " + symbol, e);
        }
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

    // ── Inner classes ─────────────────────────────────────────────────────────

    /** A single embedding of an abstract cluster type onto specific lattice sites. */
    public static class Embedding {

        private final int   clusterType;
        private final int   orbitMemberIndex;
        private final int[] siteIndices;
        private final int[] alphaIndices;

        public Embedding(int clusterType, int orbitMemberIndex, int[] siteIndices, int[] alphaIndices) {
            this.clusterType      = clusterType;
            this.orbitMemberIndex = orbitMemberIndex;
            this.siteIndices      = siteIndices;
            this.alphaIndices     = alphaIndices;
        }

        public int   getClusterType()      { return clusterType; }
        public int   getOrbitMemberIndex() { return orbitMemberIndex; }
        public int[] getSiteIndices()      { return siteIndices; }
        public int   size()                { return siteIndices.length; }
        public int[] getAlphaIndices()     { return alphaIndices; }

        @Override
        public String toString() {
            return "Embedding{type=" + clusterType + ", orbit=" + orbitMemberIndex
                 + ", sites=" + Arrays.toString(siteIndices) + "}";
        }
    }

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

    /** Intermediate representation of one orbit member, used only during template building. */
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
}
