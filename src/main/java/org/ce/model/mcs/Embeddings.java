package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.debug.MCSDebug;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
import org.ce.model.cvm.CvCfBasis;

import java.util.*;
import java.util.logging.Logger;

/**
 * Unified Hamiltonian layer for MCS.
 *
 * <p>Holds all cluster embeddings for a supercell and provides static utilities
 * for every energy/observable computation step of the MC pipeline:
 *
 * <pre>
 *  Section A — Inner data types       Embedding, CsrSiteToCfIndex, FlatEmbData,
 *                                     DeltaScratch, Vector3DKey, ClusterTemplate
 *
 *  Section B — Instance container     allEmbeddings, siteToEmbeddings,
 *                                     query methods, multiSiteEmbedCountsPerType
 *
 *  Section C — Step 1: Generation     generate(), generateCfEmbeddings()
 *                                     + helpers: buildTemplates, alphaFromSymbol,
 *                                       toVector3D, reduceMod, buildBasisValues
 *
 *  Section D — Step 2: Initial E(σ)   totalEnergyCvcf(), clusterProduct()
 *
 *  Section E — Step 5: ΔE (hot path)  deltaEExchangeCvcf(FlatEmbData)   ← primary
 *                                     deltaEExchangeCvcf(List-based)     ← fallback
 *
 *  Section F — Step 7: Observables    measureCVsFromConfig(), applyTinvTransform()
 *
 *  Section G — Index builders         buildSiteToCfIndex(), FlatEmbData.build(),
 *                                     maxEmbPerCfColumn(), totalCfEmbeddingCount()
 *
 *  Section H — Debug counters         resetDebugCounters() + static counter fields
 * </pre>
 */
public class Embeddings {

    private static final Logger LOG = Logger.getLogger(Embeddings.class.getName());

    // =========================================================================
    // Section H — Debug counters (reset at start of each run)
    // =========================================================================

    private static int     dbgClusterProdCalls    = 0;
    private static int     dbgMeasureCVCalls      = 0;
    private static int     dbgApplyTinvCalls      = 0;
    private static boolean dbgTotalEnergyTraced   = false;

    /** Reset all debug sample counters (call at start of a new run). */
    public static void resetDebugCounters() {
        dbgClusterProdCalls  = 0;
        dbgMeasureCVCalls    = 0;
        dbgApplyTinvCalls    = 0;
        dbgTotalEnergyTraced = false;
    }

    // =========================================================================
    // Section B — Instance container (holds pre-generated embedding lists)
    // =========================================================================

    private final List<Embedding>   allEmbeddings;
    private final List<Embedding>[] siteToEmbeddings;

    public Embeddings(List<Embedding> allEmbeddings, List<Embedding>[] siteToEmbeddings) {
        this.allEmbeddings    = allEmbeddings;
        this.siteToEmbeddings = siteToEmbeddings;
    }

    // ── Container queries ─────────────────────────────────────────────────────

    public List<Embedding>   getAllEmbeddings()    { return allEmbeddings; }
    public List<Embedding>[] getSiteToEmbeddings() { return siteToEmbeddings; }
    public int totalEmbeddingCount()               { return allEmbeddings.size(); }
    public int siteCount()                         { return siteToEmbeddings.length; }

    public int[] multiSiteEmbedCountsPerType(int tc) {
        int[] counts = new int[tc];
        for (Embedding e : allEmbeddings) {
            int t = e.getClusterType();
            if (t < tc && e.size() > 1) counts[t]++;
        }
        return counts;
    }

    // =========================================================================
    // Section C — Step 1: Generation  (called once per supercell build)
    // =========================================================================

    /**
     * Generates deduplicated cluster embeddings for the energy path.
     *
     * <p>For each physical site as anchor, tries all cluster templates. Deduplicates
     * by sorted site indices so each unordered cluster instance appears once.
     * Alpha assignments come from orbit member site symbols (anchor-first order).
     */
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

                int ttype    = template.getClusterType();
                int omIdx    = template.getOrbitMemberIndex();
                int anchorIdx = template.getAnchorIndex();
                List<Site> sites = clusterData.getOrbitList().get(ttype).get(omIdx).getAllSites();
                int[]      alphas = new int[sites.size()];

                // Anchor-first alpha assignment from orbit member site symbols
                int slot = 0;
                alphas[slot++] = alphaFromSymbol(sites.get(anchorIdx).getSymbol());
                for (int k = 0; k < sites.size(); k++) {
                    if (k != anchorIdx)
                        alphas[slot++] = alphaFromSymbol(sites.get(k).getSymbol());
                }

                raw.add(new Embedding(ttype, omIdx, indices, alphas, anchorIdx));
            }

            // Deduplicate: keep one embedding per unordered cluster instance
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

    /**
     * Generates per-CF-column directed embedding lists for CVCF measurement and ΔE.
     *
     * <p>Unlike {@link #generate()}, this method produces ALL directed embeddings
     * (one per anchor choice per cluster instance) without deduplication. This is
     * required for K≥3 asymmetric CF columns: e.g., the ternary 1NN pair CF u_{12}
     * = (1/N_emb) × Σ φ_1(σ_i)·φ_2(σ_j) must sum over ordered pairs (i→j) only,
     * and including the reversed embedding (j→i with swapped alphas) correctly
     * double-counts in a way that preserves the average.
     *
     * <p>Alpha assignment per CF column ({@code cfBasisIndices[l]}) is in canonical
     * {@code getAllSites()} order. Each directed embedding maps these alphas to the
     * anchor-first physical site ordering by reordering {@code cfBasisIndices[l]}
     * according to the template's {@code anchorIndex}.
     */
    public static List<List<Embedding>> generateCfEmbeddings(
            List<Vector3D>    latticePositions,
            ClusCoordListData clusterData,
            int               L,
            int[][]           cfBasisIndices,
            int[][]           lcf) {

        if (cfBasisIndices == null || latticePositions == null) return null;

        int N = latticePositions.size();

        // Build position → lattice index lookup
        Map<Vector3DKey, Integer> posToIndex = new HashMap<>();
        for (int i = 0; i < N; i++)
            posToIndex.put(new Vector3DKey(reduceMod(latticePositions.get(i), L)), i);

        // Map each CF column to its orbit type
        int ncf = cfBasisIndices.length;
        int[] colToOrbitType = new int[ncf];
        int col = 0;
        for (int t = 0; lcf != null && t < lcf.length && col < ncf; t++)
            for (int j = 0; j < lcf[t].length && col < ncf; j++)
                for (int k = 0; k < lcf[t][j] && col < ncf; k++)
                    colToOrbitType[col++] = t;

        // Build templates (same as generate(), but we use ALL templates without dedup)
        List<ClusterTemplate> templates = buildTemplates(clusterData);

        // Group templates by orbit type for efficient lookup
        Map<Integer, List<ClusterTemplate>> templatesByType = new HashMap<>();
        for (ClusterTemplate tmpl : templates)
            templatesByType.computeIfAbsent(tmpl.getClusterType(), k -> new ArrayList<>()).add(tmpl);

        // For each CF column, generate all directed embeddings
        List<List<Embedding>> cfEmbeddings = new ArrayList<>(ncf);
        for (int l = 0; l < ncf; l++) {
            int[] cfAlphas = cfBasisIndices[l];
            List<Embedding> directed = new ArrayList<>();

            if (cfAlphas == null || cfAlphas.length < 2) {
                cfEmbeddings.add(directed);
                continue;
            }

            int t = colToOrbitType[l];
            int n = cfAlphas.length; // sites per cluster
            List<ClusterTemplate> typeTemplates = templatesByType.getOrDefault(t, Collections.emptyList());

            for (ClusterTemplate template : typeTemplates) {
                int omIdx    = template.getOrbitMemberIndex();
                int anchorIdx = template.getAnchorIndex();
                Vector3D[] rel = template.getRelativeVectors();

                // Try each physical site as the anchor
                for (int i = 0; i < N; i++) {
                    Vector3D anchor = latticePositions.get(i);
                    int[] indices = new int[n];
                    boolean valid = true;

                    for (int k = 0; k < n; k++) {
                        Vector3D target = reduceMod(anchor.add(rel[k]), L);
                        Integer  j      = posToIndex.get(new Vector3DKey(target));
                        if (j == null) { valid = false; break; }
                        indices[k] = j;
                    }
                    if (!valid) continue;

                    // Map cfBasisIndices[l] (canonical getAllSites() order) to anchor-first order:
                    //   slot 0 → canonical site anchorIdx
                    //   slots 1..n-1 → remaining canonical sites in order, skipping anchorIdx
                    int[] alphas = new int[n];
                    int slot = 0;
                    alphas[slot++] = cfAlphas[anchorIdx];
                    for (int k = 0; k < n; k++) {
                        if (k != anchorIdx)
                            alphas[slot++] = cfAlphas[k];
                    }

                    directed.add(new Embedding(t, omIdx, indices, alphas, anchorIdx));
                }
            }

            cfEmbeddings.add(directed);
        }

        if (MCSDebug.ENABLED) {
            MCSDebug.separator("generateCfEmbeddings — directed embedding counts");
            for (int l = 0; l < Math.min(ncf, 10); l++) {
                MCSDebug.log("CF-GEN", "cf[%d]: orbitType=%d, nDirected=%d, cfAlphas=%s",
                        l, colToOrbitType[l], cfEmbeddings.get(l).size(),
                        Arrays.toString(cfBasisIndices[l]));
            }
        }

        return cfEmbeddings;
    }

    // =========================================================================
    // Section D — Step 2: Initial energy E(σ)
    // Called once before equilibration: totalEnergyCvcf() + clusterProduct()
    // =========================================================================

    public static double clusterProduct(Embedding e, LatticeConfig config) {
        double prod  = 1.0;
        int[]  idx    = e.getSiteIndices();
        int[]  alphas = e.getAlphaIndices();
        for (int k = 0; k < idx.length; k++)
            prod *= config.evaluateBasis(alphas[k], config.getOccupation(idx[k]));

        // ── MCS-DBG: sampled clusterProduct trace (first 3 calls) ──
        if (MCSDebug.ENABLED && dbgClusterProdCalls < 3) {
            dbgClusterProdCalls++;
            StringBuilder sb = new StringBuilder();
            for (int k = 0; k < idx.length; k++) {
                int occ = config.getOccupation(idx[k]);
                double bv = config.evaluateBasis(alphas[k], occ);
                sb.append(String.format(" site[%d]=(α=%d,σ=%d,φ=%.4f)", idx[k], alphas[k], occ, bv));
            }
            MCSDebug.log("CPROD", "type=%d prod=%.8f ←%s", e.getClusterType(), prod, sb);
        }

        return prod;
    }

    // =========================================================================
    // Section F — Step 7: Observable measurement
    // Called once per averaging sweep: measureCVsFromConfig() + applyTinvTransform()
    // Also: buildBasisValues() used during initialization
    // =========================================================================

    /** Builds the basis value lookup table [occ][alpha-1] for direct CF measurement. */
    public static double[][] buildBasisValues(int numComp) {
        double[] sequence = org.ce.model.cluster.ClusterMath.buildBasis(numComp);
        double[][] values = new double[numComp][numComp - 1];
        for (int occ = 0; occ < numComp; occ++)
            for (int alpha = 1; alpha <= numComp - 1; alpha++)
                values[occ][alpha - 1] = Math.pow(sequence[occ], alpha);
        return values;
    }

    public static double[] measureCVsFromConfig(
            LatticeConfig config,
            List<List<Embedding>> cfEmbeddings,
            double[] flatBasisMatrix,
            int ncf,
            int numComp) {
        boolean trace = MCSDebug.ENABLED && dbgMeasureCVCalls == 0;
        if (trace) dbgMeasureCVCalls++;

        int[] occ = config.getRawOcc();
        double[] v = new double[ncf];
        for (int l = 0; l < ncf; l++) {
            List<Embedding> embs = cfEmbeddings.get(l);
            if (embs == null || embs.isEmpty()) { v[l] = 0.0; continue; }
            double sum = 0.0;
            for (int ei = 0; ei < embs.size(); ei++) {
                Embedding e = embs.get(ei);
                int[] sites  = e.getSiteIndices();
                int[] alphas = e.getAlphaIndices();
                double prod  = 1.0;
                for (int k = 0; k < sites.length; k++)
                    prod *= flatBasisMatrix[occ[sites[k]] * numComp + alphas[k]];
                sum += prod;

                if (trace && l == 0 && ei < 3) {
                    MCSDebug.log("CF-MEAS", "  cf[0] emb[%d]: prod=%.8f sites=%s alphas=%s",
                            ei, prod, Arrays.toString(sites), Arrays.toString(alphas));
                }
            }
            v[l] = sum / embs.size();

            if (trace) {
                MCSDebug.log("CF-MEAS", "cf[%d]: nEmb=%d, sum=%.8f, avg=%.8f",
                        l, embs.size(), sum, v[l]);
            }
        }

        if (trace) {
            MCSDebug.vector("CF-MEAS", "measureCVsFromConfig result (uOrth)", v);
        }
        return v;
    }

    /**
     * Transforms orthogonal CFs to CVCF basis: v_cvcf = Tinv * [uOrthNonPoint | uPoint | 1].
     * Returns uOrthNonPoint unchanged when Tinv is null.
     *
     * <p>uFull layout: [uOrth(0..ncf-1) | uPoint(φ₁..φ_{K-1}) | 1.0]
     * Point CFs in natural order φ_k = Σ_s x_s · basis[s]^k (k=1..K-1),
     * matching CvCfBasis.computeRandomStateVectors.
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
        if (tRows > nonPt + nPoint)
            uFull[tRows - 1] = 1.0;

        double[] vCvcf = new double[ncf];
        for (int i = 0; i < ncf; i++) {
            double sum = 0.0;
            for (int j = 0; j < Tinv[i].length && j < tRows; j++)
                sum += Tinv[i][j] * uFull[j];
            vCvcf[i] = sum;
        }

        if (MCSDebug.ENABLED && dbgApplyTinvCalls == 0) {
            dbgApplyTinvCalls++;
            MCSDebug.separator("TINV TRANSFORM (applyTinvTransform)");
            MCSDebug.vector("TINV", "uOrthNonPoint (input)", uOrthNonPoint);
            MCSDebug.vector("TINV", "composition", composition);
            MCSDebug.vector("TINV", "uPoint (derived)", uPoint);
            MCSDebug.vector("TINV", "uFull (assembled)", uFull);
            MCSDebug.log("TINV", "Tinv dims: %dx%d, ncf=%d", Tinv.length, tRows, ncf);
            MCSDebug.vector("TINV", "vCvcf (output)", vCvcf);
        }

        return vCvcf;
    }

    // =========================================================================
    // Section G — Index builders  (called once during MCSRunner initialization)
    // buildSiteToCfIndex(), FlatEmbData.build(), maxEmbPerCfColumn(),
    // totalCfEmbeddingCount()
    // =========================================================================

    public static final class CsrSiteToCfIndex {
        public final int[] offsets;  // length N+1
        public final int[] dataL;    // length totalEntries
        public final int[] dataEI;   // length totalEntries

        public CsrSiteToCfIndex(int[] offsets, int[] dataL, int[] dataEI) {
            this.offsets = offsets;
            this.dataL   = dataL;
            this.dataEI  = dataEI;
        }
    }

    /**
     * Builds a per-site CSR index into cfEmbeddings for efficient local ΔE computation.
     * For each site, records all (cfCol, embIdx) pairs where that site participates.
     */
    public static CsrSiteToCfIndex buildSiteToCfIndex(List<List<Embedding>> cfEmbeddings, int N) {
        int[] counts = new int[N];
        for (int l = 0; l < cfEmbeddings.size(); l++) {
            List<Embedding> embs = cfEmbeddings.get(l);
            if (embs == null) continue;
            for (int ei = 0; ei < embs.size(); ei++) {
                for (int site : embs.get(ei).getSiteIndices()) {
                    counts[site]++;
                }
            }
        }

        int totalEntries = 0;
        int[] offsets = new int[N + 1];
        for (int i = 0; i < N; i++) {
            offsets[i] = totalEntries;
            totalEntries += counts[i];
        }
        offsets[N] = totalEntries;

        int[] dataL  = new int[totalEntries];
        int[] dataEI = new int[totalEntries];

        int[] currentOffsets = offsets.clone();
        for (int l = 0; l < cfEmbeddings.size(); l++) {
            List<Embedding> embs = cfEmbeddings.get(l);
            if (embs == null) continue;
            for (int ei = 0; ei < embs.size(); ei++) {
                for (int site : embs.get(ei).getSiteIndices()) {
                    int pos = currentOffsets[site]++;
                    dataL[pos]  = l;
                    dataEI[pos] = ei;
                }
            }
        }

        if (MCSDebug.ENABLED) {
            int maxPerSite = 0;
            for (int i = 0; i < N; i++) maxPerSite = Math.max(maxPerSite, counts[i]);
            MCSDebug.log("CF-IDX", "Built CSR siteToCfIndex: N=%d, totalEntries=%d, maxPerSite=%d, avgPerSite=%.1f",
                    N, totalEntries, maxPerSite, (double) totalEntries / N);
        }

        return new CsrSiteToCfIndex(offsets, dataL, dataEI);
    }

    /**
     * Computes total energy in the CVCF basis:
     *   E = N × Σ_l eciCvcf[l] × vCvcf[l]
     */
    public static double totalEnergyCvcf(
            LatticeConfig config,
            List<List<Embedding>> cfEmbeddings, double[] flatBasisMatrix,
            int ncf, double[] eciCvcf, CvCfBasis basis, int numComp) {

        int N = config.getN();
        double[] uOrth = measureCVsFromConfig(config, cfEmbeddings, flatBasisMatrix, ncf, numComp);
        double[] vCvcf = applyTinvTransform(uOrth, config.composition(), basis);

        double E = 0.0;
        for (int l = 0; l < ncf && l < eciCvcf.length; l++)
            E += eciCvcf[l] * vCvcf[l];
        E *= N;

        boolean trace = MCSDebug.ENABLED && !dbgTotalEnergyTraced;
        if (trace) {
            dbgTotalEnergyTraced = true;
            MCSDebug.separator("TOTAL ENERGY E(σ) — CVCF BASIS");
            MCSDebug.vector("E-CVCF", "uOrth (measured)", uOrth);
            MCSDebug.vector("E-CVCF", "vCvcf (after Tinv)", vCvcf);
            MCSDebug.log("E-CVCF", "E(σ) = %.10f  (E/site = %.10f, N=%d)", E, E / N, N);
        }

        return E;
    }

    /**
     * Fully flattened representation of cfEmbeddings for zero-allocation hot paths.
     *
     * <p>Layout:
     *   cfOffsets[l]..cfOffsets[l+1]     → embedding indices for CF column l
     *   embSiteStart[e]..embSiteStart[e+1] → site/alpha indices for embedding e
     *   siteData[k], alphaData[k]         → site index and alpha for the k-th site of embedding e
     */
    public static final class FlatEmbData {
        public final int   ncf;
        public final int[] cfOffsets;      // length ncf+1
        public final int[] embSiteStart;   // length totalEmb+1
        public final int[] siteData;
        public final int[] alphaData;
        public final int[] cfEmbCount;     // length ncf

        private FlatEmbData(int ncf, int[] cfOffsets, int[] embSiteStart,
                            int[] siteData, int[] alphaData, int[] cfEmbCount) {
            this.ncf          = ncf;
            this.cfOffsets    = cfOffsets;
            this.embSiteStart = embSiteStart;
            this.siteData     = siteData;
            this.alphaData    = alphaData;
            this.cfEmbCount   = cfEmbCount;
        }

        public static FlatEmbData build(List<List<Embedding>> cfEmbeddings) {
            int ncf = cfEmbeddings.size();
            int[] cfOffsets  = new int[ncf + 1];
            int[] cfEmbCount = new int[ncf];

            int totalEmb   = 0;
            int totalSites = 0;
            for (int l = 0; l < ncf; l++) {
                List<Embedding> embs = cfEmbeddings.get(l);
                int ne = (embs != null) ? embs.size() : 0;
                cfOffsets[l]  = totalEmb;
                cfEmbCount[l] = ne;
                totalEmb     += ne;
                if (embs != null)
                    for (Embedding e : embs) totalSites += e.getSiteIndices().length;
            }
            cfOffsets[ncf] = totalEmb;

            int[] embSiteStart = new int[totalEmb + 1];
            int[] siteData     = new int[totalSites];
            int[] alphaData    = new int[totalSites];

            int embIdx  = 0;
            int siteIdx = 0;
            for (int l = 0; l < ncf; l++) {
                List<Embedding> embs = cfEmbeddings.get(l);
                if (embs == null) continue;
                for (Embedding e : embs) {
                    int[] sites  = e.getSiteIndices();
                    int[] alphas = e.getAlphaIndices();
                    embSiteStart[embIdx] = siteIdx;
                    for (int k = 0; k < sites.length; k++) {
                        siteData[siteIdx]  = sites[k];
                        alphaData[siteIdx] = alphas[k];
                        siteIdx++;
                    }
                    embIdx++;
                }
            }
            embSiteStart[totalEmb] = siteIdx;

            return new FlatEmbData(ncf, cfOffsets, embSiteStart, siteData, alphaData, cfEmbCount);
        }
    }

    /**
     * Reusable scratch space for {@link #deltaEExchangeCvcf}.
     * Allocate once per ExchangeStep, reuse across all trial moves.
     */
    public static final class DeltaScratch {
        final double[] oldSumDelta;
        final double[] newSumDelta;
        final double[] deltaUOrth;
        final double[] deltaVCvcf;
        final int[]    affectedL;
        final int[]    affectedEI;
        final boolean[] seen;
        int affectedCount;
        final int[]     affectedCols;
        final boolean[] seenCol;
        int             affectedColCount;
        final int ncf;

        public DeltaScratch(int ncf, int totalEmbeddings) {
            this.ncf          = ncf;
            this.oldSumDelta  = new double[ncf];
            this.newSumDelta  = new double[ncf];
            this.deltaUOrth   = new double[ncf];
            this.deltaVCvcf   = new double[ncf];
            this.seen         = new boolean[totalEmbeddings];
            this.affectedL    = new int[totalEmbeddings];
            this.affectedEI   = new int[totalEmbeddings];
            this.affectedCols = new int[ncf];
            this.seenCol      = new boolean[ncf];
        }

        void reset(int[] lastAffectedL, int[] lastAffectedEI, int lastCount) {
            for (int a = 0; a < lastCount; a++) {
                oldSumDelta[lastAffectedL[a]] = 0.0;
                newSumDelta[lastAffectedL[a]] = 0.0;
                deltaUOrth[lastAffectedL[a]]  = 0.0;
                deltaVCvcf[lastAffectedL[a]]  = 0.0;
            }
            affectedCount = 0;
        }

        void cleanup(int maxEmbPerCol) {
            int ac = affectedCount;
            int cc = affectedColCount;
            for (int a = 0; a < ac; a++)
                seen[affectedL[a] * maxEmbPerCol + affectedEI[a]] = false;
            for (int a = 0; a < cc; a++) {
                int l = affectedCols[a];
                seenCol[l]     = false;
                oldSumDelta[l] = 0.0;
                newSumDelta[l] = 0.0;
                deltaUOrth[l]  = 0.0;
                deltaVCvcf[l]  = 0.0;
            }
            affectedCount    = 0;
            affectedColCount = 0;
        }

        static int pairKey(int cfCol, int embIdx, int maxEmbPerCol) {
            return cfCol * maxEmbPerCol + embIdx;
        }
    }

    public static int totalCfEmbeddingCount(List<List<Embedding>> cfEmbeddings) {
        int total = 0;
        for (List<Embedding> embs : cfEmbeddings)
            if (embs != null) total += embs.size();
        return total;
    }

    public static int maxEmbPerCfColumn(List<List<Embedding>> cfEmbeddings) {
        int max = 0;
        for (List<Embedding> embs : cfEmbeddings)
            if (embs != null) max = Math.max(max, embs.size());
        return max;
    }

    // =========================================================================
    // Section E — Step 5: ΔE computation (hot path, ~N×(nEquil+nAvg) calls)
    // Primary: deltaEExchangeCvcf(FlatEmbData)  — zero object indirections
    // Fallback: deltaEExchangeCvcf(List-based)  — used when flat data unavailable
    // =========================================================================

    /**
     * List-based ΔE: used when FlatEmbData is unavailable.
     * Computes ΔE for swapping occupations at sites i and j.
     */
    public static double deltaEExchangeCvcf(
            int i, int j,
            LatticeConfig config,
            List<List<Embedding>> cfEmbeddings, double[] flatBasisMatrix,
            CsrSiteToCfIndex siteToCfIndex,
            int ncf, double[] eciCvcf, CvCfBasis basis,
            DeltaScratch scratch, int maxEmbPerCol, double[] eciOrth, int numComp) {

        int[] occ = config.getRawOcc();
        int occI = occ[i];
        int occJ = occ[j];
        if (occI == occJ) return 0.0;

        int N = config.getN();

        int ac = 0;
        int startI = siteToCfIndex.offsets[i], endI = siteToCfIndex.offsets[i+1];
        for (int idx = startI; idx < endI; idx++) {
            int l  = siteToCfIndex.dataL[idx];
            int ei = siteToCfIndex.dataEI[idx];
            int key = l * maxEmbPerCol + ei;
            if (!scratch.seen[key]) {
                scratch.seen[key] = true;
                scratch.affectedL[ac]  = l;
                scratch.affectedEI[ac] = ei;
                ac++;
            }
        }
        int startJ = siteToCfIndex.offsets[j], endJ = siteToCfIndex.offsets[j+1];
        for (int idx = startJ; idx < endJ; idx++) {
            int l  = siteToCfIndex.dataL[idx];
            int ei = siteToCfIndex.dataEI[idx];
            int key = l * maxEmbPerCol + ei;
            if (!scratch.seen[key]) {
                scratch.seen[key] = true;
                scratch.affectedL[ac]  = l;
                scratch.affectedEI[ac] = ei;
                ac++;
            }
        }
        scratch.affectedCount = ac;

        for (int a = 0; a < ac; a++) {
            int l  = scratch.affectedL[a];
            int ei = scratch.affectedEI[a];
            Embedding e  = cfEmbeddings.get(l).get(ei);
            int[] sites  = e.getSiteIndices();
            int[] alphas = e.getAlphaIndices();
            double prod  = 1.0;
            for (int k = 0; k < sites.length; k++)
                prod *= flatBasisMatrix[occ[sites[k]] * numComp + alphas[k]];
            scratch.oldSumDelta[l] += prod;
        }

        occ[i] = occJ;
        occ[j] = occI;

        for (int a = 0; a < ac; a++) {
            int l  = scratch.affectedL[a];
            int ei = scratch.affectedEI[a];
            Embedding e  = cfEmbeddings.get(l).get(ei);
            int[] sites  = e.getSiteIndices();
            int[] alphas = e.getAlphaIndices();
            double prod  = 1.0;
            for (int k = 0; k < sites.length; k++)
                prod *= flatBasisMatrix[occ[sites[k]] * numComp + alphas[k]];
            scratch.newSumDelta[l] += prod;
        }

        occ[i] = occI;
        occ[j] = occJ;

        int cc = 0;
        for (int a = 0; a < ac; a++) {
            int l = scratch.affectedL[a];
            if (!scratch.seenCol[l]) {
                scratch.seenCol[l] = true;
                scratch.affectedCols[cc++] = l;
                double diff = scratch.newSumDelta[l] - scratch.oldSumDelta[l];
                if (diff != 0.0)
                    scratch.deltaUOrth[l] = diff / cfEmbeddings.get(l).size();
            }
        }
        scratch.affectedColCount = cc;

        double dE = 0.0;
        if (eciOrth != null) {
            for (int a = 0; a < cc; a++) {
                int m = scratch.affectedCols[a];
                dE += scratch.deltaUOrth[m] * eciOrth[m];
            }
        } else if (basis != null && basis.Tinv != null) {
            double[][] Tinv = basis.Tinv;
            int tCols = Tinv[0].length;
            for (int l = 0; l < ncf && l < Tinv.length; l++) {
                double sum = 0.0;
                for (int m = 0; m < ncf && m < tCols; m++)
                    sum += Tinv[l][m] * scratch.deltaUOrth[m];
                scratch.deltaVCvcf[l] = sum;
            }
            for (int l = 0; l < ncf && l < eciCvcf.length; l++)
                dE += eciCvcf[l] * scratch.deltaVCvcf[l];
        } else {
            for (int l = 0; l < ncf && l < eciCvcf.length; l++)
                dE += eciCvcf[l] * scratch.deltaUOrth[l];
        }
        dE *= N;

        for (int a = 0; a < ac; a++)
            scratch.seen[scratch.affectedL[a] * maxEmbPerCol + scratch.affectedEI[a]] = false;
        for (int a = 0; a < cc; a++) {
            int l = scratch.affectedCols[a];
            scratch.seenCol[l]     = false;
            scratch.oldSumDelta[l] = 0.0;
            scratch.newSumDelta[l] = 0.0;
            scratch.deltaUOrth[l]  = 0.0;
            scratch.deltaVCvcf[l]  = 0.0;
        }
        scratch.affectedColCount = 0;

        return dE;
    }

    /**
     * Flat-array ΔE: primary hot path with zero object indirections.
     * Uses FlatEmbData instead of List-of-Lists-of-Embedding.
     *
     * <p>NOTE: scratch is left populated after return. The caller reads
     * oldSumDelta/newSumDelta for incremental CF updates on acceptance,
     * then calls scratch.cleanup(maxEmbPerCol).
     */
    public static double deltaEExchangeCvcf(
            int i, int j,
            LatticeConfig config,
            FlatEmbData flat, double[] flatBasisMatrix,
            CsrSiteToCfIndex siteToCfIndex,
            int ncf, double[] eciCvcf, CvCfBasis basis,
            DeltaScratch scratch, int maxEmbPerCol, double[] eciOrth, int numComp) {

        int[] occ = config.getRawOcc();
        int occI = occ[i];
        int occJ = occ[j];
        if (occI == occJ) return 0.0;

        int N = config.getN();

        int ac = 0;
        int startI = siteToCfIndex.offsets[i], endI = siteToCfIndex.offsets[i + 1];
        for (int idx = startI; idx < endI; idx++) {
            int l  = siteToCfIndex.dataL[idx];
            int ei = siteToCfIndex.dataEI[idx];
            int key = l * maxEmbPerCol + ei;
            if (!scratch.seen[key]) {
                scratch.seen[key] = true;
                scratch.affectedL[ac]  = l;
                scratch.affectedEI[ac] = ei;
                ac++;
            }
        }
        int startJ = siteToCfIndex.offsets[j], endJ = siteToCfIndex.offsets[j + 1];
        for (int idx = startJ; idx < endJ; idx++) {
            int l  = siteToCfIndex.dataL[idx];
            int ei = siteToCfIndex.dataEI[idx];
            int key = l * maxEmbPerCol + ei;
            if (!scratch.seen[key]) {
                scratch.seen[key] = true;
                scratch.affectedL[ac]  = l;
                scratch.affectedEI[ac] = ei;
                ac++;
            }
        }
        scratch.affectedCount = ac;

        for (int a = 0; a < ac; a++) {
            int l   = scratch.affectedL[a];
            int ei  = scratch.affectedEI[a];
            int flatEmbIdx = flat.cfOffsets[l] + ei;
            int sStart = flat.embSiteStart[flatEmbIdx];
            int sEnd   = flat.embSiteStart[flatEmbIdx + 1];
            double prod = 1.0;
            for (int k = sStart; k < sEnd; k++)
                prod *= flatBasisMatrix[occ[flat.siteData[k]] * numComp + flat.alphaData[k]];
            scratch.oldSumDelta[l] += prod;
        }

        occ[i] = occJ;
        occ[j] = occI;

        for (int a = 0; a < ac; a++) {
            int l   = scratch.affectedL[a];
            int ei  = scratch.affectedEI[a];
            int flatEmbIdx = flat.cfOffsets[l] + ei;
            int sStart = flat.embSiteStart[flatEmbIdx];
            int sEnd   = flat.embSiteStart[flatEmbIdx + 1];
            double prod = 1.0;
            for (int k = sStart; k < sEnd; k++)
                prod *= flatBasisMatrix[occ[flat.siteData[k]] * numComp + flat.alphaData[k]];
            scratch.newSumDelta[l] += prod;
        }

        occ[i] = occI;
        occ[j] = occJ;

        int cc = 0;
        for (int a = 0; a < ac; a++) {
            int l = scratch.affectedL[a];
            if (!scratch.seenCol[l]) {
                scratch.seenCol[l] = true;
                scratch.affectedCols[cc++] = l;
                double diff = scratch.newSumDelta[l] - scratch.oldSumDelta[l];
                if (diff != 0.0)
                    scratch.deltaUOrth[l] = diff / flat.cfEmbCount[l];
            }
        }
        scratch.affectedColCount = cc;

        double dE = 0.0;
        if (eciOrth != null) {
            for (int a = 0; a < cc; a++) {
                int m = scratch.affectedCols[a];
                dE += scratch.deltaUOrth[m] * eciOrth[m];
            }
        } else if (basis != null && basis.Tinv != null) {
            double[][] Tinv = basis.Tinv;
            int tCols = Tinv[0].length;
            for (int l = 0; l < ncf && l < Tinv.length; l++) {
                double sum = 0.0;
                for (int m = 0; m < ncf && m < tCols; m++)
                    sum += Tinv[l][m] * scratch.deltaUOrth[m];
                scratch.deltaVCvcf[l] = sum;
            }
            for (int l = 0; l < ncf && l < eciCvcf.length; l++)
                dE += eciCvcf[l] * scratch.deltaVCvcf[l];
        } else {
            for (int l = 0; l < ncf && l < eciCvcf.length; l++)
                dE += eciCvcf[l] * scratch.deltaUOrth[l];
        }
        dE *= N;

        // Leave scratch populated — caller reads deltas, then calls cleanup()
        return dE;
    }

    // ── Section C helpers ─────────────────────────────────────────────────────

    /** Parses alpha index from site symbol (e.g. "s1" → 1). 1-indexed to match flatBasisMatrix layout. */
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

    // =========================================================================
    // Section A — Inner data types
    // =========================================================================

    /** A single directed embedding of a cluster onto specific physical lattice sites. */
    public static class Embedding {

        private final int   clusterType;
        private final int   orbitMemberIndex;
        private final int[] siteIndices;
        private final int[] alphaIndices;
        // anchorIndex: canonical getAllSites() position of siteIndices[0]. -1 for CF embeddings.
        private final int   anchorIndex;

        public Embedding(int clusterType, int orbitMemberIndex,
                         int[] siteIndices, int[] alphaIndices, int anchorIndex) {
            this.clusterType      = clusterType;
            this.orbitMemberIndex = orbitMemberIndex;
            this.siteIndices      = siteIndices;
            this.alphaIndices     = alphaIndices;
            this.anchorIndex      = anchorIndex;
        }

        public int   getClusterType()      { return clusterType; }
        public int   getOrbitMemberIndex() { return orbitMemberIndex; }
        public int[] getSiteIndices()      { return siteIndices; }
        public int   size()                { return siteIndices.length; }
        public int[] getAlphaIndices()     { return alphaIndices; }
        public int   getAnchorIndex()      { return anchorIndex; }

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

    /** Intermediate representation of one (orbit member, anchor) pair for template-based generation. */
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
