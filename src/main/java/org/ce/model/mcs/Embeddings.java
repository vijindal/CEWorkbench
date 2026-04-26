package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.debug.MCSDebug;
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

    // ── MCS-DBG sampling counters (to avoid flooding) ──
    private static int dbgTotalEnergyCalls  = 0;
    private static int dbgClusterProdCalls  = 0;
    private static int dbgDeltaExchCalls    = 0;
    private static int dbgDeltaSingleCalls  = 0;
    private static int dbgMeasureCVCalls    = 0;
    private static int dbgApplyTinvCalls    = 0;

    /** Reset all debug sample counters (call at start of a new run). */
    public static void resetDebugCounters() {
        dbgTotalEnergyCalls = 0;
        dbgClusterProdCalls = 0;
        dbgDeltaExchCalls   = 0;
        dbgDeltaSingleCalls = 0;
        dbgMeasureCVCalls   = 0;
        dbgApplyTinvCalls   = 0;
    }

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

    public static double totalEnergy(LatticeConfig config,
                                     Embeddings emb,
                                     double[] eci,
                                     List<List<Cluster>> orbits) {
        double sum = 0.0;
        // ── MCS-DBG: per-orbit accumulation (first call only) ──
        boolean trace = MCSDebug.ENABLED && dbgTotalEnergyCalls == 0;
        double[] orbitSums = trace ? new double[orbits.size()] : null;
        int[]    orbitCnts = trace ? new int[orbits.size()] : null;

        for (Embedding e : emb.getAllEmbeddings()) {
            int size = e.size();
            double cp = clusterProduct(e, config);
            double contrib;
            if (size > 0)
                contrib = eci[e.getClusterType()] * cp / size;
            else
                contrib = eci[e.getClusterType()] * cp;
            sum += contrib;
            if (trace) {
                int t = e.getClusterType();
                if (t < orbitSums.length) { orbitSums[t] += contrib; orbitCnts[t]++; }
            }
        }

        if (trace) {
            dbgTotalEnergyCalls++;
            MCSDebug.separator("TOTAL ENERGY E(σ) BREAKDOWN");
            MCSDebug.log("E-TOT", "Total E(σ) = %.10f  (E/site = %.10f, N=%d)", sum, sum / config.getN(), config.getN());
            for (int t = 0; t < orbitSums.length; t++) {
                if (orbitCnts[t] > 0) {
                    MCSDebug.log("E-TOT", "  orbit[%d]: eci=%.8f, embeds=%d, Σcontrib=%.10f",
                            t, t < eci.length ? eci[t] : 0.0, orbitCnts[t], orbitSums[t]);
                }
            }
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

        boolean trace = MCSDebug.ENABLED && dbgDeltaSingleCalls < 2;
        if (trace) dbgDeltaSingleCalls++;

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
            double contrib = (size > 0) ? (energyCont / size) : energyCont;
            dE += contrib;

            // ── MCS-DBG: per-embedding ΔE detail (sampled) ──
            if (trace && Math.abs(contrib) > 1e-14) {
                MCSDebug.log("ΔE-SITE", "  site=%d t=%d eci=%.6f φold=%.4f φnew=%.4f rest=%.6f cont=%.10f",
                        i, t, eci[t], phiOld, phiNew, restProduct, contrib);
            }
        }

        if (trace) {
            MCSDebug.log("ΔE-SITE", "  site=%d %d→%d TOTAL dE=%.10f", i, oldOcc, newOcc, dE);
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

        // ── MCS-DBG: sampled deltaEExchange trace (first 5 calls) ──
        if (MCSDebug.ENABLED && dbgDeltaExchCalls < 5) {
            dbgDeltaExchCalls++;
            MCSDebug.log("ΔE-EXCH", "swap(%d↔%d) occ(%d↔%d): dEi=%.10f, dEj=%.10f, ΔE=%.10f",
                    i, j, occI, occJ, dEi, dEj, dEi + dEj);
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

                // ── MCS-DBG: first 3 embedding products for CF column 0 ──
                if (trace && l == 0 && ei < 3) {
                    MCSDebug.log("CF-MEAS", "  cf[0] emb[%d]: prod=%.8f sites=%s alphas=%s",
                            ei, prod, Arrays.toString(sites), Arrays.toString(alphas));
                }
            }
            v[l] = sum / embs.size();

            // ── MCS-DBG: per-CF summary (first sweep) ──
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
            uFull[nonPt + k] = uPoint[nPoint - 1 - k]; // Reverse the order to match cfColMap
        if (tRows > nonPt + nPoint) {
            uFull[tRows - 1] = 1.0; // The constant term
        }

        double[] vCvcf = new double[ncf];
        for (int i = 0; i < ncf; i++) {
            double sum = 0.0;
            for (int j = 0; j < Tinv[i].length && j < tRows; j++)
                sum += Tinv[i][j] * uFull[j];
            vCvcf[i] = sum;
        }

        // ── MCS-DBG: full Tinv transform trace (first call only) ──
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

    // ── CVCF-basis energy computation (replaces orthogonal engine) ────────────

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
     * Builds a per-site index into cfEmbeddings for efficient local ΔE.
     * Uses CSR (Compressed Sparse Row) format for zero-allocation, contiguous iteration.
     */
    public static CsrSiteToCfIndex buildSiteToCfIndex(List<List<Embedding>> cfEmbeddings, int N) {
        // Pass 1: Count entries per site to allocate arrays
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

        // Pass 2: Populate flat arrays
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

        // ── MCS-DBG: site-to-CF index statistics ──
        if (MCSDebug.ENABLED) {
            int maxPerSite = 0;
            for (int i = 0; i < N; i++) {
                maxPerSite = Math.max(maxPerSite, counts[i]);
            }
            MCSDebug.log("CF-IDX", "Built CSR siteToCfIndex: N=%d, totalEntries=%d, maxPerSite=%d, avgPerSite=%.1f",
                    N, totalEntries, maxPerSite, (double) totalEntries / N);
        }

        return new CsrSiteToCfIndex(offsets, dataL, dataEI);
    }

    /**
     * Computes total energy in the CVCF basis:
     *   E = N × Σ_l eciCvcf[l] × vCvcf[l]
     * where vCvcf = Tinv × uOrth and uOrth is measured from cfEmbeddings.
     */
    public static double totalEnergyCvcf(
            LatticeConfig config,
            List<List<Embedding>> cfEmbeddings, double[] flatBasisMatrix,
            int ncf, double[] eciCvcf, CvCfBasis basis, int numComp) {

        int N = config.getN();

        // Step 1: Measure orthogonal CFs from cfEmbeddings
        double[] uOrth = measureCVsFromConfig(config, cfEmbeddings, flatBasisMatrix, ncf, numComp);

        // Step 2: Transform to CVCF basis
        double[] vCvcf = applyTinvTransform(uOrth, config.composition(), basis);

        // Step 3: E = N × Σ eciCvcf[l] × vCvcf[l]
        double E = 0.0;
        for (int l = 0; l < ncf && l < eciCvcf.length; l++) {
            E += eciCvcf[l] * vCvcf[l];
        }
        E *= N;

        // ── MCS-DBG: CVCF total energy breakdown ──
        boolean trace = MCSDebug.ENABLED && dbgTotalEnergyCalls == 0;
        if (trace) {
            dbgTotalEnergyCalls++;
            MCSDebug.separator("TOTAL ENERGY E(σ) — CVCF BASIS");
            MCSDebug.vector("E-CVCF", "uOrth (measured)", uOrth);
            MCSDebug.vector("E-CVCF", "vCvcf (after Tinv)", vCvcf);
            MCSDebug.log("E-CVCF", "E(σ) = %.10f  (E/site = %.10f, N=%d)", E, E / N, N);
            for (int l = 0; l < ncf && l < eciCvcf.length; l++) {
                if (Math.abs(eciCvcf[l] * vCvcf[l]) > 1e-14) {
                    MCSDebug.log("E-CVCF", "  eci[%d]=%.8f × v[%d]=%.8f = %.10f",
                            l, eciCvcf[l], l, vCvcf[l], eciCvcf[l] * vCvcf[l]);
                }
            }
        }

        return E;
    }

    // ── Pre-allocated scratch buffers for zero-allocation ΔE computation ─────

    /**
     * Reusable scratch space for {@link #deltaEExchangeCvcf}.
     * Allocate once per ExchangeStep, reuse across all trial moves.
     */
    public static final class DeltaScratch {
        final double[] oldSumDelta;
        final double[] newSumDelta;
        final double[] deltaUOrth;
        final double[] deltaVCvcf;
        /** Flat arrays for collecting affected (cfCol, embIdx) pairs without boxing. */
        final int[] affectedL;
        final int[] affectedEI;
        /** Bit-set for deduplication (replaces HashSet<Long>). Max total embedding count. */
        final boolean[] seen;
        int affectedCount;
        /** Unique CF columns touched by this move — used for sparse dot product. */
        final int[]     affectedCols;
        final boolean[] seenCol;        // dedup flag, length = ncf
        int             affectedColCount;
        final int ncf;

        public DeltaScratch(int ncf, int totalEmbeddings) {
            this.ncf = ncf;
            this.oldSumDelta = new double[ncf];
            this.newSumDelta = new double[ncf];
            this.deltaUOrth  = new double[ncf];
            this.deltaVCvcf  = new double[ncf];
            // Max affected pairs per move: entries for site_i + entries for site_j
            // Use totalEmbeddings as the upper bound for the seen-flag array
            this.seen = new boolean[totalEmbeddings];
            // Upper bound for affected pairs: sum of all per-site entries for two sites
            // In practice much smaller; use totalEmbeddings as safe upper bound
            this.affectedL    = new int[totalEmbeddings];
            this.affectedEI   = new int[totalEmbeddings];
            this.affectedCols = new int[ncf];
            this.seenCol      = new boolean[ncf];
        }

        /** Resets scratch state for the next move. Only clears what was actually used. */
        void reset(int[] lastAffectedL, int[] lastAffectedEI, int lastCount) {
            for (int a = 0; a < lastCount; a++) {
                oldSumDelta[lastAffectedL[a]] = 0.0;
                newSumDelta[lastAffectedL[a]] = 0.0;
                deltaUOrth[lastAffectedL[a]]  = 0.0;
                deltaVCvcf[lastAffectedL[a]]  = 0.0;
            }
            affectedCount = 0;
        }

        /** Computes a unique key for (cfColumn, embeddingIndex) — used as index into seen[]. */
        static int pairKey(int cfCol, int embIdx, int maxEmbPerCol) {
            return cfCol * maxEmbPerCol + embIdx;
        }
    }

    /** Computes the total number of (cfCol, embIdx) pairs across all CF columns. */
    public static int totalCfEmbeddingCount(List<List<Embedding>> cfEmbeddings) {
        int total = 0;
        int maxPerCol = 0;
        for (List<Embedding> embs : cfEmbeddings) {
            if (embs != null) {
                total += embs.size();
                maxPerCol = Math.max(maxPerCol, embs.size());
            }
        }
        return total;
    }

    /** Returns the max number of embeddings in any single CF column. */
    public static int maxEmbPerCfColumn(List<List<Embedding>> cfEmbeddings) {
        int max = 0;
        for (List<Embedding> embs : cfEmbeddings)
            if (embs != null) max = Math.max(max, embs.size());
        return max;
    }

    /**
     * Allocation-free version of {@link #deltaEExchangeCvcf} using pre-allocated scratch.
     * This is the hot-path method called ~260K times per MCS run.
     * When eciOrth is provided, bypasses the Tinv matrix-vector multiply entirely.
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

        // Collect affected pairs using flat arrays + boolean[] seen flag (no boxing)
        int ac = 0;
        
        int startI = siteToCfIndex.offsets[i];
        int endI   = siteToCfIndex.offsets[i+1];
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
        
        int startJ = siteToCfIndex.offsets[j];
        int endJ   = siteToCfIndex.offsets[j+1];
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

        // Compute old products
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

        // Temporarily swap
        occ[i] = occJ;
        occ[j] = occI;

        // Compute new products
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

        // Undo swap
        occ[i] = occI;
        occ[j] = occJ;

        // ΔuOrth[l] = (newSum - oldSum) / nEmbeddings[l]
        // Also collect unique CF columns touched (affectedCols) for sparse dot product.
        int cc = 0;
        for (int a = 0; a < ac; a++) {
            int l = scratch.affectedL[a];
            if (!scratch.seenCol[l]) {
                scratch.seenCol[l] = true;
                scratch.affectedCols[cc++] = l;
                double diff = scratch.newSumDelta[l] - scratch.oldSumDelta[l];
                if (diff != 0.0) {
                    scratch.deltaUOrth[l] = diff / cfEmbeddings.get(l).size();
                }
            }
        }
        scratch.affectedColCount = cc;

        // Compute ΔE — use eciOrth for direct dot product (bypasses Tinv multiply)
        double dE = 0.0;
        if (eciOrth != null) {
            // Fast path: ΔE = N × Σ_m ΔuOrth[m] × eciOrth[m]
            // Only unique affected CF columns have non-zero ΔuOrth — sparse iteration
            for (int a = 0; a < cc; a++) {
                int m = scratch.affectedCols[a];
                dE += scratch.deltaUOrth[m] * eciOrth[m];
            }
        } else if (basis != null && basis.Tinv != null) {
            // Fallback: Tinv matrix-vector multiply + ECI dot product
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

        // Clean up: reset only the used entries
        for (int a = 0; a < ac; a++) {
            int key = scratch.affectedL[a] * maxEmbPerCol + scratch.affectedEI[a];
            scratch.seen[key] = false;
        }
        for (int a = 0; a < cc; a++) {
            int l = scratch.affectedCols[a];
            scratch.seenCol[l]      = false;
            scratch.oldSumDelta[l]  = 0.0;
            scratch.newSumDelta[l]  = 0.0;
            scratch.deltaUOrth[l]   = 0.0;
            scratch.deltaVCvcf[l]   = 0.0;
        }
        scratch.affectedColCount = 0;

        return dE;
    }

    /**
     * Computes ΔE for exchanging sites i↔j in the CVCF basis.
     *
     * <p>Algorithm:
     * 1. For each affected CF column l, compute oldProd for affected embeddings
     * 2. Temporarily swap i↔j
     * 3. Compute newProd for the same embeddings
     * 4. ΔuOrth[l] = (newSum - oldSum) / nEmbeddings[l]
     * 5. ΔvCvcf = Tinv × ΔuOrth
     * 6. ΔE = N × Σ eciCvcf[l] × ΔvCvcf[l]
     * 7. Undo the swap
     */
    public static double deltaEExchangeCvcf(
            int i, int j,
            LatticeConfig config,
            List<List<Embedding>> cfEmbeddings, double[][] basisMatrix,
            List<int[]>[] siteToCfIndex,
            int ncf, double[] eciCvcf, CvCfBasis basis) {

        int occI = config.getOccupation(i);
        int occJ = config.getOccupation(j);
        if (occI == occJ) return 0.0;

        int N = config.getN();

        // Collect affected (cfColumn, embeddingIndex) pairs for sites i and j
        // Use a set to avoid counting shared embeddings twice
        java.util.Set<Long> seen = new java.util.HashSet<>();
        List<int[]> affectedPairs = new ArrayList<>();
        for (int[] pair : siteToCfIndex[i]) {
            long key = ((long) pair[0] << 32) | (pair[1] & 0xFFFFFFFFL);
            if (seen.add(key)) affectedPairs.add(pair);
        }
        for (int[] pair : siteToCfIndex[j]) {
            long key = ((long) pair[0] << 32) | (pair[1] & 0xFFFFFFFFL);
            if (seen.add(key)) affectedPairs.add(pair);
        }

        // Step 1: Compute old products for affected embeddings
        double[] oldSumDelta = new double[ncf]; // accumulates (prod_old) per CF column
        for (int[] pair : affectedPairs) {
            int l  = pair[0];
            int ei = pair[1];
            Embedding e  = cfEmbeddings.get(l).get(ei);
            int[] sites  = e.getSiteIndices();
            int[] alphas = e.getAlphaIndices();
            double prod  = 1.0;
            for (int k = 0; k < sites.length; k++)
                prod *= basisMatrix[config.getOccupation(sites[k])][alphas[k] - 1];
            oldSumDelta[l] += prod;
        }

        // Step 2: Temporarily swap
        config.setOccupation(i, occJ);
        config.setOccupation(j, occI);

        // Step 3: Compute new products for same embeddings
        double[] newSumDelta = new double[ncf];
        for (int[] pair : affectedPairs) {
            int l  = pair[0];
            int ei = pair[1];
            Embedding e  = cfEmbeddings.get(l).get(ei);
            int[] sites  = e.getSiteIndices();
            int[] alphas = e.getAlphaIndices();
            double prod  = 1.0;
            for (int k = 0; k < sites.length; k++)
                prod *= basisMatrix[config.getOccupation(sites[k])][alphas[k] - 1];
            newSumDelta[l] += prod;
        }

        // Step 4: Undo the swap
        config.setOccupation(i, occI);
        config.setOccupation(j, occJ);

        // Step 5: Compute ΔuOrth[l] = (newSum - oldSum) / nEmbeddings[l]
        double[] deltaUOrth = new double[ncf];
        for (int l = 0; l < ncf; l++) {
            double diff = newSumDelta[l] - oldSumDelta[l];
            if (Math.abs(diff) > 0.0) {
                int nEmbs = cfEmbeddings.get(l).size();
                deltaUOrth[l] = diff / nEmbs;
            }
        }

        // Step 6: ΔvCvcf = Tinv × ΔuOrth  (only the non-point part changes)
        double[][] Tinv = basis.Tinv;
        double[] deltaVCvcf;
        if (Tinv != null) {
            deltaVCvcf = new double[ncf];
            int tCols = Tinv[0].length;
            for (int l = 0; l < ncf && l < Tinv.length; l++) {
                double sum = 0.0;
                // Only non-point CFs change (composition is fixed in canonical MC)
                for (int m = 0; m < ncf && m < tCols; m++) {
                    sum += Tinv[l][m] * deltaUOrth[m];
                }
                deltaVCvcf[l] = sum;
            }
        } else {
            deltaVCvcf = deltaUOrth;
        }

        // Step 7: ΔE = N × Σ eciCvcf[l] × ΔvCvcf[l]
        double dE = 0.0;
        for (int l = 0; l < ncf && l < eciCvcf.length; l++) {
            dE += eciCvcf[l] * deltaVCvcf[l];
        }
        dE *= N;

        // ── MCS-DBG: sampled CVCF ΔE trace (first 5 calls) ──
        if (MCSDebug.ENABLED && dbgDeltaExchCalls < 5) {
            dbgDeltaExchCalls++;
            MCSDebug.log("ΔE-CVCF", "swap(%d↔%d) occ(%d↔%d): ΔE=%.10f, affected=%d embeddings",
                    i, j, occI, occJ, dE, affectedPairs.size());
            if (dbgDeltaExchCalls <= 2) {
                MCSDebug.vector("ΔE-CVCF", "ΔuOrth", deltaUOrth);
                MCSDebug.vector("ΔE-CVCF", "ΔvCvcf", deltaVCvcf);
                StringBuilder contribs = new StringBuilder();
                for (int l = 0; l < ncf && l < eciCvcf.length; l++) {
                    if (Math.abs(eciCvcf[l] * deltaVCvcf[l]) > 1e-14) {
                        contribs.append(String.format(" eci[%d]×Δv[%d]=%.6f×%.6f=%.8f",
                                l, l, eciCvcf[l], deltaVCvcf[l], eciCvcf[l] * deltaVCvcf[l]));
                    }
                }
                MCSDebug.log("ΔE-CVCF", "per-CF: %s → ΔE/site=%.10f", contribs, dE / N);
            }
        }

        return dE;
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
