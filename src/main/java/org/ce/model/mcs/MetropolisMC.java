package org.ce.model.mcs;

import org.ce.debug.MCSDebug;
import org.ce.model.cluster.Cluster;
import org.ce.model.cvm.CvCfBasis;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Canonical-ensemble Metropolis Monte Carlo engine.
 *
 * <p>Implements the seven-step MC algorithm and owns the sweep loops.
 * See {@link MCSRunner} for the authoritative pipeline diagram.
 *
 * <p>Inner class layout (in pipeline order):
 * <pre>
 *  Fields + constructor       Step 1 state (algorithm parameters, ECI arrays, RNG)
 *  run()                      Steps 2–7  main entry point; calls ExchangeStep + Sampler
 *  ExchangeStep               Steps 4–6  trial move, ΔE, Metropolis accept/reject
 *    └─ BlockWorker           Step 4     parallel block execution (checkerboard)
 *  Sampler                    Step 7     sweep-level observable accumulation
 *  MCResult                             immutable output DTO
 *  MCSUpdate                            real-time progress event
 *  RollingWindow                        ΔE statistics (rolling mean/σ)
 * </pre>
 */
public class MetropolisMC {

    private static final Logger LOG = Logger.getLogger(MetropolisMC.class.getName());
    private static final int DRIFT_CHECK_SWEEPS = 100;

    // ── Step 1: INITIALIZATION (fields set by MCSRunner before calling run()) ─

    private final List<List<Embeddings.Embedding>> cfEmbeddings;  // per-CF-column embeddings
    private final double[]              flatBasisMatrix;  // flat basis values [occ * numComp + alpha]
    private final Embeddings.CsrSiteToCfIndex siteToCfIndex; // per-site CF embedding index
    private final int                   ncf;          // number of non-point CFs
    private final double[]              eciCvcf;      // effective cluster interactions (CVCF basis)
    private final double[]              eciOrth;      // pre-computed: eciOrth[m] = Σ_l eci[l] × Tinv[l][m]
    private final CvCfBasis             basis;        // CVCF basis transform (Tinv)
    private final int                   numComp;      // number of species (K)
    private final double                T;            // temperature [K]
    private final double                R;            // gas constant (units of ECI)
    private final int                   nEquil;       // equilibration sweeps
    private final int                   nAvg;         // averaging sweeps
    private final Random                rng;
    private final LatticeDecomposer.DecomposedLattice decomposedLattice;
    private boolean                     parallelEnabled = false;

    private Consumer<MCSUpdate> updateListener    = null;
    private RollingWindow       deltaEWindow      = new RollingWindow(500);
    private BooleanSupplier     cancellationCheck = () -> false;

    public MetropolisMC(List<List<Embeddings.Embedding>> cfEmbeddings,
                        double[][] basisMatrix,
                        Embeddings.CsrSiteToCfIndex siteToCfIndex,
                        int ncf,
                        double[] eciCvcf,
                        double[] eciOrth,
                        CvCfBasis basis,
                        int numComp,
                        double T, int nEquil, int nAvg, double R, Random rng,
                        LatticeDecomposer.DecomposedLattice decomposedLattice) {
        if (T      <= 0) throw new IllegalArgumentException("T must be > 0");
        if (nEquil <  0) throw new IllegalArgumentException("nEquil must be >= 0");
        if (nAvg   <  1) throw new IllegalArgumentException("nAvg must be >= 1");
        this.cfEmbeddings  = cfEmbeddings;
        
        // Flatten basisMatrix: padded so index is (occ * numComp + alpha) directly
        this.flatBasisMatrix = new double[numComp * numComp];
        if (basisMatrix != null) {
            for (int occ = 0; occ < numComp; occ++) {
                for (int a = 0; a < numComp - 1; a++) {
                    this.flatBasisMatrix[occ * numComp + (a + 1)] = basisMatrix[occ][a];
                }
            }
        }

        this.siteToCfIndex = siteToCfIndex;
        this.ncf           = ncf;
        this.eciCvcf       = eciCvcf;
        this.eciOrth       = eciOrth;
        this.basis         = basis;
        this.numComp       = numComp;
        this.T             = T;
        this.R             = R;
        this.nEquil        = nEquil;
        this.nAvg          = nAvg;
        this.rng           = rng;
        this.decomposedLattice = decomposedLattice;
        // Enable parallel execution only if we have a decomposition and enough sites
        this.parallelEnabled = (decomposedLattice != null && ncf > 0);
        if (parallelEnabled) {
            LOG.info(String.format("Parallel MCS enabled: %d blocks, %d colors", 
                    decomposedLattice.numBlocks, decomposedLattice.numColors));
        }
    }

    public void setUpdateListener(Consumer<MCSUpdate> listener) { this.updateListener = listener; }
    public void setCancellationCheck(BooleanSupplier check) {
        this.cancellationCheck = check != null ? check : () -> false;
    }

    // =========================================================================
    // Steps 2–7: Main run() loop — initial energy, equilibration, averaging
    // =========================================================================

    public MCResult run(LatticeConfig config, Sampler sampler) {

        // ── Step 1: INITIALIZATION ────────────────────────────────────────────
        // config (lattice + occupations) and sampler were prepared by MCSRunner.
        // Here we wire up the move engine and compute the starting energy.
        ExchangeStep move = new ExchangeStep(cfEmbeddings, flatBasisMatrix, siteToCfIndex,
                ncf, eciCvcf, eciOrth, basis, numComp, T, R, rng, decomposedLattice);
        int N = config.getN();
        deltaEWindow.clear();
        long startTime = System.currentTimeMillis();

        // ── Step 2: COMPUTE INITIAL ENERGY E(σ) — CVCF basis ─────────────────
        // Full evaluation via CVCF cluster expansion. Subsequent steps update it
        // incrementally using ΔE (Step 5), resyncing every DRIFT_CHECK_SWEEPS.
        Embeddings.resetDebugCounters();  // reset sampled counters for this run
        double currentEnergy = Embeddings.totalEnergyCvcf(
                config, cfEmbeddings, flatBasisMatrix, ncf, eciCvcf, basis, numComp);

        // ── MCS-DBG: initial energy + composition ──
        if (MCSDebug.ENABLED) {
            double[] comp = config.composition();
            MCSDebug.separator("INITIAL ENERGY E(σ₀) — CVCF basis");
            MCSDebug.log("E-INIT", "E(σ₀) = %.10f  (E/site = %.10f)", currentEnergy, currentEnergy / N);
            MCSDebug.log("E-INIT", "N = %d, T = %.2f K, nEquil = %d, nAvg = %d", N, T, nEquil, nAvg);
            MCSDebug.vector("E-INIT", "composition", comp);
        }

        // ── Step 3: EQUILIBRATION LOOP ────────────────────────────────────────
        // N trial moves = 1 sweep. Moves are accepted/rejected (Steps 4–6) but
        // no measurements are taken. Drives the system to thermal equilibrium.
        for (int s = 0; s < nEquil; s++) {
            if (cancellationCheck.getAsBoolean())
                throw new CancellationException("Cancelled during equilibration (sweep " + s + ")");
            
            double sweepDeltaE;
            if (parallelEnabled) {
                sweepDeltaE = move.runParallelSweep(config, false);
            } else {
                sweepDeltaE = 0.0;
                for (int m = 0; m < N; m++) {
                    double dE = move.attempt(config);
                    deltaEWindow.add(dE);
                    sweepDeltaE += dE;
                }
            }
            currentEnergy += sweepDeltaE;
            
            if ((s + 1) % 100 == 0 || s + 1 == nEquil)
                emitUpdate(s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.EQUILIBRATION,
                        move.acceptRate(), startTime, null);
        }

        // Reset counters so averaging-phase accept rate is independent.
        move.resetCounters();
        sampler.reset();

        // ── Incremental CF sums: initialize from post-equilibration config ────
        // runningCfSum[l] = Σ_e prod(e) for CF column l (unnormalized).
        // Updated per accepted move by ExchangeStep; avoids O(totalEmbeddings) scan each sweep.
        double[] runningCfSum = null;
        int[]    cfEmbCount   = null;
        if (ncf > 0 && cfEmbeddings != null) {
            cfEmbCount   = new int[ncf];
            runningCfSum = new double[ncf];
            double[] initUOrth = Embeddings.measureCVsFromConfig(
                    config, cfEmbeddings, flatBasisMatrix, ncf, numComp);
            for (int l = 0; l < ncf; l++) {
                cfEmbCount[l]   = cfEmbeddings.get(l).size();
                runningCfSum[l] = initUOrth[l] * cfEmbCount[l];  // un-normalize
            }
            move.runningCfSum = runningCfSum;
            move.cfEmbCount   = cfEmbCount;
        }

        // ── Step 7: AVERAGING LOOP ────────────────────────────────────────────
        // Same move cycle as equilibration (Steps 4–6 each sweep), but now
        // Sampler records ⟨E⟩, ⟨Hmix⟩, ⟨CF⟩ every sweep for post-processing.
        for (int s = 0; s < nAvg; s++) {
            if (cancellationCheck.getAsBoolean())
                throw new CancellationException("Cancelled during averaging (sweep " + (nEquil + s) + ")");

            double sweepDeltaE;
            if (parallelEnabled) {
                sweepDeltaE = move.runParallelSweep(config, true);
            } else {
                sweepDeltaE = 0.0;
                for (int m = 0; m < N; m++) {
                    double dE = move.attempt(config);
                    deltaEWindow.add(dE);
                    sweepDeltaE += dE;
                }
            }
            currentEnergy += sweepDeltaE;
            sampler.sample(config, null, currentEnergy, flatBasisMatrix, numComp,
                           runningCfSum, cfEmbCount);

            // Periodic full-energy and CF resync to suppress floating-point drift.
            if ((s + 1) % DRIFT_CHECK_SWEEPS == 0) {
                double Hfull = Embeddings.totalEnergyCvcf(
                        config, cfEmbeddings, flatBasisMatrix, ncf, eciCvcf, basis, numComp);
                double drift = Math.abs(currentEnergy - Hfull);
                if (drift > 1e-6 * Math.max(1.0, Math.abs(Hfull)))
                    LOG.warning(String.format(
                            "Energy drift at avg sweep %d: running=%.10f, full=%.10f, diff=%.3e",
                            nEquil + s + 1, currentEnergy, Hfull, drift));
                currentEnergy = Hfull;
                // Resync runningCfSum to suppress floating-point accumulation drift
                if (runningCfSum != null) {
                    double[] freshUOrth = Embeddings.measureCVsFromConfig(
                            config, cfEmbeddings, flatBasisMatrix, ncf, numComp);
                    for (int l = 0; l < ncf; l++)
                        runningCfSum[l] = freshUOrth[l] * cfEmbCount[l];
                }
            }
            if ((nEquil + s + 1) % 100 == 0 || s + 1 == nAvg)
                emitUpdate(nEquil + s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.AVERAGING,
                        move.acceptRate(), startTime, sampler.meanCFs());
        }

        LOG.fine("MetropolisMC.run — done: acceptRate=" + String.format("%.3f", move.acceptRate()));
        return buildResult(config, sampler, move.acceptRate(), currentEnergy);
    }

    // =========================================================================
    // Steps 4–6: ExchangeStep — trial move, ΔE, Metropolis accept/reject
    // =========================================================================

    /**
     * One canonical MC step: selects a pair of sites (i, j) with different species,
     * computes ΔE (Step 5), applies the Metropolis criterion (Step 6), and if accepted
     * swaps their occupations and updates the running energy.
     *
     * <p>Canonical ensemble — total composition is strictly conserved by every swap.
     */
    static class ExchangeStep {

        private final Embeddings.FlatEmbData         flatEmbData;  // cache-friendly flat embedding arrays
        private final double[]       flatBasisMatrix;
        private final Embeddings.CsrSiteToCfIndex siteToCfIndex;
        private final int            ncf;
        private final double[]       eciCvcf;
        private final double[]       eciOrth;      // pre-computed for Tinv bypass
        private final CvCfBasis      basis;
        private final double         beta;    // 1 / (R·T)
        private final int            numComp;
        private final Random         rng;

        // Pre-allocated scratch for zero-allocation ΔE computation (serial path)
        private final Embeddings.DeltaScratch scratch;
        private final int maxEmbPerCol;
        // Per-thread scratch for parallel path — avoids per-sweep BlockWorker allocation
        private final ThreadLocal<Embeddings.DeltaScratch> perThreadScratch;

        // Incremental CF sums: runningCfSum[l] = Σ_e prod(e) for CF column l (unnormalized)
        // Updated per accepted move; used by Sampler instead of full measureCVsFromConfig each sweep.
        double[] runningCfSum;   // set externally by MetropolisMC.run() before averaging
        int[]    cfEmbCount;     // cfEmbCount[l] = number of embeddings in CF column l (constant)

        // Block-based species cache for parallel-safe canonical moves
        private int[][][] blockSpeciesSites;  // [blockIdx][species][idx] = siteIndex
        private int[][]   blockSpeciesSizes;  // [blockIdx][species] = current count
        private int[]     siteToCacheIdx;      // [siteIndex] = idx in blockSpeciesSites
        private final int[] blockOfSite;       // site-to-block mapping
        private final int numBlocks;
        private final LatticeDecomposer.DecomposedLattice decomposedLattice;
        private boolean cacheInitialized = false;

        private final AtomicLong attempts = new AtomicLong();
        private final AtomicLong accepted = new AtomicLong();

        ExchangeStep(List<List<Embeddings.Embedding>> cfEmbeddings,
                     double[] flatBasisMatrix, Embeddings.CsrSiteToCfIndex siteToCfIndex,
                     int ncf, double[] eciCvcf, double[] eciOrth, CvCfBasis basis,
                     int numComp, double T, double R, Random rng,
                     LatticeDecomposer.DecomposedLattice decomposedLattice) {
            if (T <= 0) throw new IllegalArgumentException("T must be > 0");
            if (R <= 0) throw new IllegalArgumentException("R must be > 0");
            this.flatBasisMatrix = flatBasisMatrix;
            this.siteToCfIndex = siteToCfIndex;
            this.ncf           = ncf;
            this.eciCvcf       = eciCvcf;
            this.eciOrth       = eciOrth;
            this.basis         = basis;
            this.numComp       = numComp;
            this.beta          = 1.0 / (R * T);
            this.rng           = rng;
            this.maxEmbPerCol  = (cfEmbeddings != null) ? Embeddings.maxEmbPerCfColumn(cfEmbeddings) : 0;
            
            if (decomposedLattice != null) {
                this.numBlocks = decomposedLattice.numBlocks;
                this.blockOfSite = decomposedLattice.blockOfSite;
                this.decomposedLattice = decomposedLattice;
            } else {
                this.numBlocks = 1;
                this.blockOfSite = null;
                this.decomposedLattice = null;
            }

            // Flat embedding arrays — eliminates object pointer indirections in hot product loop
            this.flatEmbData = (cfEmbeddings != null) ? Embeddings.FlatEmbData.build(cfEmbeddings) : null;
            // Primary scratch for serial phases
            int seenSize = (cfEmbeddings != null) ? cfEmbeddings.size() * maxEmbPerCol : 0;
            this.scratch = (cfEmbeddings != null) ? new Embeddings.DeltaScratch(ncf, seenSize) : null;
            // Per-thread scratch for parallel phases — one instance per ForkJoin worker, reused across sweeps
            final int seenSizeFinal = seenSize;
            this.perThreadScratch = ThreadLocal.withInitial(
                    () -> new Embeddings.DeltaScratch(ncf, seenSizeFinal));
        }

        /**
         * Orchestrates a parallel sweep by iterating over block colors and 
         * updating same-colored blocks in parallel.
         */
        double runParallelSweep(LatticeConfig config, boolean sampling) {
            rebuildCacheIfNeeded(config);
            double totalDeltaE = 0.0;
            
            // Total sites per sweep = N. We distribute these across blocks.
            int N = config.getN();
            int sitesPerBlock = N / numBlocks;

            for (int colorIdx = 0; colorIdx < decomposedLattice.numColors; colorIdx++) {
                List<Integer> blocks = decomposedLattice.colors[colorIdx];
                if (blocks.isEmpty()) continue;
                
                // Parallel update of independent blocks
                double colorDeltaE = blocks.parallelStream().mapToDouble(bIdx -> {
                    BlockWorker worker = new BlockWorker(bIdx, sitesPerBlock);
                    return worker.run(config);
                }).sum();
                totalDeltaE += colorDeltaE;
            }
            return totalDeltaE;
        }

        /**
         * Encapsulates the execution logic for a single spatial block.
         * Maintains its own scratch and RNG for thread-safety.
         */
        private final class BlockWorker {
            private final int bIdx;
            private final int numMoves;
            private final ThreadLocalRandom localRng;

            BlockWorker(int bIdx, int numMoves) {
                this.bIdx = bIdx;
                this.numMoves = numMoves;
                this.localRng = ThreadLocalRandom.current();
            }

            double run(LatticeConfig config) {
                double blockDeltaE = 0.0;
                for (int m = 0; m < numMoves; m++) {
                    blockDeltaE += attemptBlock(config);
                }
                return blockDeltaE;
            }

            private double attemptBlock(LatticeConfig config) {
                attempts.incrementAndGet();
                
                int c1 = randomNonEmptySpeciesBlock(bIdx, -1, localRng);
                int c2 = randomNonEmptySpeciesBlock(bIdx, c1, localRng);
                if (c1 < 0 || c2 < 0) return 0.0;

                int i = blockSpeciesSites[bIdx][c1][localRng.nextInt(blockSpeciesSizes[bIdx][c1])];
                int j = blockSpeciesSites[bIdx][c2][localRng.nextInt(blockSpeciesSizes[bIdx][c2])];

                Embeddings.DeltaScratch s = perThreadScratch.get();
                double dE = Embeddings.deltaEExchangeCvcf(
                        i, j, config, flatEmbData, flatBasisMatrix, siteToCfIndex,
                        ncf, eciCvcf, basis, s, maxEmbPerCol, eciOrth, numComp);

                boolean accepted_move = (dE <= 0.0 || localRng.nextDouble() < Math.exp(-beta * dE));
                if (accepted_move) {
                    synchronized (config) { // Minimal sync: occupation swap + CF update + cache
                        updateCacheForFlipBlock(i, j, c1, c2, bIdx);
                        config.setOccupation(i, c2);
                        config.setOccupation(j, c1);
                        if (runningCfSum != null) {
                            for (int a = 0; a < s.affectedColCount; a++) {
                                int l = s.affectedCols[a];
                                runningCfSum[l] += s.newSumDelta[l] - s.oldSumDelta[l];
                            }
                        }
                        accepted.incrementAndGet();
                    }
                }
                s.cleanup(maxEmbPerCol);
                return accepted_move ? dE : 0.0;
            }
        }

        double attempt(LatticeConfig config) {
            // For serial phases, use block 0 or a global fallback if needed.
            // If numBlocks > 1, we should ideally pick a random block to maintain ergodicity.
            int bIdx = (numBlocks > 1) ? rng.nextInt(numBlocks) : 0;
            rebuildCacheIfNeeded(config);
            attempts.incrementAndGet();

            int c1 = randomNonEmptySpeciesBlock(bIdx, -1, rng);
            int c2 = randomNonEmptySpeciesBlock(bIdx, c1, rng);
            if (c1 < 0 || c2 < 0) return 0.0;

            int i = blockSpeciesSites[bIdx][c1][rng.nextInt(blockSpeciesSizes[bIdx][c1])];
            int j = blockSpeciesSites[bIdx][c2][rng.nextInt(blockSpeciesSizes[bIdx][c2])];

            double dE = Embeddings.deltaEExchangeCvcf(
                    i, j, config, flatEmbData, flatBasisMatrix, siteToCfIndex,
                    ncf, eciCvcf, basis, scratch, maxEmbPerCol, eciOrth, numComp);

            boolean accepted_move = (dE <= 0.0 || rng.nextDouble() < Math.exp(-beta * dE));
            if (accepted_move) {
                updateCacheForFlipBlock(i, j, c1, c2, bIdx);
                config.setOccupation(i, c2);
                config.setOccupation(j, c1);
                // Incremental CF sum update
                if (runningCfSum != null) {
                    for (int a = 0; a < scratch.affectedColCount; a++) {
                        int l = scratch.affectedCols[a];
                        runningCfSum[l] += scratch.newSumDelta[l] - scratch.oldSumDelta[l];
                    }
                }
                accepted.incrementAndGet();
            }
            scratch.cleanup(maxEmbPerCol);
            return accepted_move ? dE : 0.0;
        }

        double acceptRate()    { long a = attempts.get(); return a == 0 ? 0.0 : (double) accepted.get() / a; }
        void   resetCounters() { attempts.set(0); accepted.set(0); }

        private int randomNonEmptySpeciesBlock(int bIdx, int exclude, Random r) {
            int count = 0;
            for (int c = 0; c < numComp; c++)
                if (c != exclude && blockSpeciesSizes[bIdx][c] > 0) count++;
            if (count == 0) return -1;
            int pick = r.nextInt(count);
            int idx  = 0;
            for (int c = 0; c < numComp; c++) {
                if (c != exclude && blockSpeciesSizes[bIdx][c] > 0) {
                    if (idx == pick) return c;
                    idx++;
                }
            }
            return -1;
        }

        private synchronized void rebuildCacheIfNeeded(LatticeConfig config) {
            if (cacheInitialized) return;
            int N = config.getN();
            int[] occ = config.getRawOcc();
            
            blockSpeciesSites = new int[numBlocks][numComp][N]; // Upper bound N
            blockSpeciesSizes = new int[numBlocks][numComp];
            siteToCacheIdx = new int[N];
            
            for (int k = 0; k < N; k++) {
                int c = occ[k];
                int bIdx = (blockOfSite != null) ? blockOfSite[k] : 0;
                int pos = blockSpeciesSizes[bIdx][c]++;
                blockSpeciesSites[bIdx][c][pos] = k;
                siteToCacheIdx[k] = pos;
            }
            cacheInitialized = true;
        }

        private void updateCacheForFlipBlock(int i, int j, int c1, int c2, int bIdx) {
            // Remove site i from block bIdx species c1 list
            int lastIdx1 = --blockSpeciesSizes[bIdx][c1];
            int lastSite1 = blockSpeciesSites[bIdx][c1][lastIdx1];
            int posI = siteToCacheIdx[i];
            blockSpeciesSites[bIdx][c1][posI] = lastSite1;
            siteToCacheIdx[lastSite1] = posI;

            // Add site i to block bIdx species c2 list
            int posI2 = blockSpeciesSizes[bIdx][c2]++;
            blockSpeciesSites[bIdx][c2][posI2] = i;
            siteToCacheIdx[i] = posI2;

            // Remove site j from block bIdx species c2 list
            int lastIdx2 = --blockSpeciesSizes[bIdx][c2];
            int lastSite2 = blockSpeciesSites[bIdx][c2][lastIdx2];
            int posJ = siteToCacheIdx[j];
            blockSpeciesSites[bIdx][c2][posJ] = lastSite2;
            siteToCacheIdx[lastSite2] = posJ;

            // Add site j to block bIdx species c1 list
            int posJ1 = blockSpeciesSizes[bIdx][c1]++;
            blockSpeciesSites[bIdx][c1][posJ1] = j;
            siteToCacheIdx[j] = posJ1;
        }
    }

    // =========================================================================
    // Step 7: MEASUREMENT — accumulate ⟨E⟩, ⟨Hmix⟩, ⟨CF⟩, Cv during averaging
    // =========================================================================

    /**
     * Accumulates thermodynamic observables during the averaging phase.
     *
     * <p>Each call to {@link #sample} corresponds to one completed averaging sweep.
     * Records:
     * <ul>
     *   <li>Running Hmix per site and Hmix² (for Cv = (⟨H²⟩ − ⟨H⟩²) / (N·R·T²))</li>
     *   <li>Mean CVCF correlation functions ⟨CF_l⟩ via Tinv transform</li>
     *   <li>Full time series of Hmix, E, and per-CF values for block-error analysis</li>
     * </ul>
     */
    static class Sampler {

        private static final Logger SLOG = Logger.getLogger(Sampler.class.getName());

        private final int         tc;
        private final int[]       orbitSizes;
        private final int         N;
        private final double      R;
        private final double[]    eci;
        private final int[]       multiSiteEmbedCounts;
        private final CvCfBasis   basis;
        private final List<List<Embeddings.Embedding>> cfEmbeddings;
        private final double[][]  basisMatrix;
        private double[]          fixedComposition; // set once before averaging; canonical MC keeps composition constant
        private boolean           hmixWarnedOnce = false;

        // Accumulators for ⟨Hmix⟩ and ⟨Hmix²⟩
        private double   sumHmix  = 0.0;
        private double   sumHmix2 = 0.0;
        private double[] sumCF;
        private long     nSamples = 0;

        // Time series for post-processing (block averaging, jackknife Cv)
        private final List<Double>   seriesHmix = new ArrayList<>();
        private final List<Double>   seriesE    = new ArrayList<>();
        private       List<Double>[] seriesCF;

        @SuppressWarnings("unchecked")
        Sampler(int N, int[] orbitSizes, List<List<Cluster>> orbits, double R,
                double[] eci, int[] multiSiteEmbedCounts, CvCfBasis basis,
                List<List<Embeddings.Embedding>> cfEmbeddings, double[][] basisMatrix) {
            if (N <= 0) throw new IllegalArgumentException("N must be > 0");
            if (R <= 0) throw new IllegalArgumentException("R must be > 0");
            this.N                    = N;
            this.tc                   = orbitSizes.length;
            this.orbitSizes           = orbitSizes.clone();
            this.R                    = R;
            this.eci                  = eci.clone();
            this.multiSiteEmbedCounts = multiSiteEmbedCounts.clone();
            this.basis                = basis;
            this.cfEmbeddings         = cfEmbeddings;
            this.basisMatrix          = basisMatrix;

            int ncf = (basis != null) ? basis.numNonPointCfs : 0;
            this.sumCF    = new double[ncf];
            this.seriesCF = new ArrayList[ncf];
            for (int l = 0; l < ncf; l++) this.seriesCF[l] = new ArrayList<>();
        }

        /**
         * Records one sweep's worth of observables from the current configuration.
         * Called once per averaging sweep by MetropolisMC after N trial moves.
         * When runningCfSum/cfEmbCount are provided, uses incremental CFs instead of
         * a full O(totalEmbeddings) measurement scan.
         */
        void sample(LatticeConfig config, Embeddings emb, double currentEnergy, double[] flatBasisMatrix, int numComp,
                    double[] runningCfSum, int[] cfEmbCount) {
            if (cfEmbeddings == null || basis == null) {
                if (!hmixWarnedOnce) {
                    SLOG.warning("CVCF measurement unavailable: cfEmbeddings or basis is null.");
                    hmixWarnedOnce = true;
                }
                nSamples++;
                seriesE.add(currentEnergy);
                return;
            }

            int ncf = basis.numNonPointCfs;

            // Use incremental CF sums if available; otherwise full measurement scan
            double[] uOrth;
            if (runningCfSum != null && cfEmbCount != null) {
                uOrth = new double[ncf];
                for (int l = 0; l < ncf; l++)
                    uOrth[l] = (cfEmbCount[l] > 0) ? runningCfSum[l] / cfEmbCount[l] : 0.0;
            } else {
                uOrth = Embeddings.measureCVsFromConfig(config, cfEmbeddings, flatBasisMatrix, ncf, numComp);
            }

            // Transform to CVCF basis: v = Tinv · u_orth
            if (basis.Tinv == null && !hmixWarnedOnce) {
                SLOG.warning("CvCf Tinv unavailable — reporting orthogonal CFs as CVCF approximation.");
                hmixWarnedOnce = true;
            }
            double[] comp = (fixedComposition != null) ? fixedComposition : config.composition();
            double[] v = Embeddings.applyTinvTransform(uOrth, comp, basis);

            // Accumulate: Hmix/site = Σ_l eci[l] · v[l]
            double hmix_per_site = 0.0;
            for (int l = 0; l < ncf; l++) {
                hmix_per_site += eci[l] * v[l];
                sumCF[l]      += v[l];
                seriesCF[l].add(v[l]);
            }

            double Hmix = hmix_per_site * N;
            sumHmix  += Hmix;
            sumHmix2 += Hmix * Hmix;   // for Cv = (⟨H²⟩ − ⟨H⟩²) / (N·R·T²)
            nSamples++;
            seriesHmix.add(Hmix);
            seriesE.add(currentEnergy);

            // ── MCS-DBG: detailed Sampler trace (first 3 averaging sweeps) ──
            if (MCSDebug.ENABLED && nSamples <= 3) {
                MCSDebug.separator("SAMPLER.sample — sweep #" + nSamples);
                MCSDebug.vector("SAMP", "uOrth (measured)", uOrth);
                MCSDebug.vector("SAMP", "v (CVCF, after Tinv)", v);
                StringBuilder contribs = new StringBuilder();
                for (int l = 0; l < ncf; l++) {
                    contribs.append(String.format(" eci[%d]*v[%d]=%.6f*%.6f=%.8f",
                            l, l, eci[l], v[l], eci[l] * v[l]));
                }
                MCSDebug.log("SAMP", "Per-CF ECI contributions:%s", contribs);
                MCSDebug.log("SAMP", "hmix_per_site = %.10f", hmix_per_site);
                MCSDebug.log("SAMP", "Hmix (total)  = %.10f  (N=%d)", Hmix, N);
                MCSDebug.log("SAMP", "currentEnergy = %.10f  (E/site = %.10f)", currentEnergy, currentEnergy / N);
                MCSDebug.log("SAMP", "running ⟨Hmix⟩/site = %.10f, nSamples = %d",
                        (sumHmix / nSamples) / N, nSamples);
            }
        }

        /** Cache the (invariant) composition so applyTinvTransform avoids O(N) recount each sweep. */
        void setFixedComposition(double[] composition) { this.fixedComposition = composition.clone(); }

        long     getSampleCount()   { return nSamples; }
        double   meanHmixPerSite()  { return nSamples == 0 ? 0.0 : (sumHmix / nSamples) / N; }

        /** Cv per site from energy fluctuations: (⟨H²⟩ − ⟨H⟩²) / (N · R · T²). */
        double heatCapacityPerSite(double T) {
            if (nSamples < 2) return 0.0;
            double mH  = sumHmix  / nSamples;
            double mH2 = sumHmix2 / nSamples;
            return (mH2 - mH * mH) / ((double) N * R * T * T);
        }

        double[] meanCFs() {
            if (basis == null || nSamples == 0) return new double[0];
            int ncf    = basis.numNonPointCfs;
            double[] r = new double[ncf];
            for (int l = 0; l < ncf; l++) r[l] = sumCF[l] / nSamples;
            return r;
        }

        List<Double> getSeriesHmix() { return new ArrayList<>(seriesHmix); }
        List<Double> getSeriesE()    { return new ArrayList<>(seriesE); }

        @SuppressWarnings("unchecked")
        List<Double>[] getSeriesCF() {
            List<Double>[] copy = new ArrayList[seriesCF.length];
            for (int i = 0; i < seriesCF.length; i++) copy[i] = new ArrayList<>(seriesCF[i]);
            return copy;
        }

        void reset() {
            sumHmix = 0; sumHmix2 = 0;
            int ncf = (basis != null) ? basis.numNonPointCfs : 0;
            sumCF = new double[ncf];
            nSamples = 0;
            seriesHmix.clear();
            seriesE.clear();
            for (int l = 0; l < ncf; l++) seriesCF[l].clear();
            hmixWarnedOnce = false;
        }
    }

    // =========================================================================
    // Result & Event Types + Infrastructure: MCResult, MCSUpdate, RollingWindow
    // =========================================================================

    /** Real-time update event emitted periodically during equilibration and averaging. */
    public static class MCSUpdate {

        public enum Phase { EQUILIBRATION, AVERAGING }

        private final int      step;
        private final double   E_total;
        private final double   deltaE;
        private final double   sigmaDE;
        private final double   meanDE;
        private final Phase    phase;
        private final double   acceptanceRate;
        private final long     timestampMs;
        private final long     elapsedMs;
        private final double[] cfs;

        public MCSUpdate(int step, double E_total, double deltaE, double sigmaDE, double meanDE,
                         Phase phase, double acceptanceRate, long timestampMs, long elapsedMs) {
            this(step, E_total, deltaE, sigmaDE, meanDE, phase, acceptanceRate, timestampMs, elapsedMs, null);
        }

        public MCSUpdate(int step, double E_total, double deltaE, double sigmaDE, double meanDE,
                         Phase phase, double acceptanceRate, long timestampMs, long elapsedMs, double[] cfs) {
            this.step           = step;
            this.E_total        = E_total;
            this.deltaE         = deltaE;
            this.sigmaDE        = sigmaDE;
            this.meanDE         = meanDE;
            this.phase          = phase;
            this.acceptanceRate = acceptanceRate;
            this.timestampMs    = timestampMs;
            this.elapsedMs      = elapsedMs;
            this.cfs            = cfs != null ? cfs.clone() : null;
        }

        public int      getStep()           { return step; }
        public double   getE_total()        { return E_total; }
        public double   getDeltaE()         { return deltaE; }
        public double   getSigmaDE()        { return sigmaDE; }
        public double   getMeanDE()         { return meanDE; }
        public Phase    getPhase()          { return phase; }
        public double   getAcceptanceRate() { return acceptanceRate; }
        public long     getTimestampMs()    { return timestampMs; }
        public long     getElapsedMs()      { return elapsedMs; }
        public double[] getCfs()            { return cfs != null ? cfs.clone() : null; }

        @Override
        public String toString() {
            return String.format("MCSUpdate{step=%d, E=%.4f, phase=%s, acceptance=%.1f%%, elapsed=%dms}",
                    step, E_total, phase, acceptanceRate * 100, elapsedMs);
        }
    }

    /**
     * Immutable raw simulation output.
     * Contains sweep-averaged quantities and raw time series for {@code MCSStatisticsProcessor}.
     */
    public static class MCResult {

        private final double   temperature;
        private final double[] composition;
        private final long     nEquilSweeps;
        private final long     nAvgSweeps;
        private final int      supercellSize;
        private final int      nSites;
        private final double   acceptRate;
        private final double   energyPerSite;
        private final double   hmixPerSite;
        private final double[] avgCFs;
        private final double[] seriesHmix;
        private final double[] seriesE;
        private final double[][] seriesCF;

        public MCResult(double temperature, double[] composition, long nEquilSweeps, long nAvgSweeps,
                        int supercellSize, int nSites, double acceptRate, double energyPerSite,
                        double hmixPerSite, double[] avgCFs, double[] seriesHmix,
                        double[] seriesE, double[][] seriesCF) {
            this.temperature   = temperature;
            this.composition   = composition.clone();
            this.nEquilSweeps  = nEquilSweeps;
            this.nAvgSweeps    = nAvgSweeps;
            this.supercellSize = supercellSize;
            this.nSites        = nSites;
            this.acceptRate    = acceptRate;
            this.energyPerSite = energyPerSite;
            this.hmixPerSite   = hmixPerSite;
            this.avgCFs        = avgCFs   != null ? avgCFs.clone()       : null;
            this.seriesHmix    = seriesHmix != null ? seriesHmix.clone() : null;
            this.seriesE       = seriesE  != null ? seriesE.clone()      : null;
            this.seriesCF      = seriesCF != null ? deepClone(seriesCF)  : null;
        }

        private static double[][] deepClone(double[][] arr) {
            double[][] copy = new double[arr.length][];
            for (int i = 0; i < arr.length; i++)
                copy[i] = arr[i] != null ? arr[i].clone() : null;
            return copy;
        }

        public double   getTemperature()  { return temperature; }
        public double[] getComposition()  { return composition.clone(); }
        public long     getNEquilSweeps() { return nEquilSweeps; }
        public long     getNAvgSweeps()   { return nAvgSweeps; }
        public int      getSupercellSize(){ return supercellSize; }
        public int      getNSites()       { return nSites; }
        public double   getAcceptRate()   { return acceptRate; }
        public double   getEnergyPerSite(){ return energyPerSite; }
        public double   getHmixPerSite()  { return hmixPerSite; }
        public double[] getAvgCFs()       { return avgCFs != null ? avgCFs.clone() : null; }
        public double[] getSeriesHmix()   { return seriesHmix != null ? seriesHmix.clone() : null; }
        public double[] getSeriesE()      { return seriesE  != null ? seriesE.clone()  : null; }
        public double[][] getSeriesCF()   { return seriesCF != null ? deepClone(seriesCF) : null; }

        @Override
        public String toString() {
            return "MCResult{T=" + temperature + ", x=" + Arrays.toString(composition)
                 + ", <E>/site=" + String.format("%.4f", energyPerSite)
                 + ", acceptRate=" + String.format("%.3f", acceptRate)
                 + ", nAvg=" + nAvgSweeps + "}";
        }
    }

    /** Efficient rolling window for tracking recent ΔE statistics (σ, mean). */
    private static final class RollingWindow {

        private final Deque<Double> window;
        private final int maxSize;
        private double sum = 0.0;
        private double sumSquares = 0.0;

        RollingWindow(int maxSize) {
            if (maxSize <= 0) throw new IllegalArgumentException("maxSize must be positive");
            this.maxSize = maxSize;
            this.window  = new ArrayDeque<>();
        }

        void add(double value) {
            if (window.size() >= maxSize) {
                double removed = window.removeFirst();
                sum -= removed;
                sumSquares -= removed * removed;
            }
            window.addLast(value);
            sum += value;
            sumSquares += value * value;
        }

        double getMean() {
            if (window.isEmpty()) return 0.0;
            return sum / window.size();
        }

        double getStdDev() {
            if (window.size() < 2) return 0.0;
            double mean     = getMean();
            double variance = (sumSquares / window.size()) - (mean * mean);
            if (variance < 0) variance = 0;
            return Math.sqrt(variance);
        }

        void clear() {
            window.clear();
            sum = 0.0;
            sumSquares = 0.0;
        }
    }

    // ── Private helpers ───────────────────────────────────────────────────────

    private void emitUpdate(int sweepNum, double currentEnergy, double sweepDeltaE,
                            MCSUpdate.Phase phase, double acceptRate, long startTime, double[] cfs) {
        if (updateListener == null) return;
        long elapsedMs = System.currentTimeMillis() - startTime;
        updateListener.accept(new MCSUpdate(
                sweepNum, currentEnergy, sweepDeltaE,
                deltaEWindow.getStdDev(), deltaEWindow.getMean(),
                phase, acceptRate, System.currentTimeMillis(), elapsedMs, cfs));
    }

    private MCResult buildResult(LatticeConfig config, Sampler sampler,
                                 double acceptRate, double currentEnergy) {
        int      L    = (int) Math.round(Math.cbrt(config.getN() / 2.0));
        double[] cfs  = sampler.meanCFs();
        double   hmix = sampler.meanHmixPerSite();
        int      N    = config.getN();
        double energy = currentEnergy / N;

        // ── MCS-DBG: CROSS-VALIDATION — recompute energy from scratch ──
        if (MCSDebug.ENABLED) {
            double freshE = Embeddings.totalEnergyCvcf(
                    config, cfEmbeddings, flatBasisMatrix, ncf, eciCvcf, basis, numComp);
            double freshPerSite = freshE / N;
            double drift = Math.abs(energy - freshPerSite);
            MCSDebug.separator("CROSS-VALIDATION (buildResult)");
            MCSDebug.log("XVAL", "running E/site   = %.10f", energy);
            MCSDebug.log("XVAL", "fresh   E/site   = %.10f  (totalEnergyCvcf recompute)", freshPerSite);
            MCSDebug.log("XVAL", "drift            = %.3e", drift);
            MCSDebug.log("XVAL", "Sampler hmix/site= %.10f", hmix);
            MCSDebug.log("XVAL", "MATCH? running vs fresh: %s", drift < 1e-6 ? "YES ✓" : "NO ✗ — ΔE accumulation bug!");
            MCSDebug.log("XVAL", "MATCH? running vs hmix:  %s",
                    Math.abs(energy - hmix) < 1e-3 ? "YES ✓" : "NO ✗ — different energy paths!");
        }

        // ── MCS-DBG: final result assembly ──
        if (MCSDebug.ENABLED) {
            MCSDebug.separator("BUILD RESULT (MetropolisMC)");
            MCSDebug.log("RESULT", "nSamples      = %d", sampler.getSampleCount());
            MCSDebug.log("RESULT", "acceptRate     = %.4f", acceptRate);
            MCSDebug.log("RESULT", "currentEnergy = %.10f  (E/site = %.10f)", currentEnergy, energy);
            MCSDebug.log("RESULT", "⟨Hmix⟩/site   = %.10f", hmix);
            MCSDebug.log("RESULT", "Cv/site (fluct)= %.10f", sampler.heatCapacityPerSite(T));
            MCSDebug.vector("RESULT", "⟨CF⟩ (mean CFs)", cfs);
        }

        return new MCResult(
                T, config.composition(), nEquil, nAvg, L, N,
                acceptRate, energy, hmix, cfs,
                toArray(sampler.getSeriesHmix()),
                toArray(sampler.getSeriesE()),
                toArray2D(sampler.getSeriesCF()));
    }

    private static double[] toArray(List<Double> list) {
        if (list == null || list.isEmpty()) return null;
        double[] arr = new double[list.size()];
        for (int i = 0; i < list.size(); i++) arr[i] = list.get(i);
        return arr;
    }

    private static double[][] toArray2D(List<Double>[] lists) {
        if (lists == null || lists.length == 0) return null;
        double[][] arr = new double[lists.length][];
        for (int i = 0; i < lists.length; i++) arr[i] = toArray(lists[i]);
        return arr;
    }
}
