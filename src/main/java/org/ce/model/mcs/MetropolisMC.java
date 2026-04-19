package org.ce.model.mcs;

import org.ce.model.cluster.Cluster;
import org.ce.model.cvm.CvCfBasis;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CancellationException;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Canonical-ensemble Metropolis Monte Carlo algorithm for alloy simulations.
 *
 * <p>This class implements the standard MC algorithm. The seven steps are explicit:
 *
 * <pre>
 *  Step 1 — INITIALIZATION    : MCSRunner builds lattice + config, calls run()
 *  Step 2 — HAMILTONIAN       : Embeddings evaluates E(σ) via cluster expansion
 *  Step 3 — EQUILIBRATION     : sweep loop, moves accepted/rejected, no sampling
 *  Step 4 — TRIAL MOVE        : ExchangeStep selects site pair (canonical swap)
 *  Step 5 — COMPUTE ΔE        : Embeddings.deltaEExchange (local, O(neighbors))
 *  Step 6 — METROPOLIS ACCEPT : ExchangeStep.accept — min(1, exp(−ΔE/kT))
 *  Step 7 — AVERAGING + RESULT: Sampler accumulates ⟨E⟩, ⟨Hmix⟩, ⟨CF⟩, Cv series
 * </pre>
 */
public class MetropolisMC {

    private static final Logger LOG = Logger.getLogger(MetropolisMC.class.getName());
    private static final int DRIFT_CHECK_SWEEPS = 100;

    // ── Step 1: INITIALIZATION (fields set by MCSRunner before calling run()) ─

    private final Embeddings           emb;      // pre-built cluster embeddings
    private final double[]             eci;      // effective cluster interactions (orthogonal basis)
    private final List<List<Cluster>>  orbits;   // cluster orbits (for energy eval)
    private final int                  numComp;  // number of species (K)
    private final double               T;        // temperature [K]
    private final double               R;        // gas constant (units of ECI)
    private final int                  nEquil;   // equilibration sweeps
    private final int                  nAvg;     // averaging sweeps
    private final Random               rng;

    private Consumer<MCSUpdate> updateListener    = null;
    private RollingWindow       deltaEWindow      = new RollingWindow(500);
    private BooleanSupplier     cancellationCheck = () -> false;

    public MetropolisMC(Embeddings emb, double[] eci, List<List<Cluster>> orbits, int numComp,
                        double T, int nEquil, int nAvg, double R, Random rng) {
        if (T      <= 0) throw new IllegalArgumentException("T must be > 0");
        if (nEquil <  0) throw new IllegalArgumentException("nEquil must be >= 0");
        if (nAvg   <  1) throw new IllegalArgumentException("nAvg must be >= 1");
        this.emb     = emb;
        this.eci     = eci;
        this.orbits  = orbits;
        this.numComp = numComp;
        this.T       = T;
        this.R       = R;
        this.nEquil  = nEquil;
        this.nAvg    = nAvg;
        this.rng     = rng;
    }

    public void setUpdateListener(Consumer<MCSUpdate> listener) { this.updateListener = listener; }
    public void setCancellationCheck(BooleanSupplier check) {
        this.cancellationCheck = check != null ? check : () -> false;
    }

    // =========================================================================
    // Main algorithm entry point
    // =========================================================================

    public MCResult run(LatticeConfig config, Sampler sampler) {

        // ── Step 1: INITIALIZATION ────────────────────────────────────────────
        // config (lattice + occupations) and sampler were prepared by MCSRunner.
        // Here we wire up the move engine and compute the starting energy.
        ExchangeStep move = new ExchangeStep(emb, eci, orbits, numComp, T, R, rng);
        int N = config.getN();
        deltaEWindow.clear();
        long startTime = System.currentTimeMillis();

        // ── Step 2: COMPUTE INITIAL ENERGY E(σ) ──────────────────────────────
        // Full evaluation via cluster expansion. Subsequent steps update it
        // incrementally using ΔE (Step 5), resyncing every DRIFT_CHECK_SWEEPS.
        double currentEnergy = Embeddings.totalEnergy(config, emb, eci, orbits);

        // ── Step 3: EQUILIBRATION LOOP ────────────────────────────────────────
        // N trial moves = 1 sweep. Moves are accepted/rejected (Steps 4–6) but
        // no measurements are taken. Drives the system to thermal equilibrium.
        for (int s = 0; s < nEquil; s++) {
            if (cancellationCheck.getAsBoolean())
                throw new CancellationException("Cancelled during equilibration (sweep " + s + ")");
            double sweepDeltaE = 0.0;
            for (int m = 0; m < N; m++) {
                double dE = move.attempt(config);   // Steps 4, 5, 6 inside
                deltaEWindow.add(dE);
                currentEnergy += dE;
                sweepDeltaE   += dE;
            }
            if ((s + 1) % 100 == 0 || s + 1 == nEquil)
                emitUpdate(s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.EQUILIBRATION,
                        move.acceptRate(), startTime, null);
        }

        // Reset counters so averaging-phase accept rate is independent.
        move.resetCounters();
        sampler.reset();

        // ── Step 7: AVERAGING LOOP ────────────────────────────────────────────
        // Same move cycle as equilibration (Steps 4–6 each sweep), but now
        // Sampler records ⟨E⟩, ⟨Hmix⟩, ⟨CF⟩ every sweep for post-processing.
        for (int s = 0; s < nAvg; s++) {
            if (cancellationCheck.getAsBoolean())
                throw new CancellationException("Cancelled during averaging (sweep " + (nEquil + s) + ")");
            double sweepDeltaE = 0.0;
            for (int m = 0; m < N; m++) {
                double dE = move.attempt(config);   // Steps 4, 5, 6 inside
                deltaEWindow.add(dE);
                currentEnergy += dE;
                sweepDeltaE   += dE;
            }
            sampler.sample(config, emb, currentEnergy);

            // Periodic full-energy resync to suppress floating-point drift.
            if ((s + 1) % DRIFT_CHECK_SWEEPS == 0) {
                double Hfull = Embeddings.totalEnergy(config, emb, eci, orbits);
                double drift = Math.abs(currentEnergy - Hfull);
                if (drift > 1e-6 * Math.max(1.0, Math.abs(Hfull)))
                    LOG.warning(String.format(
                            "Energy drift at avg sweep %d: running=%.10f, full=%.10f, diff=%.3e",
                            nEquil + s + 1, currentEnergy, Hfull, drift));
                currentEnergy = Hfull;
            }
            if ((nEquil + s + 1) % 100 == 0 || s + 1 == nAvg)
                emitUpdate(nEquil + s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.AVERAGING,
                        move.acceptRate(), startTime, sampler.meanCFs());
        }

        LOG.fine("MetropolisMC.run — done: acceptRate=" + String.format("%.3f", move.acceptRate()));
        return buildResult(config, sampler, move.acceptRate(), currentEnergy);
    }

    // =========================================================================
    // Step 4 + 6: TRIAL MOVE + METROPOLIS ACCEPT/REJECT
    // =========================================================================

    /**
     * One canonical MC step: selects a pair of sites (i, j) with different species,
     * computes ΔE (Step 5), applies the Metropolis criterion (Step 6), and if accepted
     * swaps their occupations and updates the running energy.
     *
     * <p>Canonical ensemble — total composition is strictly conserved by every swap.
     */
    static class ExchangeStep {

        private final Embeddings           emb;
        private final double[]             eci;
        private final List<List<Cluster>>  orbits;
        private final double               beta;    // 1 / (R·T)
        private final int                  numComp;
        private final Random               rng;

        private ArrayList<Integer>[] speciesSites;  // per-species site index cache
        private boolean cacheInitialized = false;

        private long attempts = 0;
        private long accepted = 0;

        ExchangeStep(Embeddings emb, double[] eci, List<List<Cluster>> orbits,
                     int numComp, double T, double R, Random rng) {
            if (T <= 0) throw new IllegalArgumentException("T must be > 0");
            if (R <= 0) throw new IllegalArgumentException("R must be > 0");
            this.emb     = emb;
            this.eci     = eci;
            this.orbits  = orbits;
            this.numComp = numComp;
            this.beta    = 1.0 / (R * T);
            this.rng     = rng;
        }

        /**
         * Performs one trial move on {@code config}.
         *
         * <ol>
         *   <li>Step 4 — randomly pick species c1 ≠ c2, then site i of c1 and site j of c2</li>
         *   <li>Step 5 — compute ΔE = E(swap) − E(current) locally via cluster expansion</li>
         *   <li>Step 6 — Metropolis: accept if ΔE ≤ 0 or rand < exp(−β·ΔE)</li>
         * </ol>
         *
         * @return ΔE if accepted, 0 if rejected (for running energy tracking)
         */
        double attempt(LatticeConfig config) {
            attempts++;
            rebuildCacheIfNeeded(config);

            // Step 4: SELECT TRIAL MOVE — pick two sites of different species
            int c1 = randomNonEmptySpecies(-1);
            int c2 = randomNonEmptySpecies(c1);
            if (c1 < 0 || c2 < 0) return 0.0;

            ArrayList<Integer> list1 = speciesSites[c1];
            ArrayList<Integer> list2 = speciesSites[c2];
            int i = list1.get(rng.nextInt(list1.size()));
            int j = list2.get(rng.nextInt(list2.size()));

            // Step 5: COMPUTE ΔE — local cluster expansion around sites i and j
            double dE = Embeddings.deltaEExchange(i, j, config, emb, eci, orbits);

            // Step 6: METROPOLIS ACCEPT/REJECT — P = min(1, exp(−ΔE / kT))
            if (accept(dE)) {
                updateCacheForFlip(i, j, c1, c2);
                config.setOccupation(i, c2);
                config.setOccupation(j, c1);
                accepted++;
                return dE;
            }
            return 0.0;   // rejected — configuration unchanged
        }

        double acceptRate()    { return attempts == 0 ? 0.0 : (double) accepted / attempts; }
        void   resetCounters() { attempts = 0; accepted = 0; }

        /** Returns true if move should be accepted: always if ΔE ≤ 0, else with probability exp(−β·ΔE). */
        private boolean accept(double dE) {
            if (dE <= 0.0) return true;
            return rng.nextDouble() < Math.exp(-beta * dE);
        }

        private int randomNonEmptySpecies(int exclude) {
            int count = 0;
            for (int c = 0; c < numComp; c++)
                if (c != exclude && speciesSites[c].size() > 0) count++;
            if (count == 0) return -1;
            int pick = rng.nextInt(count);
            int idx  = 0;
            for (int c = 0; c < numComp; c++) {
                if (c != exclude && speciesSites[c].size() > 0) {
                    if (idx == pick) return c;
                    idx++;
                }
            }
            return -1;
        }

        @SuppressWarnings("unchecked")
        private void rebuildCacheIfNeeded(LatticeConfig config) {
            if (cacheInitialized) return;
            ArrayList<Integer>[] temp = new ArrayList[numComp];
            for (int c = 0; c < numComp; c++) temp[c] = new ArrayList<>(64);
            for (int k = 0; k < config.getN(); k++) temp[config.getOccupation(k)].add(k);
            speciesSites = temp;
            cacheInitialized = true;
        }

        private void updateCacheForFlip(int i, int j, int c1, int c2) {
            speciesSites[c1].remove((Integer) i);
            speciesSites[c2].add(i);
            speciesSites[c2].remove((Integer) j);
            speciesSites[c1].add(j);
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
         */
        void sample(LatticeConfig config, Embeddings emb, double currentEnergy) {
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

            // Measure orthogonal correlation functions from current config
            double[] uOrth = Embeddings.measureCVsFromConfig(config, cfEmbeddings, basisMatrix, ncf);

            // Transform to CVCF basis: v = Tinv · u_orth
            if (basis.Tinv == null && !hmixWarnedOnce) {
                SLOG.warning("CvCf Tinv unavailable — reporting orthogonal CFs as CVCF approximation.");
                hmixWarnedOnce = true;
            }
            double[] v = Embeddings.applyTinvTransform(uOrth, config.composition(), basis);

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
        }

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
    // Infrastructure: update events, result builder, rolling statistics
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
