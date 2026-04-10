package org.ce.calculation.engine.mcs;

import org.ce.model.mcs.EmbeddingData;
import org.ce.model.mcs.Embedding;
import org.ce.model.mcs.LatticeConfig;
import org.ce.model.cluster.Cluster;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CancellationException;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.logging.Logger;

/** Metropolis Monte Carlo engine: equilibration then averaging sweeps. */
public class MCEngine {

    private static final Logger LOG = Logger.getLogger(MCEngine.class.getName());

    private final EmbeddingData       emb;
    private final double[]            eci;
    private final List<List<Cluster>> orbits;
    private final int                 numComp;
    private final double              T;
    private final double              R;
    private final int                 nEquil;
    private final int                 nAvg;
    private final Random              rng;

    private Consumer<MCSUpdate> updateListener    = null;
    private RollingWindow       deltaEWindow      = new RollingWindow(500);
    private BooleanSupplier     cancellationCheck = () -> false;
    // Sampler reference kept for passing running-avg CFs to updateListener during averaging
    private MCSampler           activeSampler     = null;

    public MCEngine(EmbeddingData emb, double[] eci, List<List<Cluster>> orbits, int numComp,
                    double T, int nEquil, int nAvg, double R, Random rng) {
        if (T      <= 0) throw new IllegalArgumentException("T must be > 0");
        if (nEquil <  0) throw new IllegalArgumentException("nEquil must be >= 0");
        if (nAvg   <  1) throw new IllegalArgumentException("nAvg must be >= 1");
        this.emb    = emb;
        this.eci    = eci;
        this.orbits = orbits;
        this.numComp = numComp;
        this.T      = T;
        this.R      = R;
        this.nEquil = nEquil;
        this.nAvg   = nAvg;
        this.rng    = rng;
    }

    public void setUpdateListener(Consumer<MCSUpdate> listener) { this.updateListener = listener; }
    public void setCancellationCheck(BooleanSupplier check) {
        this.cancellationCheck = check != null ? check : () -> false;
    }

    public MCResult run(LatticeConfig config, MCSampler sampler) {
        this.activeSampler = sampler;
        ExchangeStep step = new ExchangeStep(emb, eci, orbits, numComp, T, R, rng);
        int N = config.getN();
        deltaEWindow.clear();
        long startTime = System.currentTimeMillis();
        double currentEnergy = LocalEnergyCalc.totalEnergy(config, emb, eci, orbits);

        for (int s = 0; s < nEquil; s++) {
            if (cancellationCheck.getAsBoolean())
                throw new CancellationException("Cancelled during equilibration (sweep " + s + ")");
            double sweepDeltaE = 0.0;
            for (int m = 0; m < N; m++) {
                double dE = step.attempt(config);
                deltaEWindow.add(dE);
                currentEnergy += dE;
                sweepDeltaE   += dE;
            }
            // Emit update every 100 sweeps or on final equilibration sweep
            if ((s + 1) % 100 == 0 || s + 1 == nEquil) {
                emitUpdate(s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.EQUILIBRATION,
                        step.acceptRate(), startTime, null);
            }
        }

        step.resetCounters();
        sampler.reset();

        for (int s = 0; s < nAvg; s++) {
            if (cancellationCheck.getAsBoolean())
                throw new CancellationException("Cancelled during averaging (sweep " + (nEquil + s) + ")");
            double sweepDeltaE = 0.0;
            for (int m = 0; m < N; m++) {
                double dE = step.attempt(config);
                deltaEWindow.add(dE);
                currentEnergy += dE;
                sweepDeltaE   += dE;
            }
            sampler.sample(config, emb, currentEnergy);
            // Emit update every 100 sweeps or on final averaging sweep
            if ((nEquil + s + 1) % 100 == 0 || s + 1 == nAvg) {
                emitUpdate(nEquil + s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.AVERAGING,
                        step.acceptRate(), startTime, sampler.meanCFs());
            }
        }

        // Compute statistics from stored time series
        sampler.computeStatistics(T);

        LOG.fine("MCEngine.run — done: acceptRate=" + String.format("%.3f", step.acceptRate()));
        return buildResult(config, sampler, step.acceptRate(), currentEnergy);
    }

    // -------------------------------------------------------------------------
    // Inner classes
    // -------------------------------------------------------------------------

    /** Efficient rolling window statistics calculator. */
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

    /** Real-time update event emitted each N sweeps during a simulation run. */
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

    private void emitUpdate(int sweepNum, double currentEnergy, double sweepDeltaE,
                            MCSUpdate.Phase phase, double acceptRate, long startTime, double[] cfs) {
        if (updateListener == null) return;
        long elapsedMs = System.currentTimeMillis() - startTime;
        updateListener.accept(new MCSUpdate(
                sweepNum, currentEnergy, sweepDeltaE,
                deltaEWindow.getStdDev(), deltaEWindow.getMean(),
                phase, acceptRate, System.currentTimeMillis(), elapsedMs, cfs));
    }

    private MCResult buildResult(LatticeConfig config, MCSampler sampler,
                                 double acceptRate, double currentEnergy) {
        int      L    = (int) Math.round(Math.cbrt(config.getN() / 2.0));
        double[] cfs  = sampler.meanCFs();
        double   hmix = sampler.meanHmixPerSite();
        int      N    = config.getN();

        // Use corrected energy from time average if available, otherwise fallback
        double energy = Double.isNaN(sampler.meanEnergyPerSite())
                ? currentEnergy / N               // fallback for < 4 samples
                : sampler.meanEnergyPerSite();    // true ⟨E⟩/N

        // Use jackknife Cv if available, otherwise fallback to old estimator
        double cv = Double.isNaN(sampler.cvJackknife())
                ? sampler.heatCapacityPerSite(T)  // old estimator as fallback
                : sampler.cvJackknife();

        return new MCResult(T, config.composition(), cfs, energy, hmix, cv,
                acceptRate, nEquil, nAvg, L, N,
                sampler.getTauInt(), sampler.getStatInefficiency(), sampler.getNEff(),
                sampler.getBlockSizeUsed(), sampler.getNBlocks(),
                sampler.stdEnergyPerSite(), sampler.stdHmixPerSite(),
                sampler.stdCFs(), sampler.cvJackknife(), sampler.cvStdErr());
    }
}
