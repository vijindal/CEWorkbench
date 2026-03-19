package org.ce.domain.engine.mcs;

import org.ce.domain.cluster.Cluster;

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
            emitUpdate(s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.EQUILIBRATION,
                    step.acceptRate(), startTime);
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
            sampler.sample(config, emb);
            emitUpdate(nEquil + s + 1, currentEnergy, sweepDeltaE, MCSUpdate.Phase.AVERAGING,
                    step.acceptRate(), startTime);
        }

        LOG.fine("MCEngine.run — done: acceptRate=" + String.format("%.3f", step.acceptRate()));
        return buildResult(config, sampler, step.acceptRate(), currentEnergy);
    }

    private void emitUpdate(int sweepNum, double currentEnergy, double sweepDeltaE,
                            MCSUpdate.Phase phase, double acceptRate, long startTime) {
        if (updateListener == null) return;
        long elapsedMs = System.currentTimeMillis() - startTime;
        updateListener.accept(new MCSUpdate(
                sweepNum, currentEnergy, sweepDeltaE,
                deltaEWindow.getStdDev(), deltaEWindow.getMean(),
                phase, acceptRate, System.currentTimeMillis(), elapsedMs));
    }

    private MCResult buildResult(LatticeConfig config, MCSampler sampler,
                                 double acceptRate, double currentEnergy) {
        int      L    = (int) Math.round(Math.cbrt(config.getN() / 2.0));
        double[] cfs  = sampler.meanCFs();
        double   hmix = sampler.meanHmixPerSite();
        int      N    = config.getN();
        return new MCResult(T, config.composition(), cfs, currentEnergy / N, hmix,
                sampler.heatCapacityPerSite(T), acceptRate, nEquil, nAvg, L, N);
    }
}
