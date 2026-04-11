package org.ce.model.mcs;

import org.ce.model.mcs.EmbeddingData;
import org.ce.model.mcs.LatticeConfig;
import org.ce.model.cluster.Cluster;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * One canonical MC step: selects two sites of different occupation,
 * computes deltaE, and accepts/rejects via the Metropolis criterion.
 */
public class ExchangeStep {

    private final EmbeddingData       emb;
    private final double[]            eci;
    private final List<List<Cluster>> orbits;
    private final double              beta;
    private final double              R;
    private final int                 numComp;
    private final Random              rng;

    @SuppressWarnings("unchecked")
    private ArrayList<Integer>[] speciesSites;
    private boolean cacheInitialized = false;

    private long attempts = 0;
    private long accepted = 0;

    public ExchangeStep(EmbeddingData emb, double[] eci, List<List<Cluster>> orbits,
                        int numComp, double T, double R, Random rng) {
        if (T <= 0) throw new IllegalArgumentException("T must be > 0");
        if (R <= 0) throw new IllegalArgumentException("R must be > 0");
        this.emb    = emb;
        this.eci    = eci;
        this.orbits = orbits;
        this.numComp = numComp;
        this.R      = R;
        this.beta   = 1.0 / (R * T);
        this.rng    = rng;
    }

    public double attempt(LatticeConfig config) {
        attempts++;
        rebuildCacheIfNeeded(config);

        int c1 = randomNonEmptySpecies(-1);
        int c2 = randomNonEmptySpecies(c1);
        if (c1 < 0 || c2 < 0) return 0.0;

        ArrayList<Integer> list1 = speciesSites[c1];
        ArrayList<Integer> list2 = speciesSites[c2];
        int i = list1.get(rng.nextInt(list1.size()));
        int j = list2.get(rng.nextInt(list2.size()));

        double dE = LocalEnergyCalc.deltaEExchange(i, j, config, emb, eci, orbits);

        if (accept(dE)) {
            updateCacheForFlip(i, j, c1, c2);
            config.setOccupation(i, c2);
            config.setOccupation(j, c1);
            accepted++;
            return dE;
        }
        return 0.0;
    }

    public long   getAttempts()   { return attempts; }
    public long   getAccepted()   { return accepted; }
    public double acceptRate()    { return attempts == 0 ? 0.0 : (double) accepted / attempts; }
    public void   resetCounters() { attempts = 0; accepted = 0; }
    public void   invalidateCache() { cacheInitialized = false; }

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
