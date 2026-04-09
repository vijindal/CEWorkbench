package org.ce.domain.engine.mcs;

import java.util.Arrays;
import java.util.Random;

/**
 * Flat occupation array representing the atomic configuration of a periodic supercell.
 * Occupation encoding: 0=A, 1=B, 2=C, ...
 */
public class LatticeConfig {

    private final int   numComp;
    private final int[] occ;
    private final SiteOperatorBasis basis;

    public LatticeConfig(int N, int numComp) {
        if (N < 1)       throw new IllegalArgumentException("N must be >= 1");
        if (numComp < 2) throw new IllegalArgumentException("numComp must be >= 2");
        this.numComp = numComp;
        this.occ     = new int[N];
        this.basis   = new SiteOperatorBasis(numComp);
    }

    public LatticeConfig(int[] occ, int numComp) {
        if (numComp < 2) throw new IllegalArgumentException("numComp must be >= 2");
        this.numComp = numComp;
        this.occ     = occ.clone();
        this.basis   = new SiteOperatorBasis(numComp);
        for (int i = 0; i < this.occ.length; i++) {
            if (this.occ[i] < 0 || this.occ[i] >= numComp)
                throw new IllegalArgumentException("occ[" + i + "]=" + occ[i] + " out of range");
        }
    }

    public int getOccupation(int i) { return occ[i]; }

    public void setOccupation(int i, int occ) {
        if (occ < 0 || occ >= numComp)
            throw new IllegalArgumentException("occ must be in [0," + (numComp-1) + "], got " + occ);
        this.occ[i] = occ;
    }

    public int getN()               { return occ.length; }
    public int getNumComp()         { return numComp; }
    SiteOperatorBasis getBasis() { return basis; }

    public void randomise(double[] xFrac, Random rng) {
        if (xFrac.length != numComp)
            throw new IllegalArgumentException("xFrac length " + xFrac.length + " != numComp " + numComp);
        Arrays.fill(occ, 0);
        int[] indices = new int[occ.length];
        for (int i = 0; i < occ.length; i++) indices[i] = i;
        int placed = 0;
        for (int c = 1; c < numComp; c++) {
            int count = (int) Math.round(xFrac[c] * occ.length);
            count = Math.min(count, occ.length - placed);
            for (int i = placed; i < placed + count; i++) {
                int j = i + rng.nextInt(occ.length - i);
                int tmp = indices[i]; indices[i] = indices[j]; indices[j] = tmp;
                occ[indices[i]] = c;
            }
            placed += count;
        }
    }

    public LatticeConfig copy() { return new LatticeConfig(occ, numComp); }

    public int countSpecies(int c) {
        int count = 0;
        for (int o : occ) if (o == c) count++;
        return count;
    }

    public double[] composition() {
        double[] x = new double[numComp];
        for (int o : occ) x[o]++;
        for (int c = 0; c < numComp; c++) x[c] /= occ.length;
        return x;
    }
}
