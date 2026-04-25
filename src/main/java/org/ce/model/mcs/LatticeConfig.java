package org.ce.model.mcs;

import org.ce.debug.MCSDebug;

import java.util.Arrays;
import java.util.Random;

/**
 * Flat occupation array representing the atomic configuration of a periodic supercell.
 * Occupation encoding: 0=A, 1=B, 2=C, ...
 */
public class LatticeConfig {

    private final int   numComp;
    private final int[] occ;
    private final Basis basis;

    public LatticeConfig(int N, int numComp) {
        if (N < 1)       throw new IllegalArgumentException("N must be >= 1");
        if (numComp < 2) throw new IllegalArgumentException("numComp must be >= 2");
        this.numComp = numComp;
        this.occ     = new int[N];
        this.basis   = new Basis(numComp);
    }

    public LatticeConfig(int[] occ, int numComp) {
        if (numComp < 2) throw new IllegalArgumentException("numComp must be >= 2");
        this.numComp = numComp;
        this.occ     = occ.clone();
        this.basis   = new Basis(numComp);
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

    public int[] getRawOcc() { return occ; }

    public int getN()       { return occ.length; }
    public int getNumComp() { return numComp; }

    /** Evaluates the orthogonal basis function phi_alpha at site occupation sigma. */
    double evaluateBasis(int alpha, int sigma) {
        return basis.evaluate(alpha, sigma);
    }

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

    // ── Orthogonal point-function basis (private implementation detail) ────────

    private static final class Basis {

        private final double[][] basisMatrix; // basisMatrix[alpha-1][sigma]

        Basis(int numComp) {
            if (numComp < 2) throw new IllegalArgumentException("numComp must be >= 2");
            this.basisMatrix = buildBasis(numComp);

            // ── MCS-DBG: print full basis function table once ──
            if (MCSDebug.ENABLED) {
                MCSDebug.separator("BASIS FUNCTION TABLE (numComp=" + numComp + ")");
                StringBuilder hdr = new StringBuilder("         ");
                for (int s = 0; s < numComp; s++) hdr.append(String.format("  σ=%-5d", s));
                MCSDebug.log("BASIS", hdr.toString());
                for (int a = 0; a < basisMatrix.length; a++) {
                    StringBuilder row = new StringBuilder(String.format("  α=%-3d  ", a + 1));
                    for (int s = 0; s < basisMatrix[a].length; s++)
                        row.append(String.format("  %+.4f", basisMatrix[a][s]));
                    MCSDebug.log("BASIS", row.toString());
                }
            }
        }

        double evaluate(int alpha, int sigma) {
            return basisMatrix[alpha - 1][sigma];
        }

        private static double[][] buildBasis(int n) {
            double[] sequence = org.ce.model.cluster.ClusterMath.buildBasis(n);
            double[][] basis = new double[n - 1][n];
            for (int alpha = 1; alpha <= n - 1; alpha++)
                for (int s = 0; s < n; s++)
                    basis[alpha - 1][s] = Math.pow(sequence[s], alpha);
            return basis;
        }
    }
}
