package org.ce.domain.engine.mcs;

import org.ce.domain.cluster.Cluster;

import java.util.List;
import java.util.Arrays;

/** Accumulates running averages of thermodynamic observables during the averaging phase. */
public class MCSampler {

    private final int    tc;
    private final int[]  orbitSizes;
    private final int    N;
    private final double R;
    private final double[]            hmixCoeff;
    private final int[]               multiSiteEmbedCounts;
    private final double[]            cfNumScratch;

    private double   sumHmix  = 0.0;
    private double   sumHmix2 = 0.0;
    private double[] sumCF;
    private long     nSamples = 0;

    public MCSampler(int N, int[] orbitSizes, List<List<Cluster>> orbits, double R,
                     double[] hmixCoeff, int[] multiSiteEmbedCounts) {
        if (N <= 0) throw new IllegalArgumentException("N must be > 0");
        if (R <= 0) throw new IllegalArgumentException("R must be > 0");
        this.N                    = N;
        this.tc                   = orbitSizes.length;
        this.orbitSizes           = orbitSizes.clone();
        this.sumCF                = new double[tc];
        this.R                    = R;
        this.hmixCoeff            = hmixCoeff.clone();
        this.multiSiteEmbedCounts = multiSiteEmbedCounts.clone();
        this.cfNumScratch         = new double[tc];
    }

    public void sample(LatticeConfig config, EmbeddingData emb) {
        for (int t = 0; t < tc; t++) cfNumScratch[t] = 0.0;

        for (Embedding e : emb.getAllEmbeddings()) {
            int t    = e.getClusterType();
            int size = e.size();
            if (t >= tc || size <= 1) continue;
            double phi = LocalEnergyCalc.clusterProduct(e, config);
            cfNumScratch[t] += phi;
        }

        double hmix_per_site = 0.0;
        for (int t = 0; t < tc; t++) {
            int embedCnt = multiSiteEmbedCounts[t];
            if (embedCnt > 0) {
                double u = cfNumScratch[t] / embedCnt;
                sumCF[t]      += u;
                hmix_per_site += hmixCoeff[t] * u;
            }
        }

        double Hmix = hmix_per_site * N;
        sumHmix  += Hmix;
        sumHmix2 += Hmix * Hmix;
        nSamples++;
    }

    public long getSampleCount() { return nSamples; }

    public double meanHmixPerSite() {
        return nSamples == 0 ? 0.0 : (sumHmix / nSamples) / N;
    }

    public double heatCapacityPerSite(double T) {
        if (nSamples < 2) return 0.0;
        double mH  = sumHmix  / nSamples;
        double mH2 = sumHmix2 / nSamples;
        return (mH2 - mH * mH) / ((double) N * R * T * T);
    }

    public double[] meanCFs() {
        double[] r = new double[tc];
        if (nSamples == 0) return r;
        for (int t = 0; t < tc; t++) r[t] = sumCF[t] / nSamples;
        return r;
    }

    public void reset() {
        sumHmix = 0; sumHmix2 = 0;
        sumCF = new double[tc];
        nSamples = 0;
    }
}
