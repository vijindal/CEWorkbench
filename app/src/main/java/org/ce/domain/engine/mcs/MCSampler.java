package org.ce.domain.engine.mcs;

import org.ce.domain.cluster.Cluster;

import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

/** Accumulates running averages of thermodynamic observables during the averaging phase. */
public class MCSampler {

    private final int    tc;
    private final int[]  orbitSizes;
    private final int    N;
    private final double R;
    private final double[]            hmixCoeff;
    private final int[]               multiSiteEmbedCounts;
    private final double[]            cfNumScratch;

    // Running sums (backward compatibility)
    private double   sumHmix  = 0.0;
    private double   sumHmix2 = 0.0;
    private double[] sumCF;
    private long     nSamples = 0;

    // Full per-sweep time series stored during averaging (reset on reset())
    private final List<Double>   seriesHmix = new ArrayList<>();
    private final List<Double>   seriesE    = new ArrayList<>();
    private       List<Double>[] seriesCF;

    // Results populated by computeStatistics() after all averaging sweeps finish
    private double   tauInt            = Double.NaN;
    private double   statInefficiency  = Double.NaN;  // s = 1 + 2·τ_int
    private int      nEff              = 0;
    private int      blockSizeUsed     = 0;
    private int      nBlocks           = 0;
    private double   meanE             = Double.NaN;  // true ⟨E⟩/N
    private double   stdE              = Double.NaN;  // SEM of ⟨E⟩/N
    private double   stdHmix           = Double.NaN;  // SEM of ⟨H⟩/N
    private double[] stdCF;                           // [tc] SEM of each CF
    private double   cvJackknife       = Double.NaN;  // Cv from jackknife
    private double   cvStdErr          = Double.NaN;  // jackknife SEM of Cv

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
        this.seriesCF             = new ArrayList[tc];
        for (int t = 0; t < tc; t++) {
            this.seriesCF[t] = new ArrayList<>();
        }
    }

    public void sample(LatticeConfig config, EmbeddingData emb, double currentEnergy) {
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
            double u = (embedCnt > 0) ? cfNumScratch[t] / embedCnt : 0.0;
            if (embedCnt > 0) {
                sumCF[t]      += u;
                hmix_per_site += hmixCoeff[t] * u;
            }
            seriesCF[t].add(u);  // store normalized CF (same units as meanCFs())
        }

        double Hmix = hmix_per_site * N;
        sumHmix  += Hmix;
        sumHmix2 += Hmix * Hmix;
        nSamples++;

        // Store time series for statistics computation
        seriesHmix.add(Hmix);
        seriesE.add(currentEnergy);
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

    /**
     * Computes statistics from the stored time series after all averaging sweeps.
     * Implements Sokal's automatic windowing estimator for τ_int, automatic block
     * averaging, and jackknife resampling for unbiased Cv.
     */
    public void computeStatistics(double T) {
        int n = seriesHmix.size();
        if (n < 4) return;

        // 1. τ_int for Hmix (primary observable for block-size determination)
        double[] hmixArr = toArray(seriesHmix);
        tauInt           = computeTauInt(hmixArr);
        statInefficiency = 1.0 + 2.0 * tauInt;
        nEff             = (int) Math.max(1, Math.round(n / statInefficiency));

        // 2. Automatic block size
        blockSizeUsed = Math.max(20, (int) Math.ceil(5.0 * tauInt));
        nBlocks       = n / blockSizeUsed;   // discard trailing partial block

        if (nBlocks < 2) {
            // Not enough blocks to compute statistics reliably
            return;
        }

        // 3. Block averages for E, Hmix, each CF
        double[] blockE    = blockMeans(toArray(seriesE), blockSizeUsed, nBlocks);
        double[] blockHmix = blockMeans(toArray(seriesHmix), blockSizeUsed, nBlocks);

        meanE   = mean(blockE) / N;
        stdE    = sem(blockE) / N;
        stdHmix = sem(blockHmix) / N;
        stdCF   = new double[tc];
        for (int t = 0; t < tc; t++) {
            if (seriesCF[t].size() == n) {
                double[] blk = blockMeans(toArray(seriesCF[t]), blockSizeUsed, nBlocks);
                stdCF[t] = sem(blk);
            }
        }

        // 4. Jackknife estimate of Cv and its error
        //    Cv_full = Var(Hmix)/(N·R·T²)  computed from block means
        double[] cvJK = new double[nBlocks];
        for (int j = 0; j < nBlocks; j++) {
            // leave-one-out: mean and variance of Hmix excluding block j
            double sumH = 0, sumH2 = 0;
            for (int k = 0; k < nBlocks; k++) {
                if (k == j) continue;
                sumH  += blockHmix[k];
                sumH2 += blockHmix[k] * blockHmix[k];
            }
            int m = nBlocks - 1;
            double mH  = sumH / m;
            double varH = (sumH2 - m * mH * mH) / (m - 1);   // unbiased sample variance
            cvJK[j] = varH / ((double) N * R * T * T);
        }
        // Full Cv using all blocks
        {
            double varH = sampleVariance(blockHmix);
            cvJackknife = varH / ((double) N * R * T * T);
        }
        // Jackknife SEM: σ² = ((n-1)/n) · Σ(Cv_j - Cv_mean)²
        double cvJKmean = mean(cvJK);
        double sumSq = 0;
        for (double v : cvJK) sumSq += (v - cvJKmean) * (v - cvJKmean);
        cvStdErr = Math.sqrt((double)(nBlocks - 1) / nBlocks * sumSq);
    }

    /**
     * Sokal automatic windowing estimator for integrated autocorrelation time τ_int.
     * Stops when the windowing condition t >= 5.0 * tau is satisfied.
     */
    private static double computeTauInt(double[] x) {
        int n = x.length;
        double mean = mean(x);
        double C0 = 0;
        for (double v : x) C0 += (v - mean) * (v - mean);
        C0 /= n;
        if (C0 < 1e-15) return 0.5;
        double tau = 0.5;
        for (int t = 1; t < n / 2; t++) {
            double Ct = 0;
            for (int s = 0; s < n - t; s++) Ct += (x[s] - mean) * (x[s + t] - mean);
            Ct /= (n - t);
            double rho = Ct / C0;
            tau += rho;
            if (t >= 5.0 * tau) break;   // Sokal windowing condition
            if (rho < 0 && t > 10) break; // early stop on decorrelation
        }
        return Math.max(0.5, tau);
    }

    /** Compute block means from series x divided into nBlocks blocks of size B. */
    private static double[] blockMeans(double[] x, int B, int nBlocks) {
        double[] blocks = new double[nBlocks];
        for (int b = 0; b < nBlocks; b++) {
            double sum = 0;
            for (int i = 0; i < B; i++) {
                sum += x[b * B + i];
            }
            blocks[b] = sum / B;
        }
        return blocks;
    }

    /** Compute arithmetic mean of array. */
    private static double mean(double[] x) {
        if (x.length == 0) return 0.0;
        double sum = 0;
        for (double v : x) sum += v;
        return sum / x.length;
    }

    /** Compute standard error of the mean (SEM) = std(x) / sqrt(n). */
    private static double sem(double[] x) {
        if (x.length < 2) return 0.0;
        double m = mean(x);
        double sumSq = 0;
        for (double v : x) sumSq += (v - m) * (v - m);
        double variance = sumSq / (x.length - 1);
        return Math.sqrt(variance / x.length);
    }

    /** Compute unbiased sample variance = Σ(x-mean)² / (n-1). */
    private static double sampleVariance(double[] x) {
        if (x.length < 2) return 0.0;
        double m = mean(x);
        double sumSq = 0;
        for (double v : x) sumSq += (v - m) * (v - m);
        return sumSq / (x.length - 1);
    }

    /** Convert List<Double> to primitive double[]. */
    private static double[] toArray(List<Double> list) {
        double[] arr = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            arr[i] = list.get(i);
        }
        return arr;
    }

    // ===== Public getters for statistics =====

    public double getTauInt()             { return tauInt; }
    public double getStatInefficiency()   { return statInefficiency; }
    public int getNEff()                  { return nEff; }
    public int getBlockSizeUsed()         { return blockSizeUsed; }
    public int getNBlocks()               { return nBlocks; }
    public double meanEnergyPerSite()     { return meanE; }
    public double stdEnergyPerSite()      { return stdE; }
    public double stdHmixPerSite()        { return stdHmix; }
    public double[] stdCFs()              { return stdCF != null ? stdCF.clone() : null; }
    public double cvJackknife()           { return cvJackknife; }
    public double cvStdErr()              { return cvStdErr; }

    public void reset() {
        sumHmix = 0; sumHmix2 = 0;
        sumCF = new double[tc];
        nSamples = 0;
        seriesHmix.clear();
        seriesE.clear();
        for (int t = 0; t < tc; t++) {
            seriesCF[t].clear();
        }
        tauInt            = Double.NaN;
        statInefficiency  = Double.NaN;
        nEff              = 0;
        blockSizeUsed     = 0;
        nBlocks           = 0;
        meanE             = Double.NaN;
        stdE              = Double.NaN;
        stdHmix           = Double.NaN;
        stdCF             = null;
        cvJackknife       = Double.NaN;
        cvStdErr          = Double.NaN;
    }
}
