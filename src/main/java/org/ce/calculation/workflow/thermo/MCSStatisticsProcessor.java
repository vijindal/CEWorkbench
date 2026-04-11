package org.ce.calculation.workflow.thermo;

import java.util.List;

/**
 * Calculation-layer post-processor for MCS time series statistics.
 *
 * Accepts raw time series data from the model layer (MCSampler) and computes
 * publication-quality statistical quantities: integrated autocorrelation time (τ_int),
 * automatic block averaging, standard error of the mean (SEM), and jackknife estimates
 * of heat capacity with error bars.
 *
 * This class owns all time-series analysis — not physics evaluation.
 * The model layer (MCSampler) accumulates raw observables; this class (MCSStatisticsProcessor)
 * processes them into scientific quantities.
 */
public class MCSStatisticsProcessor {

    private final double N;
    private final double R;
    private final double T;
    private final List<Double> seriesHmix;
    private final List<Double> seriesE;
    private final List<Double>[] seriesCF;

    // Computed statistics
    private double   tauInt            = Double.NaN;
    private double   statInefficiency  = Double.NaN;  // s = 1 + 2·τ_int
    private int      nEff              = 0;
    private int      blockSizeUsed     = 0;
    private int      nBlocks           = 0;
    private double   stdEnergyPerSite  = Double.NaN;  // SEM of ⟨E⟩/N
    private double   stdHmixPerSite    = Double.NaN;  // SEM of ⟨H⟩/N
    private double[] stdCFs;                          // [numNonPointCfs] SEM of each CVCF CF
    private double   cvJackknife       = Double.NaN;  // unbiased Cv from jackknife
    private double   cvStdErr          = Double.NaN;  // jackknife SEM of Cv

    /**
     * @param N number of lattice sites
     * @param R gas constant (J/(mol·K))
     * @param T temperature (K)
     * @param seriesHmix raw Hmix time series (in J, not J/site)
     * @param seriesE raw total energy time series (in J)
     * @param seriesCF raw CVCF correlation function time series per CF
     */
    public MCSStatisticsProcessor(double N, double R, double T,
                                   List<Double> seriesHmix, List<Double> seriesE,
                                   List<Double>[] seriesCF) {
        if (N <= 0) throw new IllegalArgumentException("N must be > 0");
        if (R <= 0) throw new IllegalArgumentException("R must be > 0");
        if (T <= 0) throw new IllegalArgumentException("T must be > 0");

        this.N          = N;
        this.R          = R;
        this.T          = T;
        this.seriesHmix = seriesHmix;
        this.seriesE    = seriesE;
        this.seriesCF   = seriesCF;
    }

    /**
     * Computes statistics from the stored time series.
     * Implements Sokal's automatic windowing estimator for τ_int, automatic block
     * averaging, and jackknife resampling for unbiased Cv.
     */
    public void computeStatistics() {
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

        // 3. Block averages for E, Hmix, each CVCF CF
        double[] blockE    = blockMeans(toArray(seriesE), blockSizeUsed, nBlocks);
        double[] blockHmix = blockMeans(toArray(seriesHmix), blockSizeUsed, nBlocks);

        stdEnergyPerSite = sem(blockE) / N;
        stdHmixPerSite   = sem(blockHmix) / N;
        int ncf = (seriesCF != null) ? seriesCF.length : 0;
        stdCFs  = new double[ncf];
        for (int l = 0; l < ncf; l++) {
            if (seriesCF[l].size() == n) {
                double[] blk = blockMeans(toArray(seriesCF[l]), blockSizeUsed, nBlocks);
                stdCFs[l] = sem(blk);
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
            cvJK[j] = varH / (N * R * T * T);
        }
        // Full Cv using all blocks
        {
            double varH = sampleVariance(blockHmix);
            cvJackknife = varH / (N * R * T * T);
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

    // ===== Public getters for computed statistics =====

    public double getTauInt()             { return tauInt; }
    public double getStatInefficiency()   { return statInefficiency; }
    public int getNEff()                  { return nEff; }
    public int getBlockSizeUsed()         { return blockSizeUsed; }
    public int getNBlocks()               { return nBlocks; }
    public double getStdEnergyPerSite()   { return stdEnergyPerSite; }
    public double getStdHmixPerSite()     { return stdHmixPerSite; }
    public double[] getStdCFs()           { return stdCFs != null ? stdCFs.clone() : null; }
    public double getCvJackknife()        { return cvJackknife; }
    public double getCvStdErr()           { return cvStdErr; }
}
