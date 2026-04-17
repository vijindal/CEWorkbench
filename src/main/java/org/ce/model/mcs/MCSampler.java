package org.ce.model.mcs;

import org.ce.model.mcs.EmbeddingData;
import org.ce.model.mcs.LatticeConfig;
import org.ce.model.mcs.CvCfEvaluator;
import org.ce.model.cluster.Cluster;
import org.ce.model.cvm.CvCfBasis;

import java.util.List;
import java.util.ArrayList;
import java.util.logging.Logger;

/** Accumulates running averages of thermodynamic observables during the averaging phase. */
public class MCSampler {

    private static final Logger LOG = Logger.getLogger(MCSampler.class.getName());

    private final int         tc;
    private final int[]       orbitSizes;
    private final int         N;
    private final double      R;
    private final double[]    eci;               // CVCF-ordered ECIs (0..numNonPointCfs-1)
    private final int[]       multiSiteEmbedCounts;
    private final CvCfBasis     basis;             // null-safe: CVCF transform unavailable if null
    private final List<List<EmbeddingData.Embedding>> cfEmbeddings;
    private final double[][]    basisMatrix;       // Inden basis values for direct measurement
    private boolean             hmixWarnedOnce       = false;

    // Running sums: indexed by CVCF CF position (0..numNonPointCfs-1), not orbit type
    private double   sumHmix  = 0.0;
    private double   sumHmix2 = 0.0;
    private double[] sumCF;           // size numNonPointCfs
    private long     nSamples = 0;

    // Full per-sweep time series (CVCF indexed)
    private final List<Double>   seriesHmix = new ArrayList<>();
    private final List<Double>   seriesE    = new ArrayList<>();
    private       List<Double>[] seriesCF;  // size numNonPointCfs


    public MCSampler(int N, int[] orbitSizes, List<List<Cluster>> orbits, double R,
                     double[] eci, int[] multiSiteEmbedCounts, CvCfBasis basis,
                     List<List<EmbeddingData.Embedding>> cfEmbeddings, double[][] basisMatrix) {
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
        this.sumCF   = new double[ncf];
        this.seriesCF = new ArrayList[ncf];
        for (int l = 0; l < ncf; l++) {
            this.seriesCF[l] = new ArrayList<>();
        }
    }

    public void sample(LatticeConfig config, EmbeddingData emb, double currentEnergy) {
        if (cfEmbeddings == null || basis == null) {
            if (!hmixWarnedOnce) {
                LOG.warning("CVCF measurement unavailable: cfEmbeddings or basis is null.");
                hmixWarnedOnce = true;
            }
            nSamples++;
            seriesE.add(currentEnergy);
            return;
        }

        int ncf = basis.numNonPointCfs;

        // 1. Measure orthogonal CFs from configuration
        double[] uOrth = CvCfEvaluator.measureCVsFromConfig(config, cfEmbeddings, basisMatrix, ncf);

        // 2. Build full orthogonal vector: [uOrth(non-point) | uPoint(point)]
        //    Point CFs: uPoint[k] = Σ_sigma x[sigma] * seq[sigma]^(k+1)
        double[] v = applyTinvTransform(uOrth, config);

        // 3. Accumulate Hmix and CFs from the CVCF vector
        double hmix_per_site = 0.0;
        for (int l = 0; l < ncf; l++) {
            hmix_per_site += eci[l] * v[l];
            sumCF[l]      += v[l];
            seriesCF[l].add(v[l]);
        }

        double Hmix = hmix_per_site * N;
        sumHmix  += Hmix;
        sumHmix2 += Hmix * Hmix;
        nSamples++;

        seriesHmix.add(Hmix);
        seriesE.add(currentEnergy);
    }

    /**
     * Transforms orthogonal CFs to CVCF basis: v_cvcf = Tinv * u_orth_full.
     * Builds the full orthogonal vector [non-point CFs | point CFs] then multiplies by Tinv.
     * Returns only the non-point CVCF CFs (length = basis.numNonPointCfs).
     */
    private double[] applyTinvTransform(double[] uOrthNonPoint, LatticeConfig config) {
        double[][] Tinv = basis.Tinv;
        int ncf = basis.numNonPointCfs;

        if (Tinv == null) {
            // No transform available: return orthogonal CFs as-is (fallback, logged once)
            if (!hmixWarnedOnce) {
                LOG.warning("CvCf Tinv unavailable — reporting orthogonal CFs as CVCF approximation.");
                hmixWarnedOnce = true;
            }
            return uOrthNonPoint.clone();
        }

        // Point CFs: u_point[k] = Σ_sigma x[sigma] * basisSeq[sigma]^(k+1)
        int K = basis.numComponents;
        double[] x = config.composition();
        double[] basisSeq = org.ce.model.cluster.ClusterMath.buildBasis(K);
        int nPoint = K - 1;
        double[] uPoint = new double[nPoint];
        for (int k = 0; k < nPoint; k++) {
            for (int s = 0; s < K; s++) {
                uPoint[k] += x[s] * Math.pow(basisSeq[s], k + 1);
            }
        }

        // Assemble full orthogonal vector: [uOrthNonPoint | uPoint]
        int tRows = Tinv[0].length;
        double[] uFull = new double[tRows];
        int nonPt = Math.min(uOrthNonPoint.length, tRows);
        System.arraycopy(uOrthNonPoint, 0, uFull, 0, nonPt);
        for (int k = 0; k < nPoint && nonPt + k < tRows; k++) {
            uFull[nonPt + k] = uPoint[k];
        }

        // v_cvcf[i] = Σ_j Tinv[i][j] * uFull[j]  (only non-point rows needed)
        double[] vCvcf = new double[ncf];
        for (int i = 0; i < ncf; i++) {
            double sum = 0.0;
            for (int j = 0; j < Tinv[i].length && j < tRows; j++) {
                sum += Tinv[i][j] * uFull[j];
            }
            vCvcf[i] = sum;
        }
        return vCvcf;
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

    /**
     * Returns the mean CVCF correlation functions v[0..numNonPointCfs-1].
     * Returns an empty array when Tinv is unavailable (ternary/quaternary).
     */
    public double[] meanCFs() {
        if (basis == null || nSamples == 0) return new double[0];
        int ncf    = basis.numNonPointCfs;
        double[] r = new double[ncf];
        for (int l = 0; l < ncf; l++) r[l] = sumCF[l] / nSamples;
        return r;
    }

    /**
     * Returns the raw Hmix time series (in J, not J/site).
     * Access after all averaging sweeps for post-processing.
     */
    public List<Double> getSeriesHmix() {
        return new ArrayList<>(seriesHmix);
    }

    /**
     * Returns the raw total energy time series (in J).
     * Access after all averaging sweeps for post-processing.
     */
    public List<Double> getSeriesE() {
        return new ArrayList<>(seriesE);
    }

    /**
     * Returns the raw CVCF correlation function time series.
     * Access after all averaging sweeps for post-processing.
     */
    public List<Double>[] getSeriesCF() {
        @SuppressWarnings("unchecked")
        List<Double>[] copy = new ArrayList[seriesCF.length];
        for (int i = 0; i < seriesCF.length; i++) {
            copy[i] = new ArrayList<>(seriesCF[i]);
        }
        return copy;
    }

    public void reset() {
        sumHmix = 0; sumHmix2 = 0;
        int ncf = (basis != null) ? basis.numNonPointCfs : 0;
        sumCF = new double[ncf];
        nSamples = 0;
        seriesHmix.clear();
        seriesE.clear();
        for (int l = 0; l < ncf; l++) {
            seriesCF[l].clear();
        }
        hmixWarnedOnce        = false;
    }
}
