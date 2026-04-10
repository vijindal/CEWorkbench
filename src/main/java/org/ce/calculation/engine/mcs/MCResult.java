package org.ce.calculation.engine.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.cluster.cvcf.CvCfBasis;
import java.util.Arrays;

/** Immutable result from a completed Monte Carlo simulation. */
public class MCResult {

    private final double   temperature;
    private final double[] composition;
    private final double[] avgCFs;
    private final double   energyPerSite;
    private final double   hmixPerSite;
    private final double   heatCapacityPerSite;
    private final double   acceptRate;
    private final long     nEquilSweeps;
    private final long     nAvgSweeps;
    private final int      supercellSize;
    private final int      nSites;

    // Statistics fields (publication-quality)
    private final double   tauInt;            // integrated autocorrelation time (sweeps)
    private final double   statInefficiency;  // s = 1 + 2·τ_int
    private final int      nEff;              // effective independent samples
    private final int      blockSizeUsed;     // auto-determined block size
    private final int      nBlocks;
    private final double   stdEnergyPerSite;
    private final double   stdHmixPerSite;
    private final double[] stdCFs;
    private final double   cvJackknife;       // unbiased Cv from jackknife
    private final double   cvStdErr;          // jackknife SEM of Cv

    // Backward-compatible 11-arg constructor — delegates to new 21-arg with NaN/0 defaults
    MCResult(double temperature, double[] composition, double[] avgCFs,
             double energyPerSite, double hmixPerSite, double heatCapacityPerSite,
             double acceptRate, long nEquilSweeps, long nAvgSweeps,
             int supercellSize, int nSites) {
        this(temperature, composition, avgCFs,
             energyPerSite, hmixPerSite, heatCapacityPerSite,
             acceptRate, nEquilSweeps, nAvgSweeps,
             supercellSize, nSites,
             Double.NaN, Double.NaN, 0, 0, 0,
             Double.NaN, Double.NaN, null, Double.NaN, Double.NaN);
    }

    // New 21-arg constructor with all statistics fields
    MCResult(double temperature, double[] composition, double[] avgCFs,
             double energyPerSite, double hmixPerSite, double heatCapacityPerSite,
             double acceptRate, long nEquilSweeps, long nAvgSweeps,
             int supercellSize, int nSites,
             double tauInt, double statInefficiency, int nEff,
             int blockSizeUsed, int nBlocks,
             double stdEnergyPerSite, double stdHmixPerSite, double[] stdCFs,
             double cvJackknife, double cvStdErr) {
        this.temperature         = temperature;
        this.composition         = composition.clone();
        this.avgCFs              = avgCFs.clone();
        this.energyPerSite       = energyPerSite;
        this.hmixPerSite         = hmixPerSite;
        this.heatCapacityPerSite = heatCapacityPerSite;
        this.acceptRate          = acceptRate;
        this.nEquilSweeps        = nEquilSweeps;
        this.nAvgSweeps          = nAvgSweeps;
        this.supercellSize       = supercellSize;
        this.nSites              = nSites;
        this.tauInt              = tauInt;
        this.statInefficiency    = statInefficiency;
        this.nEff                = nEff;
        this.blockSizeUsed       = blockSizeUsed;
        this.nBlocks             = nBlocks;
        this.stdEnergyPerSite    = stdEnergyPerSite;
        this.stdHmixPerSite      = stdHmixPerSite;
        this.stdCFs              = stdCFs != null ? stdCFs.clone() : null;
        this.cvJackknife         = cvJackknife;
        this.cvStdErr            = cvStdErr;
    }

    public double   getTemperature()         { return temperature; }
    public double[] getComposition()         { return composition.clone(); }
    public double[] getAvgCFs()              { return avgCFs.clone(); }
    public double   getEnergyPerSite()       { return energyPerSite; }
    public double   getHmixPerSite()         { return hmixPerSite; }
    public double   getHeatCapacityPerSite() { return heatCapacityPerSite; }
    public double   getAcceptRate()          { return acceptRate; }
    public long     getNEquilSweeps()        { return nEquilSweeps; }
    public long     getNAvgSweeps()          { return nAvgSweeps; }
    public int      getSupercellSize()       { return supercellSize; }
    public int      getNSites()              { return nSites; }

    // Statistics getters
    public double   getTauInt()              { return tauInt; }
    public double   getStatInefficiency()    { return statInefficiency; }
    public int      getNEff()                { return nEff; }
    public int      getBlockSizeUsed()       { return blockSizeUsed; }
    public int      getNBlocks()             { return nBlocks; }
    public double   getStdEnergyPerSite()    { return stdEnergyPerSite; }
    public double   getStdHmixPerSite()      { return stdHmixPerSite; }
    public double[] getStdCFs()              { return stdCFs != null ? stdCFs.clone() : null; }
    public double   getCvJackknife()         { return cvJackknife; }
    public double   getCvStdErr()            { return cvStdErr; }

    @Override
    public String toString() {
        return "MCResult{T=" + temperature + ", x=" + Arrays.toString(composition)
             + ", <E>/site=" + String.format("%.4f", energyPerSite)
             + ", acceptRate=" + String.format("%.3f", acceptRate) + "}";
    }

    /** Diagnostic printing utilities for internal verification. */
    public static final class Debug {

        private Debug() {}

        public static void printMcsSummary(MCResult result) {
            System.out.println("============================================================");
            System.out.println("  MCS CALCULATION DEBUG SUMMARY");
            System.out.println("============================================================");
            System.out.println(String.format("  Temperature   : %.2f K", result.getTemperature()));
            System.out.println(String.format("  Composition   : %s", Arrays.toString(result.getComposition())));
            System.out.println(String.format("  Acceptance    : %.2f%%", result.getAcceptRate() * 100));
            System.out.println(String.format("  Avg Sweeps    : %d", result.getNAvgSweeps()));
            System.out.println(String.format("  Tau_int       : %.2f sweeps", result.getTauInt()));
            System.out.println(String.format("  N_eff         : %d", result.getNEff()));
            System.out.println("------------------------------------------------------------");
            System.out.println(String.format("  <E>/site      : %.6f ± %.6f J/mol",
                    result.getEnergyPerSite(), result.getStdEnergyPerSite()));
            System.out.println(String.format("  <H>/site      : %.6f ± %.6f J/mol",
                    result.getHmixPerSite(), result.getStdHmixPerSite()));
            System.out.println(String.format("  Cv (Jackknife): %.6f ± %.6f J/(mol·K)",
                    result.getCvJackknife(), result.getCvStdErr()));
            System.out.println("------------------------------------------------------------");

            double[] cfs = result.getAvgCFs();
            double[] std = result.getStdCFs();
            if (cfs != null && cfs.length > 0) {
                System.out.println("  Mean Correlation Functions (CVCF basis):");
                for (int i = 0; i < cfs.length; i++) {
                    double sigma = (std != null && i < std.length) ? std[i] : 0.0;
                    System.out.println(String.format("    CF[%2d] = %+.6f ± %.6f", i, cfs[i], sigma));
                }
            } else {
                System.out.println("  No CVCF correlation functions measured.");
            }
            System.out.println("============================================================");
        }

        public static void printEciInfo(double[] eciCvcf, double[] eciOrth, CvCfBasis basis,
                                        int N, int[] orbitSizes) {
            System.out.println("============================================================");
            System.out.println("  MCS ECI & ORBIT DIAGNOSTIC");
            System.out.println("============================================================");
            if (basis == null) {
                System.out.println("  [ERROR] CvCfBasis is null; ECI names unknown.");
            } else {
                int ncf = basis.numNonPointCfs;
                System.out.println(String.format("  CVCF Basis: %s (K=%d), numNonPointCfs=%d",
                        basis.structurePhase, basis.numComponents, ncf));
                System.out.println("  Non-zero CVCF ECIs:");
                boolean found = false;
                for (int l = 0; l < ncf; l++) {
                    if (Math.abs(eciCvcf[l]) > 1e-12) {
                        System.out.println(String.format("    idx %2d: %-8s (= %-8s) : %+.6f J/mol",
                                l, basis.eciNames.get(l), basis.cfNames.get(l), eciCvcf[l]));
                        found = true;
                    }
                }
                if (!found) System.out.println("    [WARNING] All CVCF ECIs are ZERO!");
            }

            if (orbitSizes != null && N > 0) {
                System.out.println("  Orbit Multiplicities (from ClusterData):");
                for (int t = 0; t < orbitSizes.length; t++) {
                    double mult = (double) orbitSizes[t] / N;
                    System.out.println(String.format("    Orbit %2d: counts=%-6d mult=%6.2f", t, orbitSizes[t], mult));
                }
            }

            if (eciOrth != null) {
                System.out.println("  Transformed ECIs (Metropolis ECIs, orthogonal basis):");
                for (int t = 0; t < eciOrth.length; t++) {
                    if (Math.abs(eciOrth[t]) > 1e-12) {
                        System.out.println(String.format("    Orbit %2d : %+.6f", t, eciOrth[t]));
                    }
                }
            }
            System.out.println("============================================================");
        }
    }
}

