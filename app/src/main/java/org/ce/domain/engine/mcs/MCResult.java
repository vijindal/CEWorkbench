package org.ce.domain.engine.mcs;

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
}
