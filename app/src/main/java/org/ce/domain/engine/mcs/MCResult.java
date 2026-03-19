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

    MCResult(double temperature, double[] composition, double[] avgCFs,
             double energyPerSite, double hmixPerSite, double heatCapacityPerSite,
             double acceptRate, long nEquilSweeps, long nAvgSweeps,
             int supercellSize, int nSites) {
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

    @Override
    public String toString() {
        return "MCResult{T=" + temperature + ", x=" + Arrays.toString(composition)
             + ", <E>/site=" + String.format("%.4f", energyPerSite)
             + ", acceptRate=" + String.format("%.3f", acceptRate) + "}";
    }
}
