package org.ce.model.result;

/**
 * Represents the thermodynamic equilibrium state
 * produced by CVM or MCS engines.
 */
public class EquilibriumState {

    private final double temperature;
    private final double[] composition;
    private final double energy;
    private final double freeEnergy;
    private final double entropy;         // NaN if not from CVM
    private final double enthalpyStdErr;   // NaN if not from MCS
    private final double heatCapacity;     // NaN if not from MCS
    private final double[] optimizedCFs;   // null if not from CVM
    private final double[] avgCFs;         // null if not from MCS
    private final double[] stdCFs;         // null if not from MCS

    // Backward-compatible 4-arg constructor
    public EquilibriumState(
            double temperature,
            double[] composition,
            double energy,
            double freeEnergy) {
        this(temperature, composition, energy, freeEnergy, Double.NaN, Double.NaN, Double.NaN, null, null, null);
    }

    // Full 10-arg constructor
    public EquilibriumState(
            double temperature,
            double[] composition,
            double energy,
            double freeEnergy,
            double entropy,
            double enthalpyStdErr,
            double heatCapacity,
            double[] optimizedCFs,
            double[] avgCFs,
            double[] stdCFs) {

        this.temperature    = temperature;
        this.composition    = composition.clone();
        this.energy         = energy;
        this.freeEnergy     = freeEnergy;
        this.entropy        = entropy;
        this.enthalpyStdErr = enthalpyStdErr;
        this.heatCapacity   = heatCapacity;
        this.optimizedCFs   = optimizedCFs != null ? optimizedCFs.clone() : null;
        this.avgCFs         = avgCFs != null ? avgCFs.clone() : null;
        this.stdCFs         = stdCFs != null ? stdCFs.clone() : null;
    }

    public double getTemperature() { return temperature; }

    public double[] getComposition() {
        return composition.clone();
    }

    public double getEnergy() {
        return energy;
    }

    public double getFreeEnergy() {
        return freeEnergy;
    }

    /** Returns true when {@link #getFreeEnergy()} holds a physically valid G value. */
    public boolean isFreeEnergyValid() { return !Double.isNaN(freeEnergy); }

    public double getEnthalpyStdErr() { return enthalpyStdErr; }
    public double getHeatCapacity()   { return heatCapacity; }
    public double getEntropy()        { return entropy; }
    public double[] getOptimizedCFs() { return optimizedCFs != null ? optimizedCFs.clone() : null; }
    public double[] getAvgCFs()       { return avgCFs != null ? avgCFs.clone() : null; }
    public double[] getStdCFs()       { return stdCFs != null ? stdCFs.clone() : null; }
}
