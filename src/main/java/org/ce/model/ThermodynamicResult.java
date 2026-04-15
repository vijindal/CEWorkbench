package org.ce.model;

/**
 * Thermodynamic equilibrium state produced by CVM or MCS engines.
 *
 * <p>Used by workflow and UI layers. Provides a stable, clean API with standard
 * thermodynamic field names (gibbsEnergy, enthalpy, entropy, stdEnthalpy, heatCapacity).</p>
 *
 * <p>All fields are public for direct, ergonomic access by consumers.</p>
 */
public class ThermodynamicResult {

    public final double   temperature;
    public final double[] composition;

    public final double   gibbsEnergy;
    public final double   enthalpy;
    public final double   entropy;       // NaN if not from CVM
    public final double   stdEnthalpy;   // NaN if not from MCS
    public final double   heatCapacity;  // NaN if not from MCS
    public final double[] optimizedCFs;  // null if not from CVM
    public final double[] avgCFs;        // null if not from MCS
    public final double[] stdCFs;        // null if not from MCS
    /** The specific property that was requested for this calculation (optional). */
    public final org.ce.calculation.CalculationDescriptor.Property requestedProp;

    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy,
            double entropy,
            double stdEnthalpy,
            double heatCapacity,
            double[] optimizedCFs,
            double[] avgCFs,
            double[] stdCFs,
            org.ce.calculation.CalculationDescriptor.Property requestedProp) {

        this.temperature  = temperature;
        this.composition  = composition.clone();
        this.gibbsEnergy  = gibbsEnergy;
        this.enthalpy     = enthalpy;
        this.entropy      = entropy;
        this.stdEnthalpy  = stdEnthalpy;
        this.heatCapacity = heatCapacity;
        this.optimizedCFs = optimizedCFs != null ? optimizedCFs.clone() : null;
        this.avgCFs       = avgCFs != null ? avgCFs.clone() : null;
        this.stdCFs       = stdCFs != null ? stdCFs.clone() : null;
        this.requestedProp = requestedProp;
    }

    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy,
            double entropy,
            double stdEnthalpy,
            double heatCapacity,
            double[] optimizedCFs,
            double[] avgCFs,
            double[] stdCFs) {
        this(temperature, composition, gibbsEnergy, enthalpy, entropy, stdEnthalpy, heatCapacity, optimizedCFs, avgCFs, stdCFs, null);
    }

    /** Returns true when {@link #gibbsEnergy} holds a physically valid G value. */
    public boolean isFreeEnergyValid() { return !Double.isNaN(gibbsEnergy); }
}
