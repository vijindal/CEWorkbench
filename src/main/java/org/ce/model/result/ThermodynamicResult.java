package org.ce.domain.result;

/**
 * Engine-independent thermodynamic result.
 *
 * <p>Used by workflow and UI layers. Provides a stable, clean API that decouples
 * clients from engine-specific models like EquilibriumState.</p>
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

    // Backward-compatible 4-arg constructor
    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy) {
        this(temperature, composition, gibbsEnergy, enthalpy, Double.NaN, Double.NaN, Double.NaN, null, null, null);
    }

    // Full 10-arg constructor
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
    }

    /** Returns true when {@link #gibbsEnergy} holds a physically valid G value. */
    public boolean isFreeEnergyValid() { return !Double.isNaN(gibbsEnergy); }

    /**
     * Converts engine-specific EquilibriumState to UI-ready result.
     */
    public static ThermodynamicResult from(EquilibriumState state) {
        return new ThermodynamicResult(
                state.getTemperature(),
                state.getComposition(),
                state.getFreeEnergy(),
                state.getEnergy(),
                state.getEntropy(),
                state.getEnthalpyStdErr(),
                state.getHeatCapacity(),
                state.getOptimizedCFs(),
                state.getAvgCFs(),
                state.getStdCFs()
        );
    }
}
