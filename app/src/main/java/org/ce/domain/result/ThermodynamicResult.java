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
    public final double   stdEnthalpy;   // NaN if not from MCS
    public final double   heatCapacity;  // NaN if not from MCS
    public final double[] avgCFs;        // null if not from MCS
    public final double[] stdCFs;        // null if not from MCS

    // Backward-compatible 4-arg constructor
    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy) {
        this(temperature, composition, gibbsEnergy, enthalpy, Double.NaN, Double.NaN, null, null);
    }

    // 6-arg constructor — delegates to full 8-arg with null CFs
    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy,
            double stdEnthalpy,
            double heatCapacity) {
        this(temperature, composition, gibbsEnergy, enthalpy, stdEnthalpy, heatCapacity, null, null);
    }

    // Full 8-arg constructor
    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy,
            double stdEnthalpy,
            double heatCapacity,
            double[] avgCFs,
            double[] stdCFs) {

        this.temperature  = temperature;
        this.composition  = composition;
        this.gibbsEnergy  = gibbsEnergy;
        this.enthalpy     = enthalpy;
        this.stdEnthalpy  = stdEnthalpy;
        this.heatCapacity = heatCapacity;
        this.avgCFs       = avgCFs != null ? avgCFs.clone() : null;
        this.stdCFs       = stdCFs != null ? stdCFs.clone() : null;
    }

    /**
     * Converts engine-specific EquilibriumState to UI-ready result.
     */
    public static ThermodynamicResult from(EquilibriumState state) {
        return new ThermodynamicResult(
                state.getTemperature(),
                state.getComposition(),
                state.getFreeEnergy(),
                state.getEnergy(),
                state.getEnthalpyStdErr(),
                state.getHeatCapacity(),
                state.getAvgCFs(),
                state.getStdCFs()
        );
    }
}
