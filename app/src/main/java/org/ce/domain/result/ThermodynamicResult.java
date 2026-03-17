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

    public final double temperature;
    public final double[] composition;

    public final double gibbsEnergy;
    public final double enthalpy;

    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy) {

        this.temperature = temperature;
        this.composition = composition;
        this.gibbsEnergy = gibbsEnergy;
        this.enthalpy = enthalpy;
    }

    /**
     * Converts engine-specific EquilibriumState to UI-ready result.
     *
     * @param state the engine result
     * @return UI-friendly thermodynamic result
     */
    public static ThermodynamicResult from(EquilibriumState state) {
        return new ThermodynamicResult(
                state.getTemperature(),
                state.getComposition(),
                state.getFreeEnergy(),
                state.getEnergy()
        );
    }
}
