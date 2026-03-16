package org.ce.domain.result;

/**
 * Represents the thermodynamic equilibrium state
 * produced by CVM or MCS engines.
 */
public class EquilibriumState {

    private final double temperature;
    private final double[] composition;
    private final double energy;
    private final double freeEnergy;

    public EquilibriumState(
            double temperature,
            double[] composition,
            double energy,
            double freeEnergy) {

        this.temperature = temperature;
        this.composition = composition.clone();
        this.energy = energy;
        this.freeEnergy = freeEnergy;
    }

    public double getTemperature() {
        return temperature;
    }

    public double[] getComposition() {
        return composition.clone();
    }

    public double getEnergy() {
        return energy;
    }

    public double getFreeEnergy() {
        return freeEnergy;
    }
}
