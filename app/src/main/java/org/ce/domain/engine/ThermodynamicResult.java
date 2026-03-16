package org.ce.domain.engine;

/**
 * Result of a thermodynamic calculation.
 */
public class ThermodynamicResult {

    public final double freeEnergy;

    public final double[] equilibriumVariables;

    public ThermodynamicResult(double freeEnergy, double[] equilibriumVariables) {
        this.freeEnergy = freeEnergy;
        this.equilibriumVariables = equilibriumVariables;
    }
}
