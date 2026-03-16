package org.ce.domain.engine;

import org.ce.domain.result.EquilibriumState;

/**
 * Thermodynamic engine implementing the Cluster Variation Method (CVM).
 *
 * This class will later delegate to the CVM solver implementation.
 */
public class CVMEngine implements ThermodynamicEngine {

    private final CVMInput input;

    public CVMEngine(CVMInput input) {
        this.input = input;
    }

    @Override
    public EquilibriumState solve() {

        // Placeholder implementation.
        // Real CVM solver will be integrated later.

        double energy = 0.0;
        double freeEnergy = 0.0;

        return new EquilibriumState(
                input.getTemperature(),
                input.getComposition(),
                energy,
                freeEnergy
        );
    }
}
