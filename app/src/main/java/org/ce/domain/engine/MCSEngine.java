package org.ce.domain.engine;

import org.ce.domain.result.EquilibriumState;

/**
 * Thermodynamic engine implementing Monte Carlo simulation.
 *
 * This class will later delegate to the Monte Carlo solver.
 */
public class MCSEngine implements ThermodynamicEngine {

    private final MCSInput input;

    public MCSEngine(MCSInput input) {
        this.input = input;
    }

    @Override
    public EquilibriumState solve() {

        // Placeholder implementation.
        // Real Monte Carlo solver will be integrated later.

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
