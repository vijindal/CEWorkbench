package org.ce.domain.engine;

import org.ce.domain.result.EquilibriumState;

/**
 * Thermodynamic engine implementing Monte Carlo simulation.
 *
 * Placeholder — real Monte Carlo solver will be integrated later.
 */
public class MCSEngine implements ThermodynamicEngine {

    @Override
    public EquilibriumState compute(ThermodynamicInput input) throws Exception {
        // Placeholder implementation.
        return new EquilibriumState(
                input.temperature,
                input.composition,
                0.0,
                0.0
        );
    }
}
