package org.ce.domain.engine;

import org.ce.domain.result.EquilibriumState;

/**
 * Generic interface for thermodynamic engines.
 *
 * Implementations:
 *  - CVMEngine
 *  - MCSEngine
 *
 * Each engine performs a thermodynamic calculation
 * and produces an equilibrium state.
 */
public interface ThermodynamicEngine {

    /**
     * Executes the thermodynamic calculation.
     *
     * @return equilibrium thermodynamic state
     */
    EquilibriumState solve();

}
