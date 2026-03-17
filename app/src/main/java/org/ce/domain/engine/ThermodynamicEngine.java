package org.ce.domain.engine;

import org.ce.domain.result.EquilibriumState;

/**
 * Interface for thermodynamic engines (CVM, Monte Carlo, etc.).
 */
public interface ThermodynamicEngine {

    /**
     * Runs a thermodynamic calculation.
     *
     * @param input bundled thermodynamic calculation input
     * @return equilibrium state with G, H, S, and thermodynamic properties
     * @throws Exception if calculation fails
     */
    EquilibriumState compute(ThermodynamicInput input) throws Exception;
}
