package org.ce.workflow;

import org.ce.domain.engine.CVMEngine;
import org.ce.domain.engine.CVMInput;
import org.ce.domain.engine.MCSEngine;
import org.ce.domain.engine.MCSInput;
import org.ce.domain.result.EquilibriumState;

/**
 * Central service for thermodynamic calculations.
 *
 * This class orchestrates thermodynamic engines
 * without implementing any physics itself.
 */
public class CalculationService {

    /**
     * Runs a CVM thermodynamic calculation.
     */
    public EquilibriumState runCVM(CVMInput input) {

        CVMEngine engine = new CVMEngine(input);

        return engine.solve();
    }

    /**
     * Runs a Monte Carlo thermodynamic calculation.
     */
    public EquilibriumState runMCS(MCSInput input) {

        MCSEngine engine = new MCSEngine(input);

        return engine.solve();
    }
}
