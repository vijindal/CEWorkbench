package org.ce.domain.engine;

import org.ce.domain.cluster.AllClusterData;
import org.ce.storage.cec.CECEntry;

/**
 * Interface for thermodynamic engines (CVM, Monte Carlo, etc.).
 */
public interface ThermodynamicEngine {

    /**
     * Runs a thermodynamic calculation.
     *
     * @param clusterData cluster topology and correlation functions
     * @param cec CEC database
     * @param temperature temperature
     * @param composition composition vector
     *
     * @return calculation result
     */
    ThermodynamicResult compute(
            AllClusterData clusterData,
            CECEntry cec,
            double temperature,
            double[] composition
    );
}
