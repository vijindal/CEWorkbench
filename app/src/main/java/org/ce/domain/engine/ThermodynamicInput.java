package org.ce.domain.engine;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.hamiltonian.CECEntry;

/**
 * Bundles all inputs required for thermodynamic calculations.
 *
 * <p>This is the unified input to any ThermodynamicEngine implementation
 * (CVM, Monte Carlo, hybrid, etc.).</p>
 */
public class ThermodynamicInput {

    public final AllClusterData clusterData;
    public final CECEntry cec;
    public final double temperature;
    public final double[] composition;
    public final String systemId;
    public final String systemName;

    public ThermodynamicInput(
            AllClusterData clusterData,
            CECEntry cec,
            double temperature,
            double[] composition,
            String systemId,
            String systemName) {

        this.clusterData = clusterData;
        this.cec = cec;
        this.temperature = temperature;
        this.composition = composition;
        this.systemId = systemId;
        this.systemName = systemName;
    }
}
