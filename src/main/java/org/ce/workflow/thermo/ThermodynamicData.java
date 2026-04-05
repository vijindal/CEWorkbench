package org.ce.workflow.thermo;

import static org.ce.domain.cluster.AllClusterData.ClusterData;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.hamiltonian.CECEntry;

/**
 * Holds the scientific data required for thermodynamic calculations.
 */
public class ThermodynamicData {

    public final AllClusterData clusterData;
    public final CECEntry cec;
    public final String systemId;
    public final String systemName;

    public ThermodynamicData(
            AllClusterData clusterData,
            CECEntry cec,
            String systemId,
            String systemName) {

        this.clusterData = clusterData;
        this.cec = cec;
        this.systemId = systemId;
        this.systemName = systemName;
    }
}
