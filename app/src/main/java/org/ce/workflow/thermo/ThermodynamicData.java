package org.ce.workflow.thermo;

import org.ce.domain.cluster.AllClusterData;
import org.ce.storage.cec.CECEntry;

/**
 * Holds the scientific data required for thermodynamic calculations.
 */
public class ThermodynamicData {

    public final AllClusterData clusterData;

    public final CECEntry cec;

    public ThermodynamicData(AllClusterData clusterData, CECEntry cec) {
        this.clusterData = clusterData;
        this.cec = cec;
    }
}
