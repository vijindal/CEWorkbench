package org.ce.workflow.thermo;

import org.ce.domain.hamiltonian.CECEntry;

/**
 * Holds the scientific data required for thermodynamic calculations.
 */
public class ThermodynamicData {

    public final CECEntry cec;
    public final String systemId;
    public final String systemName;

    public ThermodynamicData(
            CECEntry cec,
            String systemId,
            String systemName) {

        this.cec = cec;
        this.systemId = systemId;
        this.systemName = systemName;
    }
}
