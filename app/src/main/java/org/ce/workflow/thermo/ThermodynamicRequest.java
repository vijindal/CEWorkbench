package org.ce.workflow.thermo;

/**
 * Request describing a thermodynamic calculation.
 */
public class ThermodynamicRequest {

    public final String systemId;

    public final double temperature;

    public final double[] composition;

    public final String engineType;

    public ThermodynamicRequest(
            String systemId,
            double temperature,
            double[] composition,
            String engineType) {

        this.systemId = systemId;
        this.temperature = temperature;
        this.composition = composition;
        this.engineType = engineType;
    }
}
