package org.ce.workflow.thermo;

/**
 * Request describing a thermodynamic calculation.
 *
 * <p>Two IDs are required because cluster geometry and the Hamiltonian
 * are stored separately with different naming conventions:</p>
 * <ul>
 *   <li>{@code clusterId} -- element-agnostic, e.g. {@code BCC_A2_T_bin}</li>
 *   <li>{@code hamiltonianId} -- element-specific, e.g. {@code Nb-Ti_BCC_A2_T}</li>
 * </ul>
 */
import java.util.function.Consumer;

public class ThermodynamicRequest {

    public final String clusterId;

    public final String hamiltonianId;

    public final double temperature;

    public final double[] composition;

    public final String engineType;

    /** Optional sink for real-time progress messages. May be null. */
    public final Consumer<String> progressSink;

    public ThermodynamicRequest(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double[] composition,
            String engineType,
            Consumer<String> progressSink) {

        this.clusterId     = clusterId;
        this.hamiltonianId = hamiltonianId;
        this.temperature   = temperature;
        this.composition   = composition;
        this.engineType    = engineType;
        this.progressSink  = progressSink;
    }

    public ThermodynamicRequest(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double[] composition,
            String engineType) {
        this(clusterId, hamiltonianId, temperature, composition, engineType, null);
    }
}
