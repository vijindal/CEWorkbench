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
import org.ce.domain.engine.ProgressEvent;

import java.util.function.Consumer;

public class ThermodynamicRequest {

    public final String clusterId;
    public final String hamiltonianId;
    public final double temperature;
    public final double[] composition;
    public final String engineType;
    /** Optional sink for real-time text progress messages. May be null. */
    public final Consumer<String> progressSink;
    /** Optional sink for structured engine events (chart data). May be null. */
    public final Consumer<ProgressEvent> eventSink;
    /** MCS lattice size (default 4) */
    public final int mcsL;
    /** MCS equilibration sweeps (default 1000) */
    public final int mcsNEquil;
    /** MCS averaging sweeps (default 2000) */
    public final int mcsNAvg;

    public ThermodynamicRequest(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double[] composition,
            String engineType,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg) {

        this.clusterId     = clusterId;
        this.hamiltonianId = hamiltonianId;
        this.temperature   = temperature;
        this.composition   = composition;
        this.engineType    = engineType;
        this.progressSink  = progressSink;
        this.eventSink     = eventSink;
        this.mcsL          = mcsL;
        this.mcsNEquil     = mcsNEquil;
        this.mcsNAvg       = mcsNAvg;
    }

    public ThermodynamicRequest(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double[] composition,
            String engineType,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink) {
        this(clusterId, hamiltonianId, temperature, composition, engineType, progressSink, eventSink, 4, 1000, 2000);
    }

    public ThermodynamicRequest(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double[] composition,
            String engineType,
            Consumer<String> progressSink) {
        this(clusterId, hamiltonianId, temperature, composition, engineType, progressSink, null);
    }

    public ThermodynamicRequest(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double[] composition,
            String engineType) {
        this(clusterId, hamiltonianId, temperature, composition, engineType, null, null);
    }
}
