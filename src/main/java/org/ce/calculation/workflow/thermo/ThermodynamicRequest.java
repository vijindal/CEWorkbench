package org.ce.calculation.workflow.thermo;

import org.ce.calculation.engine.ProgressEvent;

import java.util.function.Consumer;

/**
 * Request describing a thermodynamic calculation.
 *
 * <p>Contains only calculation parameters (T, composition, MCS params, progress sinks).
 * System identity, cluster data, Hamiltonian, and engine type are all held by the
 * {@link org.ce.model.ModelSession} passed separately to the workflow.</p>
 */
public class ThermodynamicRequest {

    public final double temperature;
    public final double[] composition;
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
            double temperature,
            double[] composition,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg) {

        this.temperature  = temperature;
        this.composition  = composition;
        this.progressSink = progressSink;
        this.eventSink    = eventSink;
        this.mcsL         = mcsL;
        this.mcsNEquil    = mcsNEquil;
        this.mcsNAvg      = mcsNAvg;
    }

    public ThermodynamicRequest(
            double temperature,
            double[] composition,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink) {
        this(temperature, composition, progressSink, eventSink, 4, 1000, 2000);
    }

    public ThermodynamicRequest(
            double temperature,
            double[] composition,
            Consumer<String> progressSink) {
        this(temperature, composition, progressSink, null, 4, 1000, 2000);
    }

    public ThermodynamicRequest(
            double temperature,
            double[] composition) {
        this(temperature, composition, null, null, 4, 1000, 2000);
    }
}
