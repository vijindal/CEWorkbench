package org.ce.calculation.engine;

import org.ce.model.cluster.AllClusterData;
import org.ce.model.hamiltonian.CECEntry;

import java.util.function.Consumer;

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

    /**
     * Full constructor — {@code clusterData} must be non-null (pre-built by
     * {@link org.ce.model.ModelSession}).
     */
    public ThermodynamicInput(
            AllClusterData clusterData,
            CECEntry cec,
            double temperature,
            double[] composition,
            String systemId,
            String systemName,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg) {

        java.util.Objects.requireNonNull(clusterData,
                "clusterData must be pre-built by ModelSession — not null");
        this.clusterData  = clusterData;
        this.cec          = cec;
        this.temperature  = temperature;
        this.composition  = composition;
        this.systemId     = systemId;
        this.systemName   = systemName;
        this.progressSink = progressSink;
        this.eventSink    = eventSink;
        this.mcsL         = mcsL;
        this.mcsNEquil    = mcsNEquil;
        this.mcsNAvg      = mcsNAvg;
    }

}

