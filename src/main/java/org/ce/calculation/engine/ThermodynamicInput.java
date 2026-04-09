package org.ce.domain.engine;

import static org.ce.domain.cluster.AllClusterData.ClusterData;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.hamiltonian.CECEntry;

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

    /** Full constructor for specialists (e.g. MCS or pre-identified CVM). */
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

    /** Primary constructor for thermodynamic calculations with progress updates. */
    public ThermodynamicInput(
            CECEntry cec, double temperature, double[] composition,
            String systemId, String systemName, Consumer<String> progressSink) {
        this(null, cec, temperature, composition, systemId, systemName, progressSink, null, 4, 1000, 2000);
    }

    /** Minimal constructor for thermodynamic calculations. */
    public ThermodynamicInput(
            CECEntry cec, double temperature, double[] composition, String systemId, String systemName) {
        this(null, cec, temperature, composition, systemId, systemName, null, null, 4, 1000, 2000);
    }
}
