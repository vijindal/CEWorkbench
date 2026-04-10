package org.ce.calculation.engine;

/**
 * Structured progress events emitted by thermodynamic engines during computation.
 *
 * <p>Consumed by the UI layer to drive live chart updates. Text log messages continue
 * to flow through the separate {@code Consumer<String>} progressSink.</p>
 */
public abstract class ProgressEvent {

    /**
     * Signals the start of a new calculation — listeners should clear prior chart data.
     */
    public static final class EngineStart extends ProgressEvent {
        /** "MCS" or "CVM". */
        public final String engineType;
        /** Total sweep count for MCS; 0 for CVM (unknown upfront). */
        public final int totalSteps;

        public EngineStart(String engineType, int totalSteps) {
            this.engineType = engineType;
            this.totalSteps = totalSteps;
        }
    }

    /**
     * One MC sweep result — emitted on every sweep during MCS.
     *
     * <p>{@code cfs} contains the running-average correlation functions computed by
     * {@code MCSampler.meanCFs()} during the AVERAGING phase; it is {@code null}
     * during EQUILIBRATION (sampler has not yet accumulated data).</p>
     */
    public static final class McSweep extends ProgressEvent {
        public final int     step;
        public final int     totalSteps;
        public final double  energyPerSite;
        public final double  acceptanceRate;
        /** True during equilibration; false during averaging. */
        public final boolean equilibration;
        /** Running-average CFs; null during equilibration. */
        public final double[] cfs;

        public McSweep(int step, int totalSteps, double energyPerSite,
                       double acceptanceRate, boolean equilibration, double[] cfs) {
            this.step           = step;
            this.totalSteps     = totalSteps;
            this.energyPerSite  = energyPerSite;
            this.acceptanceRate = acceptanceRate;
            this.equilibration  = equilibration;
            this.cfs            = cfs != null ? cfs.clone() : null;
        }
    }

    /**
     * One Newton-Raphson iteration — emitted by CVM post-hoc from the stored iteration
     * trace after convergence.
     *
     * <p>Note: on solver failure {@code CVMPhaseModel.create()} throws before the trace
     * is readable, so chart events are only emitted on successful convergence.</p>
     */
    public static final class CvmIteration extends ProgressEvent {
        public final int      iteration;
        public final double   gibbsEnergy;
        public final double   gradientNorm;
        public final double   enthalpy;     // H in J/mol
        public final double   entropy;      // S in J/(mol·K)
        /** Non-point correlation functions at this iteration (length = ncf); may be null. */
        public final double[] cfs;

        public CvmIteration(int iteration, double gibbsEnergy, double gradientNorm,
                           double enthalpy, double entropy, double[] cfs) {
            this.iteration    = iteration;
            this.gibbsEnergy  = gibbsEnergy;
            this.gradientNorm = gradientNorm;
            this.enthalpy     = enthalpy;
            this.entropy      = entropy;
            this.cfs          = cfs != null ? cfs.clone() : null;
        }
    }
}
