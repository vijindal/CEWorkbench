package org.ce.domain.engine.mcs;

/** Real-time update event for MCS simulations. */
public class MCSUpdate {

    public enum Phase { EQUILIBRATION, AVERAGING }

    private final int    step;
    private final double E_total;
    private final double deltaE;
    private final double sigmaDE;
    private final double meanDE;
    private final Phase  phase;
    private final double acceptanceRate;
    private final long   timestampMs;
    private final long   elapsedMs;

    public MCSUpdate(int step, double E_total, double deltaE, double sigmaDE, double meanDE,
                     Phase phase, double acceptanceRate, long timestampMs, long elapsedMs) {
        this.step          = step;
        this.E_total       = E_total;
        this.deltaE        = deltaE;
        this.sigmaDE       = sigmaDE;
        this.meanDE        = meanDE;
        this.phase         = phase;
        this.acceptanceRate = acceptanceRate;
        this.timestampMs   = timestampMs;
        this.elapsedMs     = elapsedMs;
    }

    public int    getStep()           { return step; }
    public double getE_total()        { return E_total; }
    public double getDeltaE()         { return deltaE; }
    public double getSigmaDE()        { return sigmaDE; }
    public double getMeanDE()         { return meanDE; }
    public Phase  getPhase()          { return phase; }
    public double getAcceptanceRate() { return acceptanceRate; }
    public long   getTimestampMs()    { return timestampMs; }
    public long   getElapsedMs()      { return elapsedMs; }

    @Override
    public String toString() {
        return String.format("MCSUpdate{step=%d, E=%.4f, phase=%s, acceptance=%.1f%%, elapsed=%dms}",
                step, E_total, phase, acceptanceRate * 100, elapsedMs);
    }
}
