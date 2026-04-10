package org.ce.model;

/**
 * Immutable value object capturing the engine configuration for a calculation session.
 *
 * <p>Resolved once at session-creation time so the calculation layer never
 * needs to re-derive which engine or basis to use. CVM always uses the CVCF basis.</p>
 */
public final class EngineConfig {

    /** Engine type: {@code "CVM"} or {@code "MCS"}. */
    public final String engineType;

    public EngineConfig(String engineType) {
        this.engineType = engineType != null ? engineType.toUpperCase() : "CVM";
    }

    /** Convenience factory for CVM (always uses CVCF basis). */
    public static EngineConfig cvm() {
        return new EngineConfig("CVM");
    }

    /** Convenience factory for MCS. */
    public static EngineConfig mcs() {
        return new EngineConfig("MCS");
    }

    public boolean isCvm() { return "CVM".equals(engineType); }
    public boolean isMcs() { return "MCS".equals(engineType); }

    @Override
    public String toString() {
        return engineType;
    }
}
