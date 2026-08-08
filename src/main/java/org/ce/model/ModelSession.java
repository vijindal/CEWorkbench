package org.ce.model;

import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace.SystemId;
import org.ce.model.storage.DataStore;

import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Immutable session object holding all pre-computed state for a system.
 *
 * <p>Created once when the user confirms system identity (elements + structure + model).
 * Invalidated and rebuilt when any of those three change. Shared read-only across
 * all GUI panels and the Calculation Layer.</p>
 *
 * <p>Use {@link Builder} to construct an instance. Building performs cluster
 * identification, Hamiltonian loading, and CVCF basis resolution exactly once.
 * All calculations within the session reuse this pre-computed state, eliminating
 * redundant disk I/O and computation per calculation call.</p>
 */
public final class ModelSession {

    /**
     * Value of {@link #resolvedHamiltonianId} when the ECIs were supplied by the
     * caller rather than loaded from the Hamiltonian store.
     */
    public static final String INLINE_HAMILTONIAN_ID = "<inline>";

    /** Engine type for a calculation session. CVM always uses the CVCF basis. */
    public enum EngineConfig {
        CVM, MCS;

        public boolean isCvm() { return this == CVM; }
        public boolean isMcs() { return this == MCS; }
    }

    // =========================================================================
    // Fields
    // =========================================================================

    /** The three-part system identity that produced this session. */
    public final SystemId systemId;

    /** Loaded and validated Hamiltonian (ECI parameters). */
    public final CECEntry cecEntry;

    /**
     * The Hamiltonian ID actually used (strictly forced to end in {@code _CVCF}).
     */
    public final String resolvedHamiltonianId;

    /** Engine type decided at session-creation time. */
    public final EngineConfig engineConfig;

    // =========================================================================
    // Constructor (package-private — use Builder)
    // =========================================================================

    ModelSession(
            SystemId systemId,
            CECEntry cecEntry,
            String resolvedHamiltonianId,
            EngineConfig engineConfig) {

        this.systemId              = systemId;
        this.cecEntry              = cecEntry;
        this.resolvedHamiltonianId = resolvedHamiltonianId;
        this.engineConfig          = engineConfig;
    }

    // =========================================================================
    // Derived properties
    // =========================================================================

    /** Number of chemical components derived from {@code systemId.elements}. */
    public int numComponents() {
        return systemId.numComponents();
    }

    /** Canonical element order for this system; see {@link SystemId#elementList()}. */
    public List<String> elements() {
        return systemId.elementList();
    }

    /** Short human-readable label, e.g. {@code "Nb-Ti / BCC_A2 / T [CVM]"}. */
    public String label() {
        return systemId.elements() + " / " + systemId.structure() + " / " + systemId.model()
                + " [" + engineConfig + "]";
    }

    // =========================================================================
    // Builder
    // =========================================================================

    /**
     * Builds a {@link ModelSession} from a {@link SystemId} and {@link EngineConfig}.
     *
     * <p>Consolidates everything that previously happened on every calculation call:</p>
     * <ol>
     *   <li>Hamiltonian ID resolution (strictly enforced _CVCF suffix)</li>
     *   <li>Hamiltonian loading ({@link CECEntry} from HamiltonianStore)</li>
     * </ol>
     */
    public static class Builder {

        private static final Logger LOG = Logger.getLogger(Builder.class.getName());

        private final DataStore.HamiltonianStore hamiltonianStore;

        public Builder(DataStore.HamiltonianStore hamiltonianStore) {
            this.hamiltonianStore = hamiltonianStore;
        }

        /**
         * Builds a {@link ModelSession} synchronously.
         */
        public ModelSession build(
                SystemId systemId,
                EngineConfig engineConfig,
                Consumer<String> progressSink) throws Exception {

            emit(progressSink, "\n[Workflow] Stage 1: Loading Specifications...");

            // ── Stage 1b: Resolve Hamiltonian ID ──────────────────────────────
            String baseId = systemId.hamiltonianId();
            String resolvedId = baseId.endsWith("_CVCF") ? baseId : baseId + "_CVCF";
            
            emit(progressSink, "  [Session] Stage 1b: Resolving Hamiltonian ID...");
            if (!hamiltonianStore.exists(resolvedId)) {
                throw new IllegalStateException("Required CVCF Hamiltonian file not found: '" + resolvedId + 
                    "'. Please ensure the Hamiltonian is stored with the _CVCF suffix.");
            }
            emit(progressSink, "  [Session] ✓ Resolved to '" + resolvedId + "'");

            // ── Stage 1c: Load and validate Hamiltonian ────────────────────────
            emit(progressSink, "  [Session] Stage 1c: Loading Hamiltonian...");
            CECEntry cecEntry = hamiltonianStore.load(resolvedId);

            return assemble(systemId, engineConfig, cecEntry, resolvedId, progressSink);
        }

        /**
         * Builds a {@link ModelSession} from a caller-supplied {@link CECEntry},
         * bypassing the Hamiltonian store entirely.
         *
         * <p>Used by the JSON API when an external caller (e.g. a fitting tool)
         * supplies its own ECIs rather than referencing the stored CEC database.
         * The ECIs must already be in the CVCF basis — no transformation is applied.</p>
         *
         * <p>Note: sessions built this way must <em>not</em> be routed through
         * {@code CalculationService.getOrBuildSession}, whose cache is keyed only on
         * (systemId, engineConfig) and would collide with a stored-Hamiltonian session
         * for the same system.</p>
         */
        public ModelSession build(
                SystemId systemId,
                EngineConfig engineConfig,
                CECEntry providedEntry,
                Consumer<String> progressSink) throws Exception {

            if (providedEntry == null)
                throw new IllegalArgumentException("providedEntry must not be null");

            emit(progressSink, "\n[Workflow] Stage 1: Loading Specifications...");
            emit(progressSink, "  [Session] Using caller-supplied Hamiltonian (store bypassed)");

            return assemble(systemId, engineConfig, providedEntry, INLINE_HAMILTONIAN_ID, progressSink);
        }

        /** Shared validation + construction for both the stored and inline paths. */
        private ModelSession assemble(
                SystemId systemId,
                EngineConfig engineConfig,
                CECEntry cecEntry,
                String resolvedId,
                Consumer<String> progressSink) {

            // Validate CEC: term count must be > 0
            int termCount = cecEntry.cecTerms == null ? 0 : cecEntry.cecTerms.length;
            if (termCount <= 0) {
                throw new IllegalStateException("CEC term count (" + termCount
                        + ") is invalid (must be > 0)");
            }
            for (CECEntry.CECTerm term : cecEntry.cecTerms) {
                term.validate();
            }

            // Report the actual term count, not the declared `ncf` — an inline entry
            // may omit `ncf` (an int, defaulting to 0) and would otherwise report "0 terms".
            emit(progressSink, "  [Session] ✓ Hamiltonian loaded (" + termCount + " terms)");

            emit(progressSink, "  [Session] ✓ Session ready — " + systemId.elements()
                    + " / " + systemId.structure() + " / " + systemId.model());

            return new ModelSession(
                    systemId,
                    cecEntry,
                    resolvedId,
                    engineConfig);
        }

        private static void emit(Consumer<String> sink, String msg) {
            if (sink != null) sink.accept(msg);
        }
    }
}
