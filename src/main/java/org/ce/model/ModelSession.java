package org.ce.model;

import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace.SystemId;
import org.ce.model.storage.DataStore;

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

    /** Engine type for a calculation session. CVM always uses the CVCF basis. */
    public enum EngineConfig {
        CVM, MCS;

        /** Convenience factory — kept for source compatibility. */
        public static EngineConfig cvm() { return CVM; }

        /** Convenience factory — kept for source compatibility. */
        public static EngineConfig mcs() { return MCS; }

        public boolean isCvm() { return this == CVM; }
        public boolean isMcs() { return this == MCS; }
    }

    // =========================================================================
    // Fields
    // =========================================================================

    /** The three-part system identity that produced this session. */
    public final SystemId systemId;

    /**
     * Full cluster identification result (clusters, symmetry matrices, C-matrix,
     * CF list). Pre-computed by {@link Builder#build}.
     */
    public final PipelineResult clusterData;

    /** Loaded and validated Hamiltonian (ECI parameters). */
    public final CECEntry cecEntry;

    /**
     * The Hamiltonian ID actually used (may have a {@code _CVCF} suffix for CVM mode,
     * or equal {@link SystemId#hamiltonianId()} for MCS).
     */
    public final String resolvedHamiltonianId;

    /**
     * CVCF basis for this (structure, model, numComponents) combination.
     * Pre-resolved from {@link CvCfBasis.Registry} at session-creation time.
     */
    public final CvCfBasis cvcfBasis;

    /** Engine type decided at session-creation time. */
    public final EngineConfig engineConfig;

    // =========================================================================
    // Constructor (package-private — use Builder)
    // =========================================================================

    ModelSession(
            SystemId systemId,
            PipelineResult clusterData,
            CECEntry cecEntry,
            String resolvedHamiltonianId,
            CvCfBasis cvcfBasis,
            EngineConfig engineConfig) {

        this.systemId              = systemId;
        this.clusterData           = clusterData;
        this.cecEntry              = cecEntry;
        this.resolvedHamiltonianId = resolvedHamiltonianId;
        this.cvcfBasis             = cvcfBasis;
        this.engineConfig          = engineConfig;
    }

    // =========================================================================
    // Derived properties
    // =========================================================================

    /** Number of chemical components derived from {@code systemId.elements}. */
    public int numComponents() {
        return systemId.elements.split("-").length;
    }

    /** Short human-readable label, e.g. {@code "Nb-Ti / BCC_A2 / T [CVM]"}. */
    public String label() {
        return systemId.elements + " / " + systemId.structure + " / " + systemId.model
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
     *   <li>Cluster identification ({@link ClusterIdentificationWorkflow})</li>
     *   <li>Hamiltonian ID resolution (CVCF preference for CVM mode)</li>
     *   <li>Hamiltonian loading ({@link CECManagementWorkflow#loadAndValidateCEC})</li>
     *   <li>CVCF basis lookup ({@link CvCfBasis.Registry})</li>
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
         *
         * <p>Callers running on the EDT <em>must</em> wrap this in a
         * {@link javax.swing.SwingWorker} to avoid blocking the UI thread.</p>
         *
         * @param systemId     three-part system identity (elements, structure, model)
         * @param engineConfig engine type for this session
         * @param progressSink optional sink for progress text; may be {@code null}
         * @return fully constructed, immutable {@link ModelSession}
         * @throws Exception if cluster identification or Hamiltonian loading fails
         */
        public ModelSession build(
                SystemId systemId,
                EngineConfig engineConfig,
                Consumer<String> progressSink) throws Exception {

            emit(progressSink, "Building session for " + systemId.elements
                    + " / " + systemId.structure + " / " + systemId.model
                    + " [" + engineConfig + "] ...");

            // ── Stage 1: Cluster identification ──────────────────────────────
            emit(progressSink, "  [Session] Stage 1: Cluster identification...");
            int numComponents = systemId.elements.split("-").length;

            // Derive file names from structure/model (same convention as
            // ClusterIdentificationRequest(ThermodynamicInput)):
            //   clus/<structure>-<model>.txt  and  <structure>-SG
            String clusterFile = "clus/" + systemId.structure + "-" + systemId.model + ".txt";
            String symGroup    = systemId.structure + "-SG";

            ClusterIdentificationRequest clusterReq = ClusterIdentificationRequest.builder()
                    .structurePhase(systemId.structure)
                    .model(systemId.model)
                    .numComponents(numComponents)
                    .disorderedClusterFile(clusterFile)
                    .orderedClusterFile(clusterFile)
                    .disorderedSymmetryGroup(symGroup)
                    .orderedSymmetryGroup(symGroup)
                    .build();

            PipelineResult clusterData = ClusterCFIdentificationPipeline.runFullWorkflow(clusterReq, progressSink);
            clusterData.printSummary(System.out::println);
            emit(progressSink, "  [Session] ✓ Cluster identification complete");

            // ── Stage 2: Resolve Hamiltonian ID ──────────────────────────────
            String baseHamiltonianId = systemId.hamiltonianId();
            String resolvedHamiltonianId = resolveHamiltonianId(
                    baseHamiltonianId, engineConfig, progressSink);

            // ── Stage 3: Load and validate Hamiltonian ────────────────────────
            emit(progressSink, "  [Session] Stage 3: Loading Hamiltonian '" + resolvedHamiltonianId + "'...");
            CECEntry cecEntry = hamiltonianStore.load(resolvedHamiltonianId);

            // Validate CEC: term count must be > 0 and <= basis.numNonPointCfs
            int ncf = CvCfBasis.getNumNonPointCfs(systemId.structure, systemId.model, numComponents);
            int termCount = cecEntry.cecTerms == null ? 0 : cecEntry.cecTerms.length;
            if (termCount <= 0 || termCount > ncf) {
                throw new IllegalStateException("CEC term count (" + termCount
                        + ") is invalid (must be > 0 and <= " + ncf + ")");
            }
            for (CECEntry.CECTerm term : cecEntry.cecTerms) {
                term.validate();
            }

            emit(progressSink, "  [Session] ✓ Hamiltonian loaded (" + cecEntry.ncf + " terms)");

            // ── Stage 4: CVCF basis explicit generation ─────────────
            emit(progressSink, "  [Session] Stage 4: Generating CVCF basis...");
            CvCfBasis cvcfBasis = CvCfBasis.generate(
                    systemId.structure,
                    clusterData,
                    clusterData.getMatrixData(),
                    systemId.model,
                    progressSink
            );
            emit(progressSink, "  [Session] ✓ Basis resolved: "
                    + cvcfBasis.numNonPointCfs + " non-point CFs, "
                    + cvcfBasis.numComponents + " point variables");
            emit(progressSink, "  [Session] ✓ Session ready — " + systemId.elements
                    + " / " + systemId.structure + " / " + systemId.model);

            return new ModelSession(
                    systemId,
                    clusterData,
                    cecEntry,
                    resolvedHamiltonianId,
                    cvcfBasis,
                    engineConfig);
        }

        /**
         * Resolves the Hamiltonian ID to use.
         *
         * <p>For CVM: prefers the {@code _CVCF} suffixed Hamiltonian when available.
         * For MCS: uses the base ID directly (logic moved from
         * {@code ThermodynamicWorkflow.resolveHamiltonianIdForEngine}).</p>
         */
        private String resolveHamiltonianId(
                String baseId,
                EngineConfig engineConfig,
                Consumer<String> progressSink) {

            if (!engineConfig.isCvm()) {
                return baseId;
            }
            if (baseId == null || baseId.isBlank() || baseId.endsWith("_CVCF")) {
                return baseId;
            }

            // Preferred: exact _CVCF suffix
            String preferredId = baseId + "_CVCF";
            if (hamiltonianStore.exists(preferredId)) {
                emit(progressSink, "  [Session] CVM mode: using CVCF Hamiltonian '"
                        + preferredId + "'");
                LOG.info("CVM mode: using CVCF Hamiltonian '" + preferredId + "'");
                return preferredId;
            }

            // Legacy fallback: strip last segment, append _CVCF
            int lastUnderscore = baseId.lastIndexOf('_');
            if (lastUnderscore > 0) {
                String legacyId = baseId.substring(0, lastUnderscore) + "_CVCF";
                if (hamiltonianStore.exists(legacyId)) {
                    emit(progressSink, "  [Session] CVM mode: using CVCF Hamiltonian (legacy) '"
                            + legacyId + "'");
                    LOG.info("CVM mode: using CVCF Hamiltonian (legacy) '" + legacyId + "'");
                    return legacyId;
                }
            }

            emit(progressSink, "  [Session] CVM mode: CVCF Hamiltonian not found, "
                    + "falling back to '" + baseId + "'");
            LOG.warning("CVM mode: CVCF Hamiltonian not found, falling back to '" + baseId + "'");
            return baseId;
        }

        private static void emit(Consumer<String> sink, String msg) {
            if (sink != null) sink.accept(msg);
        }
    }
}
