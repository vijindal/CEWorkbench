package org.ce;

import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;
import org.ce.calculation.workflow.CalculationService;
import org.ce.calculation.workflow.CECManagementWorkflow;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.ModelSession;

import java.io.IOException;
import java.util.function.Consumer;

/**
 * Unified context for CE Workbench, providing shared infrastructure and services.
 *
 * <p>This class centralizes the initialization of all core services (storage,
 * workflows) using a shared {@link Workspace}. Both the GUI and CLI
 * use this context to ensure consistent behavior and API usage.</p>
 */
public class CEWorkbenchContext {

    private final Workspace workspace;
    private final DataStore.HamiltonianStore hamiltonianStore;
    private final CECManagementWorkflow cecWorkflow;
    private final CalculationService calculationService;
    private final ModelSession.Builder sessionBuilder;
    private Consumer<String> logSink = System.out::println;
    private Runnable logClearer = null;

    /**
     * Initializes a new context with the default workspace location.
     */
    public CEWorkbenchContext() {
        this(new Workspace());
    }

    /**
     * Initializes a new context with a specific workspace.
     */
    public CEWorkbenchContext(Workspace workspace) {
        this.workspace = workspace;
        this.hamiltonianStore = new DataStore.HamiltonianStore(workspace);
        this.cecWorkflow = new CECManagementWorkflow(hamiltonianStore);

        ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow();
        this.sessionBuilder = new ModelSession.Builder(hamiltonianStore);
        this.calculationService = new CalculationService(thermoWorkflow, sessionBuilder);
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    public Workspace getWorkspace() {
        return workspace;
    }

    public DataStore.HamiltonianStore getHamiltonianStore() {
        return hamiltonianStore;
    }

    public CECManagementWorkflow getCecWorkflow() {
        return cecWorkflow;
    }

    public CalculationService getCalculationService() {
        return calculationService;
    }

    public ModelSession.Builder getSessionBuilder() {
        return sessionBuilder;
    }

    /**
     * Sets the sink for log messages. This allows the GUI to route output
     * to its output panel while still allowing CLI output.
     */
    public void setLogSink(Consumer<String> logSink) {
        this.logSink = logSink;
    }

    /**
     * Unified print API for the entire application.
     * Routes message to the registered log sink.
     */
    public void log(String message) {
        if (logSink != null) {
            logSink.accept(message);
        }
    }

    /**
     * Sets the provider for clearing the log area.
     */
    public void setLogClearer(Runnable logClearer) {
        this.logClearer = logClearer;
    }

    /**
     * Clears the registered log area.
     */
    public void clearLog() {
        if (logClearer != null) {
            logClearer.run();
        }
    }

    // =========================================================================
    // High-level API (Shared between CLI and GUI)
    // =========================================================================

    /**
     * Runs Type-1a: Cluster identification.
     * Does NOT save results automatically.
     */
    public AllClusterData identifyClusters(ClusterIdentificationRequest request, Consumer<String> progressSink) throws IOException {
        return AllClusterData.identify(request, progressSink);
    }

}
