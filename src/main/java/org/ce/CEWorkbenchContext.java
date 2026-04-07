package org.ce;

import org.ce.domain.engine.cvm.CVMEngine;
import org.ce.domain.engine.mcs.MCSEngine;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.Workspace;
import org.ce.workflow.CalculationService;
import org.ce.workflow.ClusterIdentificationWorkflow;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.domain.cluster.AllClusterData;
import org.ce.workflow.cec.CECManagementWorkflow;
import org.ce.workflow.thermo.ThermodynamicWorkflow;

import java.io.IOException;
import java.util.function.Consumer;

/**
 * Unified context for CE Workbench, providing shared infrastructure and services.
 * 
 * <p>This class centralizes the initialization of all core services (storage, 
 * engines, workflows) using a shared {@link Workspace}. Both the GUI and CLI
 * use this context to ensure consistent behavior and API usage.</p>
 */
public class CEWorkbenchContext {

    private final Workspace workspace;
    private final ClusterDataStore clusterStore;
    private final HamiltonianStore hamiltonianStore;
    private final CECManagementWorkflow cecWorkflow;
    private final CalculationService calculationService;
    private final CVMEngine cvmEngine;
    private final MCSEngine mcsEngine;
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
        this.clusterStore = new ClusterDataStore(workspace);
        this.hamiltonianStore = new HamiltonianStore(workspace);
        this.cecWorkflow = new CECManagementWorkflow(hamiltonianStore, clusterStore);
        this.cvmEngine = new CVMEngine();
        this.mcsEngine = new MCSEngine();
        
        ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow(
                clusterStore, cecWorkflow, cvmEngine, mcsEngine
        );
        this.calculationService = new CalculationService(thermoWorkflow);
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    public Workspace getWorkspace() {
        return workspace;
    }

    public ClusterDataStore getClusterStore() {
        return clusterStore;
    }

    public HamiltonianStore getHamiltonianStore() {
        return hamiltonianStore;
    }

    public CECManagementWorkflow getCecWorkflow() {
        return cecWorkflow;
    }

    public CalculationService getCalculationService() {
        return calculationService;
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
        return ClusterIdentificationWorkflow.identify(request, progressSink);
    }

    /**
     * Saves AllClusterData to the storage.
     */
    public void saveClusterData(String clusterId, AllClusterData data) throws IOException {
        clusterStore.save(clusterId, data);
    }
}

