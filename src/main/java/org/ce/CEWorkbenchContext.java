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

    // =========================================================================
    // High-level API (Shared between CLI and GUI)
    // =========================================================================

    /**
     * Runs Type-1a: Cluster identification and saves the results.
     */
    public AllClusterData runType1a(String clusterId, ClusterIdentificationRequest request, Consumer<String> progressSink) throws IOException {
        AllClusterData data = ClusterIdentificationWorkflow.identify(request, progressSink);
        clusterStore.save(clusterId, data);
        return data;
    }

    /**
     * Runs Type-1b: Scaffolds a Hamiltonian from existing cluster data.
     */
    public void runType1b(String hamiltonianId, String elements, String structure, String model) throws Exception {
        cecWorkflow.scaffoldFromClusterData(hamiltonianId, elements, structure, model);
    }
}
