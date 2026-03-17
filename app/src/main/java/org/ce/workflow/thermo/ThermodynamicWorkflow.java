package org.ce.workflow.thermo;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.result.EquilibriumState;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.cec.CECEntry;
import org.ce.workflow.cec.CECManagementWorkflow;
import org.ce.workflow.thermo.ThermodynamicRequest;

/**
 * Workflow responsible for thermodynamic calculations.
 *
 * It loads the required scientific data and delegates
 * the computation to thermodynamic engines.
 */
public class ThermodynamicWorkflow {

    private final ClusterDataStore clusterStore;

    private final CECManagementWorkflow cecWorkflow;

    private final ThermodynamicEngine engine;

    public ThermodynamicWorkflow(
            ClusterDataStore clusterStore,
            CECManagementWorkflow cecWorkflow,
            ThermodynamicEngine engine) {

        this.clusterStore = clusterStore;
        this.cecWorkflow = cecWorkflow;
        this.engine = engine;
    }

    /**
     * Loads the scientific data required for thermodynamic calculations.
     */
    public ThermodynamicData loadThermodynamicData(String systemId) throws Exception {

        AllClusterData clusterData = clusterStore.load(systemId);

        CECEntry cec = cecWorkflow.loadAndValidateCEC(systemId);

        // System name: can be derived from CECEntry (e.g., elements-structurePhase)
        // For now, use systemId. Future: systemName = cec.elements + "-" + cec.structurePhase
        String systemName = systemId;

        return new ThermodynamicData(clusterData, cec, systemId, systemName);
    }

    /**
     * Prepares the data required for a thermodynamic calculation.
     */
    public ThermodynamicData prepareCalculation(ThermodynamicRequest request) throws Exception {

        return loadThermodynamicData(request.systemId);

    }

    /**
     * Runs a thermodynamic calculation.
     */
    public EquilibriumState runCalculation(ThermodynamicRequest request) throws Exception {

        ThermodynamicData data = loadThermodynamicData(request.systemId);

        ThermodynamicInput input = new ThermodynamicInput(
                data.clusterData,
                data.cec,
                request.temperature,
                request.composition,
                data.systemId,
                data.systemName
        );

        return engine.compute(input);
    }

}
