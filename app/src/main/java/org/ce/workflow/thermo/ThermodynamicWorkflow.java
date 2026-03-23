package org.ce.workflow.thermo;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.result.EquilibriumState;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.ClusterDataStore;
import org.ce.domain.hamiltonian.CECEntry;
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

    private final ThermodynamicEngine cvmEngine;

    private final ThermodynamicEngine mcsEngine;

    public ThermodynamicWorkflow(
            ClusterDataStore clusterStore,
            CECManagementWorkflow cecWorkflow,
            ThermodynamicEngine cvmEngine,
            ThermodynamicEngine mcsEngine) {

        this.clusterStore = clusterStore;
        this.cecWorkflow = cecWorkflow;
        this.cvmEngine = cvmEngine;
        this.mcsEngine = mcsEngine;
    }

    /**
     * Loads the scientific data required for thermodynamic calculations.
     *
     * @param clusterId     element-agnostic cluster dataset ID, e.g. BCC_A2_T_bin
     * @param hamiltonianId element-specific Hamiltonian ID, e.g. Nb-Ti_BCC_A2_T
     */
    public ThermodynamicData loadThermodynamicData(
            String clusterId, String hamiltonianId) throws Exception {

        AllClusterData clusterData = clusterStore.load(clusterId);

        CECEntry cec = cecWorkflow.loadAndValidateCEC(clusterId, hamiltonianId);

        String systemName = cec.elements + "_" + cec.structurePhase;

        return new ThermodynamicData(clusterData, cec, clusterId, systemName);
    }

    /**
     * Prepares the data required for a thermodynamic calculation.
     */
    public ThermodynamicData prepareCalculation(ThermodynamicRequest request) throws Exception {

        return loadThermodynamicData(request.clusterId, request.hamiltonianId);
    }

    /**
     * Runs a thermodynamic calculation.
     */
    public ThermodynamicResult runCalculation(ThermodynamicRequest request) throws Exception {

        ThermodynamicData data = loadThermodynamicData(request.clusterId, request.hamiltonianId);

        ThermodynamicInput input = new ThermodynamicInput(
                data.clusterData,
                data.cec,
                request.temperature,
                request.composition,
                data.systemId,
                data.systemName,
                request.progressSink,
                request.eventSink,
                request.mcsL,
                request.mcsNEquil,
                request.mcsNAvg
        );

        ThermodynamicEngine engine;

        switch (request.engineType) {
            case "CVM":
                engine = cvmEngine;
                break;
            case "MCS":
                engine = mcsEngine;
                break;
            default:
                throw new IllegalArgumentException("Unknown engine: " + request.engineType);
        }

        EquilibriumState state = engine.compute(input);
        return ThermodynamicResult.from(state);
    }

}
