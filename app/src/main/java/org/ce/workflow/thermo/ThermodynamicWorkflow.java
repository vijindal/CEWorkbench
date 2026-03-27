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

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Workflow responsible for thermodynamic calculations.
 *
 * It loads the required scientific data and delegates
 * the computation to thermodynamic engines.
 */
public class ThermodynamicWorkflow {

    private static final Logger LOG = Logger.getLogger(ThermodynamicWorkflow.class.getName());

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

        LOG.fine("STAGE 3a: Load cluster data (Stages 1-3 topology)...");
        AllClusterData clusterData = clusterStore.load(clusterId);
        LOG.fine("  ✓ Loaded " + clusterId);

        LOG.fine("STAGE 3b: Load Hamiltonian (ECI coefficients)...");
        CECEntry cec = cecWorkflow.loadAndValidateCEC(clusterId, hamiltonianId);
        LOG.fine("  ✓ Loaded " + hamiltonianId);

        // Extract and log CEC data
        if (cec != null) {
            LOG.fine("    elements: " + cec.elements);
            LOG.fine("    structurePhase: " + cec.structurePhase);
            if (cec.cecTerms != null && cec.cecTerms.length > 0) {
                LOG.fine("    ECI terms: " + cec.cecTerms.length + " (format: a + b*T)");
                // Show first few terms as example
                for (int i = 0; i < Math.min(3, cec.cecTerms.length); i++) {
                    LOG.finer("      [" + i + "] " + cec.cecTerms[i].name
                            + ": a=" + cec.cecTerms[i].a + ", b=" + cec.cecTerms[i].b);
                }
                if (cec.cecTerms.length > 3) {
                    LOG.finer("      ... (" + (cec.cecTerms.length - 3) + " more terms)");
                }
            } else {
                LOG.fine("    ECI terms: none");
            }
        }

        String systemName = cec.elements + "_" + cec.structurePhase;
        LOG.fine("STAGE 3c: Create ThermodynamicData bundle");
        LOG.fine("  systemName: " + systemName);

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

        LOG.info("ThermodynamicWorkflow.runCalculation — ENTER");
        LOG.info("  clusterId: " + request.clusterId);
        LOG.info("  hamiltonianId: " + request.hamiltonianId);
        LOG.info("  temperature: " + request.temperature + " K");
        LOG.info("  composition: " + Arrays.toString(request.composition));
        LOG.info("  engineType: " + request.engineType);

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
        ThermodynamicResult result = ThermodynamicResult.from(state);

        LOG.info("ThermodynamicWorkflow.runCalculation — EXIT: G=" + String.format("%.4e", state.getFreeEnergy())
                + " J/mol");
        return result;
    }

}
