package org.ce.workflow.thermo;

import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.result.EquilibriumState;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.workflow.cec.CECManagementWorkflow;

import java.util.Arrays;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Workflow responsible for thermodynamic calculations.
 *
 * It loads the required scientific data and delegates
 * the computation to thermodynamic engines.
 */
public class ThermodynamicWorkflow {

    private static final Logger LOG = Logger.getLogger(ThermodynamicWorkflow.class.getName());

    private final CECManagementWorkflow cecWorkflow;
    private final ThermodynamicEngine cvmEngine;
    private final ThermodynamicEngine mcsEngine;

    public ThermodynamicWorkflow(
            CECManagementWorkflow cecWorkflow,
            ThermodynamicEngine cvmEngine,
            ThermodynamicEngine mcsEngine) {

        this.cecWorkflow = cecWorkflow;
        this.cvmEngine = cvmEngine;
        this.mcsEngine = mcsEngine;
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
        LOG.info("  cvmBasisMode: " + request.cvmBasisMode);

        String resolvedHamiltonianId = resolveHamiltonianIdForEngine(
                request.hamiltonianId, request.engineType, request.cvmBasisMode, request.progressSink);
        
        LOG.fine("STAGE 3: Load Hamiltonian (ECI coefficients)...");
        emit(request.progressSink, "STAGE 3: Load Hamiltonian (ECI coefficients)");
        CECEntry cec = cecWorkflow.loadAndValidateCEC(request.clusterId, resolvedHamiltonianId);
        emit(request.progressSink, "  Loaded Hamiltonian: " + resolvedHamiltonianId);
        LOG.fine("  ✓ Loaded " + resolvedHamiltonianId);

        String systemName = (cec != null) ? cec.elements + "_" + cec.structurePhase : "unknown";

        ThermodynamicInput input = new ThermodynamicInput(
                null, // clusterData resolved on-the-fly by engine
                cec,
                request.temperature,
                request.composition,
                resolvedHamiltonianId,
                systemName,
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

        if (request.progressSink != null) {
            emit(request.progressSink, "");
            emit(request.progressSink, "  RESULTS AT " + request.temperature + " K (composition: " + Arrays.toString(request.composition) + ")");
            emit(request.progressSink, "  " + "-".repeat(60));
            emit(request.progressSink, String.format("  Gibbs Energy (G): %15.6f J/mol", result.gibbsEnergy));
            emit(request.progressSink, String.format("  Enthalpy (H):     %15.6f J/mol", result.enthalpy));
            if (!Double.isNaN(result.entropy)) {
                emit(request.progressSink, String.format("  Entropy (S):      %15.6f J/(mol·K)", result.entropy));
            }
            if (result.optimizedCFs != null) {
                emit(request.progressSink, "  Equilibrium CFs (non-point):");
                for (int i = 0; i < result.optimizedCFs.length; i++) {
                    emit(request.progressSink, String.format("    CF[%d] = %15.10f", i, result.optimizedCFs[i]));
                }
            }
            emit(request.progressSink, "  " + "-".repeat(60));
            emit(request.progressSink, "  ✓ Calculation successful");
            emit(request.progressSink, "");
        }

        LOG.info("ThermodynamicWorkflow.runCalculation — EXIT: G=" + String.format("%.4e", state.getFreeEnergy())
                + " J/mol");
        return result;
    }

    /**
     * For CVM runs prefer CVCF Hamiltonians when available.
     * Examples:
     *   Nb-Ti_BCC_A2_T      -> Nb-Ti_BCC_A2_T_CVCF (preferred)
     *   Nb-Ti_BCC_A2_T      -> Nb-Ti_BCC_A2_CVCF   (legacy fallback)
     */
    private String resolveHamiltonianIdForEngine(
            String requestedHamiltonianId,
            String engineType,
            String cvmBasisMode,
            Consumer<String> progressSink) {
        if (!"CVM".equalsIgnoreCase(engineType)) {
            return requestedHamiltonianId;
        }
        if (requestedHamiltonianId == null || requestedHamiltonianId.isBlank()) {
            return requestedHamiltonianId;
        }
        if ("ORTHO".equalsIgnoreCase(cvmBasisMode)) {
            emit(progressSink, "CVM mode: using ORTHO Hamiltonian '" + requestedHamiltonianId + "'");
            LOG.info("CVM mode: using ORTHO Hamiltonian '" + requestedHamiltonianId + "'");
            return requestedHamiltonianId;
        }
        if (requestedHamiltonianId.endsWith("_CVCF")) {
            return requestedHamiltonianId;
        }

        String preferredId = requestedHamiltonianId + "_CVCF";
        if (cecWorkflow.hamiltonianExists(preferredId)) {
            emit(progressSink, "CVM mode: using CVCF Hamiltonian '" + preferredId
                    + "' instead of '" + requestedHamiltonianId + "'");
            LOG.info("CVM mode: using CVCF Hamiltonian '" + preferredId
                    + "' instead of '" + requestedHamiltonianId + "'");
            return preferredId;
        }

        int lastUnderscore = requestedHamiltonianId.lastIndexOf('_');
        if (lastUnderscore > 0) {
            String legacyId = requestedHamiltonianId.substring(0, lastUnderscore) + "_CVCF";
            if (cecWorkflow.hamiltonianExists(legacyId)) {
                emit(progressSink, "CVM mode: using CVCF Hamiltonian '" + legacyId
                        + "' instead of '" + requestedHamiltonianId + "'");
                LOG.info("CVM mode: using CVCF Hamiltonian '" + legacyId
                        + "' instead of '" + requestedHamiltonianId + "'");
                return legacyId;
            }
        }

        emit(progressSink, "CVM mode: CVCF Hamiltonian not found (tried '" + preferredId
                + "' and legacy pattern), falling back to '" + requestedHamiltonianId + "'");
        LOG.warning("CVM mode: CVCF Hamiltonian not found (tried '" + preferredId
                + "' and legacy pattern), falling back to '" + requestedHamiltonianId + "'");
        return requestedHamiltonianId;
    }

    private static void emit(Consumer<String> progressSink, String line) {
        if (progressSink != null) {
            progressSink.accept(line);
        }
    }

}
