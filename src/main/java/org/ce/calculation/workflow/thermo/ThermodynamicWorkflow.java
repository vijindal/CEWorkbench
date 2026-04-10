package org.ce.calculation.workflow.thermo;

import org.ce.calculation.engine.ThermodynamicEngine;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Calculation-layer workflow for thermodynamic calculations.
 *
 * <p>Accepts a pre-built {@link ModelSession} (which already holds cluster data,
 * Hamiltonian, and resolved basis) and a {@link ThermodynamicRequest} (which
 * carries only calculation parameters: T, composition, MCS params, progress sinks).
 * No cluster identification or Hamiltonian loading occurs here.</p>
 */
public class ThermodynamicWorkflow {

    private static final Logger LOG = Logger.getLogger(ThermodynamicWorkflow.class.getName());

    private final ThermodynamicEngine cvmEngine;
    private final ThermodynamicEngine mcsEngine;

    public ThermodynamicWorkflow(
            ThermodynamicEngine cvmEngine,
            ThermodynamicEngine mcsEngine) {

        this.cvmEngine = cvmEngine;
        this.mcsEngine = mcsEngine;
    }

    /**
     * Runs a thermodynamic calculation using pre-built session state.
     *
     * <p>No cluster identification or Hamiltonian loading occurs here — both are
     * already available on {@code session}.</p>
     *
     * @param session pre-built model session holding cluster data, ECI, and basis
     * @param request calculation parameters (T, composition, MCS params, sinks)
     */
    public ThermodynamicResult runCalculation(
            ModelSession session,
            ThermodynamicRequest request) throws Exception {

        LOG.info("ThermodynamicWorkflow.runCalculation — ENTER");
        LOG.info("  session: " + session.label());
        LOG.info("  temperature: " + request.temperature + " K");
        LOG.info("  composition: " + Arrays.toString(request.composition));
        LOG.info("  engineType: " + session.engineConfig.engineType);

        emit(request.progressSink, "STAGE 3: Using pre-loaded Hamiltonian '"
                + session.resolvedHamiltonianId + "'");

        ThermodynamicEngine.Input input = new ThermodynamicEngine.Input(
                session.clusterData,
                session.cecEntry,
                request.temperature,
                request.composition,
                session.resolvedHamiltonianId,
                session.label(),
                request.progressSink,
                request.eventSink,
                request.mcsL,
                request.mcsNEquil,
                request.mcsNAvg
        );

        ThermodynamicEngine engine = selectEngine(session.engineConfig.engineType);
        ThermodynamicResult result = engine.compute(input);

        if (request.progressSink != null) {
            emit(request.progressSink, "");
            emit(request.progressSink, "  RESULTS AT " + request.temperature
                    + " K (composition: " + Arrays.toString(request.composition) + ")");
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

        LOG.info("ThermodynamicWorkflow.runCalculation — EXIT: G="
                + String.format("%.4e", result.gibbsEnergy) + " J/mol");
        return result;
    }

    private ThermodynamicEngine selectEngine(String engineType) {
        return switch (engineType) {
            case "CVM" -> cvmEngine;
            case "MCS" -> mcsEngine;
            default -> throw new IllegalArgumentException("Unknown engine: " + engineType);
        };
    }

    private static void emit(java.util.function.Consumer<String> sink, String line) {
        if (sink != null) sink.accept(line);
    }
}
