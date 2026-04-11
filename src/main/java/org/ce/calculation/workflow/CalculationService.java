package org.ce.calculation.workflow;

import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;
import org.ce.calculation.workflow.thermo.ScanWorkflows;
import org.ce.calculation.workflow.thermo.ThermodynamicRequest;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationSpecifications;
import org.ce.model.storage.Workspace.SystemId;

import java.util.List;
import java.util.function.Consumer;

/**
 * Central service for thermodynamic calculations.
 *
 * <p>Implements two primary Roles:
 * <ol>
 *   <li><b>Model Construction:</b> Orchestrates the creation of the persistent model
 *       representation ({@link ModelSession}) using {@link ModelSession.Builder}.</li>
 *   <li><b>Execution:</b> Executes {@link CalculationSpecifications} by dispatching
 *       to the appropriate internal workflow based on Property and Mode.</li>
 * </ol>
 * </p>
 */
public class CalculationService {

    private final ThermodynamicWorkflow  thermoWorkflow;
    private final ScanWorkflows.LineScan   lineScanWorkflow;
    private final ScanWorkflows.GridScan   gridScanWorkflow;
    private final ScanWorkflows.FiniteSizeScan fssWorkflow;
    private final ModelSession.Builder     sessionBuilder;

    private ModelSession cachedSession = null;

    public CalculationService(ThermodynamicWorkflow thermoWorkflow, ModelSession.Builder sessionBuilder) {
        this.thermoWorkflow   = thermoWorkflow;
        this.sessionBuilder   = sessionBuilder;
        this.lineScanWorkflow = new ScanWorkflows.LineScan(thermoWorkflow);
        this.gridScanWorkflow = new ScanWorkflows.GridScan(thermoWorkflow);
        this.fssWorkflow      = new ScanWorkflows.FiniteSizeScan(thermoWorkflow);
    }

    // =========================================================================
    // Unified Entry Point
    // =========================================================================

    /**
     * Unified execution entry point for all calculations.
     *
     * @param modelSpecs global system identity (for model construction)
     * @param calcSpecs  calculation details (property, mode, parameters)
     * @param textSink   optional sink for progress messages
     * @return single point result or the primary point of a scan
     * @throws Exception if model construction or execution fails
     */
    public ThermodynamicResult execute(
            ModelSpecifications modelSpecs,
            CalculationSpecifications calcSpecs,
            Consumer<String> textSink) throws Exception {

        // 1. Model Construction Role
        ModelSession session = getOrBuildSession(modelSpecs, textSink);

        // 2. Execution Role
        return switch (calcSpecs.getMode()) {
            case SINGLE_POINT -> runSinglePoint(session, mapToRequest(calcSpecs, textSink));

            case LINE_SCAN -> {
                List<ThermodynamicResult> rs = executeScan(modelSpecs, calcSpecs, textSink);
                yield rs.get(0);
            }

            case FINITE_SIZE_SCALING -> runFiniteSizeScan(session,
                    calcSpecs.getOrDefault(Parameter.TEMPERATURE),
                    calcSpecs.getOrDefault(Parameter.COMPOSITION),
                    calcSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                    calcSpecs.getOrDefault(Parameter.MCS_NAVG),
                    textSink, null);

            default -> throw new UnsupportedOperationException("Mode not yet implemented in unified API: " + calcSpecs.getMode());
        };
    }

    /**
     * Unified execution entry point for scans that return multiple results.
     */
    public List<ThermodynamicResult> executeScan(
            ModelSpecifications modelSpecs,
            CalculationSpecifications calcSpecs,
            Consumer<String> textSink) throws Exception {

        ModelSession session = getOrBuildSession(modelSpecs, textSink);

        if (calcSpecs.getMode() == Mode.LINE_SCAN) {
            if (calcSpecs.get(Parameter.T_START).isPresent()) {
                return runLineScanTemperature(session,
                        calcSpecs.getOrDefault(Parameter.COMPOSITION),
                        calcSpecs.getOrDefault(Parameter.T_START),
                        calcSpecs.getOrDefault(Parameter.T_END),
                        calcSpecs.getOrDefault(Parameter.T_STEP),
                        calcSpecs.getOrDefault(Parameter.MCS_L),
                        calcSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                        calcSpecs.getOrDefault(Parameter.MCS_NAVG));
            } else {
                return runLineScanComposition(session,
                        calcSpecs.getOrDefault(Parameter.TEMPERATURE),
                        calcSpecs.getOrDefault(Parameter.X_START),
                        calcSpecs.getOrDefault(Parameter.X_END),
                        calcSpecs.getOrDefault(Parameter.X_STEP),
                        calcSpecs.getOrDefault(Parameter.MCS_L),
                        calcSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                        calcSpecs.getOrDefault(Parameter.MCS_NAVG));
            }
        }
        
        throw new UnsupportedOperationException("Scan mode not implemented: " + calcSpecs.getMode());
    }

    /**
     * Unified execution entry point for 2D grid scans.
     */
    public List<List<ThermodynamicResult>> executeGridScan(
            ModelSpecifications modelSpecs,
            CalculationSpecifications calcSpecs,
            Consumer<String> textSink) throws Exception {

        ModelSession session = getOrBuildSession(modelSpecs, textSink);
        return runGridScan(session,
                calcSpecs.getOrDefault(Parameter.T_START),
                calcSpecs.getOrDefault(Parameter.T_END),
                calcSpecs.getOrDefault(Parameter.T_STEP),
                calcSpecs.getOrDefault(Parameter.X_START),
                calcSpecs.getOrDefault(Parameter.X_END),
                calcSpecs.getOrDefault(Parameter.X_STEP));
    }

    /**
     * Role: Model Construction.
     * Builds the persistent model representation if it doesn't match the current specifications.
     */
    public ModelSession getOrBuildSession(ModelSpecifications specs, Consumer<String> sink) throws Exception {
        SystemId systemId = new SystemId(specs.elements, specs.structure, specs.modelName);
        if (cachedSession != null && cachedSession.systemId.equals(systemId)
                && cachedSession.engineConfig.equals(specs.engineConfig)) {
            return cachedSession;
        }
        cachedSession = sessionBuilder.build(systemId, specs.engineConfig, sink);
        return cachedSession;
    }

    private ThermodynamicRequest mapToRequest(CalculationSpecifications specs, Consumer<String> textSink) {
        return new ThermodynamicRequest(
                specs.getOrDefault(Parameter.TEMPERATURE),
                specs.getOrDefault(Parameter.COMPOSITION),
                textSink, null,
                specs.getOrDefault(Parameter.MCS_L),
                specs.getOrDefault(Parameter.MCS_NEQUIL),
                specs.getOrDefault(Parameter.MCS_NAVG),
                specs.getOrDefault(Parameter.FIXED_CORRELATIONS)
        );
    }

    // =========================================================================
    // Internal Workflow Dispatchers
    // =========================================================================

    private ThermodynamicResult runSinglePoint(
            ModelSession session,
            ThermodynamicRequest request) throws Exception {
        return thermoWorkflow.runCalculation(session, request);
    }

    private ThermodynamicResult runFiniteSizeScan(
            ModelSession session,
            double temperature,
            double[] composition,
            int nEquil,
            int nAvg,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink) throws Exception {

        return fssWorkflow.run(session, temperature, composition,
                nEquil, nAvg, progressSink, eventSink);
    }

    private List<ThermodynamicResult> runLineScanTemperature(
            ModelSession session,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg) throws Exception {

        return lineScanWorkflow.scanTemperature(session, composition, tStart, tEnd, tStep,
                mcsL, mcsNEquil, mcsNAvg);
    }

    private List<ThermodynamicResult> runLineScanComposition(
            ModelSession session,
            double temperature,
            double xStart,
            double xEnd,
            double xStep,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg) throws Exception {

        return lineScanWorkflow.scanComposition(session, temperature, xStart, xEnd, xStep,
                mcsL, mcsNEquil, mcsNAvg);
    }

    private List<List<ThermodynamicResult>> runGridScan(
            ModelSession session,
            double tStart,
            double tEnd,
            double tStep,
            double xStart,
            double xEnd,
            double xStep) throws Exception {

        return gridScanWorkflow.scanTX(session, tStart, tEnd, tStep, xStart, xEnd, xStep);
    }
}
