package org.ce.calculation.workflow;

import org.ce.model.ProgressEvent;
import org.ce.model.ModelSession;
import org.ce.calculation.workflow.thermo.ScanWorkflows;
import org.ce.calculation.workflow.thermo.ThermodynamicRequest;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationSpecifications;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Central service for thermodynamic calculations.
 * Unified for multicomponent support and simplified workflow dispatch.
 */
public class CalculationService {

    private final ModelSession.Builder  sessionBuilder;
    private final ThermodynamicWorkflow thermoWorkflow;
    private final ScanWorkflows.LineScan lineScan;
    private final ScanWorkflows.GridScan gridScan;
    private final ScanWorkflows.FiniteSizeScan fssWorkflow;

    private ModelSession cachedSession = null;

    public CalculationService(ModelSession.Builder sessionBuilder,
                              ThermodynamicWorkflow thermoWorkflow) {
        this.sessionBuilder = sessionBuilder;
        this.thermoWorkflow = thermoWorkflow;
        this.lineScan       = new ScanWorkflows.LineScan(thermoWorkflow);
        this.gridScan       = new ScanWorkflows.GridScan(thermoWorkflow);
        this.fssWorkflow    = new ScanWorkflows.FiniteSizeScan(thermoWorkflow);
    }

    /**
     * Role: Execution Entry Point.
     * Decides which workflow to run based on CalcSpecs.
     */
    public Object execute(
            ModelSpecifications modelSpecs,
            CalculationSpecifications calcSpecs,
            Consumer<String> textSink,
            Consumer<ProgressEvent> eventSink) throws Exception {

        return switch (calcSpecs.getMode()) {
            case ANALYSIS -> runAnalysis(modelSpecs, calcSpecs, textSink, eventSink);
            case FINITE_SIZE_SCALING -> {
                ModelSession session = getOrBuildSession(modelSpecs, textSink);
                // Note: FSS currently uses legacy single-composition parameter mapping
                yield fssWorkflow.run(session,
                    calcSpecs.getOrDefault(Parameter.TEMPERATURE),
                    calcSpecs.getOrDefault(Parameter.COMPOSITION),
                    calcSpecs.getProperty(),
                    calcSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                    calcSpecs.getOrDefault(Parameter.MCS_NAVG),
                    textSink, eventSink);
            }
        };
    }

    private Object runAnalysis(
            ModelSpecifications modelSpecs,
            CalculationSpecifications calcSpecs,
            Consumer<String> textSink,
            Consumer<ProgressEvent> eventSink) throws Exception {

        ModelSession session = getOrBuildSession(modelSpecs, textSink);

        double tStart = calcSpecs.getOrDefault(Parameter.T_START);
        double tEnd   = calcSpecs.getOrDefault(Parameter.T_END);
        double tStep  = calcSpecs.getOrDefault(Parameter.T_STEP);
        boolean tVaries = Math.abs(tStart - tEnd) > 1e-6;

        double[] xStarts = calcSpecs.getOrDefault(Parameter.X_STARTS);
        double[] xEnds   = calcSpecs.getOrDefault(Parameter.X_ENDS);
        double[] xSteps  = calcSpecs.getOrDefault(Parameter.X_STEPS);

        List<ScanWorkflows.Varying> vars = new ArrayList<>();
        if (tVaries) {
            vars.add(new ScanWorkflows.Varying("T", true, -1, tStart, tEnd, tStep));
        }

        if (xStarts != null && xEnds != null) {
            for (int i = 0; i < xStarts.length; i++) {
                if (Math.abs(xStarts[i] - xEnds[i]) > 1e-6) {
                    vars.add(new ScanWorkflows.Varying("X" + (i + 2), false, i, xStarts[i], xEnds[i], xSteps[i]));
                }
            }
        }

        int mcsL = calcSpecs.getOrDefault(Parameter.MCS_L);
        int mcsNEquil = calcSpecs.getOrDefault(Parameter.MCS_NEQUIL);
        int mcsNAvg = calcSpecs.getOrDefault(Parameter.MCS_NAVG);

        if (vars.isEmpty()) {
            // 0D: Single Point
            double[] xIndep = (xStarts != null) ? xStarts : new double[session.numComponents() - 1];
            double[] fullComp = ScanWorkflows.deriveComposition(xIndep, session);
            return thermoWorkflow.runCalculation(session, new ThermodynamicRequest(
                    tStart, fullComp, calcSpecs.getProperty(), textSink, eventSink,
                    mcsL, mcsNEquil, mcsNAvg,
                    calcSpecs.getOrDefault(Parameter.FIXED_CORRELATIONS)
            ));
        } else if (vars.size() == 1) {
            // 1D: Generic Scan
            return lineScan.scan1D(session, vars.get(0), xStarts, tStart, calcSpecs.getProperty(), mcsL, mcsNEquil, mcsNAvg, textSink, eventSink);
        } else {
            // 2D: Generic Grid Scan
            if (vars.size() > 2) {
                textSink.accept("Warning: More than 2 variables vary. Truncating to first two.");
            }
            return gridScan.scan2D(session, vars.get(0), vars.get(1), xStarts, tStart, calcSpecs.getProperty(), mcsL, mcsNEquil, mcsNAvg, textSink, eventSink);
        }
    }

    /**
     * Role: Model Construction.
     */
    public ModelSession getOrBuildSession(ModelSpecifications specs, Consumer<String> sink) throws Exception {
        SystemId systemId = new SystemId(specs.elements, specs.structure, specs.modelName);
        if (cachedSession != null && cachedSession.systemId.equals(systemId)
                && cachedSession.engineConfig.equals(specs.engineConfig)) {
            if (sink != null) {
                sink.accept("[Session] Reusing cached session: " + systemId.elements + " / " + systemId.structure + " / " + systemId.model);
                sink.accept("  [Session] ✓ Basis resolved: " + cachedSession.cvcfBasis.numNonPointCfs + " non-point CFs, " 
                                                           + cachedSession.cvcfBasis.numComponents + " point variables");
            }
            return cachedSession;
        }
        cachedSession = sessionBuilder.build(systemId, specs.engineConfig, sink);
        return cachedSession;
    }

}
