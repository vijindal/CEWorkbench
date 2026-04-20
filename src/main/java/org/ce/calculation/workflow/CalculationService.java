package org.ce.calculation.workflow;

import org.ce.calculation.CalculationResult;
import org.ce.model.ProgressEvent;
import org.ce.model.ModelSession;
import org.ce.calculation.workflow.thermo.ScanWorkflows;
import org.ce.calculation.workflow.thermo.ScanWorkflows.Varying;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationSpecifications;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Central service for thermodynamic calculations.
 *
 * <p>All ANALYSIS-mode calculations (single-point, line scan, 2D grid) are unified
 * through {@link ScanWorkflows#scan}: a single-point is a 1×1 grid, a line scan is
 * N×1 or 1×M. {@link CalculationResult.Grid} is always returned for ANALYSIS mode.</p>
 */
public class CalculationService {

    private final ModelSession.Builder  sessionBuilder;
    private final ThermodynamicWorkflow thermoWorkflow;

    private ModelSession cachedSession = null;

    public CalculationService(ModelSession.Builder sessionBuilder,
                              ThermodynamicWorkflow thermoWorkflow) {
        this.sessionBuilder = sessionBuilder;
        this.thermoWorkflow = thermoWorkflow;
    }

    /**
     * Executes a calculation and returns a typed {@link CalculationResult}.
     *
     * <ul>
     *   <li>ANALYSIS → {@link CalculationResult.Grid} (1×1 for single-point,
     *       N×1 for T-scan, 1×M for X-scan, N×M for 2D grid)</li>
     *   <li>FINITE_SIZE_SCALING → {@link CalculationResult.Single}</li>
     * </ul>
     */
    public CalculationResult execute(
            ModelSpecifications modelSpecs,
            CalculationSpecifications calcSpecs,
            Consumer<String> textSink,
            Consumer<ProgressEvent> eventSink) throws Exception {

        return switch (calcSpecs.getMode()) {
            case ANALYSIS -> runAnalysis(modelSpecs, calcSpecs, textSink, eventSink);
            case FINITE_SIZE_SCALING -> {
                ModelSession session = getOrBuildSession(modelSpecs, textSink);
                yield new CalculationResult.Single(ScanWorkflows.finiteSizeScan(
                        thermoWorkflow, session,
                        calcSpecs.getOrDefault(Parameter.TEMPERATURE),
                        calcSpecs.getOrDefault(Parameter.COMPOSITION),
                        calcSpecs.getProperty(),
                        calcSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                        calcSpecs.getOrDefault(Parameter.MCS_NAVG),
                        textSink, eventSink));
            }
        };
    }

    private CalculationResult.Grid runAnalysis(
            ModelSpecifications modelSpecs,
            CalculationSpecifications calcSpecs,
            Consumer<String> textSink,
            Consumer<ProgressEvent> eventSink) throws Exception {

        ModelSession session = getOrBuildSession(modelSpecs, textSink);

        double tStart = calcSpecs.getOrDefault(Parameter.T_START);
        double tEnd   = calcSpecs.getOrDefault(Parameter.T_END);
        double tStep  = calcSpecs.getOrDefault(Parameter.T_STEP);

        double[] xStarts = calcSpecs.getOrDefault(Parameter.X_STARTS);
        double[] xEnds   = calcSpecs.getOrDefault(Parameter.X_ENDS);
        double[] xSteps  = calcSpecs.getOrDefault(Parameter.X_STEPS);

        int mcsL      = calcSpecs.getOrDefault(Parameter.MCS_L);
        int mcsNEquil = calcSpecs.getOrDefault(Parameter.MCS_NEQUIL);
        int mcsNAvg   = calcSpecs.getOrDefault(Parameter.MCS_NAVG);

        // Identify which dimensions vary
        boolean tVaries = Math.abs(tStart - tEnd) > 1e-6;

        List<Varying> xVars = new ArrayList<>();
        if (xStarts != null && xEnds != null) {
            for (int i = 0; i < xStarts.length; i++) {
                if (Math.abs(xStarts[i] - xEnds[i]) > 1e-6) {
                    xVars.add(new Varying("X" + (i + 2), false, i, xStarts[i], xEnds[i], xSteps[i]));
                }
            }
        }

        if (xVars.size() > 2 && textSink != null) {
            textSink.accept("Warning: More than 2 variables vary. Using first two X dimensions.");
        }

        // Build var1 (outer / T dimension) and var2 (inner / X dimension).
        // A sentinel with start==end and step=0 makes the loop run exactly once.
        double baseX0 = (xStarts != null && xStarts.length > 0) ? xStarts[0] : 0.0;

        Varying var1 = tVaries
                ? new Varying("T", true, -1, tStart, tEnd, tStep)
                : new Varying("T", true, -1, tStart, tStart, 0.0);

        Varying var2 = !xVars.isEmpty()
                ? xVars.get(0)
                : new Varying("X2", false, 0, baseX0, baseX0, 0.0);

        // If both T and one X vary, they become var1/var2 above.
        // If two X dims vary, use them as var1/var2 (T is fixed).
        if (!tVaries && xVars.size() >= 2) {
            var1 = xVars.get(0);
            var2 = xVars.get(1);
        }

        // Base independent composition for non-varying X dims
        double[] baseIndepComp = (xStarts != null) ? xStarts.clone()
                : new double[session.numComponents() - 1];

        return new CalculationResult.Grid(ScanWorkflows.scan(
                thermoWorkflow, session, var1, var2,
                baseIndepComp, tStart,
                calcSpecs.getProperty(),
                mcsL, mcsNEquil, mcsNAvg,
                textSink, eventSink));
    }

    /**
     * Returns a cached session when identity (elements, structure, model, engine) is unchanged,
     * otherwise builds a new one via {@link ModelSession.Builder}.
     */
    public ModelSession getOrBuildSession(ModelSpecifications specs, Consumer<String> sink) throws Exception {
        SystemId systemId = new SystemId(specs.elements, specs.structure, specs.modelName);
        if (cachedSession != null && cachedSession.systemId.equals(systemId)
                && cachedSession.engineConfig.equals(specs.engineConfig)) {
            if (sink != null) {
                sink.accept("[Session] Reusing cached session: " + systemId.elements
                        + " / " + systemId.structure + " / " + systemId.model);
                sink.accept("  [Session] ✓ Basis resolved: " + cachedSession.cvcfBasis.numNonPointCfs
                        + " non-point CFs, " + cachedSession.cvcfBasis.numComponents + " point variables");
            }
            return cachedSession;
        }
        cachedSession = sessionBuilder.build(systemId, specs.engineConfig, sink);
        return cachedSession;
    }

    /** Returns the most recently built or reused session, or null. */
    public ModelSession getLastCachedSession() { return cachedSession; }
}
