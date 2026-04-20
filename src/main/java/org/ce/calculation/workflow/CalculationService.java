package org.ce.calculation.workflow;

import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationResult;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.model.ModelSession;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Central service for thermodynamic calculations and parameter scans.
 * Absorbed orchestration logic previously in ScanWorkflows.
 */
public class CalculationService {

    private final ModelSession.Builder  sessionBuilder;
    private final ThermodynamicWorkflow thermoWorkflow;
    private ModelSession cachedSession = null;

    public CalculationService(ModelSession.Builder sessionBuilder, ThermodynamicWorkflow thermoWorkflow) {
        this.sessionBuilder = sessionBuilder;
        this.thermoWorkflow = thermoWorkflow;
    }

    public CalculationResult execute(ModelSpecifications modelSpecs, JobSpecifications jobSpecs,
                                    Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {
        ModelSession session = getOrBuildSession(modelSpecs, textSink);

        return switch (jobSpecs.getMode()) {
            case ANALYSIS -> runAnalysis(session, jobSpecs, textSink, eventSink);
            case FINITE_SIZE_SCALING -> new CalculationResult.Single(runFiniteSizeScan(session, jobSpecs, textSink, eventSink));
        };
    }

    private CalculationResult.Grid runAnalysis(ModelSession session, JobSpecifications jobSpecs,
                                              Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {
        
        double tStart = jobSpecs.getOrDefault(Parameter.T_START);
        double tEnd   = jobSpecs.getOrDefault(Parameter.T_END);
        double tStep  = jobSpecs.getOrDefault(Parameter.T_STEP);
        double[] xStarts = jobSpecs.getOrDefault(Parameter.X_STARTS);
        double[] xEnds   = jobSpecs.getOrDefault(Parameter.X_ENDS);
        double[] xSteps  = jobSpecs.getOrDefault(Parameter.X_STEPS);

        // Grid scan logic (formerly ScanWorkflows.scan)
        boolean tVaries = Math.abs(tStart - tEnd) > 1e-6;
        Varying v1 = tVaries ? new Varying("T", true, -1, tStart, tEnd, tStep)
                             : new Varying("T", true, -1, tStart, tStart, 0.0);
        
        double[] baseIndepComp = (xStarts != null) ? xStarts.clone() : new double[session.numComponents() - 1];
        Varying v2 = new Varying("X2", false, 0, baseIndepComp[0], baseIndepComp[0], 0.0);
        if (xStarts != null && xEnds != null && xSteps != null) {
            for (int i = 0; i < xStarts.length; i++) {
                 if (Math.abs(xStarts[i] - xEnds[i]) > 1e-6) {
                     v2 = new Varying("X" + (i + 2), false, i, xStarts[i], xEnds[i], xSteps[i]);
                     break;
                 }
            }
        }

        List<List<ThermodynamicResult>> grid = new ArrayList<>();
        for (double val1 = v1.start(); v1.step() > 0 ? val1 <= v1.end() + 1e-9 : val1 >= v1.end() - 1e-9; val1 += (v1.step() != 0 ? v1.step() : 1.0)) {
            List<ThermodynamicResult> row = new ArrayList<>();
            for (double val2 = v2.start(); v2.step() > 0 ? val2 <= v2.end() + 1e-9 : val2 >= v2.end() - 1e-9; val2 += (v2.step() != 0 ? v2.step() : 1.0)) {
                
                double T = tVaries ? val1 : tStart;
                double[] xIndep = baseIndepComp.clone();
                if (!v1.isTemp() && v1.compIndex() >= 0) xIndep[v1.compIndex()] = val1;
                if (!v2.isTemp() && v2.compIndex() >= 0) xIndep[v2.compIndex()] = val2;

                row.add(thermoWorkflow.runCalculation(session, new ThermodynamicWorkflow.Request(
                        T, deriveComposition(xIndep, session), jobSpecs.getProperty(), textSink, eventSink,
                        jobSpecs.getOrDefault(Parameter.MCS_L), jobSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                        jobSpecs.getOrDefault(Parameter.MCS_NAVG), null)));
                if (v2.step() == 0) break;
            }
            grid.add(row);
            if (v1.step() == 0) break;
        }
        return new CalculationResult.Grid(grid);
    }

    private ThermodynamicResult runFiniteSizeScan(ModelSession session, JobSpecifications jobSpecs,
                                                 Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {
        int[] Ls = {12, 16, 24};
        ThermodynamicResult last = null;
        for (int L : Ls) {
            last = thermoWorkflow.runCalculation(session, new ThermodynamicWorkflow.Request(
                jobSpecs.getOrDefault(Parameter.TEMPERATURE), jobSpecs.getOrDefault(Parameter.COMPOSITION),
                jobSpecs.getProperty(), textSink, eventSink, L, jobSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                jobSpecs.getOrDefault(Parameter.MCS_NAVG), null));
        }
        return last;
    }

    private double[] deriveComposition(double[] xIndep, ModelSession session) {
        double[] full = new double[session.numComponents()];
        double sum = 0;
        for (int i = 0; i < xIndep.length; i++) {
            full[i + 1] = xIndep[i];
            sum += xIndep[i];
        }
        full[0] = 1.0 - sum;
        return full;
    }

    private record Varying(String label, boolean isTemp, int compIndex, double start, double end, double step) {}

    public ModelSession getOrBuildSession(ModelSpecifications specs, Consumer<String> sink) throws Exception {
        SystemId sid = new SystemId(specs.elements(), specs.structure(), specs.modelName());
        if (cachedSession != null && cachedSession.systemId.equals(sid) && cachedSession.engineConfig.equals(specs.engineConfig())) {
            if (sink != null) sink.accept("[Session] Reusing cached session: " + sid);
            return cachedSession;
        }
        cachedSession = sessionBuilder.build(sid, specs.engineConfig(), sink);
        return cachedSession;
    }

    public ModelSession getLastCachedSession() { return cachedSession; }
}
