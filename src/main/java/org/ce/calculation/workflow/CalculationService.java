package org.ce.calculation.workflow;

import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationResult;
import org.ce.calculation.Conditions;
import org.ce.calculation.ConditionsScan;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Central service for thermodynamic calculations and parameter scans.
 *
 * <p>{@link #calculate} / {@link #calculateScan} are the single named entry point —
 * the one call CLI, GUI, and any external Java caller all route through, mirroring
 * pycalphad's {@code equilibrium(dbf, comps, phases, conditions)}. {@link #execute}
 * remains as the metadata-driven entry point used by the GUI's dynamic form
 * ({@link Registry}); it delegates to {@link #calculateScan} internally, so there is
 * exactly one calculation code path underneath both.</p>
 */
public class CalculationService {

    private final ModelSession.Builder  sessionBuilder;
    private final ThermodynamicWorkflow thermoWorkflow;
    private ModelSession cachedSession = null;

    public CalculationService(ModelSession.Builder sessionBuilder, ThermodynamicWorkflow thermoWorkflow) {
        this.sessionBuilder = sessionBuilder;
        this.thermoWorkflow = thermoWorkflow;
    }

    /** Algorithm parameters for the MCS engine. Not physical conditions — CVM callers never need these. */
    public record McsParams(int L, int nEquil, int nAvg) {
        public static final McsParams DEFAULT = new McsParams(4, 1000, 2000);
    }

    /** Single calculation point — the named entry point external callers use directly. */
    public ThermodynamicResult calculate(ModelSession session, Conditions conditions, Property property,
            McsParams mcsParams, Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {
        double[] x = conditions.resolveComposition(session.elements());
        return thermoWorkflow.runCalculation(session, new ThermodynamicWorkflow.Request(
                conditions.temperature(), x, property, textSink, eventSink,
                mcsParams.L(), mcsParams.nEquil(), mcsParams.nAvg(), null));
    }

    public ThermodynamicResult calculate(ModelSession session, Conditions conditions, Property property,
            Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {
        return calculate(session, conditions, property, McsParams.DEFAULT, textSink, eventSink);
    }

    /** Conditions scan — one call, returns every point. */
    public List<ThermodynamicResult> calculateScan(ModelSession session, ConditionsScan scan, Property property,
            McsParams mcsParams, Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {
        scan.validateAgainst(session.elements()); // fail-fast before running anything
        List<ThermodynamicResult> out = new ArrayList<>();
        for (int i = 0; i < scan.pointCount(); i++) {
            out.add(calculate(session, scan.pointAt(i), property, mcsParams, textSink, eventSink));
        }
        return out;
    }

    public List<ThermodynamicResult> calculateScan(ModelSession session, ConditionsScan scan, Property property,
            Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {
        return calculateScan(session, scan, property, McsParams.DEFAULT, textSink, eventSink);
    }

    public CalculationResult execute(ModelSpecifications modelSpecs, JobSpecifications jobSpecs,
                                    Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {

        if (textSink != null) {
            textSink.accept("\n  [Workflow] Stage 1: Loading Specifications...");
            textSink.accept(String.format("    > Model: %s / %s / %s [%s]",
                    modelSpecs.elements(), modelSpecs.structure(), modelSpecs.modelName(), modelSpecs.engineConfig()));
            textSink.accept(String.format("    > Job:   %s in %s mode",
                    jobSpecs.getProperty().displayName, jobSpecs.getMode().displayName));

            ConditionsScan scanPreview = jobSpecs.get(Parameter.CONDITIONS_SCAN)
                    .map(ConditionsScan.class::cast)
                    .orElse(null);
            if (scanPreview != null) {
                textSink.accept(String.format("      - T: %s", scanPreview.temperature()));
                textSink.accept(String.format("      - x: %s", scanPreview.composition()));
            }

            if (modelSpecs.engineConfig().isMcs()) {
                textSink.accept(String.format("      - L: %s", (Integer) jobSpecs.getOrDefault(Parameter.MCS_L)));
                textSink.accept(String.format("      - Equil: %s sweeps", (Integer) jobSpecs.getOrDefault(Parameter.MCS_NEQUIL)));
                textSink.accept(String.format("      - Avg: %s sweeps", (Integer) jobSpecs.getOrDefault(Parameter.MCS_NAVG)));
            }
        }

        ModelSession session = getOrBuildSession(modelSpecs, textSink);

        return switch (jobSpecs.getMode()) {
            case ANALYSIS -> runAnalysis(session, jobSpecs, textSink, eventSink);
        };
    }

    private CalculationResult.Grid runAnalysis(ModelSession session, JobSpecifications jobSpecs,
                                              Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {

        ConditionsScan scan = jobSpecs.get(Parameter.CONDITIONS_SCAN)
                .map(ConditionsScan.class::cast)
                .orElseGet(() -> {
                    Conditions c = jobSpecs.get(Parameter.COMPOSITION)
                            .map(Conditions.class::cast)
                            .orElseThrow(() -> new IllegalStateException(
                                    "No conditions specified (neither COMPOSITION nor CONDITIONS_SCAN set)."));
                    return ConditionsScan.fixedAt(c);
                });

        McsParams mcsParams = new McsParams(
                jobSpecs.getOrDefault(Parameter.MCS_L),
                jobSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                jobSpecs.getOrDefault(Parameter.MCS_NAVG));

        List<ThermodynamicResult> points =
                calculateScan(session, scan, jobSpecs.getProperty(), mcsParams, textSink, eventSink);

        // Grid is 1-D today (single varying axis); wrap as a single row.
        List<List<ThermodynamicResult>> grid = new ArrayList<>();
        grid.add(points);
        return new CalculationResult.Grid(grid);
    }

    public ModelSession getOrBuildSession(ModelSpecifications modelSpecs, Consumer<String> textSink) throws Exception {
        SystemId systemId = new SystemId(modelSpecs.elements(), modelSpecs.structure(), modelSpecs.modelName());
        EngineConfig engineConfig = modelSpecs.engineConfig();

        if (cachedSession != null && cachedSession.systemId.equals(systemId) && cachedSession.engineConfig == engineConfig) {
            if (textSink != null) {
                textSink.accept("  [Session] Reusing cached session: " + cachedSession.label());
            }
            return cachedSession;
        }

        if (textSink != null) {
            textSink.accept("  [Session] Building new session for " + modelSpecs.elements() + " / " + modelSpecs.structure() + " / " + modelSpecs.modelName());
        }
        cachedSession = sessionBuilder.build(systemId, engineConfig, textSink);
        return cachedSession;
    }

    public ModelSession getLastCachedSession() { return cachedSession; }
}
