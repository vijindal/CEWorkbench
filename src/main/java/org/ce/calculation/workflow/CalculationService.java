package org.ce.calculation.workflow;

import org.ce.calculation.engine.ProgressEvent;
import org.ce.model.result.ThermodynamicResult;
import org.ce.model.ModelSession;
import org.ce.calculation.workflow.thermo.ScanWorkflows;
import org.ce.calculation.workflow.thermo.ThermodynamicRequest;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;

import java.util.List;
import java.util.function.Consumer;

/**
 * Central service for thermodynamic calculations.
 *
 * <p>Acts as a dispatcher between the UI layer and the calculation workflows.
 * All methods accept a pre-built {@link ModelSession} which carries the
 * pre-computed cluster data, Hamiltonian, and engine configuration. The session
 * is created once per system identity and reused across all calculations.</p>
 */
public class CalculationService {

    private final ThermodynamicWorkflow  thermoWorkflow;
    private final ScanWorkflows.LineScan   lineScanWorkflow;
    private final ScanWorkflows.GridScan   gridScanWorkflow;
    private final ScanWorkflows.FiniteSizeScan fssWorkflow;

    public CalculationService(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow   = thermoWorkflow;
        this.lineScanWorkflow = new ScanWorkflows.LineScan(thermoWorkflow);
        this.gridScanWorkflow = new ScanWorkflows.GridScan(thermoWorkflow);
        this.fssWorkflow      = new ScanWorkflows.FiniteSizeScan(thermoWorkflow);
    }

    // =========================================================================
    // Single Point Calculation
    // =========================================================================

    /** Primary entry point — accepts a pre-built session and a fully constructed request. */
    public ThermodynamicResult runSinglePoint(
            ModelSession session,
            ThermodynamicRequest request) throws Exception {
        return thermoWorkflow.runCalculation(session, request);
    }

    // =========================================================================
    // Finite-Size Scaling Scan (production-grade MCS)
    // =========================================================================

    public ThermodynamicResult runFiniteSizeScan(
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

    // =========================================================================
    // Line Scan (1D)
    // =========================================================================

    public List<ThermodynamicResult> runLineScanTemperature(
            ModelSession session,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep) throws Exception {

        return lineScanWorkflow.scanTemperature(session, composition, tStart, tEnd, tStep);
    }

    public List<ThermodynamicResult> runLineScanTemperature(
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

    public List<ThermodynamicResult> runLineScanComposition(
            ModelSession session,
            double temperature,
            double xStart,
            double xEnd,
            double xStep) throws Exception {

        return lineScanWorkflow.scanComposition(session, temperature, xStart, xEnd, xStep);
    }

    public List<ThermodynamicResult> runLineScanComposition(
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

    // =========================================================================
    // Grid Scan (2D)
    // =========================================================================

    public List<List<ThermodynamicResult>> runGridScan(
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
