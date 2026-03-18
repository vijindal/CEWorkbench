package org.ce.workflow;

import org.ce.domain.result.ThermodynamicResult;
import org.ce.workflow.thermo.GridScanWorkflow;
import org.ce.workflow.thermo.LineScanWorkflow;
import org.ce.workflow.thermo.ThermodynamicRequest;
import org.ce.workflow.thermo.ThermodynamicWorkflow;

import java.util.List;

/**
 * Central service for thermodynamic calculations.
 *
 * <p>Acts as a dispatcher between UI and calculation workflows.
 * All methods take two IDs:</p>
 * <ul>
 *   <li>{@code clusterId} -- element-agnostic cluster data, e.g. {@code BCC_A2_T_bin}</li>
 *   <li>{@code hamiltonianId} -- element-specific ECI parameters, e.g. {@code Nb-Ti_BCC_A2_T}</li>
 * </ul>
 */
public class CalculationService {

    private final ThermodynamicWorkflow thermoWorkflow;
    private final LineScanWorkflow lineScanWorkflow;
    private final GridScanWorkflow gridScanWorkflow;

    public CalculationService(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
        this.lineScanWorkflow = new LineScanWorkflow(thermoWorkflow);
        this.gridScanWorkflow = new GridScanWorkflow(thermoWorkflow);
    }

    // =========================================================================
    // Single Point Calculation
    // =========================================================================

    public ThermodynamicResult runSinglePoint(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double[] composition,
            String engineType) throws Exception {

        return thermoWorkflow.runCalculation(
                new ThermodynamicRequest(clusterId, hamiltonianId, temperature, composition, engineType)
        );
    }

    // =========================================================================
    // Line Scan (1D)
    // =========================================================================

    public List<ThermodynamicResult> runLineScanTemperature(
            String clusterId,
            String hamiltonianId,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep,
            String engineType) throws Exception {

        return lineScanWorkflow.scanTemperature(
                clusterId, hamiltonianId, composition, tStart, tEnd, tStep, engineType
        );
    }

    public List<ThermodynamicResult> runLineScanComposition(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double xStart,
            double xEnd,
            double xStep,
            String engineType) throws Exception {

        return lineScanWorkflow.scanComposition(
                clusterId, hamiltonianId, temperature, xStart, xEnd, xStep, engineType
        );
    }

    // =========================================================================
    // Grid Scan (2D)
    // =========================================================================

    public List<List<ThermodynamicResult>> runGridScan(
            String clusterId,
            String hamiltonianId,
            double tStart,
            double tEnd,
            double tStep,
            double xStart,
            double xEnd,
            double xStep,
            String engineType) throws Exception {

        return gridScanWorkflow.scanTX(
                clusterId, hamiltonianId, tStart, tEnd, tStep, xStart, xEnd, xStep, engineType
        );
    }
}
