package org.ce.workflow;

import org.ce.domain.result.ThermodynamicResult;
import org.ce.workflow.thermo.GridScanWorkflow;
import org.ce.workflow.thermo.LineScanWorkflow;
import org.ce.workflow.thermo.ThermodynamicWorkflow;

import java.util.List;

/**
 * Central service for thermodynamic calculations.
 *
 * <p>Acts as a dispatcher between UI and calculation workflows.
 * Provides a clean, stable API for all calculation modes:
 * single point, line scan, and grid scan.</p>
 *
 * <p>UI talks ONLY to this service; all complexity is hidden.</p>
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

    /**
     * Runs a single-point thermodynamic calculation.
     *
     * @param systemId system identifier
     * @param temperature calculation temperature
     * @param composition mole fractions
     * @param engineType "CVM" or "MCS"
     * @return thermodynamic result
     */
    public ThermodynamicResult runSinglePoint(
            String systemId,
            double temperature,
            double[] composition,
            String engineType) throws Exception {

        return thermoWorkflow.runCalculation(
                new org.ce.workflow.thermo.ThermodynamicRequest(
                        systemId, temperature, composition, engineType
                )
        );
    }

    // =========================================================================
    // Line Scan (1D)
    // =========================================================================

    /**
     * Scans temperature at fixed composition.
     *
     * @param systemId system identifier
     * @param composition fixed mole fractions
     * @param tStart starting temperature
     * @param tEnd ending temperature
     * @param tStep temperature step
     * @param engineType "CVM" or "MCS"
     * @return results for each temperature point
     */
    public List<ThermodynamicResult> runLineScanTemperature(
            String systemId,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep,
            String engineType) throws Exception {

        return lineScanWorkflow.scanTemperature(
                systemId, composition, tStart, tEnd, tStep, engineType
        );
    }

    /**
     * Scans composition at fixed temperature (binary system).
     *
     * @param systemId system identifier
     * @param temperature fixed temperature
     * @param xStart starting mole fraction of component B
     * @param xEnd ending mole fraction of component B
     * @param xStep mole fraction step size
     * @param engineType "CVM" or "MCS"
     * @return results for each composition point
     */
    public List<ThermodynamicResult> runLineScanComposition(
            String systemId,
            double temperature,
            double xStart,
            double xEnd,
            double xStep,
            String engineType) throws Exception {

        return lineScanWorkflow.scanComposition(
                systemId, temperature, xStart, xEnd, xStep, engineType
        );
    }

    // =========================================================================
    // Grid Scan (2D)
    // =========================================================================

    /**
     * Scans temperature and composition (2D grid).
     *
     * <p>Returns a 2D list where grid[i][j] corresponds to
     * temperature T_i and composition x_j. Useful for generating
     * T–x heatmaps and phase diagram data.</p>
     *
     * @param systemId system identifier
     * @param tStart starting temperature
     * @param tEnd ending temperature
     * @param tStep temperature step size
     * @param xStart starting mole fraction of component B
     * @param xEnd ending mole fraction of component B
     * @param xStep mole fraction step size
     * @param engineType "CVM" or "MCS"
     * @return 2D grid of results [temperature][composition]
     */
    public List<List<ThermodynamicResult>> runGridScan(
            String systemId,
            double tStart,
            double tEnd,
            double tStep,
            double xStart,
            double xEnd,
            double xStep,
            String engineType) throws Exception {

        return gridScanWorkflow.scanTX(
                systemId,
                tStart, tEnd, tStep,
                xStart, xEnd, xStep,
                engineType
        );
    }
}
