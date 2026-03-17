package org.ce.workflow.thermo;

import org.ce.domain.result.ThermodynamicResult;

import java.util.ArrayList;
import java.util.List;

/**
 * Performs 2D scan (e.g., Temperature vs Composition).
 *
 * <p>Orchestrates a grid of single-point calculations, varying two parameters
 * simultaneously to build 2D thermodynamic maps.</p>
 */
public class GridScanWorkflow {

    private final ThermodynamicWorkflow thermoWorkflow;

    public GridScanWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }

    /**
     * Scans temperature and composition (2D grid).
     *
     * <p>Returns a 2D list where grid[i][j] corresponds to
     * temperature T_i and composition x_j.</p>
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
    public List<List<ThermodynamicResult>> scanTX(
            String systemId,
            double tStart,
            double tEnd,
            double tStep,
            double xStart,
            double xEnd,
            double xStep,
            String engineType) throws Exception {

        List<List<ThermodynamicResult>> grid = new ArrayList<>();

        for (double T = tStart; T <= tEnd; T += tStep) {

            List<ThermodynamicResult> row = new ArrayList<>();

            for (double x = xStart; x <= xEnd; x += xStep) {

                double[] composition = new double[]{1.0 - x, x};

                ThermodynamicRequest request = new ThermodynamicRequest(
                        systemId,
                        T,
                        composition,
                        engineType
                );

                row.add(thermoWorkflow.runCalculation(request));
            }

            grid.add(row);
        }

        return grid;
    }
}
