package org.ce.calculation.workflow.thermo;

import org.ce.model.result.ThermodynamicResult;
import org.ce.model.ModelSession;

import java.util.ArrayList;
import java.util.List;

/**
 * Performs 2D scan (Temperature × Composition).
 *
 * <p>Uses a pre-built {@link ModelSession} so cluster identification and
 * Hamiltonian loading happen only once for the entire grid.</p>
 */
public class GridScanWorkflow {

    private final ThermodynamicWorkflow thermoWorkflow;

    public GridScanWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }

    /**
     * Scans temperature and composition (2D grid).
     *
     * @return 2D list where {@code grid.get(i).get(j)} corresponds to T_i and x_j
     */
    public List<List<ThermodynamicResult>> scanTX(
            ModelSession session,
            double tStart,
            double tEnd,
            double tStep,
            double xStart,
            double xEnd,
            double xStep) throws Exception {

        List<List<ThermodynamicResult>> grid = new ArrayList<>();

        for (double T = tStart; T <= tEnd; T += tStep) {
            List<ThermodynamicResult> row = new ArrayList<>();
            for (double x = xStart; x <= xEnd; x += xStep) {
                double[] comp = new double[]{1.0 - x, x};
                row.add(thermoWorkflow.runCalculation(session,
                        new ThermodynamicRequest(T, comp)));
            }
            grid.add(row);
        }

        return grid;
    }
}
