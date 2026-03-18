package org.ce.workflow.thermo;

import org.ce.domain.result.ThermodynamicResult;

import java.util.ArrayList;
import java.util.List;

/**
 * Performs 1D parameter scan (temperature OR composition).
 *
 * <p>Orchestrates multiple single-point calculations, varying one parameter
 * while holding others constant.</p>
 */
public class LineScanWorkflow {

    private final ThermodynamicWorkflow thermoWorkflow;

    public LineScanWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }

    /**
     * Scans temperature at fixed composition.
     *
     * @param systemId system identifier
     * @param composition fixed mole fractions
     * @param tStart starting temperature
     * @param tEnd ending temperature
     * @param tStep temperature step size
     * @param engineType "CVM" or "MCS"
     * @return results for each temperature point
     */
    public List<ThermodynamicResult> scanTemperature(
            String clusterId,
            String hamiltonianId,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep,
            String engineType) throws Exception {

        List<ThermodynamicResult> results = new ArrayList<>();

        for (double T = tStart; T <= tEnd; T += tStep) {

            ThermodynamicRequest request = new ThermodynamicRequest(
                    clusterId,
                    hamiltonianId,
                    T,
                    composition,
                    engineType
            );

            results.add(thermoWorkflow.runCalculation(request));
        }

        return results;
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
    public List<ThermodynamicResult> scanComposition(
            String clusterId,
            String hamiltonianId,
            double temperature,
            double xStart,
            double xEnd,
            double xStep,
            String engineType) throws Exception {

        List<ThermodynamicResult> results = new ArrayList<>();

        for (double x = xStart; x <= xEnd; x += xStep) {

            double[] composition = new double[]{1.0 - x, x};

            ThermodynamicRequest request = new ThermodynamicRequest(
                    clusterId,
                    hamiltonianId,
                    temperature,
                    composition,
                    engineType
            );

            results.add(thermoWorkflow.runCalculation(request));
        }

        return results;
    }
}
