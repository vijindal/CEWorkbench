package org.ce.workflow.thermo;

import org.ce.domain.result.ThermodynamicResult;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Performs 1D parameter scan (temperature OR composition).
 *
 * <p>Orchestrates multiple single-point calculations, varying one parameter
 * while holding others constant.</p>
 */
public class LineScanWorkflow {

    private static final Logger LOG = Logger.getLogger(LineScanWorkflow.class.getName());

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

        LOG.info("LineScanWorkflow.scanTemperature — ENTER");
        LOG.info("  clusterId: " + clusterId);
        LOG.info("  hamiltonianId: " + hamiltonianId);
        LOG.info("  composition: " + java.util.Arrays.toString(composition));
        LOG.info("  T range: " + tStart + " to " + tEnd + " K, step " + tStep + " K");
        LOG.info("  engineType: " + engineType);

        List<ThermodynamicResult> results = new ArrayList<>();
        int pointNum = 0;
        int totalPoints = (int) Math.ceil((tEnd - tStart) / tStep) + 1;

        for (double T = tStart; T <= tEnd; T += tStep) {
            pointNum++;
            LOG.info("  Point " + pointNum + "/" + totalPoints + ": T=" + T + " K");

            ThermodynamicRequest request = new ThermodynamicRequest(
                    clusterId,
                    hamiltonianId,
                    T,
                    composition,
                    engineType
            );

            results.add(thermoWorkflow.runCalculation(request));
        }

        LOG.info("LineScanWorkflow.scanTemperature — EXIT: completed " + results.size() + " points");
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

        LOG.info("LineScanWorkflow.scanComposition — ENTER");
        LOG.info("  clusterId: " + clusterId);
        LOG.info("  hamiltonianId: " + hamiltonianId);
        LOG.info("  temperature: " + temperature + " K");
        LOG.info("  x range: " + xStart + " to " + xEnd + ", step " + xStep);
        LOG.info("  engineType: " + engineType);

        List<ThermodynamicResult> results = new ArrayList<>();
        int pointNum = 0;
        int totalPoints = (int) Math.ceil((xEnd - xStart) / xStep) + 1;

        for (double x = xStart; x <= xEnd; x += xStep) {
            pointNum++;
            LOG.info("  Point " + pointNum + "/" + totalPoints + ": x_B=" + String.format("%.3f", x));

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

        LOG.info("LineScanWorkflow.scanComposition — EXIT: completed " + results.size() + " points");
        return results;
    }
}
