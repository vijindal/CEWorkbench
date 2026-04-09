package org.ce.calculation.workflow.thermo;

import org.ce.model.result.ThermodynamicResult;
import org.ce.model.ModelSession;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Performs 1D parameter scan (temperature OR composition).
 *
 * <p>Orchestrates multiple single-point calculations, varying one parameter
 * while holding others constant. Uses a pre-built {@link ModelSession} so
 * cluster identification and Hamiltonian loading happen only once.</p>
 */
public class LineScanWorkflow {

    private static final Logger LOG = Logger.getLogger(LineScanWorkflow.class.getName());

    private final ThermodynamicWorkflow thermoWorkflow;

    public LineScanWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }

    /**
     * Scans temperature at fixed composition.
     */
    public List<ThermodynamicResult> scanTemperature(
            ModelSession session,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep) throws Exception {
        return scanTemperature(session, composition, tStart, tEnd, tStep, 4, 1000, 2000);
    }

    /**
     * Scans temperature at fixed composition with explicit MCS parameters.
     * MCS params are ignored when the session uses CVM.
     */
    public List<ThermodynamicResult> scanTemperature(
            ModelSession session,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg) throws Exception {

        LOG.info("LineScanWorkflow.scanTemperature — ENTER: " + session.label());
        LOG.info("  T range: " + tStart + " to " + tEnd + " K, step " + tStep + " K");

        List<ThermodynamicResult> results = new ArrayList<>();
        for (double T = tStart; T <= tEnd; T += tStep) {
            results.add(thermoWorkflow.runCalculation(session,
                    new ThermodynamicRequest(T, composition, null, null, mcsL, mcsNEquil, mcsNAvg)));
        }

        LOG.info("LineScanWorkflow.scanTemperature — EXIT: " + results.size() + " points");
        return results;
    }

    /**
     * Scans composition at fixed temperature (binary system).
     */
    public List<ThermodynamicResult> scanComposition(
            ModelSession session,
            double temperature,
            double xStart,
            double xEnd,
            double xStep) throws Exception {
        return scanComposition(session, temperature, xStart, xEnd, xStep, 4, 1000, 2000);
    }

    /**
     * Scans composition at fixed temperature with explicit MCS parameters.
     * MCS params are ignored when the session uses CVM.
     */
    public List<ThermodynamicResult> scanComposition(
            ModelSession session,
            double temperature,
            double xStart,
            double xEnd,
            double xStep,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg) throws Exception {

        LOG.info("LineScanWorkflow.scanComposition — ENTER: " + session.label());
        LOG.info("  T=" + temperature + " K, x range: " + xStart + " to " + xEnd);

        List<ThermodynamicResult> results = new ArrayList<>();
        for (double x = xStart; x <= xEnd; x += xStep) {
            double[] comp = new double[]{1.0 - x, x};
            results.add(thermoWorkflow.runCalculation(session,
                    new ThermodynamicRequest(temperature, comp, null, null, mcsL, mcsNEquil, mcsNAvg)));
        }

        LOG.info("LineScanWorkflow.scanComposition — EXIT: " + results.size() + " points");
        return results;
    }
}
