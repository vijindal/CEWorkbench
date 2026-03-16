package org.ce.domain.engine.cvm;

import org.ce.domain.engine.cvm.CVMPhaseModel;
import org.ce.domain.model.result.EngineMetrics;
import org.ce.domain.model.result.EquilibriumState;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Executor for CVM Phase Model calculations.
 *
 * <p>Bridges between CVMPhaseModel and CalculationProgressListener.
 * Provides parameter scanning and query execution with progress updates.</p>
 *
 * <p>CVMPhaseModelExecutor manages a model that can be queried multiple times
 * with different parameters.</p>
 */
public class CVMPhaseModelExecutor {

    private static final Logger LOG = Logger.getLogger(CVMPhaseModelExecutor.class.getName());

    /**
     * Executes initial thermodynamic evaluation of the CVM phase model.
     *
     * <p>Logs model information and performs first minimization if not already done.</p>
     *
     * @param model The CVM phase model (already initialized by create)
     * @return true if execution succeeds, false otherwise
     */
    public static boolean initializeModel(CVMPhaseModel model) {
        if (model == null) {
            LOG.warning("CVMPhaseModelExecutor: Model is null");
            return false;
        }

        LOG.fine("CVMPhaseModelExecutor.initializeModel — ENTER: T=" + model.getTemperature() + " K");

        try {
            // Model is already minimized in CVMPhaseModel.create()
            // Just log and return success

            EquilibriumState state = model.getEquilibriumState();

            // Extract CVM-specific diagnostics via pattern match on metrics
            EngineMetrics.CvmMetrics cvm = (EngineMetrics.CvmMetrics) state.metrics();

            LOG.fine("CVMPhaseModelExecutor.initializeModel — EXIT: T=" + model.getTemperature()
                    + " K, G=" + String.format("%.8e", state.gibbsEnergy().getAsDouble())
                    + ", iterations=" + cvm.iterations()
                    + ", convergence=" + String.format("%.2e", cvm.gradientNorm()));
            return true;

        } catch (Exception ex) {
            LOG.log(Level.WARNING, "CVMPhaseModelExecutor.initializeModel — FAILED: "
                    + ex.getClass().getSimpleName() + ": " + ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * Executes a parameter query on the CVM phase model.
     *
     * <p>Typically used for parameter scanning (T-scan, x-scan, etc).
     * Automatically triggers re-minimization if parameters have changed.</p>
     *
     * @param model The CVM phase model
     * @return unified EquilibriumState with G, H, S, CFs and CvmMetrics diagnostics
     */
    public static EquilibriumState queryModel(CVMPhaseModel model) {

        try {
            EquilibriumState state = model.getEquilibriumState();

            LOG.fine("CVMPhaseModelExecutor.queryModel — G=" + String.format("%.8e", state.gibbsEnergy().getAsDouble()));

            return state;

        } catch (Exception ex) {
            LOG.log(Level.WARNING, "CVMPhaseModelExecutor.queryModel — FAILED: " + ex.getMessage(), ex);
            return null;
        }
    }
}