package org.ce.calculation.workflow.thermo;

import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Container for all parameter scan workflows.
 * Unified for multicomponent support.
 */
public class ScanWorkflows {

    private ScanWorkflows() {}

    /**
     * Helper to define a varying parameter in a scan.
     */
    public record Varying(String label, boolean isTemp, int compIndex, double start, double end, double step) {}

    /**
     * Derives the full n-component composition array from n-1 independent components.
     * @param xIndep Independent components [x2, x3, ..., xn]
     * @param session Active session providing component count
     * @return Full normalized array [x1, x2, ..., xn]
     */
    public static double[] deriveComposition(double[] xIndep, ModelSession session) {
        int nTotal = session.numComponents();
        double[] full = new double[nTotal];
        double sumIndep = 0;
        for (int i = 0; i < xIndep.length && i + 1 < nTotal; i++) {
            full[i + 1] = xIndep[i];
            sumIndep += xIndep[i];
        }
        full[0] = 1.0 - sumIndep;
        return full;
    }

    // =========================================================================
    // Grid Scan: 2D (T x Xi, or Xi x Xj)
    // =========================================================================

    public static class GridScan {
        private final ThermodynamicWorkflow thermoWorkflow;

        public GridScan(ThermodynamicWorkflow thermoWorkflow) {
            this.thermoWorkflow = thermoWorkflow;
        }

        /**
         * Generic 2D scan over any two variables (T or Xi).
         */
        public List<List<ThermodynamicResult>> scan2D(
                ModelSession session,
                Varying var1,
                Varying var2,
                double[] baseIndepComp,
                double baseT,
                int mcsL, int mcsNEquil, int mcsNAvg,
                Consumer<String> strSink,
                Consumer<ProgressEvent> eventSink) throws Exception {

            List<List<ThermodynamicResult>> grid = new ArrayList<>();

            for (double v1 = var1.start; (var1.step > 0 ? v1 <= var1.end + 1e-9 : v1 >= var1.end - 1e-9); v1 += var1.step) {
                if (strSink != null) strSink.accept(String.format("Scanning row: %s = %.2f", var1.label, v1));
                List<ThermodynamicResult> row = new ArrayList<>();
                for (double v2 = var2.start; (var2.step > 0 ? v2 <= var2.end + 1e-9 : v2 >= var2.end - 1e-9); v2 += var2.step) {

                    double T = baseT;
                    double[] xIndep = (baseIndepComp != null) ? baseIndepComp.clone() : new double[session.numComponents()-1];

                    if (var1.isTemp) T = v1; else xIndep[var1.compIndex] = v1;
                    if (var2.isTemp) T = v2; else xIndep[var2.compIndex] = v2;

                    row.add(thermoWorkflow.runCalculation(session,
                            new ThermodynamicRequest(T, deriveComposition(xIndep, session), strSink, eventSink, mcsL, mcsNEquil, mcsNAvg)));
                    if (Math.abs(var2.step) < 1e-10) break;
                }
                grid.add(row);
                if (Math.abs(var1.step) < 1e-10) break;
            }
            return grid;
        }

    }

    // =========================================================================
    // Line Scan: 1D (T or Xi)
    // =========================================================================

    public static class LineScan {
        private static final Logger LOG = Logger.getLogger(LineScan.class.getName());
        private final ThermodynamicWorkflow thermoWorkflow;

        public LineScan(ThermodynamicWorkflow thermoWorkflow) {
            this.thermoWorkflow = thermoWorkflow;
        }

        public List<ThermodynamicResult> scan1D(
                ModelSession session,
                Varying var,
                double[] baseIndepComp,
                double baseT,
                int mcsL, int mcsNEquil, int mcsNAvg,
                Consumer<String> progressSink,
                Consumer<ProgressEvent> eventSink) throws Exception {

            List<ThermodynamicResult> results = new ArrayList<>();
            for (double v = var.start; (var.step > 0 ? v <= var.end + 1e-9 : v >= var.end - 1e-9); v += var.step) {
                double T = baseT;
                double[] xIndep = (baseIndepComp != null) ? baseIndepComp.clone() : new double[session.numComponents()-1];
                
                if (var.isTemp) T = v; else xIndep[var.compIndex] = v;

                results.add(thermoWorkflow.runCalculation(session,
                        new ThermodynamicRequest(T, deriveComposition(xIndep, session), progressSink, eventSink, mcsL, mcsNEquil, mcsNAvg)));
                
                if (Math.abs(var.step) < 1e-10) break;
            }
            return results;
        }

    }

    // =========================================================================
    // Finite-Size Scaling Scan
    // =========================================================================

    public static class FiniteSizeScan {
        private final ThermodynamicWorkflow thermoWorkflow;

        public FiniteSizeScan(ThermodynamicWorkflow thermoWorkflow) {
            this.thermoWorkflow = thermoWorkflow;
        }

        public ThermodynamicResult run(
                ModelSession session,
                double temperature,
                double[] composition,
                int nEquil,
                int nAvg,
                Consumer<String> progressSink,
                Consumer<ProgressEvent> eventSink) throws Exception {

            int[] Ls = {12, 16, 24};
            ThermodynamicResult[] results = new ThermodynamicResult[Ls.length];

            for (int i = 0; i < Ls.length; i++) {
                if (progressSink != null) progressSink.accept("FSS: Running L=" + Ls[i]);
                results[i] = thermoWorkflow.runCalculation(session,
                        new ThermodynamicRequest(temperature, composition, progressSink, eventSink, Ls[i], nEquil, nAvg));
            }

            // Simple extrapolation logic (placeholder for actual FSS physics)
            double hInf = results[results.length - 1].enthalpy;
            double sigHInf = results[results.length - 1].stdEnthalpy;
            double cvInf = results[results.length - 1].heatCapacity;
            double[] extrapCFs = results[results.length - 1].avgCFs;
            double[] extrapStdCFs = results[results.length - 1].stdCFs;

            return new ThermodynamicResult(
                    results[0].temperature,
                    results[0].composition,
                    results[0].gibbsEnergy, // Placeholder
                    hInf,
                    Double.NaN,
                    sigHInf,
                    cvInf,
                    null,
                    extrapCFs,
                    extrapStdCFs
            );
        }
    }
}
