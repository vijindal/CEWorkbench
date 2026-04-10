package org.ce.calculation.workflow.thermo;

import org.ce.calculation.engine.ProgressEvent;
import org.ce.model.result.ThermodynamicResult;
import org.ce.model.ModelSession;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Container for all parameter scan workflows.
 *
 * <p>Groups three related scan strategies that share the same job (parameter scan)
 * but differ in scan dimension. All are instantiated by {@link CalculationService}
 * with a shared {@link ThermodynamicWorkflow}.</p>
 */
public class ScanWorkflows {

    private ScanWorkflows() {}

    // =========================================================================
    // Grid Scan: 2D (Temperature × Composition)
    // =========================================================================

    /**
     * Performs 2D scan (Temperature × Composition).
     *
     * <p>Uses a pre-built {@link ModelSession} so cluster identification and
     * Hamiltonian loading happen only once for the entire grid.</p>
     */
    public static class GridScan {

        private final ThermodynamicWorkflow thermoWorkflow;

        public GridScan(ThermodynamicWorkflow thermoWorkflow) {
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

    // =========================================================================
    // Line Scan: 1D (Temperature OR Composition)
    // =========================================================================

    /**
     * Performs 1D parameter scan (temperature OR composition).
     *
     * <p>Orchestrates multiple single-point calculations, varying one parameter
     * while holding others constant. Uses a pre-built {@link ModelSession} so
     * cluster identification and Hamiltonian loading happen only once.</p>
     */
    public static class LineScan {

        private static final Logger LOG = Logger.getLogger(LineScan.class.getName());

        private final ThermodynamicWorkflow thermoWorkflow;

        public LineScan(ThermodynamicWorkflow thermoWorkflow) {
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

            LOG.info("LineScan.scanTemperature — ENTER: " + session.label());
            LOG.info("  T range: " + tStart + " to " + tEnd + " K, step " + tStep + " K");

            List<ThermodynamicResult> results = new ArrayList<>();
            for (double T = tStart; T <= tEnd; T += tStep) {
                results.add(thermoWorkflow.runCalculation(session,
                        new ThermodynamicRequest(T, composition, null, null, mcsL, mcsNEquil, mcsNAvg)));
            }

            LOG.info("LineScan.scanTemperature — EXIT: " + results.size() + " points");
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

            LOG.info("LineScan.scanComposition — ENTER: " + session.label());
            LOG.info("  T=" + temperature + " K, x range: " + xStart + " to " + xEnd);

            List<ThermodynamicResult> results = new ArrayList<>();
            for (double x = xStart; x <= xEnd; x += xStep) {
                double[] comp = new double[]{1.0 - x, x};
                results.add(thermoWorkflow.runCalculation(session,
                        new ThermodynamicRequest(temperature, comp, null, null, mcsL, mcsNEquil, mcsNAvg)));
            }

            LOG.info("LineScan.scanComposition — EXIT: " + results.size() + " points");
            return results;
        }
    }

    // =========================================================================
    // Finite-Size Scaling Scan
    // =========================================================================

    /**
     * Finite-size scaling workflow for production-grade MCS results.
     *
     * <p>Runs three MCS simulations at L=12, L=16, L=24 (N=3456, 8192, 27648 BCC sites),
     * fits each observable to the 1/N finite-size scaling model:</p>
     *
     * <pre>  O(N) = O(∞) + A/N</pre>
     *
     * <p>Uses a pre-built {@link ModelSession} so cluster identification and Hamiltonian
     * loading happen only once across all three L-value runs.</p>
     */
    public static class FiniteSizeScan {

        /** Default BCC supercell sizes: even, B2-compatible, 8× ratio in N. */
        public static final int[] DEFAULT_L_VALUES = {12, 16, 24};

        private final ThermodynamicWorkflow thermoWorkflow;

        public FiniteSizeScan(ThermodynamicWorkflow thermoWorkflow) {
            this.thermoWorkflow = thermoWorkflow;
        }

        /**
         * Runs MCS at three L values, extrapolates to thermodynamic limit.
         *
         * @return ThermodynamicResult whose enthalpy/stdEnthalpy/heatCapacity/avgCFs/stdCFs
         *         hold the extrapolated (∞-limit) values and their uncertainties.
         */
        public ThermodynamicResult run(
                ModelSession session,
                double temperature,
                double[] composition,
                int nEquil,
                int nAvg,
                Consumer<String> strSink,
                Consumer<ProgressEvent> evtSink) throws Exception {

            int[] lValues = DEFAULT_L_VALUES;
            int nRuns = lValues.length;

            double[]               N       = new double[nRuns];
            double[]               invN    = new double[nRuns];
            ThermodynamicResult[]  results = new ThermodynamicResult[nRuns];

            for (int i = 0; i < nRuns; i++) {
                int L = lValues[i];
                N[i]    = 2.0 * L * L * L;
                invN[i] = 1.0 / N[i];

                if (strSink != null) strSink.accept(String.format(
                    "\n  ══ FSS run %d/%d: L=%d  (%d BCC sites) ══════════════════════════",
                    i + 1, nRuns, L, (int) N[i]));

                ThermodynamicRequest req = new ThermodynamicRequest(
                    temperature, composition, strSink, evtSink, L, nEquil, nAvg);

                results[i] = thermoWorkflow.runCalculation(session, req);
            }

            return extrapolateAndReport(lValues, N, invN, results, strSink);
        }

        // =====================================================================
        // Extrapolation
        // =====================================================================

        private ThermodynamicResult extrapolateAndReport(
                int[] lValues, double[] N, double[] invN,
                ThermodynamicResult[] results,
                Consumer<String> sink) {

            int nRuns = results.length;

            double[] H    = new double[nRuns];
            double[] sigH = new double[nRuns];
            double[] Cv   = new double[nRuns];
            double[] sigCv = new double[nRuns];

            for (int i = 0; i < nRuns; i++) {
                ThermodynamicResult r = results[i];
                H[i]    = r.enthalpy;
                sigH[i] = Double.isNaN(r.stdEnthalpy) || r.stdEnthalpy < 1e-12
                            ? Math.abs(r.enthalpy) * 0.005 : r.stdEnthalpy;
                Cv[i]   = Double.isNaN(r.heatCapacity) ? 0.0 : r.heatCapacity;
                sigCv[i] = 0.1 * Math.max(1e-8, Cv[i]);
            }

            double[] fitH  = weightedLinearFit(invN, H, sigH);
            double hInf    = fitH[0];
            double sigHInf = fitH[1];
            double hSlope  = fitH[2];

            double[] fitCv = weightedLinearFit(invN, Cv, sigCv);
            double cvInf   = fitCv[0];

            if (sink != null) {
                sink.accept("\n  ══ Finite-Size Scaling: O(N) = O(∞) + A/N ══════════════════════");
                sink.accept(String.format("  %-5s  %-8s  %-25s  %-20s",
                    "L", "N", "⟨H⟩/site (J/mol)", "Cv (J/mol·K)"));
                sink.accept("  " + "─".repeat(65));
                for (int i = 0; i < nRuns; i++) {
                    sink.accept(String.format("  %-5d  %-8d  %+.6f ± %-12.6f  %.6f",
                        lValues[i], (int) N[i], H[i], sigH[i], Cv[i]));
                }
                sink.accept("  " + "─".repeat(65));
                sink.accept(String.format("  %-5s  %-8s  %+.6f ± %-12.6f  %.6f",
                    "∞", "∞", hInf, sigHInf, cvInf));
                sink.accept(String.format("  Finite-size coefficient A = %.2f J/mol  " +
                    "(convergence scale N* ≈ %.0f sites)", hSlope, Math.abs(hSlope / hInf)));
            }

            double[] extrapCFs    = null;
            double[] extrapStdCFs = null;

            if (results[0].avgCFs != null) {
                int nCF = results[0].avgCFs.length;
                extrapCFs    = new double[nCF];
                extrapStdCFs = new double[nCF];

                if (sink != null) {
                    sink.accept("\n  ══ Extrapolated Correlation Functions ═══════════════════════════");
                    sink.accept(String.format("  %-8s  %-14s  %-14s  %-14s  %s",
                        "CF", "L=" + lValues[0], "L=" + lValues[1], "L=" + lValues[2], "∞ (extrap.)"));
                    sink.accept("  " + "─".repeat(75));
                }

                for (int t = 0; t < nCF; t++) {
                    double[] cfVals = new double[nRuns];
                    double[] cfSigs = new double[nRuns];
                    boolean nonzero = false;

                    for (int i = 0; i < nRuns; i++) {
                        cfVals[i] = results[i].avgCFs != null ? results[i].avgCFs[t] : 0.0;
                        double s  = (results[i].stdCFs != null && t < results[i].stdCFs.length)
                                    ? results[i].stdCFs[t] : 0.0;
                        cfSigs[i] = s > 1e-12 ? s : 1e-6;
                        if (Math.abs(cfVals[i]) > 1e-8) nonzero = true;
                    }

                    double[] fitCF   = weightedLinearFit(invN, cfVals, cfSigs);
                    extrapCFs[t]     = fitCF[0];
                    extrapStdCFs[t]  = fitCF[1];

                    if (sink != null && nonzero) {
                        sink.accept(String.format("  CF[%2d]    %+.6f      %+.6f      %+.6f      %+.6f ± %.6f",
                            t, cfVals[0], cfVals[1], cfVals[2], fitCF[0], fitCF[1]));
                    }
                }
            }

            if (sink != null) {
                sink.accept("\n  ✓ Finite-size scan complete. Extrapolated values are thermodynamic-limit estimates.");
            }

            return new ThermodynamicResult(
                    results[0].temperature,
                    results[0].composition,
                    hInf,
                    Double.NaN,
                    hInf,
                    sigHInf,
                    cvInf,
                    null,
                    extrapCFs,
                    extrapStdCFs
            );
        }

        /**
         * Weighted linear fit: y = a + b*x, with weights w_i = 1/σ_i².
         *
         * @return {a, σ_a, b, σ_b}
         */
        static double[] weightedLinearFit(double[] x, double[] y, double[] sigma) {
            int n = x.length;
            double W = 0, Wx = 0, Wy = 0, Wxx = 0, Wxy = 0;
            for (int i = 0; i < n; i++) {
                double w = 1.0 / (sigma[i] * sigma[i]);
                W   += w;
                Wx  += w * x[i];
                Wy  += w * y[i];
                Wxx += w * x[i] * x[i];
                Wxy += w * x[i] * y[i];
            }
            double D = W * Wxx - Wx * Wx;
            if (Math.abs(D) < 1e-100) return new double[]{y[n - 1], 0.0, 0.0, 0.0};
            double a    = (Wxx * Wy  - Wx  * Wxy) / D;
            double b    = (W   * Wxy - Wx  * Wy ) / D;
            double sigA = Math.sqrt(Wxx / D);
            double sigB = Math.sqrt(W   / D);
            return new double[]{a, sigA, b, sigB};
        }
    }
}
