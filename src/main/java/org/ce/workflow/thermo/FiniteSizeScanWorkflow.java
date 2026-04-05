package org.ce.workflow.thermo;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import org.ce.domain.engine.ProgressEvent;
import org.ce.domain.result.ThermodynamicResult;

import java.util.function.Consumer;

/**
 * Finite-size scaling workflow for production-grade MCS results.
 *
 * <p>Runs three MCS simulations at L=12, L=16, L=24 (N=3456, 8192, 27648 BCC sites),
 * fits each observable to the 1/N finite-size scaling model:</p>
 *
 * <pre>  O(N) = O(∞) + A/N</pre>
 *
 * <p>and reports the extrapolated thermodynamic-limit values with uncertainties
 * for ⟨H⟩/site, ⟨E⟩/site, Cv, and all correlation functions.</p>
 *
 * <p><b>BCC L-value rationale:</b> L must be even for B2 supercell compatibility
 * (2-atom primitive cell). L=12/16/24 gives an 8× range in N, ensuring the
 * 1/N fit is well-conditioned. L=24 (27648 sites) alone has sub-percent
 * finite-size error for most BCC alloys away from phase transitions.</p>
 */
public class FiniteSizeScanWorkflow {

    /** Default BCC supercell sizes: even, B2-compatible, 8× ratio in N. */
    public static final int[] DEFAULT_L_VALUES = {12, 16, 24};

    private final ThermodynamicWorkflow thermoWorkflow;

    public FiniteSizeScanWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }

    /**
     * Runs MCS at three L values, extrapolates to thermodynamic limit.
     *
     * @return ThermodynamicResult whose enthalpy/stdEnthalpy/heatCapacity/avgCFs/stdCFs
     *         hold the extrapolated (∞-limit) values and their uncertainties.
     */
    public ThermodynamicResult run(
            String clusterId,
            String hamiltonianId,
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
                clusterId, hamiltonianId, temperature, composition,
                "MCS", strSink, evtSink, L, nEquil, nAvg);

            results[i] = thermoWorkflow.runCalculation(req);
        }

        return extrapolateAndReport(lValues, N, invN, results, strSink);
    }

    // =========================================================================
    // Extrapolation
    // =========================================================================

    private ThermodynamicResult extrapolateAndReport(
            int[] lValues, double[] N, double[] invN,
            ThermodynamicResult[] results,
            Consumer<String> sink) {

        int nRuns = results.length;

        // ── Collect per-L observables ─────────────────────────────────────────
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
            sigCv[i] = 0.1 * Math.max(1e-8, Cv[i]);   // rough 10% if stdCv not available
        }

        // ── Fit ⟨H⟩ vs 1/N ───────────────────────────────────────────────────
        double[] fitH  = weightedLinearFit(invN, H, sigH);
        double hInf    = fitH[0];
        double sigHInf = fitH[1];
        double hSlope  = fitH[2];

        // ── Fit Cv vs 1/N ────────────────────────────────────────────────────
        double[] fitCv = weightedLinearFit(invN, Cv, sigCv);
        double cvInf   = fitCv[0];

        // ── Summary table ─────────────────────────────────────────────────────
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

        // ── CFs ──────────────────────────────────────────────────────────────
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
                hInf,           // gibbsEnergy (H proxy — G not directly available from MCS)
                Double.NaN,     // entropy
                hInf,           // enthalpy = extrapolated ⟨H⟩/site
                sigHInf,        // stdEnthalpy = σ from 1/N fit
                cvInf,          // heatCapacity = extrapolated Cv
                null,           // optimizedCFs
                extrapCFs,
                extrapStdCFs
        );
    }

    // =========================================================================
    // Statistics
    // =========================================================================

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
