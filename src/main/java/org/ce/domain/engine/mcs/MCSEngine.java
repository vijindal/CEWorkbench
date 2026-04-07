package org.ce.domain.engine.mcs;

import static org.ce.domain.cluster.AllClusterData.ClusterData;

import static org.ce.domain.cluster.ClusterResults.*;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import org.ce.domain.cluster.CMatrix;
import org.ce.domain.engine.ProgressEvent;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.engine.mcs.MCResult;
import org.ce.domain.engine.mcs.MCSRunner;
import org.ce.domain.engine.mcs.MCSUpdate;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.domain.hamiltonian.CECEvaluator;
import org.ce.domain.result.EquilibriumState;

import java.util.logging.Logger;

import java.util.function.Consumer;

/**
 * Thermodynamic engine implementing Monte Carlo simulation (canonical ensemble).
 *
 * <p>Default parameters: L=4 supercell (128 BCC sites), 1000 equilibration sweeps,
 * 2000 averaging sweeps. ECI = a + b*T from CECEntry.</p>
 */
public class MCSEngine implements ThermodynamicEngine {

    private static final double GAS_CONSTANT      = 8.314;  // J/(mol·K)

    private static final Logger LOG = Logger.getLogger(MCSEngine.class.getName());

    @Override
    public EquilibriumState compute(ThermodynamicInput input) throws Exception {

        ClusCoordListResult clusterData =
                input.clusterData.getDisorderedClusterResult().getDisClusterData();

        // Validate C-matrix dimensions match basis (same check as CVMEngine)
        String structurePhase = input.cec.structurePhase;
        String model = input.cec.model;
        int numComponents = input.composition.length;
        var basis = CvCfBasis.Registry.INSTANCE.get(structurePhase, model, numComponents);

        CMatrix.Result cmatResult = input.clusterData.getCMatrixResult();
        cmatResult.validateCols(
                basis.totalCfs(),
                "C-matrix dimension mismatch (basis.numNonPointCfs=" + basis.numNonPointCfs
                + " + " + basis.numComponents + " point variables)"
        );
        LOG.fine("✓ C-matrix dimensions valid: " + basis.totalCfs() + " columns");

        // 4. Evaluate ECI at temperature using unified CECEvaluator.
        // Extraction is logged at INFO level to assist in debugging scaling/basis issues.
        double[] eci = CECEvaluator.evaluate(input.cec, input.temperature, basis, "MCS");

        int L           = input.mcsL;
        int nEquil      = input.mcsNEquil;
        int nAvg        = input.mcsNAvg;
        int totalSweeps = nEquil + nAvg;
        int N           = 2 * L * L * L;  // BCC sites

        MCSRunner.Builder builder = MCSRunner.builder()
                .clusterData(clusterData)
                .eci(eci)
                .numComp(input.composition.length)
                .T(input.temperature)
                .composition(input.composition)
                .nEquil(nEquil)
                .nAvg(nAvg)
                .L(L)
                .seed(System.currentTimeMillis())
                .R(GAS_CONSTANT)
                .basis(basis)
                .cmatResult(cmatResult)
                .cancellationCheck(Thread.currentThread()::isInterrupted);

        Consumer<String> strSink = input.progressSink;
        Consumer<ProgressEvent> evtSink = input.eventSink;

        if (strSink != null || evtSink != null) {
            if (strSink != null) {
                strSink.accept(String.format("  MCS started: L=%d (%d sites), %d equil + %d avg sweeps",
                        L, N, nEquil, nAvg));
            }
            // EngineStart emitted here — before build().run() — so the chart clears before sweeps arrive
            if (evtSink != null) {
                evtSink.accept(new ProgressEvent.EngineStart("MCS", totalSweeps));
            }
            final int sitesCount = N;
            final int sweepCount = totalSweeps;
            builder.updateListener(mcUpdate -> {
                if (strSink != null && (mcUpdate.getStep() % 100 == 0 || mcUpdate.getStep() == sweepCount)) {
                    strSink.accept(String.format("  [%-5s] sweep %4d/%-4d  <E>/site=%9.5f  accept=%5.1f%%",
                            mcUpdate.getPhase(),
                            mcUpdate.getStep(), sweepCount,
                            mcUpdate.getE_total() / sitesCount,
                            mcUpdate.getAcceptanceRate() * 100));
                }
                if (evtSink != null) {
                    evtSink.accept(new ProgressEvent.McSweep(
                            mcUpdate.getStep(), sweepCount,
                            mcUpdate.getE_total() / sitesCount,
                            mcUpdate.getAcceptanceRate(),
                            mcUpdate.getPhase() == MCSUpdate.Phase.EQUILIBRATION,
                            null));  // CFs no longer charted
                }
            });
        }

        MCResult result = builder.build().run();

        // [DEBUG] Temporary post-run diagnostic print
        mcsDebugData.printMcsSummary(result);

        // Post-run statistics summary
        if (strSink != null) {
            strSink.accept(
                "  ── MCS Statistics ─────────────────────────────────────────");
            strSink.accept(String.format(
                "  τ_int = %.1f sweeps  s = %.3f  n_eff = %d  (block size = %d, %d blocks)",
                result.getTauInt(), result.getStatInefficiency(),
                result.getNEff(), result.getBlockSizeUsed(), result.getNBlocks()));
            strSink.accept(String.format(
                "  ⟨H⟩/site = %.6f ± %.6f  J/mol",
                result.getHmixPerSite(), nanZero(result.getStdHmixPerSite())));
            strSink.accept(String.format(
                "  ⟨E⟩/site = %.6f ± %.6f  J/mol",
                result.getEnergyPerSite(), nanZero(result.getStdEnergyPerSite())));
            strSink.accept(String.format(
                "  Cv       = %.6f ± %.6f  J/(mol·K)  [jackknife]",
                result.getCvJackknife(), nanZero(result.getCvStdErr())));

            // Correlation functions
            double[] cfs    = result.getAvgCFs();
            double[] stdCFs = result.getStdCFs();
            if (cfs != null && cfs.length > 0) {
                strSink.accept(
                    "  ── Correlation Functions ────────────────────────────────");
                for (int i = 0; i < cfs.length; i++) {
                    double sigma = (stdCFs != null && i < stdCFs.length) ? stdCFs[i] : Double.NaN;
                    strSink.accept(String.format(
                        "  CF[%2d] = %+.6f ± %.6f",
                        i, cfs[i], nanZero(sigma)));
                }
            }

            // Parameter recommendations
            emitParameterRecommendations(strSink, result, L, nEquil, nAvg);
        }

        return new EquilibriumState(
                result.getTemperature(),
                result.getComposition(),
                result.getHmixPerSite(),
                Double.NaN,                 // freeEnergy — not available from canonical MCS
                Double.NaN,                 // entropy — not available from canonical MCS
                result.getStdHmixPerSite(), // enthalpy std error from MCS
                result.getCvJackknife(),    // heat capacity from jackknife
                null,                       // optimizedCFs (CVM only)
                result.getAvgCFs(),         // mean correlation functions
                result.getStdCFs()          // CF standard errors from block averaging
        );
    }

    /**
     * Emits post-run parameter recommendations to the log sink.
     *
     * <p>Three tiers of advice:</p>
     * <ul>
     *   <li><b>nAvg</b> — exact: derived from τ_int to reach n_eff ≥ 100</li>
     *   <li><b>nEquil</b> — heuristic: 20·τ_int ensures decorrelation before averaging</li>
     *   <li><b>L</b> — heuristic: fixed floor + τ_int proxy for correlation length</li>
     * </ul>
     */
    private static void emitParameterRecommendations(
            Consumer<String> sink, MCResult result, int L, int nEquil, int nAvg) {

        double tauInt = result.getTauInt();
        double s      = result.getStatInefficiency();
        int    nEff   = result.getNEff();

        boolean anyIssue = false;

        // ── nAvg ────────────────────────────────────────────────────────────
        if (nEff < 100) {
            int nAvgMin = (int) Math.ceil(100.0 * s);
            sink.accept(String.format(
                "  *** WARNING: n_eff = %d < 100 — error bars are unreliable", nEff));
            sink.accept(String.format(
                "      → set nAvg ≥ %d  (= ⌈100 × s⌉, guarantees n_eff ≥ 100)", nAvgMin));
            anyIssue = true;
        } else if (nEff < 200) {
            sink.accept(String.format(
                "  ⚠  n_eff = %d — marginal; consider nAvg ≥ %d for publication",
                nEff, (int) Math.ceil(200.0 * s)));
            anyIssue = true;
        }

        // ── nEquil ──────────────────────────────────────────────────────────
        int nEquilMin = Math.max(200, (int) Math.ceil(20.0 * tauInt));
        if (nEquil < nEquilMin) {
            sink.accept(String.format(
                "  *** WARNING: nEquil = %d may be too short (τ_int = %.1f sweeps)",
                nEquil, tauInt));
            sink.accept(String.format(
                "      → set nEquil ≥ %d  (= max(200, ⌈20 × τ_int⌉))", nEquilMin));
            anyIssue = true;
        }

        // ── L ───────────────────────────────────────────────────────────────
        if (L < 6) {
            sink.accept(String.format(
                "  ⚠  L = %d (%d sites) — minimum for any meaningful result is L = 6 (432 sites)",
                L, 2 * L * L * L));
            anyIssue = true;
        } else if (L < 10) {
            sink.accept(String.format(
                "  ⚠  L = %d (%d sites) — small supercell; thermodynamic finite-size effects are likely",
                L, 2 * L * L * L));
            sink.accept(
                "      → verify by comparing ⟨H⟩/site at L and L+4 (e.g. L=" + (L + 4) + ", " + (2*(L+4)*(L+4)*(L+4)) + " sites)");
            anyIssue = true;
        }
        if (tauInt > 50) {
            sink.accept(String.format(
                "  ⚠  τ_int = %.1f sweeps is large — correlation length may approach box size (L = %d)",
                tauInt, L));
            anyIssue = true;
        }

        // Size convergence check is always recommended — finite-size effects are
        // thermodynamic and cannot be detected from τ_int or n_eff alone
        sink.accept(String.format(
            "  ℹ  Size check: compare ⟨H⟩/site at L=%d and L=%d to confirm thermodynamic limit",
            L, L + 4));

        if (!anyIssue) {
            sink.accept("  ✓  Statistical parameters look sufficient (n_eff ≥ 100, nEquil ≥ 20·τ_int)");
        }
    }

    /** Helper: returns 0 if value is NaN, otherwise returns the value. */
    private static double nanZero(double value) {
        return Double.isNaN(value) ? 0.0 : value;
    }
}
