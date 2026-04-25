package org.ce.model.mcs;

import org.ce.debug.MCSDebug;
import org.ce.model.ModelSession;
import org.ce.model.PhysicsConstants;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.mcs.MetropolisMC.MCSUpdate;
import org.ce.model.mcs.MetropolisMC.MCResult;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Persistent MCS model for a fixed (system, temperature) point.
 *
 * <p>Now split into {@link MCSGeometry} (lattice data) and {@link MCSRunner} (ECIs at T).
 * This allows reusing the expensive geometry across different temperature points.</p>
 */
public class MCSRunner {

    private static final Logger LOG = Logger.getLogger(MCSRunner.class.getName());

    private final MCSGeometry geo;
    private final double[]    eciCvcf;      // ECIs in CVCF basis
    private final double[]    eciOrth;      // Pre-computed: eciOrth[m] = Σ_l eci[l] × Tinv[l][m]
    private final double      T;
    private final double      R;
    private final Embeddings.CsrSiteToCfIndex siteToCfIndex;  // per-site CF embedding index (CSR format)

    private MCSRunner(MCSGeometry geo, double[] eciCvcf, double[] eciOrth,
                      Embeddings.CsrSiteToCfIndex siteToCfIndex, double T, double R) {
        this.geo = geo;
        this.eciCvcf = eciCvcf.clone();
        this.eciOrth = eciOrth;
        this.siteToCfIndex = siteToCfIndex;
        this.T = T;
        this.R = R;

        LOG.info(String.format("MCSRunner ready — T=%.1f K, numComp=%d", T, geo.numComp));
    }

    /**
     * Factory method to build a runner for a specific temperature using an existing geometry.
     */
    public static MCSRunner forTemperature(MCSGeometry geo, ModelSession session, double T, Consumer<String> progressSink) {
        // Evaluate Hamiltonian at T — ECI in CVCF basis
        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "MCS", progressSink);

        // Build site-to-cfEmbedding index for efficient ΔE computation
        Embeddings.CsrSiteToCfIndex siteToCfIndex = Embeddings.buildSiteToCfIndex(geo.cfEmbeddings, geo.nSites());

        // Pre-compute eciOrth: eciOrth[m] = Σ_l eci[l] × Tinv[l][m]
        // This collapses the Tinv matrix-vector multiply + ECI dot product
        // into a single dot product in the hot loop.
        double[] eciOrth = computeEciOrth(eciCvcf, geo.basis);

        // ── MCS-DBG: ECI vectors after Hamiltonian evaluation ──
        if (MCSDebug.ENABLED) {
            MCSDebug.separator("ECI EVALUATION at T=" + T + " K");
            if (geo.basis != null) {
                MCSDebug.log("ECI", "Basis: %s, K=%d, numNonPointCfs=%d",
                        geo.basis.structurePhase, geo.basis.numComponents, geo.basis.numNonPointCfs);
                MCSDebug.log("ECI", "Tinv: %s",
                        geo.basis.Tinv != null
                                ? (geo.basis.Tinv.length + "x" + geo.basis.Tinv[0].length)
                                : "NULL");
            }
            MCSDebug.vector("ECI", "eciCvcf (CVCF basis — used for E and ΔE)", eciCvcf);
            if (eciOrth != null) MCSDebug.vector("ECI", "eciOrth (pre-computed for fast ΔE)", eciOrth);
            MCSDebug.log("ECI", "NOTE: buildEciByOrbitType REMOVED — engine uses CVCF ECIs directly");
        }

        return new MCSRunner(geo, eciCvcf, eciOrth, siteToCfIndex, T, PhysicsConstants.R_GAS);
    }

    /**
     * Pre-computes eciOrth[m] = Σ_l eciCvcf[l] × Tinv[l][m].
     * This allows ΔE = N × Σ_m ΔuOrth[m] × eciOrth[m] — a single dot product
     * instead of Tinv matrix-vector multiply followed by ECI dot product.
     */
    private static double[] computeEciOrth(double[] eciCvcf, CvCfBasis basis) {
        if (basis == null || basis.Tinv == null) return null;
        double[][] Tinv = basis.Tinv;
        int ncf = basis.numNonPointCfs;
        int tCols = Tinv[0].length;
        double[] eciOrth = new double[tCols];
        for (int m = 0; m < tCols; m++) {
            double sum = 0.0;
            for (int l = 0; l < ncf && l < eciCvcf.length && l < Tinv.length; l++) {
                sum += eciCvcf[l] * Tinv[l][m];
            }
            eciOrth[m] = sum;
        }
        return eciOrth;
    }

    /** Holds the MCResult and Sampler from one run. */
    public static final class MCSRunResult {
        public final MCResult result;
        public final MetropolisMC.Sampler sampler;

        public MCSRunResult(MCResult result, MetropolisMC.Sampler sampler) {
            this.result  = result;
            this.sampler = sampler;
        }
    }

    public MCSRunResult run(
            double[] xFrac,
            int nEquil,
            int nAvg,
            long seed,
            Consumer<String> progressSink,
            Consumer<MCSUpdate> updateListener,
            BooleanSupplier cancellationCheck) {

        int N = geo.nSites();
        LOG.fine(String.format("MCSRunner.run — N=%d, T=%.1f K, nEquil=%d, nAvg=%d", N, T, nEquil, nAvg));

        emit(progressSink, String.format(
                "  Step 1 [INIT]   lattice: L=%d, N=%d sites, K=%d components, T=%.1f K, x=%s",
                geo.L, N, geo.numComp, T, Arrays.toString(xFrac)));

        Random rng = new Random(seed);
        LatticeConfig config = new LatticeConfig(N, geo.numComp);
        config.randomise(xFrac, rng);
        emit(progressSink, String.format("  Step 1 [INIT]   random occupation set (seed=%d)", seed));

        emit(progressSink, "  Step 2 [HAMILTONIAN]  computing E(σ₀) via CVCF cluster expansion...");

        int ncf = geo.basis != null ? geo.basis.numNonPointCfs : eciCvcf.length;

        MetropolisMC.Sampler sampler = new MetropolisMC.Sampler(
                N, geo.orbitSizes, geo.orbits, R, eciCvcf,
                geo.multiSiteEmbedCounts, geo.basis, geo.cfEmbeddings, geo.basisMatrix);

        MetropolisMC engine = new MetropolisMC(
                geo.cfEmbeddings, geo.basisMatrix, siteToCfIndex,
                ncf, eciCvcf, eciOrth, geo.basis,
                geo.numComp, T, nEquil, nAvg, R, rng);
        if (cancellationCheck != null) engine.setCancellationCheck(cancellationCheck);

        emit(progressSink, String.format(
                "  Step 3 [EQUILIBRATION]  %d sweeps × %d trial moves = %,d Metropolis steps",
                nEquil, N, (long) nEquil * N));

        final boolean[] avgStarted = {false};
        Consumer<MCSUpdate> wrappedListener = mcUpdate -> {
            if (!avgStarted[0] && mcUpdate.getPhase() == MCSUpdate.Phase.AVERAGING) {
                avgStarted[0] = true;
                emit(progressSink, "  Step 4–6 [MOVE/ΔE/ACCEPT]  trial move → ΔE → Metropolis criterion (each sweep)");
                emit(progressSink, String.format("  Step 7 [AVERAGING]  %d sweeps, recording ⟨E⟩ ⟨Hmix⟩ ⟨CF⟩ per sweep", nAvg));
            }
            if (updateListener != null) updateListener.accept(mcUpdate);
        };
        engine.setUpdateListener(wrappedListener);

        MCResult result = engine.run(config, sampler);

        emit(progressSink, String.format(
                "  Step 7 [DONE]   accept=%.1f%%  ⟨E⟩/site=%.6f J/mol  ⟨Hmix⟩/site=%.6f J/mol",
                result.getAcceptRate() * 100, result.getEnergyPerSite(), result.getHmixPerSite()));
        
        double[] cfs = result.getAvgCFs();
        if (cfs != null && cfs.length > 0) {
            StringBuilder sb = new StringBuilder("  Step 7 [CFs]    ⟨CF⟩ =");
            for (double cf : cfs) sb.append(String.format("  %+.5f", cf));
            emit(progressSink, sb.toString());
        }

        return new MCSRunResult(result, sampler);
    }

    public int nSites() { return geo.nSites(); }
    public double temperature() { return T; }

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }

    public static final class Debug {
        private Debug() {}
        public static void printMcsSummary(MCResult result) {
            System.out.println("============================================================");
            System.out.println("  MCS SIMULATION OUTPUT");
            System.out.println("============================================================");
            System.out.printf("  Temperature   : %.2f K%n", result.getTemperature());
            System.out.printf("  Composition   : %s%n", Arrays.toString(result.getComposition()));
            System.out.printf("  Acceptance    : %.2f%%%n", result.getAcceptRate() * 100);
            System.out.printf("  Equilib       : %d sweeps%n", result.getNEquilSweeps());
            System.out.printf("  Averaging     : %d sweeps%n", result.getNAvgSweeps());
            System.out.println("------------------------------------------------------------");
            System.out.printf("  <E>/site      : %.6f J/mol%n", result.getEnergyPerSite());
            System.out.printf("  <H>/site      : %.6f J/mol%n", result.getHmixPerSite());
            System.out.println("============================================================");
        }
    }
}
