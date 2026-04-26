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
 * Entry point for canonical-ensemble Metropolis Monte Carlo thermodynamic calculations.
 *
 * <p>Holds pre-computed ECIs and temperature for one (geometry, T) point.
 * Geometry is provided by {@link MCSGeometry} and can be reused across temperature points.
 *
 * <h2>Pipeline</h2>
 * <pre>
 *  Step 1  INITIALIZATION      MCSGeometry.build()             — BCC lattice, cluster embeddings,
 *                                                                 CVCF basis (temperature-independent)
 *          TEMPERATURE SETUP   MCSRunner.forTemperature()      — ECI evaluation, eciOrth pre-compute,
 *                                                                 CSR site-to-CF index
 *          CONFIG INIT         MCSRunner.run()                 — random occupation, LatticeConfig
 *
 *  Step 2  INITIAL ENERGY      Embeddings.totalEnergyCvcf()    — E(σ₀) = N × Σ eciCvcf[l] × vCvcf[l]
 *
 *  Step 3  EQUILIBRATION       MetropolisMC.run() equil loop   — nEquil sweeps × N moves, no sampling
 *
 *  Step 4  TRIAL MOVE          ExchangeStep.attempt()          — select site pair (i,j) with σᵢ ≠ σⱼ
 *                              ExchangeStep.runParallelSweep() — checkerboard parallel variant
 *
 *  Step 5  COMPUTE ΔE          Embeddings.deltaEExchangeCvcf() — local O(neighbors), zero-allocation
 *                                                                 via FlatEmbData + DeltaScratch
 *
 *  Step 6  METROPOLIS ACCEPT   ExchangeStep (inline)           — min(1, exp(−ΔE / R·T))
 *                                                                 incremental CF sum update on accept
 *
 *  Step 7  AVERAGING           MetropolisMC.run() avg loop     — nAvg sweeps × N moves
 *          MEASUREMENT         Sampler.sample()                — ⟨E⟩, ⟨Hmix⟩, ⟨CF⟩ per sweep
 *                                                                 via incremental runningCfSum
 *
 *  Result  MCResult                                            — immutable DTO with all averages
 * </pre>
 *
 * <h2>Key Data Flow</h2>
 * <pre>
 *  MCSGeometry  ──►  MCSRunner  ──►  MetropolisMC.run()
 *  (geometry)        (ECIs, T)        │
 *                                     ├─ ExchangeStep  ──►  Embeddings.deltaEExchangeCvcf()
 *                                     │                      (FlatEmbData, DeltaScratch)
 *                                     └─ Sampler       ──►  Embeddings.applyTinvTransform()
 *                                                            (runningCfSum → uOrth → vCvcf)
 * </pre>
 *
 * <h2>Parallelism</h2>
 * <p>For L ≥ 4, the supercell is partitioned into a 2×2×2 (or 4×4×4 for L ≥ 12)
 * checkerboard by {@link LatticeDecomposer}. Same-color blocks are updated in parallel
 * via {@code parallelStream()}. Each block uses a {@code ThreadLocal<DeltaScratch>}
 * and {@code ThreadLocalRandom} for zero-contention hot paths.
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
     * Step 1 (TEMPERATURE SETUP) — builds a runner for a specific temperature using an existing geometry.
     * Evaluates the Hamiltonian ECIs at T and pre-computes the eciOrth fast-path vector.
     */
    public static MCSRunner forTemperature(MCSGeometry geo, ModelSession session, double T, Consumer<String> progressSink) {
        // Step 1: ECI evaluation — Hamiltonian at T in CVCF basis
        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "MCS", progressSink);

        // Step 1: CSR site-to-CF index for O(neighbors) ΔE computation
        Embeddings.CsrSiteToCfIndex siteToCfIndex = Embeddings.buildSiteToCfIndex(geo.cfEmbeddings, geo.nSites());

        // Step 1: eciOrth[m] = Σ_l eciCvcf[l] × Tinv[l][m] — collapses matrix multiply + ECI dot
        // product into a single dot product in the hot loop (Steps 5–6).
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

    /**
     * Steps 1–7 (CONFIG INIT → RESULT) — runs the full MC simulation at fixed (T, xFrac).
     * Delegates Steps 2–7 to {@link MetropolisMC#run}.
     */
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
        // Composition is invariant in canonical MC — cache it once to avoid O(N) recount per sweep
        sampler.setFixedComposition(config.composition());

        // Compute spatial decomposition for parallel execution
        double rMax = LatticeDecomposer.computeRMax(geo.cfEmbeddings, geo.positions);
        int blocksPerDim = (geo.L >= 12) ? 4 : 2;
        LatticeDecomposer.DecomposedLattice dl = LatticeDecomposer.decompose(geo.positions, geo.L, rMax, blocksPerDim);

        MetropolisMC engine = new MetropolisMC(
                geo.cfEmbeddings, geo.basisMatrix, siteToCfIndex,
                ncf, eciCvcf, eciOrth, geo.basis,
                geo.numComp, T, nEquil, nAvg, R, rng, dl);
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
