package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.ModelSession;
import org.ce.model.PhysicsConstants;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.mcs.MetropolisMC.MCSUpdate;
import org.ce.model.mcs.MetropolisMC.MCResult;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Persistent MCS model for a fixed (system, temperature) point.
 *
 * <p>Built once from a {@link ModelSession} at a given temperature. Pre-computes
 * ECIs (CVCF → orthogonal basis transform), embeddings, and orbit structure.
 * Multiple calls to {@link #run} then sweep different compositions or seeds on
 * the same pre-built geometry — no re-allocation of cluster data between runs.</p>
 *
 * <p>When temperature changes a new {@code MCSRunner} must be built (ECIs depend on T).
 * When only composition or sweep counts change, reuse the existing instance.</p>
 */
public class MCSRunner {

    private static final Logger LOG = Logger.getLogger(MCSRunner.class.getName());

    // ── Model-layer state (fixed for the lifetime of this instance) ───────────
    private final ClusCoordListData       clusterData;
    private final double[]                eciOrth;      // ECIs in orthogonal basis (Step 3 done at build)
    private final double[]                eciCvcf;      // ECIs in CVCF basis (for sampler / CF measurement)
    private final int                     numComp;
    private final double                  T;
    private final double                  R;
    private final int                     L;
    private final List<Vector3D>          positions;    // pre-built lattice positions
    private final Embeddings              emb;          // pre-built cluster embeddings
    private final List<List<Cluster>>     orbits;
    private final int[]                   orbitSizes;
    private final int[]                   multiSiteEmbedCounts;
    private final CvCfBasis               basis;
    private final List<List<Embeddings.Embedding>> cfEmbeddings;  // null if basis unavailable
    private final double[][]              basisMatrix;              // null if basis unavailable

    private MCSRunner(Builder b) {
        this.clusterData          = b.clusterData;
        this.eciCvcf              = b.eci.clone();
        this.numComp              = b.numComp;
        this.T                    = b.T;
        this.R                    = b.R;
        this.L                    = b.L;
        this.basis                = b.basis;

        // Step 3: transform ECIs to orthogonal basis (done once at build time)
        int tc = clusterData.getTc();
        this.eciOrth = buildEciByOrbitType(eciCvcf, tc, basis);

        // Step 4: build lattice geometry and embeddings (done once at build time)
        this.positions = (b.customPositions != null) ? b.customPositions : buildBCCPositions(L);
        this.emb       = Embeddings.generate(positions, clusterData, L);

        this.orbits    = clusterData.getOrbitList();
        this.orbitSizes = new int[tc];
        for (int t = 0; t < tc; t++) orbitSizes[t] = orbits.get(t).size();

        this.multiSiteEmbedCounts = emb.multiSiteEmbedCountsPerType(tc);

        if (basis != null && b.matrixData != null) {
            this.cfEmbeddings = Embeddings.generateCfEmbeddings(
                    emb.getAllEmbeddings(), clusterData, b.matrixData.getCfBasisIndices(), b.lcf);
            this.basisMatrix  = Embeddings.buildBasisValues(numComp);
        } else {
            this.cfEmbeddings = null;
            this.basisMatrix  = null;
        }

        int N = positions.size();
        LOG.info(String.format("MCSRunner built — L=%d, N=%d, T=%.1f K, numComp=%d", L, N, T, numComp));
        Debug.printEciInfo(eciCvcf, eciOrth, basis, N, orbitSizes);
    }

    /** Holds the MCResult and Sampler from one run, for post-processing by the calculation layer. */
    public static final class MCSRunResult {
        public final MCResult result;
        public final MetropolisMC.Sampler sampler;

        public MCSRunResult(MCResult result, MetropolisMC.Sampler sampler) {
            this.result  = result;
            this.sampler = sampler;
        }
    }

    /**
     * Runs one MC calculation at the model's fixed temperature for the given composition.
     *
     * <p>Safe to call multiple times on the same instance (e.g. different compositions
     * in a scan). Each call initializes a fresh {@link LatticeConfig} from {@code xFrac}
     * and runs independent equilibration + averaging sweeps.</p>
     *
     * <p>Step banners (Step 1–7) are written to {@code progressSink} at the start/end
     * of each major algorithm phase. Per-sweep progress is reported via
     * {@code updateListener} at the existing 100-sweep cadence.</p>
     *
     * @param xFrac             mole fractions (length == numComp, must sum to 1)
     * @param nEquil            equilibration sweeps
     * @param nAvg              averaging sweeps
     * @param seed              RNG seed (use {@code System.currentTimeMillis()} for random)
     * @param progressSink      optional text sink for step banners; may be {@code null}
     * @param updateListener    optional per-sweep callback; may be {@code null}
     * @param cancellationCheck optional interrupt check; may be {@code null}
     */
    public MCSRunResult run(
            double[] xFrac,
            int nEquil,
            int nAvg,
            long seed,
            Consumer<String> progressSink,
            Consumer<MCSUpdate> updateListener,
            BooleanSupplier cancellationCheck) {

        int N = positions.size();
        LOG.fine(String.format("MCSRunner.run — N=%d, T=%.1f K, nEquil=%d, nAvg=%d", N, T, nEquil, nAvg));

        // ── Step 1: INITIALIZATION ────────────────────────────────────────────
        emit(progressSink, String.format(
                "  Step 1 [INIT]   lattice: L=%d, N=%d sites, K=%d components, T=%.1f K, x=%s",
                L, N, numComp, T, Arrays.toString(xFrac)));

        Random rng = new Random(seed);
        LatticeConfig config = new LatticeConfig(N, numComp);
        config.randomise(xFrac, rng);
        emit(progressSink, String.format(
                "  Step 1 [INIT]   random occupation set (seed=%d)", seed));

        // ── Step 2: INITIAL ENERGY E(σ) — computed inside MetropolisMC.run() ──
        emit(progressSink, "  Step 2 [HAMILTONIAN]  computing E(σ₀) via cluster expansion...");

        MetropolisMC.Sampler sampler = new MetropolisMC.Sampler(
                N, orbitSizes, orbits, R, eciCvcf,
                multiSiteEmbedCounts, basis, cfEmbeddings, basisMatrix);

        MetropolisMC engine = new MetropolisMC(emb, eciOrth, orbits, numComp, T, nEquil, nAvg, R, rng);
        if (cancellationCheck != null) engine.setCancellationCheck(cancellationCheck);

        // ── Step 3: EQUILIBRATION — banner before sweep loop ─────────────────
        emit(progressSink, String.format(
                "  Step 3 [EQUILIBRATION]  %d sweeps × %d trial moves = %,d Metropolis steps",
                nEquil, N, (long) nEquil * N));

        // Wrap updateListener to inject phase-change banners at first avg sweep
        final boolean[] avgStarted = {false};
        Consumer<MCSUpdate> wrappedListener = mcUpdate -> {
            if (!avgStarted[0] && mcUpdate.getPhase() == MCSUpdate.Phase.AVERAGING) {
                avgStarted[0] = true;
                // Steps 4–6 are the per-sweep inner loop; banner at phase boundary
                emit(progressSink, String.format(
                        "  Step 4–6 [MOVE/ΔE/ACCEPT]  trial move → ΔE → Metropolis criterion (each sweep)"));
                emit(progressSink, String.format(
                        "  Step 7 [AVERAGING]  %d sweeps, recording ⟨E⟩ ⟨Hmix⟩ ⟨CF⟩ per sweep",
                        nAvg));
            }
            if (updateListener != null) updateListener.accept(mcUpdate);
        };
        engine.setUpdateListener(wrappedListener);

        MCResult result = engine.run(config, sampler);

        // ── Post-run summary ──────────────────────────────────────────────────
        emit(progressSink, String.format(
                "  Step 7 [DONE]   accept=%.1f%%  ⟨E⟩/site=%.6f J/mol  ⟨Hmix⟩/site=%.6f J/mol",
                result.getAcceptRate() * 100, result.getEnergyPerSite(), result.getHmixPerSite()));
        double[] cfs = result.getAvgCFs();
        if (cfs != null && cfs.length > 0) {
            StringBuilder sb = new StringBuilder("  Step 7 [CFs]    ⟨CF⟩ =");
            for (double cf : cfs) sb.append(String.format("  %+.5f", cf));
            emit(progressSink, sb.toString());
        }

        LOG.fine(String.format("MCSRunner.run — EXIT: acceptRate=%.3f, <E>/site=%.6f",
                result.getAcceptRate(), result.getEnergyPerSite()));
        return new MCSRunResult(result, sampler);
    }

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }

    /** Number of lattice sites (2·L³ for BCC). */
    public int nSites() { return positions.size(); }

    /** Temperature this model was built for. */
    public double temperature() { return T; }

    // ── ECI transform ─────────────────────────────────────────────────────────

    private static double[] buildEciByOrbitType(double[] eciCvcf, int tc, CvCfBasis basis) {
        if (basis == null || basis.Tinv == null) {
            LOG.warning("buildEciByOrbitType: Tinv unavailable — passing CVCF eci[] directly to MCEngine");
            double[] out = new double[tc];
            System.arraycopy(eciCvcf, 0, out, 0, Math.min(eciCvcf.length, tc));
            return out;
        }
        double[][] Tinv = basis.Tinv;
        int nCvcf = eciCvcf.length;
        double[] eciOrth = new double[tc];
        for (int t = 0; t < tc; t++) {
            double sum = 0.0;
            for (int l = 0; l < nCvcf && l < Tinv.length; l++)
                if (t < Tinv[l].length) sum += eciCvcf[l] * Tinv[l][t];
            eciOrth[t] = sum;
        }
        return eciOrth;
    }

    // ── Lattice geometry ──────────────────────────────────────────────────────

    public static List<Vector3D> buildBCCPositions(int L) {
        if (L < 1) throw new IllegalArgumentException("L must be >= 1");
        List<Vector3D> pos = new ArrayList<>(2 * L * L * L);
        for (int ix = 0; ix < L; ix++)
            for (int iy = 0; iy < L; iy++)
                for (int iz = 0; iz < L; iz++) {
                    pos.add(new Vector3D(ix,       iy,       iz      ));
                    pos.add(new Vector3D(ix + 0.5, iy + 0.5, iz + 0.5));
                }
        return pos;
    }

    // ── Diagnostics ───────────────────────────────────────────────────────────

    /** Diagnostic printing utilities. */
    public static final class Debug {

        private Debug() {}

        public static void printMcsSummary(MCResult result) {
            System.out.println("============================================================");
            System.out.println("  MCS SIMULATION OUTPUT");
            System.out.println("============================================================");
            System.out.printf("  Temperature   : %.2f K%n", result.getTemperature());
            System.out.printf("  Composition   : %s%n", java.util.Arrays.toString(result.getComposition()));
            System.out.printf("  Acceptance    : %.2f%%%n", result.getAcceptRate() * 100);
            System.out.printf("  Equilib       : %d sweeps%n", result.getNEquilSweeps());
            System.out.printf("  Averaging     : %d sweeps%n", result.getNAvgSweeps());
            System.out.println("------------------------------------------------------------");
            System.out.printf("  <E>/site      : %.6f J/mol%n", result.getEnergyPerSite());
            System.out.printf("  <H>/site      : %.6f J/mol%n", result.getHmixPerSite());
            System.out.println("  (Correlation functions and error bars computed by");
            System.out.println("   MCSStatisticsProcessor in calculation layer)");
            System.out.println("------------------------------------------------------------");

            double[] cfs = result.getAvgCFs();
            if (cfs != null && cfs.length > 0) {
                System.out.println("  Mean Correlation Functions (CVCF basis):");
                for (int i = 0; i < cfs.length; i++)
                    System.out.printf("    CF[%2d] = %+.6f%n", i, cfs[i]);
            } else {
                System.out.println("  No CVCF correlation functions measured.");
            }

            double[] seriesH = result.getSeriesHmix();
            if (seriesH != null)
                System.out.printf("  Raw time series available: Hmix (%d points), E (%d points)%n",
                        seriesH.length, result.getSeriesE() != null ? result.getSeriesE().length : 0);

            System.out.println("============================================================");
        }

        public static void printEciInfo(double[] eciCvcf, double[] eciOrth,
                                        CvCfBasis basis, int N, int[] orbitSizes) {
            System.out.println("============================================================");
            System.out.println("  MCS ECI & ORBIT DIAGNOSTIC");
            System.out.println("============================================================");
            if (basis == null) {
                System.out.println("  [ERROR] CvCfBasis is null; ECI names unknown.");
            } else {
                int ncf = basis.numNonPointCfs;
                System.out.printf("  CVCF Basis: %s (K=%d), numNonPointCfs=%d%n",
                        basis.structurePhase, basis.numComponents, ncf);
                System.out.println("  Non-zero CVCF ECIs:");
                boolean found = false;
                for (int l = 0; l < ncf; l++) {
                    if (Math.abs(eciCvcf[l]) > 1e-12) {
                        System.out.printf("    idx %2d: %-8s (= %-8s) : %+.6f J/mol%n",
                                l, basis.eciNames.get(l), basis.cfNames.get(l), eciCvcf[l]);
                        found = true;
                    }
                }
                if (!found) System.out.println("    [WARNING] All CVCF ECIs are ZERO!");
            }

            if (orbitSizes != null && N > 0) {
                System.out.println("  Orbit Multiplicities (from ClusterData):");
                for (int t = 0; t < orbitSizes.length; t++)
                    System.out.printf("    Orbit %2d: counts=%-6d mult=%6.2f%n",
                            t, orbitSizes[t], (double) orbitSizes[t] / N);
            }

            if (eciOrth != null) {
                System.out.println("  Transformed ECIs (Metropolis ECIs, orthogonal basis):");
                for (int t = 0; t < eciOrth.length; t++)
                    if (Math.abs(eciOrth[t]) > 1e-12)
                        System.out.printf("    Orbit %2d : %+.6f%n", t, eciOrth[t]);
            }
            System.out.println("============================================================");
        }
    }

    // ── Builder ───────────────────────────────────────────────────────────────

    public static Builder builder() { return new Builder(); }

    public static class Builder {
        private ClusCoordListData           clusterData;
        private double[]                    eci;
        private int                         numComp         = 2;
        private double                      T;
        private int                         L               = 12;
        private List<Vector3D>              customPositions = null;
        private double                      R               = PhysicsConstants.R_GAS;
        private CvCfBasis                   basis           = null;
        private CMatrixPipeline.CMatrixData matrixData      = null;
        private int[][]                     lcf             = null;

        private Builder() {}

        /**
         * Initialises all model fields from a {@link ModelSession} at the given temperature.
         *
         * <p>Steps performed here (model-layer concerns):</p>
         * <ul>
         *   <li>Step 2: {@link CECEvaluator#evaluate} → ECI array in CVCF basis</li>
         *   <li>Step 2b: C-matrix dimension validation</li>
         * </ul>
         * Steps 3 (ECI transform) and 4 (embeddings) are executed in the {@link MCSRunner}
         * constructor so they are shared across all subsequent {@link #run} calls.
         */
        public Builder session(ModelSession s, double temperature) {
            CMatrixPipeline.CMatrixData matData = s.clusterData.getMatrixData();

            // Step 2b: validate C-matrix columns match basis
            int cmatCols = (matData.getCmat().isEmpty()
                            || matData.getCmat().get(0).isEmpty()
                            || matData.getCmat().get(0).get(0).length == 0)
                           ? 0 : matData.getCmat().get(0).get(0)[0].length;
            if (cmatCols != s.cvcfBasis.totalCfs())
                throw new IllegalStateException("C-matrix dimension mismatch (cmatCols="
                        + cmatCols + ", basis.totalCfs=" + s.cvcfBasis.totalCfs() + ")");

            // Step 2: evaluate Hamiltonian at T — ECI in CVCF basis
            this.eci        = CECEvaluator.evaluate(s.cecEntry, temperature, s.cvcfBasis, "MCS");
            this.T          = temperature;
            this.clusterData = s.clusterData.getDisorderedClusterResult().getDisClusterData();
            this.numComp    = s.numComponents();
            this.basis      = s.cvcfBasis;
            this.matrixData = matData;
            this.lcf        = s.clusterData.getLcf();
            return this;
        }

        // Fine-grained setters (kept for tests that build without a ModelSession)
        public Builder clusterData(ClusCoordListData d)         { this.clusterData = d;       return this; }
        public Builder eci(double[] e)                           { this.eci = e;               return this; }
        public Builder numComp(int n)                            { this.numComp = n;           return this; }
        public Builder T(double t)                               { this.T = t;                 return this; }
        public Builder L(int l)                                  { this.L = l;                 return this; }
        public Builder latticePositions(List<Vector3D> pos)      { this.customPositions = pos; return this; }
        public Builder R(double r)                               { this.R = r;                 return this; }
        public Builder basis(CvCfBasis b)                        { this.basis = b;             return this; }
        public Builder matrixData(CMatrixPipeline.CMatrixData d) { this.matrixData = d;        return this; }
        public Builder lcf(int[][] l)                            { this.lcf = l;               return this; }

        public MCSRunner build() {
            if (clusterData == null) throw new IllegalStateException("clusterData required");
            if (eci == null)         throw new IllegalStateException("eci required");
            if (T <= 0)              throw new IllegalStateException("T must be > 0");
            if (basis == null)
                LOG.warning("MCSRunner built without a CvCfBasis — CVCF CF measurement unavailable");
            return new MCSRunner(this);
        }

        private static final Logger LOG = Logger.getLogger(Builder.class.getName());
    }
}
