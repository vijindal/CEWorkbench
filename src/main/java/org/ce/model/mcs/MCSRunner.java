package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.mcs.EmbeddingData;
import org.ce.model.mcs.EmbeddingGenerator;
import org.ce.model.mcs.LatticeConfig;
import org.ce.model.mcs.CvCfEvaluator;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.CMatrix;
import org.ce.model.cluster.ClusterResults.ClusCoordListResult;
import org.ce.model.cluster.cvcf.CvCfBasis;
import org.ce.model.mcs.MCSUpdate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.logging.Logger;

/** Top-level orchestrator for the MCS engine path. */
public class MCSRunner {

    private static final Logger LOG = Logger.getLogger(MCSRunner.class.getName());

    private final ClusCoordListResult clusterData;
    private final double[]            eci;
    private final int                 numComp;
    private final double              T;
    private final double[]            xFrac;
    private final double              R;
    private final int                 nEquil;
    private final int                 nAvg;
    private final int                 L;
    private final List<Vector3D>      customPositions;
    private final long                seed;
    private final Consumer<MCSUpdate> updateListener;
    private final BooleanSupplier     cancellationCheck;
    private final CvCfBasis           basis;
    private final CMatrix.Result       cmatResult;

    private MCSRunner(Builder b) {
        this.clusterData       = b.clusterData;
        this.eci               = b.eci.clone();
        this.numComp           = b.numComp;
        this.T                 = b.T;
        this.xFrac             = b.xFrac.clone();
        this.nEquil            = b.nEquil;
        this.nAvg              = b.nAvg;
        this.L                 = b.L;
        this.customPositions   = b.customPositions;
        this.seed              = b.seed;
        this.R                 = b.R;
        this.updateListener    = b.updateListener;
        this.cancellationCheck = b.cancellationCheck;
        this.basis             = b.basis;
        this.cmatResult        = b.cmatResult;
    }

    /**
     * Holds the result and sampler from an MCS run.
     * The sampler contains raw time series for post-processing by the calculation layer.
     */
    public static final class MCSRunResult {
        public final MCResult result;
        public final MCSampler sampler;

        public MCSRunResult(MCResult result, MCSampler sampler) {
            this.result = result;
            this.sampler = sampler;
        }
    }

    public MCSRunResult run() {
        Random rng = new Random(seed);

        List<Vector3D> positions = (customPositions != null) ? customPositions : buildBCCPositions(L);
        int N = positions.size();
        LOG.fine(String.format("MCSRunner.run — L=%d, N=%d, T=%.1f K, nEquil=%d, nAvg=%d", L, N, T, nEquil, nAvg));

        EmbeddingData emb = EmbeddingGenerator.generateEmbeddings(positions, clusterData, L);

        LatticeConfig config = new LatticeConfig(N, numComp);
        config.randomise(xFrac, rng);

        int tc = clusterData.getTc();
        int[] orbitSizes = new int[tc];
        List<List<Cluster>> orbits = clusterData.getOrbitList();
        for (int t = 0; t < tc; t++) orbitSizes[t] = orbits.get(t).size();

        int[]    multiSiteEmbedCounts = emb.multiSiteEmbedCountsPerType(tc);

        // New Direct Measurement Structures
        List<List<EmbeddingData.Embedding>> cfEmbeddings = null;
        double[][] basisMatrix = null;
        if (basis != null && cmatResult != null) {
            cfEmbeddings = EmbeddingGenerator.generateCfEmbeddings(emb.getAllEmbeddings(), clusterData, cmatResult.getCfBasisIndices());
            basisMatrix  = CvCfEvaluator.buildBasisValues(numComp);
        }

        MCSampler sampler = new MCSampler(N, orbitSizes, orbits, R, eci, multiSiteEmbedCounts, basis,
                                         cfEmbeddings, basisMatrix);

        // Build orbit-type-indexed ECI array for MCEngine/LocalEnergyCalc.
        // eci[] is CVCF-ordered (eci[l] = ECI for basis.cfNames[l]).
        // MCEngine uses eci[clusterType] — indexed by orbit type t.
        // Transform: eciOrth[t] = Σ_l  eci_cvcf[l] * Tinv[l][t]
        // This is the effective orthogonal-basis ECI for orbit type t, assuming
        // orbit type t corresponds to T-matrix row t (confirmed for binary BCC_A2).
        double[] eciOrth = buildEciByOrbitType(eci, tc, basis);

        // [DEBUG] Diagnostic print of Hamiltonian terms and orbit multiplicities
        MCResult.Debug.printEciInfo(eci, eciOrth, basis, N, orbitSizes);

        MCEngine engine = new MCEngine(emb, eciOrth, orbits, numComp, T, nEquil, nAvg, R, rng);
        if (updateListener    != null) engine.setUpdateListener(updateListener);
        if (cancellationCheck != null) engine.setCancellationCheck(cancellationCheck);

        MCResult result = engine.run(config, sampler);
        LOG.fine(String.format("MCSRunner.run — EXIT: acceptRate=%.3f, <E>/site=%.6f",
                result.getAcceptRate(), result.getEnergyPerSite()));
        return new MCSRunResult(result, sampler);
    }

    /**
     * Converts a CVCF-ordered ECI array (indexed by CF position l) into an
     * orbit-type-indexed ECI array (indexed by orbit type t) suitable for
     * {@link LocalEnergyCalc}.
     *
     * <p>Uses the Tinv transform: {@code eciOrth[t] = Σ_l eci_cvcf[l] * Tinv[l][t]}.
     * This is the effective orthogonal-basis ECI for orbit type t, given that orbit
     * type t corresponds to T-matrix row t.</p>
     *
     * <p>When {@code basis} or {@code basis.Tinv} is null (ternary/quaternary without
     * Tinv), falls back to passing the CVCF eci[] directly with a WARNING — this will
     * be wrong for ternary/quaternary but avoids a crash.</p>
     */
    private static double[] buildEciByOrbitType(double[] eciCvcf, int tc, CvCfBasis basis) {
        if (basis == null || basis.Tinv == null) {
            LOG.warning("buildEciByOrbitType: Tinv unavailable — passing CVCF eci[] directly to MCEngine "
                    + "(ECI indexing will be wrong for multi-component systems)");
            // Pad or trim to tc length
            double[] out = new double[tc];
            System.arraycopy(eciCvcf, 0, out, 0, Math.min(eciCvcf.length, tc));
            return out;
        }
        double[][] Tinv = basis.Tinv;   // [numCFs][numRows(T)] = [totalCfs][tc]
        int nCvcf = eciCvcf.length;     // = numNonPointCfs
        double[] eciOrth = new double[tc];
        for (int t = 0; t < tc; t++) {
            double sum = 0.0;
            for (int l = 0; l < nCvcf && l < Tinv.length; l++) {
                if (t < Tinv[l].length) sum += eciCvcf[l] * Tinv[l][t];
            }
            eciOrth[t] = sum;
        }
        return eciOrth;
    }

    public static List<Vector3D> buildBCCPositions(int L) {
        if (L < 1) throw new IllegalArgumentException("L must be >= 1");
        List<Vector3D> pos = new ArrayList<>(2 * L * L * L);
        for (int ix = 0; ix < L; ix++)
            for (int iy = 0; iy < L; iy++)
                for (int iz = 0; iz < L; iz++) {
                    pos.add(new Vector3D(ix,         iy,         iz        ));
                    pos.add(new Vector3D(ix + 0.5,   iy + 0.5,   iz + 0.5  ));
                }
        return pos;
    }

    public static Builder builder() { return new Builder(); }

    public static class Builder {
        private ClusCoordListResult clusterData;
        private double[]            eci;
        private int                 numComp           = 2;
        private double              T;
        private double[]            xFrac;
        private int                 nEquil            = 500;
        private int                 nAvg              = 500;
        private int                 L                 = 12;
        private List<Vector3D>      customPositions   = null;
        private long                seed              = 0L;
        private double              R                 = 1.0;
        private Consumer<MCSUpdate> updateListener    = null;
        private BooleanSupplier     cancellationCheck = null;
        private CvCfBasis           basis             = null;
        private CMatrix.Result       cmatResult        = null;

        private Builder() {}

        public Builder clusterData(ClusCoordListResult d)       { this.clusterData = d;        return this; }
        public Builder eci(double[] e)                          { this.eci = e;                return this; }
        public Builder numComp(int n)                           { this.numComp = n;            return this; }
        public Builder T(double t)                              { this.T = t;                  return this; }
        public Builder composition(double[] x)                  { this.xFrac = x.clone();      return this; }
        public Builder compositionBinary(double xB)             { this.xFrac = new double[]{1.0 - xB, xB}; return this; }
        public Builder nEquil(int n)                            { this.nEquil = n;             return this; }
        public Builder nAvg(int n)                              { this.nAvg = n;               return this; }
        public Builder L(int l)                                 { this.L = l;                  return this; }
        public Builder latticePositions(List<Vector3D> pos)     { this.customPositions = pos;  return this; }
        public Builder seed(long s)                             { this.seed = s;               return this; }
        public Builder R(double r)                              { this.R = r;                  return this; }
        public Builder updateListener(Consumer<MCSUpdate> l)    { this.updateListener = l;     return this; }
        public Builder cancellationCheck(BooleanSupplier check) { this.cancellationCheck = check; return this; }
        public Builder basis(CvCfBasis b)                       { this.basis = b;              return this; }
        public Builder cmatResult(CMatrix.Result r)              { this.cmatResult = r;         return this; }

        public MCSRunner build() {
            if (clusterData == null) throw new IllegalStateException("clusterData required");
            if (eci == null)         throw new IllegalStateException("eci required");
            if (T <= 0)              throw new IllegalStateException("T must be > 0");
            if (xFrac == null) {
                xFrac = new double[numComp];
                for (int c = 0; c < numComp; c++) xFrac[c] = 1.0 / numComp;
            }
            if (basis == null) {
                LOG.warning("MCSRunner built without a CvCfBasis — "
                        + "CVCF CF measurement will be unavailable for this run");
            }
            return new MCSRunner(this);
        }
    }
}
