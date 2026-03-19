package org.ce.domain.engine.mcs;

import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.ClusCoordListResult;
import org.ce.domain.cluster.Vector3D;

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
    }

    public MCResult run() {
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

        double[] hmixCoeff           = emb.computeHmixCoeff(eci, tc);
        int[]    multiSiteEmbedCounts = emb.multiSiteEmbedCountsPerType(tc);
        MCSampler sampler = new MCSampler(N, orbitSizes, orbits, R, hmixCoeff, multiSiteEmbedCounts);

        MCEngine engine = new MCEngine(emb, eci, orbits, numComp, T, nEquil, nAvg, R, rng);
        if (updateListener    != null) engine.setUpdateListener(updateListener);
        if (cancellationCheck != null) engine.setCancellationCheck(cancellationCheck);

        MCResult result = engine.run(config, sampler);
        LOG.fine(String.format("MCSRunner.run — EXIT: acceptRate=%.3f, <E>/site=%.6f",
                result.getAcceptRate(), result.getEnergyPerSite()));
        return result;
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
        private int                 nEquil            = 1000;
        private int                 nAvg              = 2000;
        private int                 L                 = 4;
        private List<Vector3D>      customPositions   = null;
        private long                seed              = 0L;
        private double              R                 = 1.0;
        private Consumer<MCSUpdate> updateListener    = null;
        private BooleanSupplier     cancellationCheck = null;

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

        public MCSRunner build() {
            if (clusterData == null) throw new IllegalStateException("clusterData required");
            if (eci == null)         throw new IllegalStateException("eci required");
            if (T <= 0)              throw new IllegalStateException("T must be > 0");
            if (xFrac == null) {
                xFrac = new double[numComp];
                for (int c = 0; c < numComp; c++) xFrac[c] = 1.0 / numComp;
            }
            return new MCSRunner(this);
        }
    }
}
