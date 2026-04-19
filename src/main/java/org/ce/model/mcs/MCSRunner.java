package org.ce.model.mcs;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.mcs.MetropolisMC.MCSUpdate;
import org.ce.model.mcs.MetropolisMC.MCResult;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.logging.Logger;

/** Top-level orchestrator for the MCS engine path. */
public class MCSRunner {

    private static final Logger LOG = Logger.getLogger(MCSRunner.class.getName());

    private final ClusCoordListData clusterData;
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
    private final CMatrixPipeline.CMatrixData matrixData;
    private final int[][]              lcf;

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
        this.matrixData        = b.matrixData;
        this.lcf               = b.lcf;
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

    public MCSRunResult run() {
        Random rng = new Random(seed);

        List<Vector3D> positions = (customPositions != null) ? customPositions : buildBCCPositions(L);
        int N = positions.size();
        LOG.fine(String.format("MCSRunner.run — L=%d, N=%d, T=%.1f K, nEquil=%d, nAvg=%d", L, N, T, nEquil, nAvg));

        Embeddings emb = Embeddings.generate(positions, clusterData, L);

        LatticeConfig config = new LatticeConfig(N, numComp);
        config.randomise(xFrac, rng);

        int tc = clusterData.getTc();
        int[] orbitSizes = new int[tc];
        List<List<Cluster>> orbits = clusterData.getOrbitList();
        for (int t = 0; t < tc; t++) orbitSizes[t] = orbits.get(t).size();

        int[] multiSiteEmbedCounts = emb.multiSiteEmbedCountsPerType(tc);

        List<List<Embeddings.Embedding>> cfEmbeddings = null;
        double[][] basisMatrix = null;
        if (basis != null && matrixData != null) {
            cfEmbeddings = Embeddings.generateCfEmbeddings(emb.getAllEmbeddings(), clusterData, matrixData.getCfBasisIndices(), lcf);
            basisMatrix  = Embeddings.buildBasisValues(numComp);
        }

        MetropolisMC.Sampler sampler = new MetropolisMC.Sampler(N, orbitSizes, orbits, R, eci,
                multiSiteEmbedCounts, basis, cfEmbeddings, basisMatrix);

        double[] eciOrth = buildEciByOrbitType(eci, tc, basis);

        MCSRunner.Debug.printEciInfo(eci, eciOrth, basis, N, orbitSizes);

        MetropolisMC engine = new MetropolisMC(emb, eciOrth, orbits, numComp, T, nEquil, nAvg, R, rng);
        if (updateListener    != null) engine.setUpdateListener(updateListener);
        if (cancellationCheck != null) engine.setCancellationCheck(cancellationCheck);

        MCResult result = engine.run(config, sampler);
        LOG.fine(String.format("MCSRunner.run — EXIT: acceptRate=%.3f, <E>/site=%.6f",
                result.getAcceptRate(), result.getEnergyPerSite()));
        return new MCSRunResult(result, sampler);
    }

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

    public static Builder builder() { return new Builder(); }

    public static class Builder {
        private ClusCoordListData clusterData;
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
        private CMatrixPipeline.CMatrixData matrixData = null;
        private int[][]             lcf               = null;

        private Builder() {}

        public Builder clusterData(ClusCoordListData d)        { this.clusterData = d;       return this; }
        public Builder eci(double[] e)                          { this.eci = e;               return this; }
        public Builder numComp(int n)                           { this.numComp = n;           return this; }
        public Builder T(double t)                              { this.T = t;                 return this; }
        public Builder composition(double[] x)                  { this.xFrac = x.clone();     return this; }
        public Builder compositionBinary(double xB)             { this.xFrac = new double[]{1.0 - xB, xB}; return this; }
        public Builder nEquil(int n)                            { this.nEquil = n;            return this; }
        public Builder nAvg(int n)                              { this.nAvg = n;              return this; }
        public Builder L(int l)                                 { this.L = l;                 return this; }
        public Builder latticePositions(List<Vector3D> pos)     { this.customPositions = pos; return this; }
        public Builder seed(long s)                             { this.seed = s;              return this; }
        public Builder R(double r)                              { this.R = r;                 return this; }
        public Builder updateListener(Consumer<MCSUpdate> l)    { this.updateListener = l;    return this; }
        public Builder cancellationCheck(BooleanSupplier check) { this.cancellationCheck = check; return this; }
        public Builder basis(CvCfBasis b)                       { this.basis = b;             return this; }
        public Builder matrixData(CMatrixPipeline.CMatrixData d){ this.matrixData = d;        return this; }
        public Builder lcf(int[][] l)                           { this.lcf = l;               return this; }

        public MCSRunner build() {
            if (clusterData == null) throw new IllegalStateException("clusterData required");
            if (eci == null)         throw new IllegalStateException("eci required");
            if (T <= 0)              throw new IllegalStateException("T must be > 0");
            if (xFrac == null) {
                xFrac = new double[numComp];
                for (int c = 0; c < numComp; c++) xFrac[c] = 1.0 / numComp;
            }
            if (basis == null)
                LOG.warning("MCSRunner built without a CvCfBasis — CVCF CF measurement unavailable");
            return new MCSRunner(this);
        }
    }
}
