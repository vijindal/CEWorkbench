package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.PhysicsConstants;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.mcs.MetropolisMC.MCSUpdate;
import org.ce.model.mcs.MetropolisMC.MCResult;

import java.util.Arrays;
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
    private final double[]    eciOrth;      // ECIs in orthogonal basis
    private final double[]    eciCvcf;      // ECIs in CVCF basis
    private final double      T;
    private final double      R;

    private MCSRunner(MCSGeometry geo, double[] eciCvcf, double[] eciOrth, double T, double R) {
        this.geo = geo;
        this.eciCvcf = eciCvcf.clone();
        this.eciOrth = eciOrth.clone();
        this.T = T;
        this.R = R;

        LOG.info(String.format("MCSRunner ready — T=%.1f K, numComp=%d", T, geo.numComp));
        Debug.printEciInfo(eciCvcf, eciOrth, geo.basis, geo.nSites(), geo.orbitSizes);
    }

    /**
     * Factory method to build a runner for a specific temperature using an existing geometry.
     */
    public static MCSRunner forTemperature(MCSGeometry geo, ModelSession session, double T, Consumer<String> progressSink) {
        // Evaluate Hamiltonian at T — ECI in CVCF basis
        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "MCS", progressSink);
        
        // Transform ECIs to orthogonal basis
        int tc = geo.clusterData.getTc();
        double[] eciOrth = buildEciByOrbitType(eciCvcf, tc, geo.basis);
        
        return new MCSRunner(geo, eciCvcf, eciOrth, T, PhysicsConstants.R_GAS);
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

        emit(progressSink, "  Step 2 [HAMILTONIAN]  computing E(σ₀) via cluster expansion...");

        MetropolisMC.Sampler sampler = new MetropolisMC.Sampler(
                N, geo.orbitSizes, geo.orbits, R, eciCvcf,
                geo.multiSiteEmbedCounts, geo.basis, geo.cfEmbeddings, geo.basisMatrix);

        MetropolisMC engine = new MetropolisMC(geo.emb, eciOrth, geo.orbits, geo.numComp, T, nEquil, nAvg, R, rng);
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

    private static double[] buildEciByOrbitType(double[] eciCvcf, int tc, org.ce.model.cvm.CvCfBasis basis) {
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

        public static void printEciInfo(double[] eciCvcf, double[] eciOrth,
                                        org.ce.model.cvm.CvCfBasis basis, int N, int[] orbitSizes) {
            System.out.println("============================================================");
            System.out.println("  MCS ECI & ORBIT DIAGNOSTIC");
            System.out.println("============================================================");
            if (basis != null) {
                System.out.printf("  CVCF Basis: %s (K=%d), numNonPointCfs=%d%n",
                        basis.structurePhase, basis.numComponents, basis.numNonPointCfs);
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
}
