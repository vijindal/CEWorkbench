package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;

import java.util.Random;

/**
 * Wall-clock benchmark: Embeddings.deltaEExchangeCvcf (flat, primary hot path)
 * vs. deltaEExchangeCvcfV2 (ATAT-style single-pass) on the same sequence of
 * random exchange trials, same geometry, same ECIs.
 *
 * <p>Both methods are already verified correct (DeltaEVerificationTest) — this
 * measures only whether V2's single-pass-per-embedding approach is actually
 * faster than V1's evaluate-before/evaluate-after approach in practice.
 */
public class DeltaEBenchmark {

    public static void main(String[] args) throws Exception {
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);

        String elements = args.length > 0 ? args[0] : "Nb-Ti-V-Zr";
        int L = args.length > 1 ? Integer.parseInt(args[1]) : 8;
        double T = args.length > 2 ? Double.parseDouble(args[2]) : 1000.0;
        int warmupTrials = args.length > 3 ? Integer.parseInt(args[3]) : 200_000;
        int timedTrials  = args.length > 4 ? Integer.parseInt(args[4]) : 2_000_000;

        Workspace.SystemId id = new Workspace.SystemId(elements, "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);

        System.out.println("=== Delta-E Benchmark (" + elements + " BCC_A2, L=" + L + ", T=" + T + ") ===");

        MCSGeometry geo = MCSGeometry.build(session, L, null);
        int N = geo.nSites();
        int ncf = geo.basis != null ? geo.basis.numNonPointCfs : geo.ncf;
        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "MCS-DeltaE-Bench");
        double[] eciOrth = computeEciOrth(eciCvcf, geo.basis, ncf);

        System.out.println("N sites = " + N + ", ncf = " + ncf
                + ", warmup=" + warmupTrials + ", timed=" + timedTrials);

        Random rng = new Random(123);
        LatticeConfig config = new LatticeConfig(N, geo.numComp);
        double[] xFrac = new double[geo.numComp];
        java.util.Arrays.fill(xFrac, 1.0 / geo.numComp);
        config.randomise(xFrac, rng);

        int maxEmbPerCol = Embeddings.maxEmbPerCfColumn(geo.cfEmbeddings);
        Embeddings.DeltaScratch scratch = new Embeddings.DeltaScratch(ncf, ncf * maxEmbPerCol);

        // Pre-generate the exact same (i,j) trial sequence for both methods, so both
        // engines see identical work and identical branch patterns (fair comparison).
        System.out.println("Generating trial sequence...");
        long genStart = System.nanoTime();
        int totalTrials = warmupTrials + timedTrials;
        int[] iSeq = new int[totalTrials];
        int[] jSeq = new int[totalTrials];
        for (int t = 0; t < totalTrials; t++) {
            int i, j;
            do {
                i = rng.nextInt(N);
                j = rng.nextInt(N);
            } while (config.getOccupation(i) == config.getOccupation(j));
            iSeq[t] = i;
            jSeq[t] = j;
            // Apply the swap so the walk explores a random trajectory (like real MC),
            // not the same fixed pre-swap state repeated.
            int occI = config.getOccupation(i), occJ = config.getOccupation(j);
            config.setOccupation(i, occJ);
            config.setOccupation(j, occI);
            if (t > 0 && t % 20000 == 0) {
                System.out.println("  ... " + t + "/" + totalTrials + " (" + (System.nanoTime() - genStart) / 1e9 + " s)");
            }
        }
        System.out.println("Trial sequence generated in " + (System.nanoTime() - genStart) / 1e9 + " s");

        // Reset config to the same random start for a fair, identical replay for each method.
        Random rng2 = new Random(123);
        LatticeConfig configV1 = new LatticeConfig(N, geo.numComp);
        configV1.randomise(xFrac, rng2);
        LatticeConfig configV2 = new LatticeConfig(N, geo.numComp);
        configV2.randomise(xFrac, new Random(123));

        double sinkV1 = runTimed("V1 (evaluate-before/after)", true,
                configV1, iSeq, jSeq, geo, eciCvcf, eciOrth, scratch, maxEmbPerCol, ncf, warmupTrials, timedTrials);
        double sinkV2 = runTimed("V2 (ATAT single-pass)", false,
                configV2, iSeq, jSeq, geo, eciCvcf, eciOrth, scratch, maxEmbPerCol, ncf, warmupTrials, timedTrials);

        System.out.println();
        System.out.println("(energy sums below exist only to prevent JIT dead-code elimination — not a physical quantity)");
        System.out.printf("V1 sink = %.6f, V2 sink = %.6f%n", sinkV1, sinkV2);
    }

    private static double runTimed(String label, boolean useV1,
            LatticeConfig config, int[] iSeq, int[] jSeq, MCSGeometry geo,
            double[] eciCvcf, double[] eciOrth, Embeddings.DeltaScratch scratch, int maxEmbPerCol, int ncf,
            int warmupTrials, int timedTrials) {

        double sink = 0.0;
        System.out.println("  [" + label + "] warmup starting...");
        long warmupStart = System.nanoTime();

        for (int t = 0; t < warmupTrials; t++) {
            int i = iSeq[t], j = jSeq[t];
            double dE = useV1
                    ? Embeddings.deltaEExchangeCvcf(i, j, config, geo.flatEmbData, geo.flatBasisMatrix,
                            geo.siteToCfIndex, ncf, eciCvcf, geo.basis, scratch, maxEmbPerCol, eciOrth, geo.numComp)
                    : Embeddings.deltaEExchangeCvcfV2(i, j, config, geo.flatEmbData, geo.flatBasisMatrix,
                            geo.siteToCfIndex, ncf, eciCvcf, geo.basis, scratch, maxEmbPerCol, eciOrth, geo.numComp);
            scratch.cleanup(maxEmbPerCol);
            sink += dE;
            int occI = config.getOccupation(i), occJ = config.getOccupation(j);
            config.setOccupation(i, occJ);
            config.setOccupation(j, occI);
            if (t > 0 && t % 5000 == 0) {
                System.out.println("  [" + label + "] warmup " + t + "/" + warmupTrials
                        + " (" + (System.nanoTime() - warmupStart) / 1e9 + " s)");
            }
        }
        System.out.println("  [" + label + "] warmup done in " + (System.nanoTime() - warmupStart) / 1e9 + " s");

        long start = System.nanoTime();
        for (int t = warmupTrials; t < warmupTrials + timedTrials; t++) {
            int i = iSeq[t], j = jSeq[t];
            double dE = useV1
                    ? Embeddings.deltaEExchangeCvcf(i, j, config, geo.flatEmbData, geo.flatBasisMatrix,
                            geo.siteToCfIndex, ncf, eciCvcf, geo.basis, scratch, maxEmbPerCol, eciOrth, geo.numComp)
                    : Embeddings.deltaEExchangeCvcfV2(i, j, config, geo.flatEmbData, geo.flatBasisMatrix,
                            geo.siteToCfIndex, ncf, eciCvcf, geo.basis, scratch, maxEmbPerCol, eciOrth, geo.numComp);
            scratch.cleanup(maxEmbPerCol);
            sink += dE;
            int occI = config.getOccupation(i), occJ = config.getOccupation(j);
            config.setOccupation(i, occJ);
            config.setOccupation(j, occI);
            if ((t - warmupTrials) > 0 && (t - warmupTrials) % 20000 == 0) {
                System.out.println("  [" + label + "] timed " + (t - warmupTrials) + "/" + timedTrials
                        + " (" + (System.nanoTime() - start) / 1e9 + " s)");
            }
        }
        long elapsedNanos = System.nanoTime() - start;

        double nsPerCall = (double) elapsedNanos / timedTrials;
        System.out.printf("%-28s  %,d calls in %.3f s  ->  %.1f ns/call%n",
                label, timedTrials, elapsedNanos / 1e9, nsPerCall);
        return sink;
    }

    private static double[] computeEciOrth(double[] eciCvcf, org.ce.model.cvm.CvCfBasis basis, int ncf) {
        if (basis == null || basis.Tinv == null) return null;
        double[][] Tinv = basis.Tinv;
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
}
