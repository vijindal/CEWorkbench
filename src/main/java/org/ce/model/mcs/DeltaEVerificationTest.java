package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;

import java.util.Random;

/**
 * Verifies that the incremental ΔE hot path (Embeddings.deltaEExchangeCvcf) agrees
 * with the finite-difference ΔE obtained from two independent full-energy evaluations
 * (Embeddings.totalEnergyCvcf) before and after the same site exchange.
 *
 * <p>This checks Stage 4 (initial/full energy) and Stage 6 (trial-move ΔE) of the MCS
 * pipeline together, including the eciOrth pre-computation and both the flat (FlatEmbData)
 * and list-based ΔE implementations.
 */
public class DeltaEVerificationTest {

    public static void main(String[] args) throws Exception {
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);

        String elements = args.length > 0 ? args[0] : "Nb-Ti";
        int L = args.length > 1 ? Integer.parseInt(args[1]) : 3;
        double T = args.length > 2 ? Double.parseDouble(args[2]) : 1000.0;

        Workspace.SystemId id = new Workspace.SystemId(elements, "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);

        System.out.println("=== Delta-E Verification (" + elements + " BCC_A2, L=" + L + ", T=" + T + ") ===");

        MCSGeometry geo = MCSGeometry.build(session, L, null);
        int N = geo.nSites();
        int ncf = geo.basis != null ? geo.basis.numNonPointCfs : geo.ncf;

        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "MCS-DeltaE-Test");
        double[] eciOrth = computeEciOrth(eciCvcf, geo.basis, ncf);

        System.out.println("N sites = " + N + ", ncf = " + ncf);

        Random rng = new Random(42);
        LatticeConfig config = new LatticeConfig(N, geo.numComp);
        double[] xFrac = new double[geo.numComp];
        java.util.Arrays.fill(xFrac, 1.0 / geo.numComp);
        config.randomise(xFrac, rng);

        int maxEmbPerCol = Embeddings.maxEmbPerCfColumn(geo.cfEmbeddings);
        // "seen" is indexed by key = l * maxEmbPerCol + ei (dense grid over ncf x maxEmbPerCol),
        // so it must be sized ncf * maxEmbPerCol, not the sum of actual per-column embedding
        // counts (totalCfEmbeddingCount) — matches MetropolisMC's own scratch sizing.
        Embeddings.DeltaScratch scratch = new Embeddings.DeltaScratch(ncf, ncf * maxEmbPerCol);

        int trials = 200;
        int mismatches = 0;
        double maxAbsErr = 0.0;
        double maxRelErr = 0.0;

        for (int trial = 0; trial < trials; trial++) {
            int i = rng.nextInt(N);
            int j = rng.nextInt(N);
            if (config.getOccupation(i) == config.getOccupation(j)) continue;

            double eBefore = Embeddings.totalEnergyCvcf(
                    config, geo.cfEmbeddings, geo.flatBasisMatrix, ncf, eciCvcf, geo.basis, geo.numComp);

            // Flat-path hot ΔE (primary production path)
            double dEFlat = Embeddings.deltaEExchangeCvcf(
                    i, j, config, geo.flatEmbData, geo.flatBasisMatrix, geo.siteToCfIndex,
                    ncf, eciCvcf, geo.basis, scratch, maxEmbPerCol, eciOrth, geo.numComp);
            scratch.cleanup(maxEmbPerCol);

            // List-based fallback ΔE, independently computed on the same pre-swap state
            double dEList = Embeddings.deltaEExchangeCvcf(
                    i, j, config, geo.cfEmbeddings, geo.flatBasisMatrix, geo.siteToCfIndex,
                    ncf, eciCvcf, geo.basis, scratch, maxEmbPerCol, eciOrth, geo.numComp);

            // ATAT-style single-pass ΔE (V2), independently computed on the same pre-swap state
            double dEV2 = Embeddings.deltaEExchangeCvcfV2(
                    i, j, config, geo.flatEmbData, geo.flatBasisMatrix, geo.siteToCfIndex,
                    ncf, eciCvcf, geo.basis, scratch, maxEmbPerCol, eciOrth, geo.numComp);
            scratch.cleanup(maxEmbPerCol);

            // Apply the swap for real, then measure the true finite-difference ΔE.
            int occI = config.getOccupation(i), occJ = config.getOccupation(j);
            config.setOccupation(i, occJ);
            config.setOccupation(j, occI);

            double eAfter = Embeddings.totalEnergyCvcf(
                    config, geo.cfEmbeddings, geo.flatBasisMatrix, ncf, eciCvcf, geo.basis, geo.numComp);

            double dEFiniteDiff = eAfter - eBefore;

            double errFlat = Math.abs(dEFlat - dEFiniteDiff);
            double errList = Math.abs(dEList - dEFiniteDiff);
            double errV2   = Math.abs(dEV2 - dEFiniteDiff);
            double relFlat = errFlat / Math.max(1e-9, Math.abs(dEFiniteDiff));

            maxAbsErr = Math.max(maxAbsErr, Math.max(errFlat, Math.max(errList, errV2)));
            maxRelErr = Math.max(maxRelErr, relFlat);

            if (errFlat > 1e-6 || errList > 1e-6 || errV2 > 1e-6) {
                mismatches++;
                System.out.printf("  [MISMATCH] trial=%d i=%d j=%d  dE_flat=%.10f  dE_list=%.10f  dE_v2=%.10f  dE_finiteDiff=%.10f  errFlat=%.3e errList=%.3e errV2=%.3e%n",
                        trial, i, j, dEFlat, dEList, dEV2, dEFiniteDiff, errFlat, errList, errV2);
            } else if (trial < 5) {
                System.out.printf("  [OK] trial=%d i=%d j=%d  dE_flat=%.10f  dE_v2=%.10f  dE_finiteDiff=%.10f  err=%.3e%n",
                        trial, i, j, dEFlat, dEV2, dEFiniteDiff, errFlat);
            }

            // Revert the swap so successive trials start from the same reference config
            // (keeps composition fixed and avoids compounding drift across trials).
            config.setOccupation(i, occI);
            config.setOccupation(j, occJ);
        }

        System.out.println();
        System.out.println("=== Summary ===");
        System.out.println("Trials run (excluding same-species no-ops): " + trials);
        System.out.println("Mismatches (> 1e-6 abs error): " + mismatches);
        System.out.printf("Max abs error: %.3e%n", maxAbsErr);
        System.out.printf("Max rel error: %.3e%n", maxRelErr);
        System.out.println(mismatches == 0 ? "RESULT: PASS" : "RESULT: FAIL");
    }

    /** Mirrors MCSRunner.computeEciOrth (package-private there) for standalone testing. */
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
