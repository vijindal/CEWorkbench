package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;

import java.util.Random;

/**
 * Faithfully reproduces MetropolisMC.run()'s real trial-move sequence — genuine
 * Metropolis accept/reject via MetropolisMC.ExchangeStep.attempt() (not an
 * "accept everything" stress walk like DeltaEVerificationTestV2Rigorous used) —
 * starting from a freshly randomized configuration, and checks the running
 * energy against a from-scratch totalEnergyCvcf recomputation after EVERY
 * single step (not every 100 sweeps like the production drift check).
 *
 * <p>This directly answers: starting from random initialization, is the
 * incremental ΔE actually being computed and applied correctly at every real
 * MC step, including rejected moves (which must leave config/energy untouched)?
 */
public class MetropolisTrajectoryVerification {

    public static void main(String[] args) throws Exception {
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);

        String elements = args.length > 0 ? args[0] : "Nb-Ti-V";
        int L = args.length > 1 ? Integer.parseInt(args[1]) : 4;
        double T = args.length > 2 ? Double.parseDouble(args[2]) : 1000.0;
        int nSteps = args.length > 3 ? Integer.parseInt(args[3]) : 2000;
        long seed = args.length > 4 ? Long.parseLong(args[4]) : 7;

        Workspace.SystemId id = new Workspace.SystemId(elements, "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);

        MCSGeometry geo = MCSGeometry.build(session, L, null);
        int N = geo.nSites();
        int ncf = geo.basis != null ? geo.basis.numNonPointCfs : geo.ncf;
        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "MetroTraj-Verify");
        double[] eciOrth = computeEciOrth(eciCvcf, geo.basis, ncf);

        System.out.println("=== Metropolis Trajectory Verification (" + elements + " BCC_A2, L=" + L
                + ", T=" + T + ", steps=" + nSteps + ", seed=" + seed + ") ===");
        System.out.println("N=" + N + " ncf=" + ncf);

        Random rng = new Random(seed);
        LatticeConfig config = new LatticeConfig(N, geo.numComp);
        double[] xFrac = new double[geo.numComp];
        java.util.Arrays.fill(xFrac, 1.0 / geo.numComp);
        config.randomise(xFrac, rng);   // genuinely random start, exactly as MetropolisMC.run() receives it

        double R = org.ce.model.PhysicsConstants.R_GAS;

        // No spatial decomposition (serial path) — exercises the exact same
        // ExchangeStep.attempt() the serial equilibration/averaging loop calls.
        MetropolisMC.ExchangeStep move = new MetropolisMC.ExchangeStep(
                geo.cfEmbeddings, geo.flatBasisMatrix, geo.siteToCfIndex,
                ncf, eciCvcf, eciOrth, geo.basis, geo.numComp, T, R, rng, null);

        double currentEnergy = Embeddings.totalEnergyCvcf(
                config, geo.cfEmbeddings, geo.pointCfEmbeddings, geo.flatBasisMatrix, ncf, eciCvcf, geo.basis, geo.numComp);
        double initialFullEnergy = currentEnergy;
        System.out.printf("Initial E (from random config) = %.6f (E/site = %.6f)%n", currentEnergy, currentEnergy / N);

        int mismatches = 0;
        double maxAbsErr = 0.0, maxRelErr = 0.0;

        for (int step = 0; step < nSteps; step++) {
            double dEreturned = move.attempt(config);
            currentEnergy += dEreturned;

            double trueEnergy = Embeddings.totalEnergyCvcf(
                    config, geo.cfEmbeddings, geo.pointCfEmbeddings, geo.flatBasisMatrix, ncf, eciCvcf, geo.basis, geo.numComp);
            double err = Math.abs(currentEnergy - trueEnergy);
            double relErr = err / Math.max(1.0, Math.abs(trueEnergy));
            maxAbsErr = Math.max(maxAbsErr, err);
            maxRelErr = Math.max(maxRelErr, relErr);

            if (err > 1e-6) {
                mismatches++;
                System.out.printf("  [MISMATCH] step=%d dEreturned=%.10f currentEnergy=%.10f trueEnergy=%.10f err=%.3e%n",
                        step, dEreturned, currentEnergy, trueEnergy, err);
            } else if (step < 5 || step == nSteps - 1) {
                System.out.printf("  [OK] step=%d dEreturned=%.10f currentEnergy=%.10f trueEnergy=%.10f err=%.3e acceptRate=%.3f%n",
                        step, dEreturned, currentEnergy, trueEnergy, err, move.acceptRate());
            }
        }

        System.out.println();
        System.out.println("=== Summary ===");
        System.out.println("Steps: " + nSteps);
        System.out.printf("Final accept rate: %.4f%n", move.acceptRate());
        System.out.println("Mismatches (>1e-6 abs err vs from-scratch recompute, checked EVERY step): " + mismatches);
        System.out.printf("Max abs error: %.3e%n", maxAbsErr);
        System.out.printf("Max rel error: %.3e%n", maxRelErr);
        System.out.printf("Initial full energy: %.6f, final running energy: %.6f%n", initialFullEnergy, currentEnergy);
        System.out.println(mismatches == 0 ? "RESULT: PASS" : "RESULT: FAIL");
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
