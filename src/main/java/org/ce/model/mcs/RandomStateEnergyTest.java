package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;

import java.util.Random;

/**
 * Ground-truth cross-check of the MCS energy formula against CVM's independently
 * derived and regression-tested enthalpy formula, at the fully random (disordered,
 * infinite-temperature) configuration — the one physical state both engines must
 * agree on exactly, without depending on either engine's minimizer or Metropolis
 * sampler being correct.
 *
 * <p>CVM: H_random = Σ_l eci[l] · u_random[l], where u_random comes from
 * {@link CVMGibbsModel#computeRandomCFs} — a closed-form combinatorial formula
 * for cluster probabilities of independently, randomly placed atoms at a fixed
 * global composition (no minimization involved).
 *
 * <p>MCS: H_random ≈ (1/N) Σ_l eciCvcf[l] · v_cvcf[l], measured directly from a
 * genuinely random {@link LatticeConfig} (no equilibration/Metropolis sampling —
 * just {@link LatticeConfig#randomise}) via {@link Embeddings#totalEnergyCvcf}.
 * As the supercell size L grows, finite-size sampling noise in the random
 * configuration's measured CFs should shrink, so the two values should converge.
 *
 * <p>If MCS's energy formula (basis definition, ECI mapping, Tinv transform, or
 * normalization) were wrong in some way that both deltaEExchangeCvcf and
 * totalEnergyCvcf shared, the earlier DeltaEVerificationTest would not have
 * caught it (self-consistency only) — this test would.
 */
public class RandomStateEnergyTest {

    public static void main(String[] args) throws Exception {
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);

        String elements = args.length > 0 ? args[0] : "Nb-Ti";
        double T = args.length > 1 ? Double.parseDouble(args[1]) : 1000.0;
        int[] Ls = { 4, 6, 8, 12, 16 };

        Workspace.SystemId id = new Workspace.SystemId(elements, "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);
        int numComp = session.numComponents();
        double[] xFrac = new double[numComp];
        java.util.Arrays.fill(xFrac, 1.0 / numComp);

        System.out.println("=== Random-State Energy Cross-Check (" + elements + " BCC_A2, T=" + T + ", equiatomic) ===");

        // ── Ground truth: CVM's closed-form random-state H (no minimization) ──
        CVMGibbsModel cvm = new CVMGibbsModel();
        cvm.initialize(session.systemId.elements(), session.systemId.structure(), session.systemId.model(),
                session.cecEntry, null);
        double[] uRandom = cvm.computeRandomCFs(xFrac);
        double[] eciCvm = CECEvaluator.evaluate(session.cecEntry, T, cvm.getBasis(), "CVM-random-check", null);
        // eci order matches u order (both indexed by non-point CF column l); ncf = uRandom.length
        double hCvmRandom = 0.0;
        for (int l = 0; l < uRandom.length && l < eciCvm.length; l++) hCvmRandom += eciCvm[l] * uRandom[l];

        System.out.println("CVM random-state u (CVCF, no minimization): " + java.util.Arrays.toString(uRandom));
        System.out.printf("CVM random-state H = Sum eci[l]*u[l] = %.6f J/mol%n%n", hCvmRandom);

        // ── MCS: measure the same quantity from genuinely random (unequilibrated) configs ──
        Random rng = new Random(7);
        for (int L : Ls) {
            MCSGeometry geo = MCSGeometry.build(session, L, null);
            int N = geo.nSites();
            int ncf = geo.basis != null ? geo.basis.numNonPointCfs : geo.ncf;
            double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "MCS-random-check");

            // Average over several independent random configs at this L to reduce noise
            int nSamples = 20;
            double sumH = 0.0, sumH2 = 0.0;
            for (int s = 0; s < nSamples; s++) {
                LatticeConfig config = new LatticeConfig(N, geo.numComp);
                config.randomise(xFrac, rng);
                double E = Embeddings.totalEnergyCvcf(
                        config, geo.cfEmbeddings, geo.flatBasisMatrix, ncf, eciCvcf, geo.basis, geo.numComp);
                double hPerSite = E / N;
                sumH += hPerSite;
                sumH2 += hPerSite * hPerSite;
            }
            double meanH = sumH / nSamples;
            double stdH = Math.sqrt(Math.max(0, sumH2 / nSamples - meanH * meanH));
            double diff = meanH - hCvmRandom;

            System.out.printf("L=%-3d N=%-5d  MCS random H = %.6f +/- %.6f J/mol   diff-from-CVM = %+.6f%n",
                    L, N, meanH, stdH, diff);
        }

        System.out.println();
        System.out.println("Expect |diff-from-CVM| to shrink toward 0 as L grows (finite-size sampling noise).");
        System.out.println("A converged, non-shrinking diff would indicate a real formula mismatch, not noise.");
    }
}
