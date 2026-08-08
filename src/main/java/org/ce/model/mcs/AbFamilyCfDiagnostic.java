package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;
import org.ce.model.cvm.CvCfBasis;

import java.util.List;
import java.util.Random;

/**
 * Isolates the AB-family CF anomaly found by CvmVsMcsComparisonTest: dumps the
 * orthogonal (pre-Tinv) CF vector uOrth alongside the CVCF (post-Tinv) vector v,
 * so the bug can be localized to either measureCVsFromConfig (embedding structure /
 * basis product evaluation) or applyTinvTransform (Tinv matrix / point-CF assembly).
 */
public class AbFamilyCfDiagnostic {

    public static void main(String[] args) throws Exception {
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);
        Workspace.SystemId id = new Workspace.SystemId("Nb-Ti-V", "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);

        MCSGeometry geo = MCSGeometry.build(session, 4, null);
        int N = geo.nSites();
        int ncf = geo.basis != null ? geo.basis.numNonPointCfs : geo.ncf;
        List<String> cfNames = CvCfBasis.getNonPointCfNames("BCC_A2", "T", geo.numComp);

        Random rng = new Random(42);
        LatticeConfig config = new LatticeConfig(N, geo.numComp);
        double[] xFrac = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
        config.randomise(xFrac, rng);

        System.out.println("=== AB-Family CF Diagnostic (Nb-Ti-V, L=4, N=" + N + ") ===");
        System.out.println("numComp=" + geo.numComp + " ncf=" + ncf);

        // Step 1: raw orthogonal CFs, straight from measureCVsFromConfig — no Tinv involved.
        double[] uOrth = Embeddings.measureCVsFromConfig(config, geo.cfEmbeddings, geo.flatBasisMatrix, ncf, geo.numComp);

        System.out.println();
        System.out.println("uOrth (pre-Tinv, orthogonal basis) — should all be in a sane range (~[-1,1]):");
        for (int l = 0; l < ncf; l++) {
            System.out.printf("  uOrth[%2d] (%-10s) = %10.5f   nEmb=%d%n",
                    l, cfNames.get(l), uOrth[l], geo.cfEmbeddings.get(l).size());
        }

        // Step 2: apply Tinv to get CVCF basis.
        double[] v = Embeddings.applyTinvTransform(uOrth, config.composition(), geo.basis);

        System.out.println();
        System.out.println("v (post-Tinv, CVCF basis):");
        for (int l = 0; l < ncf; l++) {
            System.out.printf("  v[%2d]     (%-10s) = %10.5f%n", l, cfNames.get(l), v[l]);
        }

        // Step 3: dump the Tinv matrix rows for the AB-family columns specifically,
        // to see if the transform coefficients themselves look wrong.
        System.out.println();
        System.out.println("Tinv dimensions: " + geo.basis.Tinv.length + " x " + geo.basis.Tinv[0].length);
        System.out.println("Tinv rows for AB-family CF columns (index : name : row):");
        for (int l = 0; l < ncf; l++) {
            if (cfNames.get(l).contains("AB")) {
                System.out.printf("  row[%2d] (%s): %s%n", l, cfNames.get(l), java.util.Arrays.toString(geo.basis.Tinv[l]));
            }
        }

        // Step 4: cross-check against CvCfBasis's own random-state formula (ground truth,
        // no MCS code involved at all) at this exact equiatomic composition.
        System.out.println();
        System.out.println("=== Cross-check: composition() vs xFrac ===");
        System.out.println("config.composition() = " + java.util.Arrays.toString(config.composition()));
        System.out.println("xFrac (requested)    = " + java.util.Arrays.toString(xFrac));

        // Step 5: independent ground truth for uOrth at the random state — a closed-form
        // formula (ClusterCFIdentificationPipeline.computeRandomCFs), zero dependence on
        // MCS's Embeddings/generateCfEmbeddings code, zero dependence on Tinv/CVM's minimizer.
        // This directly tests whether the raw orthogonal measurement itself is correct,
        // independent of which basis (orthogonal vs CVCF) anyone later reports results in.
        String structure = "BCC_A2";
        org.ce.model.cluster.ClusterIdentificationRequest cfg =
                org.ce.model.cluster.ClusterIdentificationRequest.fromSystem("Nb-Ti-V", structure, "T");
        java.util.List<org.ce.model.cluster.Cluster> disClusters =
                org.ce.model.storage.InputLoader.parseClusterFile(cfg.getDisorderedClusterFile());
        disClusters.replaceAll(org.ce.model.cluster.Cluster::sorted);
        org.ce.model.cluster.SpaceGroup disSg =
                org.ce.model.storage.InputLoader.parseSpaceGroup(cfg.getDisorderedSymmetryGroup());
        java.util.List<org.ce.model.cluster.Cluster> ordClusters =
                org.ce.model.storage.InputLoader.parseClusterFile(cfg.getOrderedClusterFile());
        ordClusters.replaceAll(org.ce.model.cluster.Cluster::sorted);
        org.ce.model.cluster.SpaceGroup ordSg =
                org.ce.model.storage.InputLoader.parseSpaceGroup(cfg.getOrderedSymmetryGroup());

        var pr = org.ce.model.cluster.ClusterCFIdentificationPipeline.run(
                disClusters, disSg.getOperations(), ordClusters, ordSg.getOperations(),
                cfg.getTransformationMatrix(),
                new double[] { cfg.getTranslationVector().getX(), cfg.getTranslationVector().getY(), cfg.getTranslationVector().getZ() },
                geo.numComp, null);

        double[] uOrthGroundTruth = pr.computeRandomCFs(xFrac);

        System.out.println();
        System.out.println("=== Ground truth: closed-form uOrth at random state (no MCS code involved) ===");
        System.out.println("(non-point columns only, first " + ncf + " entries)");
        for (int l = 0; l < ncf; l++) {
            double gt = uOrthGroundTruth[l];
            double measured = uOrth[l];
            double diff = measured - gt;
            System.out.printf("  [%2d] %-10s  groundTruth=%10.5f  measured=%10.5f  diff=%10.5f%s%n",
                    l, cfNames.get(l), gt, measured, diff,
                    Math.abs(diff) > 0.05 ? "  <-- MISMATCH" : "");
        }
    }
}
