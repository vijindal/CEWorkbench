package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult;
import org.ce.model.cluster.ClusterMath;
import org.ce.model.mcs.Embeddings;
import org.ce.model.mcs.LatticeConfig;
import org.ce.model.mcs.MCSGeometry;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.Arrays;

/**
 * Verifies CF measurement for uniform all-A configurations.
 *
 * Section A: Analytic expected values from cfBasisIndices.
 *   uOrth_allA[l] = basis[0]^(sum of cfBasisIndices[l])
 *
 * Section B: Measured values from Embeddings.measureCVsFromConfig on an all-A LatticeConfig.
 *
 * Pass: |measured[l] - expected[l]| < 1e-10 for all l.
 */
public class VerifyCFs {

    static final double TOLERANCE = 1e-10;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());

        String[][] systems = {
            {"Nb-Ti",      "2"},
            {"Nb-Ti-V",    "3"},
            {"Nb-Ti-V-Zr", "4"},
        };

        boolean allPass = true;

        for (String[] sys : systems) {
            String elements = sys[0];
            int K = Integer.parseInt(sys[1]);

            System.out.println("\n" + "=".repeat(70));
            System.out.println("  " + elements + " (K=" + K + ")  all-A config");
            System.out.println("=".repeat(70));

            // ── Section A: Analytic expected values ──────────────────────────
            ClusterIdentificationRequest config =
                ClusterIdentificationRequest.fromSystem(elements, "BCC_A2", "T");
            PipelineResult pr = ClusterCFIdentificationPipeline.runFullWorkflow(config, null);
            int[][] cfBasisIndices = pr.getCfBasisIndices();
            int ncf = pr.getNcf();

            double[] basis = ClusterMath.buildBasis(K);
            double basisA = basis[0];

            double[] expected = new double[ncf];
            for (int l = 0; l < ncf; l++) {
                int alphaSum = 0;
                for (int a : cfBasisIndices[l]) alphaSum += a;
                expected[l] = Math.pow(basisA, alphaSum);
            }

            System.out.println("Section A — analytic expected values (basis[0]=" + basisA + "):");
            for (int l = 0; l < ncf; l++) {
                System.out.printf("  l=%-3d  alphas=%-18s  expected=%+.6f%n",
                        l, Arrays.toString(cfBasisIndices[l]), expected[l]);
            }

            // ── Section B: Measured via Embeddings ───────────────────────────
            SystemId id = new SystemId(elements, "BCC_A2", "T");
            ModelSession session = builder.build(id, EngineConfig.MCS, null);
            MCSGeometry geo = MCSGeometry.build(session, 4, null); // L=4: 128 sites

            int N = geo.nSites();
            int numComp = session.numComponents();

            // All-A: every site = occupation 0
            LatticeConfig latticeConfig = new LatticeConfig(N, numComp);
            Arrays.fill(latticeConfig.getRawOcc(), 0);

            double[][] basisMatrix2d = Embeddings.buildBasisValues(numComp);
            // Flatten: flatBasisMatrix[occ * numComp + alpha] = phi_alpha(occ)
            // MetropolisMC uses: flatBasisMatrix[occ * numComp + (a+1)] = basisMatrix[occ][a]
            double[] flatBasisMatrix = new double[numComp * numComp];
            for (int o = 0; o < numComp; o++)
                for (int a = 0; a < numComp - 1; a++)
                    flatBasisMatrix[o * numComp + (a + 1)] = basisMatrix2d[o][a];

            double[] measured = Embeddings.measureCVsFromConfig(
                    latticeConfig, geo.cfEmbeddings(), flatBasisMatrix, ncf, numComp);

            // ── Compare ──────────────────────────────────────────────────────
            System.out.println("\nSection B — measured vs expected:");
            System.out.printf("  %-5s  %-18s  %-14s  %-14s  %-8s%n",
                    "CF#", "alphas", "expected", "measured", "status");
            System.out.println("  " + "-".repeat(68));

            int pass = 0, fail = 0;
            for (int l = 0; l < ncf; l++) {
                double diff = Math.abs(measured[l] - expected[l]);
                boolean ok = diff < TOLERANCE;
                if (ok) pass++; else { fail++; allPass = false; }
                System.out.printf("  %-5d  %-18s  %+14.6f  %+14.6f  %s%n",
                        l, Arrays.toString(cfBasisIndices[l]),
                        expected[l], measured[l],
                        ok ? "PASS" : "FAIL  diff=" + String.format("%.3e", diff));
            }
            System.out.println("  " + "-".repeat(68));
            System.out.printf("  %s: %d/%d CFs correct%n",
                    fail == 0 ? "PASS" : "FAIL", pass, ncf);
        }

        System.out.println("\n" + "=".repeat(70));
        System.out.println(allPass ? "ALL SYSTEMS PASS" : "SOME FAILURES — see above");
        System.out.println("=".repeat(70));
    }
}
