package org.ce.scratch;

import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult;
import org.ce.model.cluster.ClusterMath;

import java.util.Arrays;

/**
 * Prints cfBasisIndices and analytic all-A CF values for binary, ternary, quaternary BCC_A2 T-model.
 *
 * For all-A (every site = species 0):
 *   uOrth[l] = basis[0]^(sum of cfBasisIndices[l])
 */
public class PrintCfBasis {
    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        String[][] systems = {
            {"Nb-Ti",      "2"},
            {"Nb-Ti-V",    "3"},
            {"Nb-Ti-V-Zr", "4"},
        };

        for (String[] sys : systems) {
            String elements = sys[0];
            int K = Integer.parseInt(sys[1]);

            ClusterIdentificationRequest config =
                ClusterIdentificationRequest.fromSystem(elements, "BCC_A2", "T");
            PipelineResult pr = ClusterCFIdentificationPipeline.runFullWorkflow(config, null);

            int[][] cfBasisIndices = pr.getCfBasisIndices();
            int ncf = pr.getNcf();  // non-point CFs only
            int tcf = pr.getTcf();  // total including point CFs

            double[] basis = ClusterMath.buildBasis(K);
            double basisA = basis[0];

            System.out.println("\n=== " + elements + " (K=" + K + ") ===");
            System.out.println("basis        = " + Arrays.toString(basis));
            System.out.println("basis[0] (A) = " + basisA);
            System.out.println("ncf (non-point) = " + ncf + ",  tcf (total) = " + tcf);
            System.out.println();
            System.out.printf("%-6s  %-24s  %-10s  %-14s%n", "CF#", "alpha pattern", "alpha sum", "uOrth(all-A)");
            System.out.println("-".repeat(62));
            for (int l = 0; l < ncf; l++) {
                int[] alphas = cfBasisIndices[l];
                int alphaSum = 0;
                for (int a : alphas) alphaSum += a;
                double expected = Math.pow(basisA, alphaSum);
                System.out.printf("%-6d  %-24s  %-10d  %+.6f%n",
                        l, Arrays.toString(alphas), alphaSum, expected);
            }
            // Also print point CFs
            System.out.println();
            System.out.println("Point CFs (l = ncf .. tcf-1):");
            for (int l = ncf; l < tcf - 1; l++) {  // tcf-1 = empty cluster
                int[] alphas = cfBasisIndices[l];
                int alphaSum = 0;
                for (int a : alphas) alphaSum += a;
                double expected = Math.pow(basisA, alphaSum);
                System.out.printf("  l=%-4d  alphas=%-10s  uPoint(all-A)=%+.6f%n",
                        l, Arrays.toString(alphas), expected);
            }
        }
    }
}
