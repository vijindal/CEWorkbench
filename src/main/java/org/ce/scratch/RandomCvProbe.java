package org.ce.scratch;

import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.SpaceGroup;
import org.ce.model.cluster.StructurePhaseRegistry;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.storage.InputLoader;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

/**
 * Prints the random-state CV verification in <em>both</em> bases for a given
 * system:
 *
 * <ol>
 *   <li><b>Orthogonal</b> -- {@link CMatrixPipeline#verifyRandomCVs}, which
 *       checks every cluster variable against the analytic random-state value
 *       {@code (1/K)^n} for an n-site cluster at equiatomic composition.</li>
 *   <li><b>CVCF</b> -- the same evaluation against the Stage 4 CVCF C-matrix,
 *       using the CVCF random CF vector.</li>
 * </ol>
 *
 * <p>At the random state every configuration of an n-site cluster is equally
 * likely, so each cluster probability must be exactly {@code (1/K)^n}
 * regardless of basis. That makes this an analytic ground-truth check on the
 * C-matrix, independent of any Hamiltonian or energy expression.</p>
 */
public final class RandomCvProbe {

    public static void main(String[] args) {
        String elements = args.length > 0 ? args[0] : "Nb-Ti";
        String structure = args.length > 1 ? args[1] : "BCC_A2";
        String model = args.length > 2 ? args[2] : "T";

        Consumer<String> sink = System.out::println;

        List<String> elementList = List.of(elements.split("-"));
        int K = elementList.size();

        String parent = StructurePhaseRegistry.parentOf(structure);
        String mod = model.replace("_CVCF", "");

        List<Cluster> disClusters = InputLoader.parseClusterFile("clus/" + parent + "-" + mod + ".txt");
        disClusters.replaceAll(Cluster::sorted);
        SpaceGroup disSG = InputLoader.parseSpaceGroup(parent + "-SG");

        List<Cluster> ordClusters = InputLoader.parseClusterFile("clus/" + structure + "-" + mod + ".txt");
        ordClusters.replaceAll(Cluster::sorted);
        SpaceGroup ordSG = InputLoader.parseSpaceGroup(structure + "-SG");

        System.out.println("################################################################");
        System.out.printf("#  RANDOM-STATE CV VERIFICATION  --  %s / %s / %s  (K=%d)%n",
                elements, structure, model, K);
        System.out.println("################################################################");

        ClusterCFIdentificationPipeline.PipelineResult pr = ClusterCFIdentificationPipeline.run(
                disClusters, disSG.getOperations(), ordClusters, ordSG.getOperations(),
                ordSG.getRotateMat(), ordSG.getTranslateMat(), K, null);

        CMatrixPipeline.CMatrixData cmatOrth = CMatrixPipeline.run(
                pr.toClusterIdentificationResult(), pr.toCFIdentificationResult(),
                ordClusters, K, null);

        double[] x = new double[K];
        Arrays.fill(x, 1.0 / K);

        // ---------------- 1. ORTHOGONAL BASIS ----------------
        System.out.println("\n\n****************************************************************");
        System.out.println("*  BASIS 1 of 2:  ORTHOGONAL");
        System.out.println("****************************************************************");
        CMatrixPipeline.verifyRandomCVs(x, pr, cmatOrth, sink);

        // ---------------- 2. CVCF BASIS ----------------
        System.out.println("\n\n****************************************************************");
        System.out.println("*  BASIS 2 of 2:  CVCF");
        System.out.println("****************************************************************");

        CvCfBasis basis = CvCfBasis.generate(structure, pr, cmatOrth, model, null);
        CMatrixPipeline.CMatrixData cmatCvcf = basis.cvcfCMatrixData;

        int ncf = basis.numNonPointCfs;
        double[] uCvcfFull = randomCvcfFull(pr, basis, x, ncf);

        System.out.println("\n  Composition: " + Arrays.toString(x));
        System.out.println("  Random CVCF CF vector uFull = [u ; x] (length " + uCvcfFull.length + "):");
        System.out.println("    " + Arrays.toString(uCvcfFull));

        double[][][] cv = CMatrixPipeline.evaluateCVs(
                uCvcfFull, cmatCvcf.getCmat(), cmatCvcf.getLcv(), pr.getTcdis(), pr.getLc());

        int bad = 0;
        for (int t = 0; t < cv.length; t++) {
            int n = pr.getDisClusData().getClusCoordList().get(t).getAllSites().size();
            double expected = Math.pow(1.0 / K, n);
            for (int j = 0; j < cv[t].length; j++) {
                System.out.printf("%n  Cluster Type t=%d, Group j=%d (n=%d sites):%n", t, j, n);
                double sum = 0;
                int[] w = cmatCvcf.getWcv().get(t).get(j);
                for (int v = 0; v < cv[t][j].length; v++) {
                    double val = cv[t][j][v];
                    double err = Math.abs(val - expected);
                    sum += w[v] * val;
                    String msg = String.format("    CV[%2d] = %12.8f (Expected: %12.8f, Diff: %.2e)",
                            v, val, expected, err);
                    if (err > 1e-9) { msg += " [!] DISCREPANCY"; bad++; }
                    System.out.println(msg);
                }
                System.out.printf("    sum(wcv*CV) = %.12f  (must be 1)%n", sum);
            }
        }

        System.out.println("\n================================================================");
        System.out.printf("CVCF RESULT: %s   (%d discrepancies)%n", bad == 0 ? "PASS" : "FAIL", bad);
        System.out.println("================================================================");
    }

    /**
     * Random-state CVCF CF vector in {@code uFull = [u ; x]} layout. Prefers
     * the pipeline's own {@code equiatomicCVCF} when present; otherwise
     * transforms the orthogonal random CFs into the CVCF basis via
     * {@code Tinv}.
     */
    private static double[] randomCvcfFull(
            ClusterCFIdentificationPipeline.PipelineResult pr, CvCfBasis basis, double[] x, int ncf) {

        double[] equi = pr.getEquiatomicCVCF();
        if (equi != null) {
            return equi;
        }
        double[] orthFull = pr.computeRandomCFs(x);
        double[][] tinv = basis.Tinv;
        int width = ncf + x.length;
        double[] out = new double[width];
        for (int i = 0; i < width && i < tinv.length; i++) {
            double s = 0;
            for (int k = 0; k < tinv[i].length && k < orthFull.length; k++) {
                s += tinv[i][k] * orthFull[k];
            }
            out[i] = s;
        }
        return out;
    }
}
