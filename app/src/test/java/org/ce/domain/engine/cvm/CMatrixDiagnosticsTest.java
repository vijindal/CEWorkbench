package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CMatrixResult;
import org.ce.domain.cluster.cvcf.BccA2CvCfTransformations;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.cluster.cvcf.CvCfBasisTransformer;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.Workspace;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;

import java.util.List;

/**
 * Fetch and compare C-matrices in orthogonal vs CVCF basis.
 */
public class CMatrixDiagnosticsTest {

    @Test
    @DisplayName("Fetch C-matrices: Orthogonal vs CVCF basis")
    void testFetchCMatrices() throws Exception {
        Workspace workspace = new Workspace();
        ClusterDataStore clusterStore = new ClusterDataStore(workspace);
        AllClusterData allData = clusterStore.load("BCC_A2_T_bin");

        System.out.println("\n" + "=" .repeat(90));
        System.out.println("C-MATRIX COMPARISON: ORTHOGONAL vs CVCF BASIS");
        System.out.println("=" .repeat(90));
        System.out.println();

        CMatrixResult oldCmat = allData.getCMatrixResult();
        CvCfBasis cvcfBasis = BccA2CvCfTransformations.binaryBasis();
        CMatrixResult newCmat = CvCfBasisTransformer.transform(oldCmat, cvcfBasis);

        // Get first type and group (NN cluster)
        List<List<double[][]>> oldCmatList = oldCmat.getCmat();
        List<List<double[][]>> newCmatList = newCmat.getCmat();

        System.out.println("ORTHOGONAL BASIS - First Type (t=0), First Group (j=0):");
        System.out.println("-" .repeat(90));
        double[][] oldCm = oldCmatList.get(0).get(0);
        int nv = oldCm.length;
        int oldCols = oldCm[0].length;
        System.out.printf("Dimensions: %d CVs × %d columns%n", nv, oldCols);
        System.out.println();
        printCMatrix(oldCm, "Orthogonal C-matrix");
        System.out.println();

        System.out.println("CVCF BASIS - First Type (t=0), First Group (j=0):");
        System.out.println("-" .repeat(90));
        double[][] newCm = newCmatList.get(0).get(0);
        int newCols = newCm[0].length;
        System.out.printf("Dimensions: %d CVs × %d columns%n", nv, newCols);
        System.out.println();
        printCMatrix(newCm, "CVCF C-matrix");
        System.out.println();

        System.out.println("TRANSFORMATION MATRIX T:");
        System.out.println("-" .repeat(90));
        double[][] T = cvcfBasis.T;
        System.out.printf("Dimensions: %d × %d%n", T.length, T[0].length);
        System.out.println();
        printMatrix(T);
        System.out.println();

        System.out.println("=" .repeat(90));
    }

    private void printCMatrix(double[][] cm, String title) {
        System.out.println(title + ":");
        System.out.printf("%-6s", "CV");
        for (int c = 0; c < cm[0].length; c++) {
            System.out.printf("%12s", "[" + c + "]");
        }
        System.out.println();
        System.out.println("-" .repeat(Math.max(80, 6 + cm[0].length * 12)));

        for (int v = 0; v < cm.length; v++) {
            System.out.printf("%-6d", v);
            for (int c = 0; c < cm[v].length; c++) {
                System.out.printf("%12.6f", cm[v][c]);
            }
            System.out.println();
        }
    }

    private void printMatrix(double[][] matrix) {
        System.out.printf("%-8s", "Row");
        for (int c = 0; c < matrix[0].length; c++) {
            System.out.printf("%12s", "[" + c + "]");
        }
        System.out.println();
        System.out.println("-" .repeat(Math.max(80, 8 + matrix[0].length * 12)));

        for (int r = 0; r < matrix.length; r++) {
            System.out.printf("%-8d", r);
            for (int c = 0; c < matrix[r].length; c++) {
                System.out.printf("%12.6f", matrix[r][c]);
            }
            System.out.println();
        }
    }
}
