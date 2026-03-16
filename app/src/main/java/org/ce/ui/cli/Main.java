package org.ce.ui.cli;

import org.ce.workflow.ClusterIdentificationWorkflow;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.domain.cluster.*;
import org.ce.domain.cluster.Vector3D;
import org.ce.domain.cluster.AllClusterData;
import org.ce.storage.input.InputLoader;

import java.util.List;

/**
 * Main entry point for running the complete CVM identification pipeline
 * including cluster, CF, and C-matrix calculations.
 */
public class Main {

    public static void main(String[] args) {
        System.out.println("================================================================================");
        System.out.println("                    COMPLETE CVM IDENTIFICATION PIPELINE");
        System.out.println("================================================================================");

        try {
            // Build configuration for A2 system (binary)
            ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                    .disorderedClusterFile("clus/A2-T.txt")      // Disordered phase
                    .orderedClusterFile("clus/A2-T.txt")         // Ordered phase (using A2 for demo)
                    .disorderedSymmetryGroup("A2-SG")            // Disordered symmetry
                    .orderedSymmetryGroup("A2-SG")               // Ordered symmetry
                    .transformationMatrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}})
                    .translationVector(new Vector3D(0,0,0))
                    .numComponents(5)                             // Binary system
                    .build();

            System.out.println("\nConfiguration:");
            System.out.println("  Disordered: clus/A2-T.txt with A2-SG");
            System.out.println("  Ordered: clus/A2-T.txt with A2-SG");
            System.out.println("  Components: 2 (binary)");

            // ========================================================================
            // STAGE 1-2: Run the pipeline (cluster + CF identification)
            // ========================================================================
            System.out.println("\nStage 1-2: Running cluster and CF identification...");
            AllClusterData result = ClusterIdentificationWorkflow.identify(config);

            System.out.println("\nCluster and CF identification completed!");
            System.out.println(result.getSummary());

            // Print detailed results
            result.printResults();

            // ========================================================================
            // STAGE 3: Build C-matrix
            // ========================================================================
            System.out.println("\nStage 3: Building C-matrix...");

            // Load maximal clusters for C-matrix builder
            List<Cluster> maxClusters = InputLoader.parseClusterFile("clus/A2-T.txt");

            // Build C-matrix
            CMatrixResult cmatrixResult = CMatrixBuilder.build(
                    result.getDisorderedClusterResult(),
                    result.getDisorderedCFResult(),
                    maxClusters,
                    config.getNumComponents()
            );

            System.out.println("C-matrix construction completed!");
            printCMatrixResults(cmatrixResult);

        } catch (Exception e) {
            System.err.println("\nError running pipeline:");
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        System.out.println("================================================================================");
    }

    /**
     * Prints C-matrix results with coefficients in a clean format.
     */
    private static void printCMatrixResults(CMatrixResult result) {
        System.out.println("\n================================================================================");
        System.out.println("                        C-MATRIX IDENTIFICATION RESULT");
        System.out.println("================================================================================");

        int[][] lcv = result.getLcv();
        List<List<int[]>> wcv = result.getWcv();
        List<List<double[][]>> cmat = result.getCmat();
        int[][] cfBasisIndices = result.getCfBasisIndices();

        System.out.println("\nC-MATRIX STRUCTURE");
        System.out.println("-----------------");
        System.out.println("Cluster types: " + lcv.length);
        System.out.println("Total CFs: " + cfBasisIndices.length);

        // Print C-matrix coefficients for each cluster type
        for (int t = 0; t < lcv.length; t++) {
            System.out.println("\n[CLUSTER TYPE t=" + t + "]");
            System.out.println("  Groups: " + lcv[t].length);

            for (int j = 0; j < lcv[t].length; j++) {
                int ncv = lcv[t][j];
                System.out.println("\n  Group j=" + j + " (" + ncv + " cluster variables):");

                // Print C-matrix coefficients
                if (j < cmat.get(t).size()) {
                    double[][] cm = cmat.get(t).get(j);
                    if (cm != null && cm.length > 0) {
                        int rows = cm.length;
                        int cols = cm[0].length;

                        System.out.println("    C-matrix (" + rows + " x " + cols + " coefficients):");

                        // Print header row with CF indices
                        System.out.print("             ");
                        for (int col = 0; col < cols; col++) {
                            System.out.printf("      CF[%d]  ", col);
                        }
                        System.out.println();

                        // Print each CV row
                        for (int i = 0; i < rows; i++) {
                            System.out.printf("    CV[%d]:  ", i);
                            for (int col = 0; col < cols; col++) {
                                System.out.printf("%10.6f  ", cm[i][col]);
                            }
                            System.out.println();
                        }
                    }
                }

                // Print weights
                if (j < wcv.get(t).size()) {
                    int[] weights = wcv.get(t).get(j);
                    System.out.print("    Weights: [");
                    for (int k = 0; k < weights.length; k++) {
                        System.out.print(weights[k]);
                        if (k < weights.length - 1) System.out.print(", ");
                    }
                    System.out.println("]");
                }
            }
        }

        // Print CF basis indices
        System.out.println("\n\nCORRELATION FUNCTION BASIS INDICES");
        System.out.println("----------------------------------");
        for (int col = 0; col < cfBasisIndices.length; col++) {
            int[] bases = cfBasisIndices[col];
            System.out.printf("CF[%d]: [", col);
            for (int k = 0; k < bases.length; k++) {
                System.out.print(bases[k]);
                if (k < bases.length - 1) System.out.print(", ");
            }
            System.out.println("]");
        }

        System.out.println("\n================================================================================\n");
    }
}
