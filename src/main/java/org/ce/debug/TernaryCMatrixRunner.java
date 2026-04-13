package org.ce.debug;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import org.ce.model.cluster.*;
import org.ce.model.storage.InputLoader;

/**
 * Diagnostic runner for Nb-Ti-V ternary C-matrix generation.
 */
public class TernaryCMatrixRunner {

    public static void main(String[] args) {
        String structure = "BCC_A2";
        String model = "T";
        int K = 4; // Ternary Nb-Ti-V

        System.out.println("=== Ternary C-Matrix Diagnostic Run: Nb-Ti-V / BCC_A2 / T ===");

        try {
            // 1. Setup Resource Paths
            Path inputsDir = Paths.get("data", "CEWorkbench", "inputs");
            String clusterFile = "clus/" + structure + "-" + model + ".txt";
            String symGroup = structure + "-SG";

            System.out.println("Loading clusters from: " + inputsDir.resolve(clusterFile));
            System.out.println("Loading symmetry from: " + inputsDir.resolve("sym").resolve(symGroup + ".txt"));

            // 2. Load basic inputs to build the Request
            List<Cluster> disorderedClusters = InputLoader.parseClusterFileFromPath(inputsDir, clusterFile);
            List<Cluster> orderedClusters = disorderedClusters; // Same for A2
            SpaceGroup sg = InputLoader.parseSpaceGroupFromPath(inputsDir, symGroup);

            ClusterIdentificationRequest request = ClusterIdentificationRequest.builder()
                    .numComponents(K)
                    .structurePhase(structure)
                    .model(model)
                    .disorderedClusterFile(clusterFile)
                    .orderedClusterFile(clusterFile)
                    .disorderedSymmetryGroup(symGroup)
                    .orderedSymmetryGroup(symGroup)
                    .transformationMatrix(sg.getRotateMat())
                    .translationVector(new org.ce.model.cluster.ClusterPrimitives.Vector3D(
                            sg.getTranslateMat()[0], sg.getTranslateMat()[1], sg.getTranslateMat()[2]))
                    .build();

            // 3. Run Identification Stages 1 and 2
            System.out.println("\n[Stage 1] Identifying Clusters...");
            ClusterIdentificationResult clusterResult = ClusterIdentifier.identify(
                    disorderedClusters,
                    sg.getOperations(),
                    orderedClusters,
                    sg.getOperations(),
                    request.getTransformationMatrix(),
                    new double[] { request.getTranslationVector().getX(),
                            request.getTranslationVector().getY(),
                            request.getTranslationVector().getZ() });

            System.out.println("[Stage 2] Identifying CF Orbits...");
            CFIdentificationResult cfResult = CFIdentifier.identify(
                    clusterResult,
                    clusterResult.getDisClusterData().getClusCoordList(),
                    sg.getOperations(),
                    orderedClusters,
                    sg.getOperations(),
                    sg.getRotateMat(),
                    sg.getTranslateMat(),
                    K);

            // 4. Run the Instrumented C-Matrix Generation (Stage 3)
            System.out.println("\n[Stage 3] Executing Instrumented C-Matrix Generation...");
            List<Cluster> maxClusters = disorderedClusters;
            CMatrixGenerator.runAndPrint(
                    clusterResult,
                    cfResult,
                    maxClusters,
                    K,
                    System.out::println);

            System.out.println("\n=== Run Complete ===");

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
