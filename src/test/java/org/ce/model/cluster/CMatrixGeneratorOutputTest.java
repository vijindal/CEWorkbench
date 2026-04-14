package org.ce.model.cluster;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Disabled;

import java.util.*;
import java.util.function.Consumer;
import java.nio.file.*;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.storage.Workspace;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Manual inspection test for C-matrix output.
 *
 * Loads real ClusterIdentificationResult and CFIdentificationResult for BCC_A2
 * and calls CMatrixGenerator.runAndPrint() to display formatted C-matrix blocks.
 */
@DisplayName("CMatrixGenerator Output (Manual Inspection)")
class CMatrixGeneratorOutputTest {

    /**
     * Builds a ModelSession for the given system, runs C-matrix output to both
     * stdout and a file, and returns the file path.
     */
    private Path runForSystemAndFile(String elements, String structure, String model, String filename)
            throws Exception {

        System.out.println("\n>>> Running: " + elements + " / " + structure + " / " + model);

        // Build system ID
        Workspace workspace = new Workspace();
        Workspace.SystemId systemId = new Workspace.SystemId(elements, structure, model);

        // Determine number of elements (binary, ternary, quaternary, etc)
        int K = elements.split("-").length;

        // Build ModelSession with a real progress sink (not null)
        org.ce.model.storage.DataStore.HamiltonianStore hamiltonianStore =
            new org.ce.model.storage.DataStore.HamiltonianStore(workspace);
        ModelSession.Builder builder = new ModelSession.Builder(hamiltonianStore);
        ModelSession session = builder.build(
            systemId,
            EngineConfig.CVM,
            line -> { }  // no-op progress sink (silent but never null)
        );

        // Extract results from session
        AllClusterData clusterData = session.clusterData;
        ClusterIdentificationResult clusterResult = clusterData.getOrderedClusterResult();
        CFIdentificationResult cfResult = clusterData.getOrderedCFResult();
        List<Cluster> maxClusters = clusterResult.getOrdClusterData().getCoordList().get(0);

        // Collect output to both list and stdout
        List<String> lines = new ArrayList<>();
        Consumer<String> tee = line -> {
            lines.add(line);
            System.out.println(line);
        };

        // Run and collect output
        CMatrixGenerator.runAndPrint(clusterResult, cfResult, maxClusters, K, tee);

        // Write to file
        Path outFile = Path.of(System.getProperty("java.io.tmpdir"), filename);
        Files.write(outFile, lines, StandardOpenOption.CREATE, StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING);
        System.out.println("  >>> Output written to: " + outFile);

        return outFile;
    }

    @Test
    @Disabled("Requires Nb-Ti-Al Hamiltonian file to be present")
    @DisplayName("Binary BCC_A2 T-model")
    void testBinaryBCCA2T() throws Exception {
        Path outFile = runForSystemAndFile("Nb-Ti", "BCC_A2", "T", "cmat_K2_binary.txt");
        assertTrue(Files.exists(outFile), "Output file should exist: " + outFile);
    }

    @Test
    @Disabled("Requires Nb-Ti-Al Hamiltonian file to be present")
    @DisplayName("Ternary BCC_A2 T-model")
    void testTernaryBCCA2T() throws Exception {
        Path outFile = runForSystemAndFile("Nb-Ti-Al", "BCC_A2", "T", "cmat_K3_ternary.txt");
        assertTrue(Files.exists(outFile), "Output file should exist: " + outFile);
    }

    @Test
    @Disabled("Requires Nb-Ti-Al-Mo Hamiltonian file to be present")
    @DisplayName("Quaternary BCC_A2 T-model")
    void testQuaternaryBCCA2T() throws Exception {
        Path outFile = runForSystemAndFile("Nb-Ti-Al-Mo", "BCC_A2", "T", "cmat_K4_quaternary.txt");
        assertTrue(Files.exists(outFile), "Output file should exist: " + outFile);
    }

    @Test
    @DisplayName("All Systems (prints file paths to console)")
    void printAllSystems() throws Exception {
        System.out.println("\n" + "=".repeat(100));
        System.out.println("RUNNING C-MATRIX GENERATOR TESTS (with available Hamiltonians)");
        System.out.println("=".repeat(100));

        Path nbTi = null;
        Path nbTiV = null;
        Path nbTiVZr = null;

        // Test binary system (Nb-Ti Hamiltonian exists)
        try {
            System.out.println("\n>>> [K=2] Attempting binary system: Nb-Ti / BCC_A2 / T ...");
            nbTi = runForSystemAndFile("Nb-Ti", "BCC_A2", "T", "cmat_K2_Nb-Ti.txt");
            System.out.println("    ✓ SUCCESS - File: " + nbTi);
        } catch (Exception e) {
            System.out.println("    ✗ FAILED: " + e.getClass().getSimpleName() + " - " + e.getMessage());
        }

        // Test ternary system (Nb-Ti-V Hamiltonian exists)
        try {
            System.out.println("\n>>> [K=3] Attempting ternary system: Nb-Ti-V / BCC_A2 / T ...");
            nbTiV = runForSystemAndFile("Nb-Ti-V", "BCC_A2", "T", "cmat_K3_Nb-Ti-V.txt");
            System.out.println("    ✓ SUCCESS - File: " + nbTiV);
        } catch (Exception e) {
            System.out.println("    ✗ FAILED: " + e.getClass().getSimpleName() + " - " + e.getMessage());
        }

        // Test quaternary system (Nb-Ti-V-Zr Hamiltonian exists)
        try {
            System.out.println("\n>>> [K=4] Attempting quaternary system: Nb-Ti-V-Zr / BCC_A2 / T ...");
            nbTiVZr = runForSystemAndFile("Nb-Ti-V-Zr", "BCC_A2", "T", "cmat_K4_Nb-Ti-V-Zr.txt");
            System.out.println("    ✓ SUCCESS - File: " + nbTiVZr);
        } catch (Exception e) {
            System.out.println("    ✗ FAILED: " + e.getClass().getSimpleName() + " - " + e.getMessage());
        }

        System.out.println("\n" + "=".repeat(100));
        System.out.println("OUTPUT FILES CREATED:");
        if (nbTi != null) System.out.println("  [K=2] Nb-Ti:        " + nbTi);
        if (nbTiV != null) System.out.println("  [K=3] Nb-Ti-V:      " + nbTiV);
        if (nbTiVZr != null) System.out.println("  [K=4] Nb-Ti-V-Zr:   " + nbTiVZr);
        System.out.println("=".repeat(100));

        // All three systems should succeed and create files
        assertTrue(nbTi != null && Files.exists(nbTi), "Binary (Nb-Ti) output file should exist");
        assertTrue(nbTiV != null && Files.exists(nbTiV), "Ternary (Nb-Ti-V) output file should exist");
        assertTrue(nbTiVZr != null && Files.exists(nbTiVZr), "Quaternary (Nb-Ti-V-Zr) output file should exist");
    }
}
