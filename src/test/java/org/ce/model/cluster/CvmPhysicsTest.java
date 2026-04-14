package org.ce.model.cluster;

import org.ce.CEWorkbenchContext;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.CvmPhysicsVerifier;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.function.Consumer;

import static org.junit.jupiter.api.Assertions.*;

/**
 * JUnit entry point for CVM Physics Verification.
 * This test runs the identification pipeline for BCC_A2 systems
 * (both Binary and Ternary) and invokes the CvmPhysicsVerifier.
 * 
 * Results are saved to 'physics_verification_report.txt'.
 */
public class CvmPhysicsTest {

    @Test
    public void testBccA2Physics() {
        String reportFile = "physics_verification_report.txt";
        System.out.println("Running testBccA2Physics (Binary + Ternary)...");
        System.out.println("Results will be saved to: " + reportFile);
        
        try (PrintWriter writer = new PrintWriter(new FileWriter(reportFile))) {
            // Compound logger: writes to both console and file
            Consumer<String> logger = msg -> {
                System.out.println(msg);
                writer.println(msg);
            };

            CvmPhysicsVerifier verifier = new CvmPhysicsVerifier(logger);

            // =================================================================
            // PART 1: BINARY BCC_A2
            // =================================================================
            logger.accept("### PART 1: BINARY BCC_A2 SYSTEM ###");
            ClusterIdentificationRequest binaryRequest = ClusterIdentificationRequest.builder()
                    .structurePhase("BCC_A2")
                    .model("T")
                    .numComponents(2)
                    .build();
            
            AllClusterData binaryData = AllClusterData.identify(binaryRequest, logger);
            assertNotNull(binaryData, "Binary identification failed");
            
            verifier.verifyAll(binaryData);
            writer.println("\n");
            writer.flush();

            // =================================================================
            // PART 2: TERNARY Nb-Ti-V (1000K)
            // =================================================================
            logger.accept("\n### PART 2: TERNARY Nb-Ti-V SYSTEM ###");
            ClusterIdentificationRequest ternaryRequest = ClusterIdentificationRequest.builder()
                    .structurePhase("BCC_A2")
                    .model("T")
                    .numComponents(3)
                    .build();
            
            AllClusterData ternaryData = AllClusterData.identify(ternaryRequest, logger);
            assertNotNull(ternaryData, "Ternary identification failed");

            // Requested composition: 0.33, 0.33, 0.34
            double[][] ternaryComps = {{0.33, 0.33, 0.34}};
            verifier.verifyAll(ternaryData, "Nb-Ti-V", 1000.0, ternaryComps);

            logger.accept("\nFull physics verification suite completed successfully.");
            writer.flush();
            
        } catch (Exception e) {
            fail("Physics verification failed with exception: " + e.getMessage());
        }
    }
}
