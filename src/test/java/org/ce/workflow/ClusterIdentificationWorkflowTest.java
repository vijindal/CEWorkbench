package org.ce.workflow;

import org.ce.workflow.ClusterIdentificationWorkflow;
import org.ce.workflow.ClusterIdentificationRequest;
import static org.ce.domain.cluster.ClusterPrimitives.*;
import org.ce.domain.cluster.AllClusterData;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertNotNull;

/**
 * Tests for IdentificationPipeline.
 *
 * This test suite validates the complete identification pipeline:
 *     disordered & ordered cluster files
 *     → cluster identification (Stages 1a & 1b)
 *     → CF identification (Stages 2a & 2b)
 */
public class ClusterIdentificationWorkflowTest {

    @Test
    public void testA2B2SystemWithBinaryBasis() {
        // Build the configuration for BCC_A2 binary system
        ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                .disorderedClusterFile("clus/BCC_A2-T.txt")    // Disordered phase (A2) clusters
                .orderedClusterFile("clus/BCC_A2-T.txt")       // Ordered phase clusters
                .disorderedSymmetryGroup("BCC_A2-SG")          // Disordered phase symmetry
                .orderedSymmetryGroup("BCC_A2-SG")             // Ordered phase symmetry
                .transformationMatrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}})  // Identity
                .translationVector(new Vector3D(0,0,0))    // No translation
                .numComponents(2)                           // Binary system
                .structurePhase("BCC_A2")                   // Use BCC_A2 structure (supported basis)
                .build();

        // Run the identification workflow
        AllClusterData result = ClusterIdentificationWorkflow.identify(config);

        // Verify results are not null
        assertNotNull(result, "AllClusterData should not be null");
        assertNotNull(result.getDisorderedClusterResult(), "Disordered cluster result should not be null");
        assertNotNull(result.getOrderedClusterResult(), "Ordered cluster result should not be null");
        assertNotNull(result.getDisorderedCFResult(), "Disordered CF result should not be null");
        assertNotNull(result.getOrderedCFResult(), "Ordered CF result should not be null");

        // Print results
        System.out.println("\n" + result.getSummary());
        result.printResults();
    }
}
