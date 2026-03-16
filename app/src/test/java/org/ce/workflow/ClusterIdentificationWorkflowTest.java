package org.ce.workflow;

import org.ce.workflow.ClusterIdentificationWorkflow;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.domain.cluster.Vector3D;
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
        // Build the configuration for A2 (disordered) to B2 (ordered) system
        ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                .disorderedClusterFile("clus/A2-T.txt")    // Disordered phase clusters
                .orderedClusterFile("clus/A2-T.txt")       // Ordered phase clusters (using A2 for testing)
                .disorderedSymmetryGroup("A2-SG")          // Disordered phase symmetry
                .orderedSymmetryGroup("A2-SG")             // Ordered phase symmetry
                .transformationMatrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}})  // Identity
                .translationVector(new Vector3D(0,0,0))    // No translation
                .numComponents(2)                           // Binary system
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
