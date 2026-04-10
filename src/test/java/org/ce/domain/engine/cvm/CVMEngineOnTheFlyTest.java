package org.ce.domain.engine.cvm;

import org.ce.calculation.engine.ThermodynamicEngine;
import org.ce.calculation.engine.cvm.CVMEngine;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cluster.AllClusterData;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.ThermodynamicResult;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class CVMEngineOnTheFlyTest {

    @Test
    public void testOnTheFlyBccA2() throws Exception {
        CVMEngine engine = new CVMEngine();
        
        // Setup a minimal CEC for BCC_A2
        CECEntry cec = new CECEntry();
        cec.structurePhase = "BCC_A2";
        cec.model = "T"; 
        cec.elements = "Nb-Ti-V";
        cec.cecTerms = new CECEntry.CECTerm[0]; 

        double[] composition = {0.33, 0.33, 0.34};
        double temperature = 1000.0;

        // Perform identification (required by new ThermodynamicInput design)
        ClusterIdentificationRequest idRequest = ClusterIdentificationRequest.builder()
                .disorderedClusterFile("clus/BCC_A2-T.txt")
                .orderedClusterFile("clus/BCC_A2-T.txt")
                .disorderedSymmetryGroup("BCC_A2-SG")
                .orderedSymmetryGroup("BCC_A2-SG")
                .numComponents(3)
                .structurePhase("BCC_A2")
                .model("T")
                .build();
        AllClusterData clusterData = AllClusterData.identify(idRequest);

        ThermodynamicEngine.Input input = new ThermodynamicEngine.Input(
                clusterData,
                cec,
                temperature,
                composition,
                "BCC_A2_test",
                "BCC_A2 Test System",
                System.out::println,
                null,   // eventSink
                4,      // mcsL
                1000,   // nEquil
                2000    // nAvg
        );

        // First run: calls identification
        long start1 = System.currentTimeMillis();
        ThermodynamicResult result1 = engine.compute(input);
        long end1 = System.currentTimeMillis();

        System.out.println("First run (with identification) took: " + (end1 - start1) + " ms");
        assertNotNull(result1);

        // Second run: should use cache
        long start2 = System.currentTimeMillis();
        ThermodynamicResult result2 = engine.compute(input);
        long end2 = System.currentTimeMillis();

        System.out.println("Second run (cached) took: " + (end2 - start2) + " ms");
        assertNotNull(result2);

        assertTrue((end2 - start2) < (end1 - start1), "Cached run should be faster");
    }
}
