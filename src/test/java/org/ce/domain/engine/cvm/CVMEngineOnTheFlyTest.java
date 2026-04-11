package org.ce.domain.engine.cvm;

import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cluster.AllClusterData;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CVMSolver;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.cluster.cvcf.CvCfBasis;
import org.ce.model.hamiltonian.CECEvaluator;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class CVMEngineOnTheFlyTest {

    @Test
    public void testOnTheFlyBccA2() throws Exception {
        // Setup a minimal CEC for BCC_A2
        CECEntry cec = new CECEntry();
        cec.structurePhase = "BCC_A2";
        cec.model = "T";
        cec.elements = "Nb-Ti-V";
        cec.cecTerms = new CECEntry.CECTerm[0];

        double[] composition = {0.33, 0.33, 0.34};
        double temperature = 1000.0;

        // Perform identification
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

        // Get basis and build model
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", "T", 3);
        CVMGibbsModel model = new CVMGibbsModel(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                clusterData.getCMatrixResult(),
                basis);

        // Evaluate ECI
        double[] eci = CECEvaluator.evaluate(cec, temperature, basis, "CVM");

        // Run solver (first time)
        long start1 = System.currentTimeMillis();
        CVMSolver solver = new CVMSolver();
        CVMSolver.EquilibriumResult result1 = solver.minimize(model, composition, temperature, eci, 1.0e-5);
        long end1 = System.currentTimeMillis();

        System.out.println("First run (with identification) took: " + (end1 - start1) + " ms");
        assertNotNull(result1);
        assertTrue(result1.converged, "Solver should converge");

        // Run solver again (should be faster due to cached identification)
        long start2 = System.currentTimeMillis();
        CVMSolver.EquilibriumResult result2 = solver.minimize(model, composition, temperature, eci, 1.0e-5);
        long end2 = System.currentTimeMillis();

        System.out.println("Second run (cached) took: " + (end2 - start2) + " ms");
        assertNotNull(result2);
        assertTrue(result2.converged, "Solver should converge");
    }
}
