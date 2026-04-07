package org.ce.domain.engine.cvm;

import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.domain.result.EquilibriumState;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class CVMEngineOnTheFlyTest {

    @Test
    public void testOnTheFlyBccA2() throws Exception {
        CVMEngine engine = new CVMEngine();
        
        // Setup a minimal CEC for BCC_A2
        CECEntry cec = new CECEntry();
        cec.structurePhase = "BCC_A2";
        cec.model = "T_CVCF"; // Suffix to test sanitization
        cec.elements = "Nb-Ti-V";
        // Minimal terms (just point terms to satisfy basis)
        cec.cecTerms = new CECTerm[0]; 

        double[] composition = {0.33, 0.33, 0.34};
        double temperature = 1000.0;

        ThermodynamicInput input = new ThermodynamicInput(
                cec,
                temperature,
                composition,
                "BCC_A2_test",
                "BCC_A2 Test System",
                System.out::println
        );

        // First run: calls identification
        long start1 = System.currentTimeMillis();
        EquilibriumState state1 = engine.compute(input);
        long end1 = System.currentTimeMillis();
        
        System.out.println("First run (with identification) took: " + (end1 - start1) + " ms");
        assertNotNull(state1);

        // Second run: should use cache
        long start2 = System.currentTimeMillis();
        EquilibriumState state2 = engine.compute(input);
        long end2 = System.currentTimeMillis();
        
        System.out.println("Second run (cached) took: " + (end2 - start2) + " ms");
        assertNotNull(state2);
        
        assertTrue((end2 - start2) < (end1 - start1), "Cached run should be faster");
    }
}
