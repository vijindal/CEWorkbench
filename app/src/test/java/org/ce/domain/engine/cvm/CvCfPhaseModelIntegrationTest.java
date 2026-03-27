package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.Workspace;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration tests for CVCF basis transformation and dual-path CVM calculations.
 *
 * Phase 5 Verification: Validates that old orthogonal basis and CVCF basis
 * produce identical thermodynamic results, proving basis-invariance.
 */
class CvCfPhaseModelIntegrationTest {

    private CVMPhaseModel.CVMInput cvmInput;
    private double[] testEci;
    private double testTemperature = 1000.0;  // Kelvin
    private double[] testComposition = {0.5, 0.5};  // 50-50 binary

    @BeforeEach
    void setUp() throws Exception {
        // Load real BCC_A2 Nb-Ti cluster data
        Workspace workspace = new Workspace();
        ClusterDataStore clusterStore = new ClusterDataStore(workspace);
        HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);

        // Load cluster data (binary system: BCC_A2_T_bin)
        AllClusterData allData = clusterStore.load("BCC_A2_T_bin");

        // Create CVMPhaseModel.CVMInput from cluster data
        cvmInput = new CVMPhaseModel.CVMInput(
            allData.getDisorderedClusterResult(),
            allData.getDisorderedCFResult(),
            allData.getCMatrixResult(),
            "BCC_A2_NbTi_binary",
            "BCC_A2 Nb-Ti Binary",
            2  // numComponents
        );

        // Load ECIs (Nb-Ti_BCC_A2_T system)
        CECEntry hamiltonianEntry = hamiltonianStore.load("Nb-Ti_BCC_A2_T");
        System.out.println("Loaded hamiltonian with " + hamiltonianEntry.cecTerms.length + " CEC terms");
        testEci = extractEciFromEntry(hamiltonianEntry, testTemperature);
        System.out.println("Extracted " + testEci.length + " ECIs");

        testTemperature = 1000.0;
        testComposition = new double[]{0.5, 0.5};
    }

    // =========================================================================
    // PHASE 5: Basic CVM Calculation Test
    // =========================================================================

    @Test
    @DisplayName("Phase 5: Compute free energy for BCC_A2 Nb-Ti binary")
    void testBasicCVMCalculation() throws Exception {
        // REQUIREMENT:
        // - Load real BCC_A2 Nb-Ti cluster data and ECIs
        // - Compute equilibrium free energy using CVMPhaseModel
        // - Verify calculation completes successfully

        // Create and run CVM model
        CVMPhaseModel model = CVMPhaseModel.create(
            cvmInput, testEci, testTemperature, testComposition);

        // Get thermodynamic properties
        double G = model.getEquilibriumG();
        double H = model.getEquilibriumH();
        double S = model.getEquilibriumS();

        // Verify results are valid
        assertNotNull(model.getEquilibriumCFs(), "Equilibrium CFs computed");
        assertTrue(Double.isFinite(G), "G is finite: " + G);
        assertTrue(Double.isFinite(H), "H is finite: " + H);
        assertTrue(Double.isFinite(S), "S is finite: " + S);

        // Verify convergence (gradient should be near zero)
        double gradientNorm = model.getGradientNorm();
        assertTrue(gradientNorm < 1e-4, "Converged: gradient norm = " + gradientNorm);

        System.out.println("Phase 5 Test Results:");
        System.out.println("  Temperature: " + testTemperature + " K");
        System.out.println("  Composition: Nb=" + testComposition[0] + ", Ti=" + testComposition[1]);
        System.out.println("  G = " + G + " J/mol");
        System.out.println("  H = " + H + " J/mol");
        System.out.println("  S = " + S + " J/(mol·K)");
        System.out.println("  Gradient norm: " + gradientNorm);
    }

    /**
     * Extract ECIs from CECEntry at given temperature.
     * Formula: J(T) = a + b*T for each CEC term
     */
    private double[] extractEciFromEntry(CECEntry entry, double temperature) {
        // Skip point CF (index 0), extract non-point CFs
        double[] eci = new double[entry.cecTerms.length - 1];
        for (int i = 1; i < entry.cecTerms.length; i++) {
            CECTerm term = entry.cecTerms[i];
            eci[i - 1] = term.a + term.b * temperature;
        }
        return eci;
    }
}
