package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.Workspace;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;

/**
 * Test to calculate G, H, S and CFs for Nb-Ti BCC_A2 at given conditions.
 */
public class CvmCalculationTest {

    @Test
    @DisplayName("Calculate G, H, S, and CFs for Nb-Ti BCC_A2 at 1000 K, 0.5 Ti")
    void testCvmCalculation() throws Exception {
        // System parameters
        String clusterId = "BCC_A2_T_bin";
        String hamiltonianId = "Nb-Ti_BCC_A2_T";
        double temperature = 1000.0;  // Kelvin
        double tiAtomicFraction = 0.5;
        double nbAtomicFraction = 1.0 - tiAtomicFraction;

        // Load data
        Workspace workspace = new Workspace();
        ClusterDataStore clusterStore = new ClusterDataStore(workspace);
        HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);

        System.out.println("\n" + "=" .repeat(80));
        System.out.println("CVM CALCULATION - Nb-Ti BCC_A2");
        System.out.println("=" .repeat(80));
        System.out.println();

        // Load cluster data
        System.out.println("Loading cluster data: " + clusterId);
        AllClusterData allData = clusterStore.load(clusterId);
        System.out.println("  ✓ Loaded successfully");
        System.out.println();

        // Load ECIs
        System.out.println("Loading ECIs: " + hamiltonianId);
        CECEntry hamiltonianEntry = hamiltonianStore.load(hamiltonianId);
        System.out.println("  ✓ Found " + hamiltonianEntry.cecTerms.length + " CEC terms");
        System.out.println("  ✓ ncf = " + hamiltonianEntry.ncf);
        double[] eci = extractEciFromEntry(hamiltonianEntry, temperature);
        System.out.println("  ✓ Extracted " + eci.length + " ECIs at " + temperature + " K");
        System.out.println();

        // Print ECI values
        System.out.println("ECI VALUES (J/mol):");
        System.out.println("-" .repeat(80));
        for (int i = 0; i < eci.length && i < hamiltonianEntry.cecTerms.length - 1; i++) {
            CECTerm term = hamiltonianEntry.cecTerms[i + 1];
            System.out.printf("  %s: %.6f%n", term.name, eci[i]);
        }
        System.out.println();

        // Create CVMPhaseModel
        System.out.println("Creating CVMPhaseModel...");
        CVMPhaseModel.CVMInput input = new CVMPhaseModel.CVMInput(
            allData.getDisorderedClusterResult(),
            allData.getDisorderedCFResult(),
            allData.getCMatrixResult(),
            "Nb-Ti_BCC_A2",
            "Nb-Ti Binary BCC_A2",
            2  // numComponents
        );

        double[] composition = {nbAtomicFraction, tiAtomicFraction};

        System.out.println("Debugging information:");
        System.out.println("  ECI length: " + eci.length);
        System.out.println("  Composition length: " + composition.length);
        System.out.println("  Temperature: " + temperature);

        CVMPhaseModel model;
        try {
            model = CVMPhaseModel.create(input, eci, temperature, composition);
            System.out.println("  ✓ Model created and minimized");
        } catch (IllegalArgumentException e) {
            System.out.println("  ✗ Error creating model: " + e.getMessage());
            throw e;
        }
        System.out.println();

        // Print results
        System.out.println("THERMODYNAMIC PROPERTIES AT " + temperature + " K:");
        System.out.println("-" .repeat(80));
        System.out.printf("Composition: Nb = %.4f, Ti = %.4f%n", nbAtomicFraction, tiAtomicFraction);
        System.out.println();

        double G = model.getEquilibriumG();
        double H = model.getEquilibriumH();
        double S = model.getEquilibriumS();
        double[] cfs = model.getEquilibriumCFs();
        double gradientNorm = model.getGradientNorm();

        System.out.printf("Gibbs Energy (G):     %20.10f J/mol%n", G);
        System.out.printf("Enthalpy (H):         %20.10f J/mol%n", H);
        System.out.printf("Entropy (S):          %20.10f J/(mol·K)%n", S);
        System.out.printf("Gradient Norm:        %20.10e%n", gradientNorm);
        System.out.println();

        // Verify thermodynamic relation
        double calculatedG = H - temperature * S;
        System.out.println("VERIFICATION:");
        System.out.println("-" .repeat(80));
        System.out.printf("H - T·S (calculated): %20.10f J/mol%n", calculatedG);
        System.out.printf("Direct G:             %20.10f J/mol%n", G);
        System.out.printf("Difference:           %20.10e J/mol%n", Math.abs(G - calculatedG));
        System.out.println();

        // Print CFs
        System.out.println("CORRELATION FUNCTIONS (CFs):");
        System.out.println("-" .repeat(80));
        if (cfs != null) {
            System.out.printf("Number of CFs: %d%n", cfs.length);
            System.out.println();
            for (int i = 0; i < cfs.length; i++) {
                System.out.printf("  CF[%d] = %+20.15f%n", i, cfs[i]);
            }
        } else {
            System.out.println("  CFs not available");
        }
        System.out.println();

        // Verify convergence
        boolean converged = gradientNorm < 1e-4;
        System.out.println("CONVERGENCE CHECK:");
        System.out.println("-" .repeat(80));
        System.out.printf("Gradient norm < 1e-4: %s (norm = %.10e)%n",
            converged ? "✓ YES" : "✗ NO", gradientNorm);
        System.out.println();

        System.out.println("=" .repeat(80));
    }

    /**
     * Extract ECIs from CECEntry at given temperature.
     * Formula: J(T) = a + b*T for each CEC term
     *
     * Note: Extract all cecTerms as-is. The ncf field in CECEntry
     * specifies the total number of non-point CFs.
     */
    private double[] extractEciFromEntry(CECEntry entry, double temperature) {
        double[] eci = new double[entry.ncf];
        for (int i = 0; i < entry.ncf && i < entry.cecTerms.length; i++) {
            CECTerm term = entry.cecTerms[i];
            eci[i] = term.a + term.b * temperature;
        }
        return eci;
    }
}
