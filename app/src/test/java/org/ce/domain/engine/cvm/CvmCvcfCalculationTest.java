package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.cvcf.CvCfIntegration;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.Workspace;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;

/**
 * Test to calculate G, H, S and CFs for Nb-Ti BCC_A2 in CVCF basis.
 */
public class CvmCvcfCalculationTest {

    @Test
    @DisplayName("Calculate G, H, S, and CFs for Nb-Ti BCC_A2 CVCF at 1000 K, 0.5 Ti")
    void testCvmCvcfCalculation() throws Exception {
        // System parameters
        String clusterId = "BCC_A2_T_bin";
        String cvcfHamiltonianId = "Nb-Ti_BCC_A2_CVCF";
        double temperature = 1000.0;  // Kelvin
        double tiAtomicFraction = 0.5;
        double nbAtomicFraction = 1.0 - tiAtomicFraction;

        // Load data
        Workspace workspace = new Workspace();
        ClusterDataStore clusterStore = new ClusterDataStore(workspace);
        HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);

        System.out.println("\n" + "=" .repeat(80));
        System.out.println("CVM CALCULATION - Nb-Ti BCC_A2 (CVCF BASIS)");
        System.out.println("=" .repeat(80));
        System.out.println();

        // Load cluster data
        System.out.println("Loading cluster data: " + clusterId);
        AllClusterData allData = clusterStore.load(clusterId);
        System.out.println("  ✓ Loaded successfully");
        System.out.println();

        // Load CVCF ECIs
        System.out.println("Loading CVCF basis ECIs: " + cvcfHamiltonianId);
        CECEntry cvcfEntry = hamiltonianStore.load(cvcfHamiltonianId);
        System.out.println("  ✓ Found " + cvcfEntry.cecTerms.length + " CEC terms");
        System.out.println("  ✓ ncf = " + cvcfEntry.ncf);
        // Debug: Compare ncf values
        int ncfCluster = allData.getDisorderedCFResult().getNcf();
        int ncfHamiltonian = cvcfEntry.ncf;
        System.out.println("DEBUG: ncf (cluster data) = " + ncfCluster);
        System.out.println("DEBUG: ncf (hamiltonian) = " + ncfHamiltonian);

        double[] cvcfEci = extractEciFromEntry(cvcfEntry, temperature);
        // Ensure ECI array is exactly ncf required by cluster data
        if (cvcfEci.length != ncfCluster) {
            double[] trimmed = new double[ncfCluster];
            for (int i = 0; i < ncfCluster && i < cvcfEci.length; i++) {
                trimmed[i] = cvcfEci[i];
            }
            cvcfEci = trimmed;
            System.out.println("  ✓ Trimmed/extended ECIs to match ncf required by cluster data: " + ncfCluster);
        }
        System.out.println("  ✓ Extracted " + cvcfEci.length + " ECIs at " + temperature + " K");
        System.out.println();

        // Print CVCF ECI values
        System.out.println("CVCF BASIS ECI VALUES (J/mol):");
        System.out.println("-" .repeat(80));
        System.out.println("DEBUG: cvcfEci.length = " + cvcfEci.length);
        System.out.println("DEBUG: cecTerms.length = " + cvcfEntry.cecTerms.length);
        System.out.print("DEBUG: cvcfEci = [");
        for (int i = 0; i < cvcfEci.length; i++) {
            System.out.print(cvcfEci[i]);
            if (i < cvcfEci.length - 1) System.out.print(", ");
        }
        System.out.println("]");
        System.out.print("DEBUG: cecTerms = [");
        for (int i = 0; i < cvcfEntry.cecTerms.length; i++) {
            System.out.print(cvcfEntry.cecTerms[i].name);
            if (i < cvcfEntry.cecTerms.length - 1) System.out.print(", ");
        }
        System.out.println("]");
        for (int i = 0; i < cvcfEci.length; i++) {
            String name = (i < cvcfEntry.cecTerms.length) ? cvcfEntry.cecTerms[i].name : ("ECI[" + i + "]");
            System.out.printf("  %s: %15.6f%n", name, cvcfEci[i]);
        }
        System.out.println();

        // Create CVCF integration
        System.out.println("Creating CVCF basis integration...");
        CvCfIntegration cvcfIntegration = CvCfIntegration.forBccA2Binary(
            allData.getCMatrixResult(),
            "Nb",
            "Ti"
        );
        System.out.println("  ✓ CVCF integration created");
        System.out.println("  ✓ Component pairs: " + cvcfIntegration.getComponentPairs().size());
        System.out.println();

        // Create CVMPhaseModel with CVCF integration
        System.out.println("Creating CVMPhaseModel with CVCF basis...");
        CVMPhaseModel.CVMInput input = new CVMPhaseModel.CVMInput(
            allData.getDisorderedClusterResult(),
            allData.getDisorderedCFResult(),
            allData.getCMatrixResult(),
            "Nb-Ti_BCC_A2_CVCF",
            "Nb-Ti Binary BCC_A2 (CVCF)",
            2  // numComponents
        );

        double[] composition = {nbAtomicFraction, tiAtomicFraction};

        System.out.println("Debug information:");
        System.out.println("  ECI length: " + cvcfEci.length);
        System.out.println("  Composition length: " + composition.length);
        System.out.println("  Temperature: " + temperature);

        CVMPhaseModel cvcfModel;
        try {
            // Create model with CVCF integration
            cvcfModel = CVMPhaseModel.create(input, cvcfEci, temperature, composition, cvcfIntegration);
            System.out.println("  ✓ Model created with CVCF basis and minimized");
        } catch (Exception e) {
            System.out.println("  ✗ Error creating CVCF model: " + e.getMessage());
            throw e;
        }
        System.out.println();

        // Print results
        System.out.println("THERMODYNAMIC PROPERTIES AT " + temperature + " K:");
        System.out.println("-" .repeat(80));
        System.out.printf("Composition: Nb = %.4f, Ti = %.4f%n", nbAtomicFraction, tiAtomicFraction);
        System.out.println();

        double G = cvcfModel.getEquilibriumG();
        double H = cvcfModel.getEquilibriumH();
        double S = cvcfModel.getEquilibriumS();
        double[] cfs = cvcfModel.getEquilibriumCFs();
        double gradientNorm = cvcfModel.getGradientNorm();

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
        System.out.println("CORRELATION FUNCTIONS (CFs) - CVCF BASIS:");
        System.out.println("-" .repeat(80));
        if (cfs != null) {
            System.out.printf("Number of CFs: %d%n", cfs.length);
            System.out.println();
            String[] cfNames = {"v2NbTi1", "v2NbTi2", "v3NbTi", "v4NbTi"};
            for (int i = 0; i < cfs.length; i++) {
                String name = (i < cfNames.length) ? cfNames[i] : "CF[" + i + "]";
                System.out.printf("  %s = %+20.15f%n", name, cfs[i]);
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
