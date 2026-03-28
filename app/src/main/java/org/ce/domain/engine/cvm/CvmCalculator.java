package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.cvcf.BccA2TModelCvCfTransformations;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.Workspace;

/**
 * Utility for calculating CVM thermodynamic properties.
 */
public class CvmCalculator {

    /**
     * Calculate G, H, S and CFs for a binary system at given conditions.
     */
    public static void main(String[] args) throws Exception {
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

        System.out.println("=" .repeat(80));
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
        double[] eci = extractEciFromEntry(hamiltonianEntry, temperature);
        System.out.println("  ✓ Extracted " + eci.length + " non-point ECIs at " + temperature + " K");
        System.out.println();

        // Create CVMPhaseModel
        System.out.println("Creating CVMPhaseModel...");
        CVMPhaseModel.CVMInput input = new CVMPhaseModel.CVMInput(
            allData.getDisorderedClusterResult(),
            allData.getDisorderedCFResult(),
            allData.getCMatrixResult(),
            "Nb-Ti_BCC_A2",
            "Nb-Ti Binary BCC_A2",
            2,  // numComponents
            BccA2TModelCvCfTransformations.binaryBasis()
        );

        double[] composition = {nbAtomicFraction, tiAtomicFraction};
        CVMPhaseModel model = CVMPhaseModel.create(input, eci, temperature, composition);
        System.out.println("  ✓ Model created successfully");
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

        System.out.printf("Gibbs Energy (G):     %15.6f J/mol%n", G);
        System.out.printf("Enthalpy (H):         %15.6f J/mol%n", H);
        System.out.printf("Entropy (S):          %15.6f J/(mol·K)%n", S);
        System.out.printf("Gradient Norm:        %15.6e%n", gradientNorm);
        System.out.println();

        // Print CFs
        System.out.println("CORRELATION FUNCTIONS (CFs):");
        System.out.println("-" .repeat(80));
        if (cfs != null) {
            System.out.printf("Number of CFs: %d%n", cfs.length);
            for (int i = 0; i < cfs.length; i++) {
                System.out.printf("  CF[%d] = %15.10f%n", i, cfs[i]);
            }
        } else {
            System.out.println("  CFs not available");
        }
        System.out.println();

        // Verify convergence
        boolean converged = gradientNorm < 1e-4;
        System.out.println("CONVERGENCE CHECK:");
        System.out.println("-" .repeat(80));
        System.out.printf("Gradient norm < 1e-4: %s%n", converged ? "✓ YES" : "✗ NO");
        System.out.println();

        System.out.println("=" .repeat(80));
    }

    /**
     * Extract ECIs from CECEntry at given temperature.
     * Formula: J(T) = a + b*T for each CEC term
     */
    private static double[] extractEciFromEntry(CECEntry entry, double temperature) {
        // Skip point CF (index 0), extract non-point CFs
        double[] eci = new double[entry.cecTerms.length - 1];
        for (int i = 1; i < entry.cecTerms.length; i++) {
            CECTerm term = entry.cecTerms[i];
            eci[i - 1] = term.a + term.b * temperature;
        }
        return eci;
    }
}
