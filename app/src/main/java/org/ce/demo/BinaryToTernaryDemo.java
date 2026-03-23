package org.ce.demo;

import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.domain.hamiltonian.NumericalCECTransformer;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.Workspace;
import org.ce.workflow.cec.CECManagementWorkflow;

/**
 * Demonstration: Transform binary CECs to ternary basis
 *
 * Example: Transform Nb-Ti binary (BCC_A2_T) to Nb-Ti-V ternary basis
 */
public class BinaryToTernaryDemo {

    public static void main(String[] args) throws Exception {
        System.out.println("═══════════════════════════════════════════════════════════");
        System.out.println("  CEC Basis Transformation: Binary → Ternary");
        System.out.println("═══════════════════════════════════════════════════════════\n");

        // Initialize storage
        Workspace workspace = new Workspace();
        HamiltonianStore store = new HamiltonianStore(workspace);
        CECManagementWorkflow workflow = new CECManagementWorkflow(store, null);

        // ───────────────────────────────────────────────────────────────────
        // Example 1: Transform Nb-Ti → Nb-Ti-V
        // ───────────────────────────────────────────────────────────────────
        System.out.println("Example 1: Nb-Ti → Nb-Ti-V (Ternary)\n");
        System.out.println("Binary system: Nb-Ti_BCC_A2_T");
        System.out.println("  Species: {Nb=0, Ti=1}");
        System.out.println("\nTernary system: Nb-Ti-V_BCC_A2_T");
        System.out.println("  Species: {Nb=0, Ti=1, V=2}");
        System.out.println("\nSpecies mapping: Nb→Nb, Ti→Ti (V is new)\n");

        // Load binary Hamiltonian
        CECEntry binaryEntry = store.load("Nb-Ti_BCC_A2_T");
        System.out.println("Loaded binary Hamiltonian:");
        System.out.println("  Elements: " + binaryEntry.elements);
        System.out.println("  Structure: " + binaryEntry.structurePhase);
        System.out.println("  NCF: " + binaryEntry.ncf);
        System.out.println("  Reference: " + binaryEntry.reference);
        System.out.println();

        // Display binary CEC values
        System.out.println("Binary ECI values (a + b*T):");
        for (int i = 0; i < binaryEntry.cecTerms.length; i++) {
            CECTerm term = binaryEntry.cecTerms[i];
            String bStr = (term.b == 0) ? "" : String.format(" + %.4f*T", term.b);
            System.out.printf("  %s: %.2f%s", term.name, term.a, bStr);
            if (term.description != null) {
                System.out.print("  [" + term.description + "]");
            }
            System.out.println();
        }
        System.out.println();

        // Create species mapping
        NumericalCECTransformer.SpeciesMapping mapping =
            new NumericalCECTransformer.SpeciesMapping();
        mapping.addMapping(0, 0);  // Nb in ternary → Nb in binary
        mapping.addMapping(1, 1);  // Ti in ternary → Ti in binary
        // V (index 2) has no mapping in binary

        // Transform
        try {
            CECEntry ternaryEntry = workflow.transformBinaryToTernary(
                "Nb-Ti_BCC_A2_T",
                "Nb-Ti-V_BCC_A2_T",
                "Nb-Ti-V",
                mapping
            );

            System.out.println("✓ Transformation successful!\n");
            System.out.println("Ternary ECI values (transformed to 3-species basis):");
            for (int i = 0; i < ternaryEntry.cecTerms.length; i++) {
                CECTerm term = ternaryEntry.cecTerms[i];
                String bStr = (term.b == 0) ? "" : String.format(" + %.4f*T", term.b);
                System.out.printf("  %s: %.2f%s", term.name, term.a, bStr);
                if (term.description != null) {
                    System.out.print("  [" + term.description + "]");
                }
                System.out.println();
            }
            System.out.println();
            System.out.println("Saved to: hamiltonians/" + ternaryEntry.elements +
                             "_" + ternaryEntry.structurePhase +
                             "_" + ternaryEntry.model + "/hamiltonian.json");

        } catch (Exception e) {
            System.out.println("✗ Transformation failed: " + e.getMessage());
            e.printStackTrace();
        }

        // ───────────────────────────────────────────────────────────────────
        // Example 2: Transform Ti-V → Nb-Ti-V
        // ───────────────────────────────────────────────────────────────────
        System.out.println("\n\n" + "═".repeat(63));
        System.out.println("Example 2: Ti-V → Nb-Ti-V (Ternary)\n");
        System.out.println("Binary system: Ti-V_BCC_A2_T");
        System.out.println("  Species: {Ti=0, V=1}");
        System.out.println("\nTernary system: Nb-Ti-V_BCC_A2_T");
        System.out.println("  Species: {Nb=0, Ti=1, V=2}");
        System.out.println("\nSpecies mapping: Ti→0 (in Ti-V), V→1 (in Ti-V)\n");

        // Load binary Hamiltonian
        CECEntry binaryEntry2 = store.load("Ti-V_BCC_A2_T");
        System.out.println("Loaded binary Hamiltonian:");
        System.out.println("  Elements: " + binaryEntry2.elements);
        System.out.println("  NCF: " + binaryEntry2.ncf);
        System.out.println();

        System.out.println("Binary ECI values (a + b*T):");
        for (int i = 0; i < binaryEntry2.cecTerms.length; i++) {
            CECTerm term = binaryEntry2.cecTerms[i];
            String bStr = (term.b == 0) ? "" : String.format(" + %.4f*T", term.b);
            System.out.printf("  %s: %.2f%s", term.name, term.a, bStr);
            if (term.description != null) {
                System.out.print("  [" + term.description + "]");
            }
            System.out.println();
        }
        System.out.println();

        // Create species mapping for Ti-V
        NumericalCECTransformer.SpeciesMapping mapping2 =
            new NumericalCECTransformer.SpeciesMapping();
        mapping2.addMapping(1, 0);  // Ti in ternary (index 1) → Ti in Ti-V (index 0)
        mapping2.addMapping(2, 1);  // V in ternary (index 2) → V in Ti-V (index 1)
        // Nb (index 0) has no mapping in Ti-V

        try {
            CECEntry ternaryEntry2 = workflow.transformBinaryToTernary(
                "Ti-V_BCC_A2_T",
                "Nb-Ti-V_BCC_A2_T_v2",
                "Nb-Ti-V",
                mapping2
            );

            System.out.println("✓ Transformation successful!\n");
            System.out.println("Ternary ECI values (transformed to 3-species basis):");
            for (int i = 0; i < ternaryEntry2.cecTerms.length; i++) {
                CECTerm term = ternaryEntry2.cecTerms[i];
                String bStr = (term.b == 0) ? "" : String.format(" + %.4f*T", term.b);
                System.out.printf("  %s: %.2f%s", term.name, term.a, bStr);
                if (term.description != null) {
                    System.out.print("  [" + term.description + "]");
                }
                System.out.println();
            }
            System.out.println();

        } catch (Exception e) {
            System.out.println("✗ Transformation failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n" + "═".repeat(63));
        System.out.println("Demonstration complete!");
    }
}
