package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.ClusterVariableEvaluator;
import org.ce.domain.cluster.cvcf.BccA2CvCfTransformations;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.cluster.cvcf.CvCfBasisTransformer;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.hamiltonian.CECTerm;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.Workspace;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;

import java.util.List;

/**
 * Detailed trace of G expression in orthogonal vs CVCF basis.
 */
public class CVMGExpressionTest {

    @Test
    @DisplayName("Trace full G expression: Orthogonal vs CVCF basis")
    void testGExpression() throws Exception {
        Workspace workspace = new Workspace();
        ClusterDataStore clusterStore = new ClusterDataStore(workspace);
        HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);

        System.out.println("\n" + "=" .repeat(100));
        System.out.println("G EXPRESSION ANALYSIS: ORTHOGONAL vs CVCF BASIS");
        System.out.println("=" .repeat(100));
        System.out.println();

        // Load data
        AllClusterData allData = clusterStore.load("BCC_A2_T_bin");
        CECEntry orthoEntry = hamiltonianStore.load("Nb-Ti_BCC_A2_T");
        CECEntry cvcfEntry = hamiltonianStore.load("Nb-Ti_BCC_A2_CVCF");

        double temperature = 1000.0;
        double[] composition = {0.5, 0.5};

        System.out.println("SYSTEM PARAMETERS:");
        System.out.println("-" .repeat(100));
        System.out.printf("Temperature: %.1f K%n", temperature);
        System.out.printf("Composition: %.2f, %.2f%n", composition[0], composition[1]);
        System.out.println();

        // Extract ECIs
        double[] eciOrtho = extractEciFromEntry(orthoEntry, temperature);
        double[] eciCvcf = extractEciFromEntry(cvcfEntry, temperature);

        System.out.println("ORTHOGONAL BASIS ECIs:");
        System.out.println("-" .repeat(100));
        for (int i = 0; i < eciOrtho.length && i < orthoEntry.cecTerms.length; i++) {
            System.out.printf("  ECI[%d] (%s): %12.2f J/mol%n", i, orthoEntry.cecTerms[i].name, eciOrtho[i]);
        }
        System.out.println();

        System.out.println("CVCF BASIS ECIs:");
        System.out.println("-" .repeat(100));
        for (int i = 0; i < eciCvcf.length && i < cvcfEntry.cecTerms.length; i++) {
            System.out.printf("  ECI[%d] (%s): %12.2f J/mol%n", i, cvcfEntry.cecTerms[i].name, eciCvcf[i]);
        }
        System.out.println();

        // Test u vector
        double[] u = {0.01, 0.01, 0.01, 0.01};  // Small test CF values
        System.out.println("TEST U VECTOR (independent CFs):");
        System.out.println("-" .repeat(100));
        for (int i = 0; i < u.length; i++) {
            System.out.printf("  u[%d] = %.4f%n", i, u[i]);
        }
        System.out.println();

        // Compute enthalpy component
        System.out.println("ENTHALPY CALCULATION:");
        System.out.println("-" .repeat(100));
        System.out.println("H = Σ mhdis[t] · Σ_l eci[l] · u[l]");
        System.out.println();

        double[] mhdis = allData.getDisorderedClusterResult().getMultiplicities();
        System.out.println("HSP multiplicities:");
        for (int i = 0; i < Math.min(3, mhdis.length); i++) {
            System.out.printf("  mhdis[%d] = %.2f%n", i, mhdis[i]);
        }
        System.out.println();

        double H_ortho = 0.0;
        for (int i = 0; i < eciOrtho.length && i < u.length; i++) {
            double contribution = eciOrtho[i] * u[i];
            H_ortho += mhdis.get(0) * contribution;
            System.out.printf("  eci_ortho[%d] * u[%d] = %.2f * %.4f = %.6f%n",
                i, i, eciOrtho[i], u[i], contribution);
        }
        System.out.printf("H (orthogonal): %.6f J/mol%n", H_ortho);
        System.out.println();

        double H_cvcf = 0.0;
        for (int i = 0; i < eciCvcf.length && i < u.length; i++) {
            double contribution = eciCvcf[i] * u[i];
            H_cvcf += mhdis.get(0) * contribution;
            System.out.printf("  eci_cvcf[%d] * u[%d] = %.2f * %.4f = %.6f%n",
                i, i, eciCvcf[i], u[i], contribution);
        }
        System.out.printf("H (CVCF): %.6f J/mol%n", H_cvcf);
        System.out.println();

        // Build full CF vectors
        System.out.println("FULL CF VECTOR CONSTRUCTION:");
        System.out.println("-" .repeat(100));
        int[][] cfBasisIndices = allData.getCMatrixResult().getCfBasisIndices();
        int tcf = allData.getDisorderedCFResult().getNcf();
        int ncf = eciOrtho.length;

        double[] uFullOrtho = ClusterVariableEvaluator.buildFullCFVector(
            u, composition, 2, cfBasisIndices, ncf, tcf);
        System.out.println("Full u (orthogonal basis):");
        for (int i = 0; i < uFullOrtho.length; i++) {
            System.out.printf("  u_full[%d] = %.6f%n", i, uFullOrtho[i]);
        }
        System.out.println();

        ncf = eciCvcf.length;
        CvCfBasis cvcfBasis = BccA2CvCfTransformations.binaryBasis();
        double[] uFullCvcf = ClusterVariableEvaluator.buildFullCFVector(
            u, composition, 2, cfBasisIndices, ncf, cvcfBasis.totalCfs());
        System.out.println("Full u (CVCF basis) - should be SAME as orthogonal:");
        for (int i = 0; i < uFullCvcf.length; i++) {
            System.out.printf("  u_full[%d] = %.6f (diff from ortho: %.2e)%n",
                i, uFullCvcf[i], Math.abs(uFullCvcf[i] - uFullOrtho[i]));
        }
        System.out.println();

        System.out.println("=" .repeat(100));
    }

    private double[] extractEciFromEntry(CECEntry entry, double temperature) {
        double[] eci = new double[entry.ncf];
        for (int i = 0; i < entry.ncf && i < entry.cecTerms.length; i++) {
            CECTerm term = entry.cecTerms[i];
            eci[i] = term.a + term.b * temperature;
        }
        return eci;
    }
}
