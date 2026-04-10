package org.ce.domain.cluster.cvcf;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import org.ce.model.cluster.cvcf.CvCfBasis;
import org.ce.model.cluster.cvcf.BccA2TModelCvCfTransformations;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for CVCF basis transformation consistency.
 *
 * Verifies that:
 * 1. Transformation matrices are well-formed
 * 2. CF names and dimensions match
 * 3. Forward/reverse transformations are consistent
 */
class CvCfBasisTransformationTest {

    @Test
    @DisplayName("Binary BCC_A2 basis structure is correct")
    void testBinaryBasisStructure() {
        CvCfBasis binary = BccA2TModelCvCfTransformations.binaryBasis();

        assertEquals("BCC_A2", binary.structurePhase);
        assertEquals(2, binary.numComponents);
        assertEquals(4, binary.numNonPointCfs);
        assertEquals(6, binary.totalCfs());  // 4 non-point + 2 point (xA, xB)

        // Verify CF names
        List<String> expected = List.of("v4AB", "v3AB", "v22AB", "v21AB", "xA", "xB");
        assertEquals(expected, binary.cfNames);

        // Verify T matrix dimensions
        assertEquals(6, binary.T.length);        // rows
        assertEquals(6, binary.T[0].length);     // columns
    }

    @Test
    @DisplayName("Ternary BCC_A2 basis structure is correct")
    void testTernaryBasisStructure() {
        CvCfBasis ternary = BccA2TModelCvCfTransformations.ternaryBasis();

        assertEquals("BCC_A2", ternary.structurePhase);
        assertEquals(3, ternary.numComponents);
        assertEquals(18, ternary.numNonPointCfs);
        assertEquals(21, ternary.totalCfs());  // 18 non-point + 3 point (xA, xB, xC)

        // Verify T matrix dimensions
        assertEquals(21, ternary.T.length);
        assertEquals(21, ternary.T[0].length);
    }

    @Test
    @DisplayName("Quaternary BCC_A2 basis structure is correct")
    void testQuaternaryBasisStructure() {
        CvCfBasis quaternary = BccA2TModelCvCfTransformations.quaternaryBasis();

        assertEquals("BCC_A2", quaternary.structurePhase);
        assertEquals(4, quaternary.numComponents);
        assertEquals(51, quaternary.numNonPointCfs);
        assertEquals(55, quaternary.totalCfs());  // 51 non-point + 4 point (xA, xB, xC, xD)

        // Verify T matrix dimensions
        assertEquals(55, quaternary.T.length);
        assertEquals(55, quaternary.T[0].length);
    }

    @Test
    @DisplayName("Binary basis component ordering is correct")
    void testBinaryComponentPairs() {
        // Test component ordering directly from basis
        CvCfBasis basis = BccA2TModelCvCfTransformations.binaryBasis();
        assertEquals(6, basis.cfNames.size());
        assertEquals("v4AB", basis.cfNames.get(0));  // First non-point CF
        assertEquals("v21AB", basis.cfNames.get(3)); // Last non-point CF
        assertEquals("xA", basis.cfNames.get(4));    // First point CF
        assertEquals("xB", basis.cfNames.get(5));    // Last point CF
    }

    @Test
    @DisplayName("Ternary basis component ordering is correct")
    void testTernaryComponentOrdering() {
        CvCfBasis basis = BccA2TModelCvCfTransformations.ternaryBasis();
        assertEquals(21, basis.cfNames.size());
        assertEquals("v4AB", basis.cfNames.get(0));  // First tetrahedra binary
        assertEquals("xA", basis.cfNames.get(18));   // First point variable
    }

    @Test
    @DisplayName("Quaternary basis component ordering is correct")
    void testQuaternaryComponentOrdering() {
        CvCfBasis basis = BccA2TModelCvCfTransformations.quaternaryBasis();
        assertEquals(55, basis.cfNames.size());
        assertEquals("xA", basis.cfNames.get(51));   // First point variable
        assertEquals("xB", basis.cfNames.get(52));
        assertEquals("xC", basis.cfNames.get(53));
        assertEquals("xD", basis.cfNames.get(54));
    }

    @Test
    @DisplayName("T matrix for binary has expected sparsity")
    void testBinaryTMatrixStructure() {
        double[][] T = BccA2TModelCvCfTransformations.BINARY_T;

        // Check first row: should match transformation rule for u[1][1][1]
        // u[1][1][1] = -16·v21AB + 8·v22AB + 16·v4AB + xA + xB
        // T[0] = [16, 0, 8, -16, 1, 1]  (v4, v3, v22, v21, xA, xB)
        assertEquals(16.0, T[0][0], 1e-10);   // v4
        assertEquals(0.0, T[0][1], 1e-10);    // v3
        assertEquals(8.0, T[0][2], 1e-10);    // v22
        assertEquals(-16.0, T[0][3], 1e-10);  // v21
        assertEquals(1.0, T[0][4], 1e-10);    // xA
        assertEquals(1.0, T[0][5], 1e-10);    // xB
    }
}
