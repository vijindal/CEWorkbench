package org.ce.domain.cluster.cvcf;

import org.junit.jupiter.api.Test;

import java.util.Optional;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for {@link CvCfBasis.Registry}.
 */
class CvCfBasisRegistryTest {

    @Test
    void testIsSupportedBccA2Binary() {
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", 2));
    }

    @Test
    void testIsSupportedBccA2Ternary() {
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", 3));
    }

    @Test
    void testIsSupportedBccA2Quaternary() {
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", 4));
    }

    @Test
    void testIsSupportedUnregisteredStructure() {
        assertFalse(CvCfBasis.Registry.INSTANCE.isSupported("FCC_L12", 2));
    }

    @Test
    void testIsSupportedUnsupportedNumComponents() {
        assertFalse(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", 5));
    }

    @Test
    void testGetBccA2Binary() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", 2);
        assertNotNull(basis);
        assertEquals("BCC_A2", basis.structurePhase);
        assertEquals(2, basis.numComponents);
        assertEquals(6, basis.totalCfs());
    }

    @Test
    void testGetBccA2Ternary() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", 3);
        assertNotNull(basis);
        assertEquals("BCC_A2", basis.structurePhase);
        assertEquals(3, basis.numComponents);
        assertEquals(21, basis.totalCfs());
    }

    @Test
    void testGetBccA2Quaternary() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", 4);
        assertNotNull(basis);
        assertEquals("BCC_A2", basis.structurePhase);
        assertEquals(4, basis.numComponents);
        assertEquals(55, basis.totalCfs());
    }

    @Test
    void testGetThrowsOnUnregisteredStructure() {
        IllegalArgumentException ex = assertThrows(
                IllegalArgumentException.class,
                () -> CvCfBasis.Registry.INSTANCE.get("FCC_L12", 2)
        );
        assertTrue(ex.getMessage().contains("FCC_L12"));
        assertTrue(ex.getMessage().contains("Supported"));
    }

    @Test
    void testGetThrowsOnUnsupportedNumComponents() {
        IllegalArgumentException ex = assertThrows(
                IllegalArgumentException.class,
                () -> CvCfBasis.Registry.INSTANCE.get("BCC_A2", 5)
        );
        assertTrue(ex.getMessage().contains("BCC_A2"));
        assertTrue(ex.getMessage().contains("5"));
    }

    @Test
    void testFindReturnsOptionalForValidKey() {
        Optional<CvCfBasis> opt = CvCfBasis.Registry.INSTANCE.find("BCC_A2", 2);
        assertTrue(opt.isPresent());
        assertEquals("BCC_A2", opt.get().structurePhase);
        assertEquals(2, opt.get().numComponents);
    }

    @Test
    void testFindReturnsEmptyForUnregisteredStructure() {
        Optional<CvCfBasis> opt = CvCfBasis.Registry.INSTANCE.find("FCC_L12", 2);
        assertTrue(opt.isEmpty());
    }

    @Test
    void testFindReturnsEmptyForUnsupportedNumComponents() {
        Optional<CvCfBasis> opt = CvCfBasis.Registry.INSTANCE.find("BCC_A2", 5);
        assertTrue(opt.isEmpty());
    }

    @Test
    void testSupportedSummaryContainsAllValidCombinations() {
        String summary = CvCfBasis.Registry.INSTANCE.supportedSummary();
        assertNotNull(summary);

        // Should list BCC_A2 with K=2, 3, 4
        assertTrue(summary.contains("BCC_A2"));
        assertTrue(summary.contains("K=2"));
        assertTrue(summary.contains("K=3"));
        assertTrue(summary.contains("K=4"));
    }

    @Test
    void testSupportedSummaryDoesNotListUnsupportedKeys() {
        String summary = CvCfBasis.Registry.INSTANCE.supportedSummary();

        // Should not list unsupported K values
        assertFalse(summary.contains("K=1"));
        assertFalse(summary.contains("K=5"));
    }

    @Test
    void testBasisTransformationMatricesNonNull() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", 2);
        assertNotNull(basis.T);
        assertTrue(basis.T.length > 0);
        assertTrue(basis.T[0].length > 0);
    }

    @Test
    void testBasisCfNamesNonEmpty() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", 3);
        assertNotNull(basis.cfNames);
        assertEquals(basis.totalCfs(), basis.cfNames.size());
        assertTrue(basis.numNonPointCfs > 0);
        assertTrue(basis.numNonPointCfs <= basis.totalCfs());
    }
}
