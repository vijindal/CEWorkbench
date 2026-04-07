package org.ce.domain.cluster.cvcf;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Optional;

/**
 * Unit tests for {@link CvCfBasis.Registry}.
 */
class CvCfBasisRegistryTest {

    @Test
    void testIsSupported() {
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "T", 2));
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "T", 3));
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "T", 4));
        
        // Test normalization
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "t", 2));
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "T_CVCF", 2));
        assertTrue(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "CVCF", 2));
    }

    @Test
    void testIsSupportedUnsupportedStructure() {
        assertFalse(CvCfBasis.Registry.INSTANCE.isSupported("UNKNOWN", "T", 2));
    }

    @Test
    void testIsSupportedUnsupportedModel() {
        assertFalse(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "UNKNOWN", 2));
    }

    @Test
    void testIsSupportedUnsupportedNumComponents() {
        assertFalse(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "T", 1));
        assertFalse(CvCfBasis.Registry.INSTANCE.isSupported("BCC_A2", "T", 5));
    }

    @Test
    void testGet() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", "T", 2);
        assertNotNull(basis);
        assertEquals("BCC_A2", basis.structurePhase);
        assertEquals(2, basis.numComponents);
    }

    @Test
    void testGetTernary() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", "T", 3);
        assertNotNull(basis);
        assertEquals(3, basis.numComponents);
    }

    @Test
    void testGetQuaternary() {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", "T", 4);
        assertNotNull(basis);
        assertEquals(4, basis.numComponents);
    }

    @Test
    void testGetThrowsOnUnsupportedStructure() {
        assertThrows(IllegalArgumentException.class, () -> 
            CvCfBasis.Registry.INSTANCE.get("UNKNOWN", "T", 2));
    }

    @Test
    void testGetThrowsOnUnsupportedModel() {
        assertThrows(IllegalArgumentException.class, () -> 
            CvCfBasis.Registry.INSTANCE.get("BCC_A2", "UNKNOWN", 2));
    }

    @Test
    void testGetThrowsOnUnsupportedNumComponents() {
        assertThrows(IllegalArgumentException.class, () -> 
            CvCfBasis.Registry.INSTANCE.get("BCC_A2", "T", 5));
    }

    @Test
    void testFind() {
        Optional<CvCfBasis> opt = CvCfBasis.Registry.INSTANCE.find("BCC_A2", "T", 2);
        assertTrue(opt.isPresent());
        assertEquals(2, opt.get().numComponents);
    }

    @Test
    void testFindReturnsEmptyForUnsupportedStructure() {
        Optional<CvCfBasis> opt = CvCfBasis.Registry.INSTANCE.find("UNKNOWN", "T", 2);
        assertFalse(opt.isPresent());
    }

    @Test
    void testFindReturnsEmptyForUnsupportedNumComponents() {
        Optional<CvCfBasis> opt = CvCfBasis.Registry.INSTANCE.find("BCC_A2", "T", 5);
        assertFalse(opt.isPresent());
    }

    @Test
    void testSupportedSummary() {
        String summary = CvCfBasis.Registry.INSTANCE.supportedSummary();
        assertNotNull(summary);
        assertTrue(summary.contains("BCC_A2"));
        assertTrue(summary.contains("model=T"));
        assertTrue(summary.contains("K=2"));
        assertTrue(summary.contains("K=3"));
        assertTrue(summary.contains("K=4"));
    }
}
