package org.ce.domain.cluster.cvcf;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Maps CVCF CF names to column indices for a given (structure, numComponents) combination.
 *
 * <p>Example usage:</p>
 * <pre>
 *   CvCfCFRegistry registry = CvCfCFRegistry.forBccA2Binary();
 *   int col = registry.indexOf("v21AB");  // returns 0
 *   int col = registry.indexOf("xA");     // returns 4
 * </pre>
 */
public final class CvCfCFRegistry {

    private final String structurePhase;
    private final int numComponents;
    private final List<String> cfNames;
    private final Map<String, Integer> nameToIndex;

    private CvCfCFRegistry(String structurePhase, int numComponents, List<String> cfNames) {
        this.structurePhase = structurePhase;
        this.numComponents = numComponents;
        this.cfNames = cfNames;
        this.nameToIndex = new HashMap<>();
        for (int i = 0; i < cfNames.size(); i++) {
            nameToIndex.put(cfNames.get(i), i);
        }
    }

    /** Returns the column index for the given CF name. */
    public int indexOf(String cfName) {
        Integer idx = nameToIndex.get(cfName);
        if (idx == null) {
            throw new IllegalArgumentException(
                "CF name '" + cfName + "' not found in " + structurePhase
                + " " + numComponents + "-component CVCF basis. "
                + "Available: " + cfNames);
        }
        return idx;
    }

    /** Returns true if the given CF name is in this registry. */
    public boolean contains(String cfName) {
        return nameToIndex.containsKey(cfName);
    }

    /** Returns all CF names in column order. */
    public List<String> getCfNames() {
        return cfNames;
    }

    // =========================================================================
    // Factory methods
    // =========================================================================

    /** Returns a registry for binary BCC_A2. */
    public static CvCfCFRegistry forBccA2Binary() {
        return new CvCfCFRegistry("BCC_A2", 2, BccA2CvCfTransformations.BINARY_CF_NAMES);
    }

    /** Returns a registry for ternary BCC_A2. */
    public static CvCfCFRegistry forBccA2Ternary() {
        return new CvCfCFRegistry("BCC_A2", 3, BccA2CvCfTransformations.TERNARY_CF_NAMES);
    }

    /** Returns a registry for quaternary BCC_A2. */
    public static CvCfCFRegistry forBccA2Quaternary() {
        return new CvCfCFRegistry("BCC_A2", 4, BccA2CvCfTransformations.QUATERNARY_CF_NAMES);
    }
}
