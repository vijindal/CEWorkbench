package org.ce.domain.cluster.cvcf;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Optional;
import java.util.function.IntFunction;

/**
 * Single registry mapping (structurePhase, numComponents) → CvCfBasis.
 *
 * <p>Replaces the duplicate inline maps previously scattered across
 * {@code CVMEngine}, {@code Main}, and {@code DataPreparationPanel}.</p>
 *
 * <p>Usage:</p>
 * <pre>
 *   CvCfBasis basis = CvCfBasisRegistry.INSTANCE.get("BCC_A2", 3);
 *   boolean ok     = CvCfBasisRegistry.INSTANCE.isSupported("BCC_A2", 2);
 * </pre>
 */
public final class CvCfBasisRegistry {

    /** Singleton. */
    public static final CvCfBasisRegistry INSTANCE = new CvCfBasisRegistry();

    /** Key: "STRUCTURE|numComponents", Value: factory producing the basis. */
    private final Map<String, IntFunction<CvCfBasis>> factories = new LinkedHashMap<>();

    private CvCfBasisRegistry() {
        register("BCC_A2", BccA2TModelCvCfTransformations::basisForNumComponents);
    }

    /**
     * Registers a basis factory for a structure phase.
     * The factory receives {@code numComponents} and returns the matching {@link CvCfBasis}.
     * Supported numComponents values are determined by the factory itself.
     *
     * @param structurePhase structure identifier, e.g. "BCC_A2"
     * @param factory        function: numComponents → CvCfBasis
     */
    public void register(String structurePhase, IntFunction<CvCfBasis> factory) {
        factories.put(structurePhase, factory);
    }

    /**
     * Returns true if (structurePhase, numComponents) is supported.
     *
     * <p>A combination is supported if a factory is registered for the structure
     * phase AND the factory does not throw for the given numComponents.</p>
     */
    public boolean isSupported(String structurePhase, int numComponents) {
        IntFunction<CvCfBasis> factory = factories.get(structurePhase);
        if (factory == null) return false;
        try {
            return factory.apply(numComponents) != null;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Returns an Optional with the basis, or empty if unsupported.
     */
    public Optional<CvCfBasis> find(String structurePhase, int numComponents) {
        IntFunction<CvCfBasis> factory = factories.get(structurePhase);
        if (factory == null) return Optional.empty();
        try {
            CvCfBasis basis = factory.apply(numComponents);
            return Optional.ofNullable(basis);
        } catch (Exception e) {
            return Optional.empty();
        }
    }

    /**
     * Returns the basis for (structurePhase, numComponents), or throws with a
     * clear message listing all supported keys.
     *
     * @throws IllegalArgumentException if the combination is not supported
     */
    public CvCfBasis get(String structurePhase, int numComponents) {
        return find(structurePhase, numComponents)
                .orElseThrow(() -> unsupportedError(structurePhase, numComponents));
    }

    /**
     * Returns a human-readable summary of all registered structure phases.
     * Supported numComponents values are discovered by probing 2, 3, 4.
     */
    public String supportedSummary() {
        StringBuilder sb = new StringBuilder("Supported (structurePhase, numComponents) keys:\n");
        for (String structure : factories.keySet()) {
            for (int k = 2; k <= 4; k++) {
                if (isSupported(structure, k)) {
                    sb.append("  ").append(structure).append(", K=").append(k).append("\n");
                }
            }
        }
        return sb.toString().stripTrailing();
    }

    // -------------------------------------------------------------------------

    private IllegalArgumentException unsupportedError(String structurePhase, int numComponents) {
        return new IllegalArgumentException(
                "No CVCF basis registered for (structurePhase=" + structurePhase
                + ", numComponents=" + numComponents + ").\n"
                + supportedSummary());
    }
}
