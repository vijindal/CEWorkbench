package org.ce.domain.cluster.cvcf;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.IntFunction;

/**
 * Holds all data for one (structure, numComponents) CVCF basis combination.
 *
 * <p>The CVCF basis uses cluster variables (CV probabilities) as correlation
 * functions (CFs), making cluster expansion coefficients (CECs) invariant
 * across system orders. Binary CECs can be used directly in ternary/quaternary
 * systems without re-fitting.</p>
 *
 * <p>The transformation from old orthogonal CFs to new CVCF CFs is simply:</p>
 * <pre>
 *   u_old[i] = Σ_j T[i][j] · v_new[j]
 * </pre>
 * where {@code v_new} includes both non-point CFs (v21AB, v22AB, ...) and
 * <em>all</em> point variables (xA, xB, ... for each component). Listing every
 * point variable explicitly means expressions like {@code xA + xB} or
 * {@code -xA + xB} that appear in the old CF rules are captured by T alone —
 * no separate constant vector is needed.
 *
 * <p>The corresponding cmatrix transformation is:</p>
 * <pre>
 *   cmat_new = cmat_old · T
 * </pre>
 * applied per CV row per cluster group (t, j).
 */
public final class CvCfBasis {

    /** Structure identifier, e.g. "BCC_A2". */
    public final String structurePhase;

    /** Number of chemical components (2 = binary, 3 = ternary, 4 = quaternary). */
    public final int numComponents;

    /**
     * Ordered list of all CVCF CF names.
     *
     * <p>The first {@code numNonPointCfs} entries are non-point CFs (the
     * optimization variables). The remaining entries are point variables
     * (xA, xB, ...) — one per component, all listed explicitly so that
     * the T matrix captures all linear combinations without a constant term.</p>
     *
     * <p>Example for binary BCC_A2: [v21AB, v22AB, v3AB, v4AB, xA, xB]</p>
     */
    public final List<String> cfNames;

    /**
     * Number of non-point CFs (the optimization variables).
     * Point variables count = {@code cfNames.size() - numNonPointCfs}.
     */
    public final int numNonPointCfs;

    /**
     * Ordered list of ECI (Effective Cluster Interaction) names for CVCF.
     * Corresponds to the first {@code numNonPointCfs} entries of {@code cfNames}.
     */
    public final List<String> eciNames;

    /**
     * Transformation matrix: {@code T[i_old][j_new]}.
     *
     * <p>Rows: old orthogonal CFs (in Mathematica order u[type][group][1]).</p>
     * <p>Columns: new CVCF CFs in the order of {@code cfNames}.</p>
     */
    public final double[][] T;

    /**
     * Inverse of T: {@code Tinv[j_new][i_old]}.
     *
     * <p>Used to transform an orthogonal-basis vector back to CVCF basis:
     * {@code v[j] = Σ_i Tinv[j][i] * u_orth[i]}.</p>
     *
     * <p>Must not be null; computed as the matrix inverse of T.</p>
     */
    public final double[][] Tinv;

    /** O(1) CF-name to index map. Built once in the constructor. */
    private final Map<String, Integer> cfNameIndex;

    public CvCfBasis(String structurePhase,
                     int numComponents,
                     List<String> cfNames,
                     List<String> eciNames,
                     int numNonPointCfs,
                     double[][] T) {
        this(structurePhase, numComponents, cfNames, eciNames, numNonPointCfs, T, null);
    }

    public CvCfBasis(String structurePhase,
                     int numComponents,
                     List<String> cfNames,
                     List<String> eciNames,
                     int numNonPointCfs,
                     double[][] T,
                     double[][] Tinv) {
        this.structurePhase = structurePhase;
        this.numComponents  = numComponents;
        this.cfNames        = cfNames;
        this.eciNames       = eciNames;
        this.numNonPointCfs = numNonPointCfs;
        this.T              = T;
        this.Tinv           = Tinv;
        Map<String, Integer> idx = new HashMap<>(cfNames.size() * 2);
        for (int i = 0; i < cfNames.size(); i++) idx.put(cfNames.get(i), i);
        this.cfNameIndex = idx;
    }

    /**
     * Returns the index of the given CF name in {@code cfNames}, or -1 if not found.
     * O(1) via internal HashMap.
     */
    public int indexOfCf(String name) {
        return cfNameIndex.getOrDefault(name, -1);
    }

    /** Total number of CFs (non-point + point). */
    public int totalCfs() {
        return cfNames.size();
    }

    // =========================================================================
    // Registry
    // =========================================================================

    /**
     * Single registry mapping (structurePhase, numComponents) → CvCfBasis.
     *
     * <p>Replaces the standalone CvCfBasisRegistry. Usage:</p>
     * <pre>
     *   CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", 3);
     * </pre>
     */
    public static final class Registry {

        /** Singleton. */
        public static final Registry INSTANCE = new Registry();

        /** Key: structurePhase, Value: factory producing the basis. */
        private final Map<String, IntFunction<CvCfBasis>> factories = new LinkedHashMap<>();

        private Registry() {
            register("BCC_A2", BccA2TModelCvCfTransformations::basisForNumComponents);
        }

        /**
         * Registers a basis factory for a structure phase.
         */
        public void register(String structurePhase, IntFunction<CvCfBasis> factory) {
            factories.put(structurePhase, factory);
        }

        /**
         * Returns true if (structurePhase, numComponents) is supported.
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
         * Returns the basis for (structurePhase, numComponents), or throws.
         */
        public CvCfBasis get(String structurePhase, int numComponents) {
            return find(structurePhase, numComponents)
                    .orElseThrow(() -> unsupportedError(structurePhase, numComponents));
        }

        /**
         * Returns a human-readable summary of all registered structure phases.
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

        private IllegalArgumentException unsupportedError(String structurePhase, int numComponents) {
            return new IllegalArgumentException(
                    "No CVCF basis registered for (structurePhase=" + structurePhase
                    + ", numComponents=" + numComponents + ").\n"
                    + supportedSummary());
        }
    }
}
