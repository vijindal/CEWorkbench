package org.ce.domain.cluster.cvcf;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.IntFunction;
import org.ce.domain.cluster.LinearAlgebra;
import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.CMatrix;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.ClusterMath;

/**
 * Holds all data for one (structure, model, numComponents) CVCF basis combination.
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

    // =========================================================================
    // Fields
    // =========================================================================

    /** Structure identifier, e.g. "BCC_A2". */
    public final String structurePhase;

    /**
     * CVM model identifier, e.g. "T", "TO".
     * Normalized to upper-case on construction.
     */
    public final String model;

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
     * <p>Example for binary BCC_A2: [v4AB, v3AB, v22AB, v21AB, xA, xB]</p>
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
     */
    public final double[][] Tinv;

    /** O(1) CF-name to index map. Built once in the constructor. */
    private final Map<String, Integer> cfNameIndex;

    // =========================================================================
    // Constructors
    // =========================================================================

    /**
     * Legacy constructor (no model). {@code model} defaults to empty string.
     * Tinv is computed on demand from T if not supplied.
     */
    public CvCfBasis(String structurePhase,
                     int numComponents,
                     List<String> cfNames,
                     List<String> eciNames,
                     int numNonPointCfs,
                     double[][] T) {
        this(structurePhase, "", numComponents, cfNames, eciNames, numNonPointCfs, T, null);
    }

    /**
     * Legacy constructor with Tinv but no model. {@code model} defaults to empty string.
     */
    public CvCfBasis(String structurePhase,
                     int numComponents,
                     List<String> cfNames,
                     List<String> eciNames,
                     int numNonPointCfs,
                     double[][] T,
                     double[][] Tinv) {
        this(structurePhase, "", numComponents, cfNames, eciNames, numNonPointCfs, T, Tinv);
    }

    /**
     * Full constructor.
     *
     * @param structurePhase structure identifier, e.g. "BCC_A2"
     * @param model          CVM model identifier, e.g. "T", "TO" (normalized to upper-case)
     * @param numComponents  number of chemical components
     * @param cfNames        ordered CF names (non-point first, then point variables xA, xB, ...)
     * @param eciNames       ECI names corresponding to each non-point CF
     * @param numNonPointCfs number of non-point (optimizer) CFs
     * @param T              orthogonal→CVCF transformation matrix
     * @param Tinv           inverse of T, or null (computed on demand via LinearAlgebra.invert)
     */
    public CvCfBasis(String structurePhase,
                     String model,
                     int numComponents,
                     List<String> cfNames,
                     List<String> eciNames,
                     int numNonPointCfs,
                     double[][] T,
                     double[][] Tinv) {
        this.structurePhase  = structurePhase;
        this.model           = (model == null ? "" : model.toUpperCase());
        this.numComponents   = numComponents;
        this.cfNames         = cfNames;
        this.eciNames        = eciNames;
        this.numNonPointCfs  = numNonPointCfs;
        this.T               = T;
        this.Tinv            = Tinv;
        Map<String, Integer> idx = new HashMap<>(cfNames.size() * 2);
        for (int i = 0; i < cfNames.size(); i++) idx.put(cfNames.get(i), i);
        this.cfNameIndex = idx;
    }

    /**
     * Factory method to create a basis dynamically using CvCfMatrixGenerator.
     */
    public static CvCfBasis dynamic(
            String structurePhase,
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            CMatrix.Result orthMatrix,
            String model,
            java.util.function.Consumer<String> sink) {
        return CvCfMatrixGenerator.generate(structurePhase, clusterResult, cfResult, orthMatrix, model, sink);
    }

    /** Overload without sink for backward compatibility. */
    public static CvCfBasis dynamic(
            String structurePhase,
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            CMatrix.Result orthMatrix,
            String model) {
        return dynamic(structurePhase, clusterResult, cfResult, orthMatrix, model, null);
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
    // Random-state (disordered state) computation
    // =========================================================================

    /**
     * Computes the full disordered-state (random) CVCF vector for N-R solver initialization.
     *
     * <h3>Algorithm (derived from Stage 2 CF definitions and Stage 3 data)</h3>
     * <ol>
     *   <li><b>Point CFs (orthogonal)</b>: compute the K-1 orthogonal point CFs from
     *       composition using the Inden (1992) basis:
     *       {@code ⟨σ^k⟩ = Σ_i x_i · t_i^k}  for k=1…K-1,
     *       where {@code t_i} comes from {@link ClusterMath#buildBasis(int)}.</li>
     *   <li><b>Non-point CFs (orthogonal)</b>: for each CF column, its random value is
     *       the product of the point CFs at the basis indices recorded by Stage 2:
     *       {@code u_rand[col] = Π_{b ∈ cfBasisIndices[col]} ⟨σ^b⟩}
     *       (indices are 1-based powers of σ).</li>
     *   <li><b>Full orthogonal vector</b>: assemble [u_non-point | u_point], then
     *       append the empty-cluster value {@code 1.0} so the vector length matches
     *       the row count of the T matrix (which has one row per orthogonal CF
     *       including the empty cluster).</li>
     *   <li><b>CVCF transform</b>: {@code v_full = T_inv · u_orth_full},
     *       where the first {@code numNonPointCfs} entries are the optimizer variables
     *       and the last K entries are the mole fractions (point variables in CVCF basis).</li>
     * </ol>
     *
     * @param moleFractions      mole fractions x_i (length K, Σ = 1)
     * @param orthCfBasisIndices Stage-3 decoration metadata: {@code orthCfBasisIndices[col]}
     *                           contains the 1-based σ-powers used in the orthogonal CF
     *                           at column {@code col}. Must not be null.
     * @return full CVCF vector of length {@link #totalCfs()} = numNonPointCfs + K
     * @throws IllegalArgumentException if {@code orthCfBasisIndices} is null or empty
     */
    public double[] computeRandomState(double[] moleFractions, int[][] orthCfBasisIndices) {

        if (orthCfBasisIndices == null || orthCfBasisIndices.length == 0) {
            throw new IllegalArgumentException(
                "orthCfBasisIndices must not be null or empty. "
                + "Ensure Stage 3 (CMatrix.buildOrthogonal) has been run and its result saved.");
        }

        int K    = moleFractions.length;
        int nxcf = K - 1;                         // # orthogonal point CFs
        int orthTcf = orthCfBasisIndices.length;   // total orthogonal CFs (non-point + point)
        int orthNcf = orthTcf - nxcf;             // # orthogonal non-point CFs

        // ── Step 1: K-1 orthogonal point CFs from composition ───────────────
        // pointCF[k] = ⟨σ^(k+1)⟩ = Σ_i x_i · basis_i^(k+1)   k = 0 … nxcf-1
        double[] basisVec = ClusterMath.buildBasis(K);
        double[] pointCF  = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            for (int i = 0; i < K; i++) {
                pointCF[k] += moleFractions[i] * Math.pow(basisVec[i], k + 1);
            }
        }

        // ── Step 2: orthogonal non-point random CFs ──────────────────────────
        // u_rand[col] = Π_{b ∈ cfBasisIndices[col]} pointCF[b-1]   (b is 1-based power)
        double[] uNonPoint = new double[orthNcf];
        for (int col = 0; col < orthNcf; col++) {
            double val = 1.0;
            for (int b : orthCfBasisIndices[col]) val *= pointCF[b - 1];
            uNonPoint[col] = val;
        }

        // ── Step 3: full orthogonal vector [u_non-point | u_point | 1.0] ────
        // Point CFs are placed at their correct columns using cfBasisIndices.
        // The last entry is the empty-cluster = 1.0, matching T's row count.
        double[] uPoint = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            int col = orthNcf + k;
            int power = orthCfBasisIndices[col][0]; // single σ-power decoration
            uPoint[k] = pointCF[power - 1];
        }

        double[][] tInv = resolvedTinv();
        int tRows = tInv[0].length;                // columns of Tinv = rows of T

        // Build u_orth_full of length tRows
        double[] uOrthFull = new double[tRows];
        System.arraycopy(uNonPoint, 0, uOrthFull, 0, orthNcf);
        System.arraycopy(uPoint,    0, uOrthFull, orthNcf, nxcf);
        if (tRows == orthTcf + 1) {
            uOrthFull[tRows - 1] = 1.0;            // empty-cluster row = 1
        } else if (tRows != orthTcf) {
            throw new IllegalStateException(
                "T row count mismatch for (" + structurePhase + ", " + model
                + ", K=" + K + "): T.rows=" + tRows
                + ", orthTcf=" + orthTcf
                + "  (expected T.rows == orthTcf or orthTcf+1)");
        }

        // ── Step 4: transform to CVCF ────────────────────────────────────────
        // v_full[j] = Σ_i Tinv[j][i] · uOrthFull[i]   for j = 0 … totalCfs()-1
        // The last K entries (point variables in CVCF basis) equal the mole fractions.
        int ncf  = numNonPointCfs;
        int tcf  = totalCfs();
        double[] vFull = new double[tcf];
        for (int j = 0; j < tcf; j++) {
            double sum = 0.0;
            double[] row = tInv[j];
            for (int i = 0; i < row.length; i++) sum += row[i] * uOrthFull[i];
            vFull[j] = sum;
        }

        // Override the point-variable entries with exact mole fractions
        // (the T_inv rows for point variables should already recover x_i exactly,
        // but we set them explicitly to avoid any floating-point drift)
        for (int i = 0; i < K; i++) vFull[ncf + i] = moleFractions[i];

        return vFull;
    }

    /**
     * Returns Tinv: uses the stored value if available, otherwise inverts T on demand.
     */
    private double[][] resolvedTinv() {
        if (Tinv != null) return Tinv;
        if (T    == null) throw new IllegalStateException(
            "CvCfBasis (" + structurePhase + ", " + model + ", K=" + numComponents
            + ") has neither Tinv nor T.");
        return LinearAlgebra.invert(T);
    }

    // =========================================================================
    // Registry
    // =========================================================================

    /**
     * Single registry mapping (structurePhase, model, numComponents) → CvCfBasis.
     *
     * <p>Usage:</p>
     * <pre>
     *   CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get("BCC_A2", "T", 3);
     * </pre>
     */
    public static final class Registry {

        /** Singleton. */
        public static final Registry INSTANCE = new Registry();

        /** Key: structurePhase, Value: Map of (model -> factory). */
        private final Map<String, Map<String, IntFunction<CvCfBasis>>> factories = new LinkedHashMap<>();

        private Registry() {
            register("BCC_A2", "T", BccA2TModelCvCfTransformations::basisForNumComponents);
        }

        /**
         * Registers a basis factory for a structure phase and model.
         */
        public void register(String structurePhase, String model, IntFunction<CvCfBasis> factory) {
            factories.computeIfAbsent(structurePhase, k -> new LinkedHashMap<>())
                     .put(model.toUpperCase(), factory);
        }

        /**
         * Returns true if (structurePhase, model, numComponents) is supported.
         */
        public boolean isSupported(String structurePhase, String model, int numComponents) {
            String baseModel = normalizeModel(model);
            Map<String, IntFunction<CvCfBasis>> modelMap = factories.get(structurePhase);
            if (modelMap == null) return false;
            
            IntFunction<CvCfBasis> factory = modelMap.get(baseModel);
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
        public Optional<CvCfBasis> find(String structurePhase, String model, int numComponents) {
            String baseModel = normalizeModel(model);
            Map<String, IntFunction<CvCfBasis>> modelMap = factories.get(structurePhase);
            if (modelMap == null) return Optional.empty();

            IntFunction<CvCfBasis> factory = modelMap.get(baseModel);
            if (factory == null) return Optional.empty();

            try {
                CvCfBasis basis = factory.apply(numComponents);
                return Optional.ofNullable(basis);
            } catch (Exception e) {
                return Optional.empty();
            }
        }

        /**
         * Returns the basis for (structurePhase, model, numComponents), or throws.
         */
        public CvCfBasis get(String structurePhase, String model, int numComponents) {
            return find(structurePhase, model, numComponents)
                    .orElseThrow(() -> unsupportedError(structurePhase, model, numComponents));
        }

        private String normalizeModel(String model) {
            if (model == null) return "";
            String m = model.toUpperCase();
            if (m.endsWith("_CVCF")) {
                return m.substring(0, m.length() - 5);
            }
            if (m.equals("CVCF")) {
                return "T"; // Default to T model if only CVCF specified
            }
            return m;
        }

        /**
         * Returns a human-readable summary of all registered structure phases and models.
         */
        public String supportedSummary() {
            StringBuilder sb = new StringBuilder("Supported (structurePhase, model, numComponents) keys:\n");
            for (String structure : factories.keySet()) {
                Map<String, IntFunction<CvCfBasis>> modelMap = factories.get(structure);
                for (String model : modelMap.keySet()) {
                    for (int k = 2; k <= 4; k++) {
                        if (isSupported(structure, model, k)) {
                            sb.append("  ").append(structure)
                              .append(", model=").append(model)
                              .append(", K=").append(k).append("\n");
                        }
                    }
                }
            }
            return sb.toString().stripTrailing();
        }

        private IllegalArgumentException unsupportedError(String structurePhase, String model, int numComponents) {
            return new IllegalArgumentException(
                    "No CVCF basis registered for (structurePhase=" + structurePhase
                    + ", model=" + model
                    + ", numComponents=" + numComponents + ").\n"
                    + supportedSummary());
        }
    }
}
