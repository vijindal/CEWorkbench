package org.ce.domain.cluster.cvcf;

import java.util.List;

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
     * <p>May be {@code null} if not provided — in that case CVCF random
     * initialization will fall back to the zero vector.</p>
     */
    public final double[][] Tinv;

    public CvCfBasis(String structurePhase,
                     int numComponents,
                     List<String> cfNames,
                     int numNonPointCfs,
                     double[][] T) {
        this(structurePhase, numComponents, cfNames, numNonPointCfs, T, null);
    }

    public CvCfBasis(String structurePhase,
                     int numComponents,
                     List<String> cfNames,
                     int numNonPointCfs,
                     double[][] T,
                     double[][] Tinv) {
        this.structurePhase = structurePhase;
        this.numComponents = numComponents;
        this.cfNames = cfNames;
        this.numNonPointCfs = numNonPointCfs;
        this.T = T;
        this.Tinv = Tinv;
    }

    /** Total number of CFs (non-point + point). */
    public int totalCfs() {
        return cfNames.size();
    }
}
