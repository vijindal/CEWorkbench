package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import org.ce.domain.cluster.cvcf.CvCfBasis;

import java.util.Arrays;
import java.util.List;

/**
 * Evaluates cluster variables (CVs) from correlation functions (CFs)
 * using the C-matrix.
 *
 * <p>For each HSP cluster type {@code t} and ordered-phase group {@code j},
 * the cluster variables are computed as:</p>
 * <pre>
 *   CV[t][j][v] = Î£_{k=0}^{tcf-1} cmat[t][j][v][k] Â· u_full[k] + cmat[t][j][v][tcf]
 * </pre>
 * where {@code u_full} contains all CF values (ncf non-point CFs followed by
 * nxcf point CFs, with point CFs derived from composition).
 *
 * <p>Supports arbitrary K-component systems.  For a K-component system,
 * the Kâˆ’1 point CFs are computed from the mole fractions using:</p>
 * <pre>
 *   âŸ¨s^kâŸ© = Î£_{i=0}^{K-1} x_i Â· t_i^k     for k = 1, â€¦, Kâˆ’1
 * </pre>
 * where {@code t_i} is the {@link ClusterMath#buildBasis(int)} value
 * for component i.</p>
 *
 * <p>Corresponds to Mathematica {@code cvRules} evaluation.</p>
 */
public final class ClusterVariableEvaluator {

    private ClusterVariableEvaluator() { /* utility class */ }

    // =========================================================================
    // Random-state initial guess (Mathematica: uRandRules)
    // =========================================================================

    /**
     * Computes the random-state CF values for every CF column.
     *
     * <p>At the random (completely disordered) state, multi-site CFs factor as
     * products of individual point CFs:</p>
     * <pre>
     *   CF_random[col] = Î _{b âˆˆ basisIndices[col]} pointCF[b âˆ’ 1]
     * </pre>
     * where {@code pointCF[k] = âŸ¨Ïƒ^{k+1}âŸ© = Î£_i x_i Â· t_i^{k+1}}.
     *
     * <p>This is the Java equivalent of Mathematica's {@code uRandRules}.
     * For K=2 binary at equimolar, Ïƒ=0, so all random non-point CFs = 0
     * (matching the old u=0 initial guess).  For Kâ‰¥3, random CFs involving
     * Ïƒâ‚‚ are nonzero, which is critical for starting with all positive CVs.</p>
     *
     * @param moleFractions   mole fractions (length K, Î£ = 1)
     * @param numElements     number of components K
     * @param cfBasisIndices  decoration patterns from CMatrix.Result:
     *                        cfBasisIndices[col] = basis indices for CF col
     * @param ncf             number of non-point (independent) CFs
     * @param tcf             total number of CFs
     * @return non-point CFs at the random state (length ncf)
     */
    public static double[] computeRandomCFs(
            double[] moleFractions,
            int numElements,
            int[][] cfBasisIndices,
            int ncf,
            int tcf) {

        // Compute the Kâˆ’1 point CFs from composition
        double[] basis = ClusterMath.buildBasis(numElements);
        int nxcf = tcf - ncf;
        double[] pointCFs = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            int power = k + 1;
            for (int i = 0; i < numElements; i++) {
                pointCFs[k] += moleFractions[i] * Math.pow(basis[i], power);
            }
        }

        // For each non-point CF, random value = product of point CFs
        // for each site's basis-index decoration
        double[] uRandom = new double[ncf];
        for (int l = 0; l < ncf; l++) {
            int[] indices = cfBasisIndices[l];
            double val = 1.0;
            for (int b : indices) {
                val *= pointCFs[b - 1]; // basisIndex is 1-based
            }
            uRandom[l] = val;
        }

        return uRandom;
    }

    /**
     * Builds the full CF vector for the C-matrix multiplication.
     *
     * <p>{@code u_full} has length {@code tcf}:</p>
     * <ul>
     *   <li>Indices {@code 0..ncf-1}: the non-point CFs being optimised</li>
     *   <li>Indices {@code ncf..tcf-1}: point CFs determined by composition</li>
     * </ul>
     *
     * <p>The Kâˆ’1 point CFs are placed in the order determined by the CF
     * identification pipeline (not necessarily ascending power order).
     * Each point CF column has a single basis-index decoration (from
     * {@code cfBasisIndices}) that specifies the power: Ïƒ^{basisIndex}.</p>
     *
     * @param u              non-point CF values (length ncf)
     * @param moleFractions  mole fractions of all K components (length K, Î£ = 1)
     * @param numElements    number of chemical components K (â‰¥ 2)
     * @param cfBasisIndices per-CF basis-index decorations from CMatrix.Result
     * @param ncf            number of non-point CFs
     * @param tcf            total number of CFs
     * @return full CF vector (length tcf)
     */
    public static double[] buildFullCFVector(double[] u, double[] moleFractions,
                                             int numElements, int[][] cfBasisIndices,
                                             int ncf, int tcf) {
        double[] uFull = new double[tcf];
        System.arraycopy(u, 0, uFull, 0, ncf);

        // Pre-compute all Kâˆ’1 point CF values: pointCF[k] = âŸ¨Ïƒ^{k+1}âŸ©
        double[] basis = ClusterMath.buildBasis(numElements);
        int nxcf = tcf - ncf;
        double[] pointCFValues = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            int power = k + 1;
            for (int i = 0; i < numElements; i++) {
                pointCFValues[k] += moleFractions[i] * Math.pow(basis[i], power);
            }
        }

        // Place each point CF in correct column using cfBasisIndices
        for (int k = 0; k < nxcf; k++) {
            int col = ncf + k;
            int basisIndex = cfBasisIndices[col][0]; // single decoration for point CF
            // basisIndex is 1-based â†’ power = basisIndex â†’ pointCFValues index = basisIndex-1
            uFull[col] = pointCFValues[basisIndex - 1];
        }
        return uFull;
    }

    /**
     * Convenience overload for binary systems.
     *
     * @param u           non-point CF values (length ncf)
     * @param composition mole fraction of component B (0 â‰¤ xB â‰¤ 1)
     * @param cfBasisIndices per-CF basis-index decorations from CMatrix.Result
     * @param ncf         number of non-point CFs
     * @param tcf         total number of CFs
     * @return full CF vector (length tcf)
     */
    public static double[] buildFullCFVector(double[] u, double composition,
                                             int[][] cfBasisIndices, int ncf, int tcf) {
        return buildFullCFVector(u, new double[]{1.0 - composition, composition},
                2, cfBasisIndices, ncf, tcf);
    }

    // =========================================================================
    // CVCF-basis methods (no cfBasisIndices needed)
    // =========================================================================

    /**
     * Builds the full CVCF vector for C-matrix multiplication.
     *
     * <p>In the CVCF basis the point variables are simply the mole fractions
     * themselves — one entry per component — so no {@code cfBasisIndices} is needed:</p>
     * <pre>
     *   vFull[0..ncf-1]     = v  (non-point optimization variables)
     *   vFull[ncf + i]      = x[i]   for i = 0..K-1
     * </pre>
     *
     * @param v    non-point CVCF values (length ncf)
     * @param x    mole fractions of all K components (length K = tcf - ncf)
     * @param ncf  number of non-point CFs
     * @param tcf  total number of CFs (ncf + K)
     * @return full CVCF vector (length tcf)
     */
    public static double[] buildFullCVCFVector(double[] v, double[] x, int ncf, int tcf) {
        double[] vFull = Arrays.copyOf(v, tcf);
        for (int i = 0; i < x.length; i++) {
            vFull[ncf + i] = x[i];
        }
        return vFull;
    }

    /**
     * Computes exact random-state non-point CVCF values for solver initialization.
     *
     * <p>Algorithm:
     * <ol>
     *   <li>Compute exact orthogonal random CFs from composition and orthogonal
     *       {@code cfBasisIndices} (Mathematica uRandRules equivalent).</li>
     *   <li>Transform orthogonal full CF vector to CVCF full vector using
     *       {@code v = T^{-1} u}.</li>
     *   <li>Return the first {@code ncf} non-point CVCF values.</li>
     * </ol></p>
     *
     * <p>Handles both transformation-table conventions:
     * {@code T.rows == oldTcf} and {@code T.rows == oldTcf + 1} (extra empty-cluster row).</p>
     *
     * @param x mole fractions (length K)
     * @param basis CVCF basis description
     * @param orthCfBasisIndices orthogonal-basis decoration metadata from Stage-3
     * @return non-point CVCF random-state values (length ncf)
     */
    /**
     * Generates the full disordered-state vector (v_nonpoint + x) directly.
     * This avoids potential basis-index mismatches by using the target probabilities directly.
     *
     * @param x mole fractions (length K)
     * @param basis CVCF basis
     * @param orthCfBasisIndices orthogonal site-decoration metadata
     * @return full CF vector (length tcf = ncf + K)
     */
    public static double[] generateFullRandomCVCFs(double[] x, CvCfBasis basis, int[][] orthCfBasisIndices) {
        int ncf = basis.numNonPointCfs;
        int tcf = basis.totalCfs();

        // Specific algebraic fallback for ternary BCC_A2 (VJ rules 2026-04-04)
        if (x.length == 3 && ncf == 18 && tcf == 21) {
            System.out.println("DEBUG: Using Ternary BCC_A2 algebraic fallback for random state (VJ rules)");
            return generateTernaryBccA2RandomCFs(x);
        }

        // 1. Get the non-point v from the transformation
        double[] vNonPoint = computeRandomCVCFs(x, basis, orthCfBasisIndices);

        // 2. Build the full vector (ncf + K)
        double[] vFull = new double[tcf];
        System.arraycopy(vNonPoint, 0, vFull, 0, ncf);

        // 3. Last K entries are exactly the mole fractions
        for (int i = 0; i < x.length; i++) {
            if (ncf + i < tcf) {
                vFull[ncf + i] = x[i];
            }
        }
        return vFull;
    }

    /**
     * Algebraic implementation of ternary BCC_A2 disordered state rules provided by user.
     * Orders: tetra (bin, tern), triangle (bin, tern), pair (II, I), points.
     */
    private static double[] generateTernaryBccA2RandomCFs(double[] x) {
        double xA = x[0], xB = x[1], xC = x[2];
        double[] v = new double[21];

        // Tetrahedra (bin, tern)
        v[0] = xA * xA * xB * xB; // v4AB
        v[1] = xA * xA * xC * xC; // v4AC
        v[2] = xB * xB * xC * xC; // v4BC
        v[3] = xA * xA * xB * xC; // v4ABC1
        v[4] = xA * xB * xB * xC; // v4ABC2
        v[5] = xA * xB * xC * xC; // v4ABC3

        // Triangles (bin, tern)
        v[6] = -xA * xA * xB + xA * xB * xB; // v3AB
        v[7] = -xA * xA * xC + xA * xC * xC; // v3AC
        v[8] = -xB * xB * xC + xB * xC * xC; // v3BC
        v[9]  = xA * xB * xC; // v3ABC1
        v[10] = xA * xB * xC; // v3ABC2
        v[11] = xA * xB * xC; // v3ABC3

        // Pairs (II, I)
        v[12] = xA * xB; // v22AB
        v[13] = xA * xC; // v22AC
        v[14] = xB * xC; // v22BC
        v[15] = xA * xB; // v21AB
        v[16] = xA * xC; // v21AC
        v[17] = xB * xC; // v21BC

        // Point variables
        v[18] = xA;
        v[19] = xB;
        v[20] = xC;

        return v;
    }

    public static double[] computeRandomCVCFs(double[] x, CvCfBasis basis, int[][] orthCfBasisIndices) {
        int ncf = basis.numNonPointCfs;
        int tcf = basis.totalCfs();
        if (orthCfBasisIndices == null || orthCfBasisIndices.length == 0) {
            // Legacy cluster_data files may not carry orthogonal decoration metadata.
            // Fallback: use composition-only orthogonal vector (non-point = 0, point from composition),
            // then transform via T^{-1}. This keeps the solver runnable, though less accurate
            // than the exact uRandRules path.
            return computeRandomCVCFsFallback(x, basis, ncf, tcf);
        }

        int oldTcf = orthCfBasisIndices.length;
        int oldNxcf = x.length - 1;
        int oldNcf = oldTcf - oldNxcf;
        if (oldNcf < 0) {
            throw new IllegalArgumentException("Invalid orthogonal CF dimensions: oldTcf=" + oldTcf
                    + ", numComponents=" + x.length);
        }

        // Step 1: exact orthogonal random full CF vector (non-point + point).
        double[] uOldNonPoint = computeRandomCFs(x, x.length, orthCfBasisIndices, oldNcf, oldTcf);
        double[] uOldFull = buildFullCFVector(uOldNonPoint, x, x.length, orthCfBasisIndices, oldNcf, oldTcf);

        System.out.println("DEBUG: Orthogonal CF Decorations and Values:");
        for (int i = 0; i < Math.min(uOldFull.length, 25); i++) {
            System.out.println("  u[" + i + "]: val=" + String.format("%.8e", uOldFull[i]) 
                + " indices=" + Arrays.toString(orthCfBasisIndices[i]));
        }

        // Align with T row convention: either exact rows or extra empty-cluster row.
        int tRows = basis.T.length;
        double[] uForTransform;
        if (tRows == oldTcf) {
            uForTransform = uOldFull;
        } else if (tRows == oldTcf + 1) {
            uForTransform = Arrays.copyOf(uOldFull, tRows);
            uForTransform[tRows - 1] = 1.0; // empty-cluster / composition-sum row
        } else {
            throw new IllegalArgumentException("Basis/orthogonal size mismatch for random init: T rows="
                    + tRows + ", orthogonal tcf=" + oldTcf);
        }

        System.out.println("DEBUG: Orthogonal Random state (uOldFull): " + Arrays.toString(uForTransform));

        // Step 2: v_nonpoint = Tinv[:ncf, :] * u_for_transform
        double[][] tInv = basis.Tinv != null ? basis.Tinv : LinearAlgebra.invert(basis.T);
        double[] vRand = new double[ncf];
        for (int j = 0; j < ncf; j++) {
            for (int i = 0; i < tInv[0].length; i++) {
                vRand[j] += tInv[j][i] * uForTransform[i];
            }
        }
        return vRand;
    }

    private static double[] computeRandomCVCFsFallback(double[] x, CvCfBasis basis, int ncf, int tcf) {
        int tRows = basis.T.length;
        int nxcf = x.length - 1;
        double[] basisVec = ClusterMath.buildBasis(x.length);
        double[] pointCFs = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            int power = k + 1;
            for (int i = 0; i < x.length; i++) {
                pointCFs[k] += x[i] * Math.pow(basisVec[i], power);
            }
        }

        double[] uApprox;
        if (tRows == tcf || tRows == tcf + 1) {
            uApprox = new double[tRows];
            // Fix for ternary: properly map K-1 point CFs to uApprox
            // In disordered state, u_linear (basisVec^1) and u_quadratic (basisVec^2) ...
            int start = tcf - (x.length - 1);
            if (tRows == tcf + 1) start = tcf - (x.length - 1); // skip empty sum row? No.

            // Placement depends on T row convention.
            // For BCC_A2_T_tern: Row 18 is quadratic, 19 is linear, 20 is constant.
            // Fallback placement: assume linear at T-3, quadratic at T-2, const at T-1?
            // Actually, we use the pointCFs calculated from basis powers.
            for (int k = 0; k < nxcf && (tcf - nxcf + k) < uApprox.length; k++) {
                uApprox[tcf - nxcf + k] = pointCFs[k];
            }
            if (tRows == tcf + 1) {
                uApprox[tRows - 1] = 1.0;
            }
        } else {
            throw new IllegalArgumentException("Cannot build fallback random init: T rows="
                    + tRows + ", expected " + tcf + " or " + (tcf + 1));
        }

        System.out.println("DEBUG: Orthogonal Random state (fallback uApprox): " + Arrays.toString(uApprox));

        double[][] tInv = basis.Tinv != null ? basis.Tinv : LinearAlgebra.invert(basis.T);
        double[] vRand = new double[ncf];
        for (int j = 0; j < ncf; j++) {
            for (int i = 0; i < tInv[0].length; i++) {
                vRand[j] += tInv[j][i] * uApprox[i];
            }
        }
        return vRand;
    }

    /**
     * Evaluates all cluster variables from the full CF vector.
     *
     * @param uFull   full CF vector (length tcf)
     * @param cmat    C-matrix data: {@code cmat.get(t).get(j)[v][k]},
     *                last column (index tcf) is the constant term
     * @param lcv     CV counts: {@code lcv[t][j]}
     * @param tcdis   number of HSP cluster types
     * @param lc      ordered clusters per HSP type: {@code lc[t]}
     * @return cluster variables: {@code cv[t][j][v]}
     */
    public static double[][][] evaluate(
            double[] uFull,
            List<List<double[][]>> cmat,
            int[][] lcv,
            int tcdis,
            int[] lc) {

        int tcf = uFull.length;
        double[][][] cv = new double[tcdis][][];

        for (int t = 0; t < tcdis; t++) {
            cv[t] = new double[lc[t]][];
            for (int j = 0; j < lc[t]; j++) {
                double[][] cm = cmat.get(t).get(j);  // [lcv[t][j]] x [tcf] or [tcf+1]
                int nv = lcv[t][j];
                cv[t][j] = new double[nv];
                if (nv == 0) {
                    continue;
                }
                int cols = cm[0].length;
                boolean hasConstant = cols == (tcf + 1);
                if (!(hasConstant || cols == tcf)) {
                    throw new IllegalArgumentException(
                            "Invalid cmat width at (t=" + t + ", j=" + j + "): got " + cols
                                    + ", expected " + tcf + " (no constant) or " + (tcf + 1)
                                    + " (with constant)");
                }
                for (int v = 0; v < nv; v++) {
                    double val = hasConstant ? cm[v][tcf] : 0.0;
                    for (int k = 0; k < tcf; k++) {
                        val += cm[v][k] * uFull[k];
                    }
                    cv[t][j][v] = val;
                }
            }
        }
        return cv;
    }
}

