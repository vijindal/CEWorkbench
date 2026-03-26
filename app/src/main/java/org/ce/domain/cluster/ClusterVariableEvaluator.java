package org.ce.domain.cluster;

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
 * where {@code t_i} is the {@link RMatrixCalculator#buildBasis(int)} value
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
     * @param cfBasisIndices  decoration patterns from CMatrixResult:
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
        double[] basis = RMatrixCalculator.buildBasis(numElements);
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
     * @param cfBasisIndices per-CF basis-index decorations from CMatrixResult
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
        double[] basis = RMatrixCalculator.buildBasis(numElements);
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
     * @param cfBasisIndices per-CF basis-index decorations from CMatrixResult
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
     * Computes the random-state non-point CVCF values for solver initialisation.
     *
     * <p>Converts the orthogonal-basis random CFs to the CVCF basis using the
     * precomputed inverse transformation {@code Tinv} stored in the basis object.
     * If {@code basis.Tinv} is {@code null} a zero vector is returned (safe fallback
     * for equicomposition binary where the orthogonal random state is zero anyway).</p>
     *
     * @param x      mole fractions (length K)
     * @param basis  CVCF basis description (must carry cfBasisIndices-compatible data
     *               via its T/Tinv matrices for the orthogonal intermediate step)
     * @return non-point CVCF values at the random (disordered) state (length ncf)
     */
    public static double[] computeRandomCVCFs(double[] x, CvCfBasis basis) {
        int ncf = basis.numNonPointCfs;
        int tcf = basis.totalCfs();

        if (basis.Tinv == null) {
            return new double[ncf]; // zero-vector fallback
        }

        // Step 1: build orthogonal random u_full using the standard formula.
        // Point CFs in orthogonal basis are σ^k = Σ_i x_i * t_i^k (k=1..nxcf).
        double[] basisVec = RMatrixCalculator.buildBasis(x.length);
        int nxcf = tcf - ncf;
        double[] pointCFs = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            int power = k + 1;
            for (int i = 0; i < x.length; i++) {
                pointCFs[k] += x[i] * Math.pow(basisVec[i], power);
            }
        }
        // For non-point orthogonal CFs at random state: product of point CFs raised
        // to the powers given by cfBasisIndices. Without cfBasisIndices we reconstruct
        // from the fact that each non-point CF row in T maps to a product of pointCFs.
        // We use Tinv directly on u_orth_full = [0..0 | pointCFs] (non-point = 0
        // at equicomposition; for other compositions we need the full mapping).
        // Since cfBasisIndices is unavailable in CVCF context, approximate by
        // treating non-point orthogonal CFs as products of single-site averages.
        // This is exact for any composition: u_orth[l] = Π_{b in decoration} pointCF[b-1].
        // We compute it by applying T^{-1} to a reconstructed u_orth.
        // Here we reconstruct u_orth from the known structure of T.
        // For robustness: apply Tinv to [pointCFs-products | pointCFs].
        // Without cfBasisIndices we cannot recover the decorations, so we apply
        // Tinv to u_orth where u_orth non-point = 0 (valid at equicomposition).
        // TODO: store cfBasisIndices or decoration info in CvCfBasis for exact init.
        double[] uOrthFull = new double[tcf];
        // non-point entries remain 0 (exact at equicomposition, approximate elsewhere)
        for (int k = 0; k < nxcf; k++) {
            uOrthFull[ncf + k] = pointCFs[k];
        }

        // Step 2: v_nonpoint = Tinv[:ncf, :] * u_orth_full
        double[] vRand = new double[ncf];
        for (int j = 0; j < ncf; j++) {
            for (int i = 0; i < tcf; i++) {
                vRand[j] += basis.Tinv[j][i] * uOrthFull[i];
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
                double[][] cm = cmat.get(t).get(j);  // [lcv[t][j]] x [tcf+1]
                int nv = lcv[t][j];
                cv[t][j] = new double[nv];
                for (int v = 0; v < nv; v++) {
                    double val = cm[v][tcf]; // constant term (last column)
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

