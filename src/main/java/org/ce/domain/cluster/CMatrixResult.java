package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;

/**
 * Holds C-matrix data: coefficients, counts, weights, and random CF generation.
 *
 * <p>All information needed to evaluate cluster variables and generate random
 * correlation function (CF) values at the random (disordered) state.</p>
 */
public final class CMatrixResult {

    private final List<List<double[][]>> cmat;
    private final int[][] lcv;
    private final List<List<int[]>> wcv;

    /**
     * Per-CF basis-index decoration patterns.
     *
     * <p>{@code cfBasisIndices[col]} = array of basis indices (1-based) for
     * the CF at column {@code col} in the C-matrix.  For a multi-site CF
     * decorated with site operators Ïƒ^{kâ‚}, Ïƒ^{kâ‚‚}, â€¦, this array is
     * {@code [kâ‚, kâ‚‚, â€¦]}.  Point CFs (single-site) have a single entry.</p>
     *
     * <p>This is the Java equivalent of Mathematica's {@code cfSiteOpList}
     * decoration info, used to compute the random-state CF values
     * ({@code uRandRules}) as products of point CFs.</p>
     *
     * <p>For a K-component system: 
     * <ul>
     *   <li>cfBasisIndices tells us which site operators (1..K-1) each CF uses</li>
     *   <li>Site operators are computed from composition via basis vectors</li>
     *   <li>Random CF = product of site operators for that CF's decoration</li>
     * </ul></p>
     */
    private final int[][] cfBasisIndices;

    /**
     * Names of the CVCF CF corresponding to each column index 0..totalCfs-1.
     *
     * <p>Populated by CvCfBasisTransformer after the transform;
     * null for orthogonal results.</p>
     */
    private final List<String> cmatCfNames;

    @JsonCreator
    public CMatrixResult(
            @JsonProperty("cmat")            List<List<double[][]>> cmat,
            @JsonProperty("lcv")             int[][]                lcv,
            @JsonProperty("wcv")             List<List<int[]>>      wcv,
            @JsonProperty("cfBasisIndices")  int[][]                cfBasisIndices,
            @JsonProperty("cmatCfNames")     List<String>           cmatCfNames) {
        this.cmat = cmat;
        this.lcv = lcv;
        this.wcv = wcv;
        this.cfBasisIndices = cfBasisIndices;
        this.cmatCfNames = cmatCfNames;
    }

    /**
     * Default no-arg constructor for deserialization (e.g., JSON loading).
     */
    public CMatrixResult() {
        this.cmat = null;
        this.lcv = null;
        this.wcv = null;
        this.cfBasisIndices = null;
        this.cmatCfNames = null;
    }

    public List<List<double[][]>> getCmat() {
        return cmat;
    }

    public int[][] getLcv() {
        return lcv;
    }

    public List<List<int[]>> getWcv() {
        return wcv;
    }

    /**
     * Returns the per-CF basis-index decoration patterns.
     * See field documentation for details.
     *
     * @return cfBasisIndices[col] = basis indices for CF at column col
     */
    public int[][] getCfBasisIndices() {
        return cfBasisIndices;
    }

    /**
     * Returns the names of the CVCF CF corresponding to each column index.
     * Null for C-matrices in the orthogonal basis.
     *
     * @return list of CF names, length = basis.totalCfs()
     */
    public List<String> getCmatCfNames() {
        return cmatCfNames;
    }

    // =========================================================================
    // Random-state CF generation (Mathematica: uRandRules)
    // =========================================================================

    /**
     * Evaluates random correlation function (CF) values at the disordered state.
     *
     * <p>For a K-component system at the random state, multi-site CFs factor as
     * products of site-operator expectation values:</p>
     * <pre>
     *   u_random[icf] = Î _j  âŸ¨Ïƒ^{basisIndices[icf][j]}âŸ©
     * </pre>
     * where the site operators are computed from composition using basis vectors.
     *
     * <p>This method:</p>
     * <ol>
     *   <li>Computes basis vectors for the system using RMatrixCalculator</li>
     *   <li>Evaluates Kâˆ’1 site operators (point CFs) from composition</li>
     *   <li>For each CF, multiplies the appropriate site operators based on cfBasisIndices</li>
     *   <li>Returns the ncf non-point CF values at random state</li>
     * </ol>
     *
     * <p><b>Design:</b> This method consolidates random CF generation into the data model,
     * making it self-contained and generalized for any K-component system. All necessary
     * information (cfBasisIndices, composition, K) is provided at call time.</p>
     *
     * @param moleFractions mole fractions of all K components (length K, Î£ = 1)
     * @param numElements   number of chemical components K (â‰¥ 2)
     * @return random-state CF values for non-point CFs (length = number of ncf)
     * @throws IllegalArgumentException if moleFractions length != numElements
     *
     * @see RMatrixCalculator#buildBasis(int)
     * @see ClusterVariableEvaluator#computeRandomCFs(double[], int, int[][], int, int)
     */
    /**
     * Validates that the first cmat block has exactly {@code expectedCols} columns.
     *
     * @throws IllegalStateException if the actual column count does not match
     */
    public void validateCols(int expectedCols, String contextMsg) {
        if (cmat == null || cmat.isEmpty()) return;
        var firstType = cmat.get(0);
        if (firstType == null || firstType.isEmpty()) return;
        double[][] firstBlock = firstType.get(0);
        if (firstBlock == null || firstBlock.length == 0) return;
        int actual = firstBlock[0].length;
        if (actual != expectedCols) {
            throw new IllegalStateException(contextMsg
                    + ": expected " + expectedCols + " columns, got " + actual);
        }
    }

    public double[] evaluateRandomCFs(double[] moleFractions, int numElements) {
        if (moleFractions == null) {
            throw new IllegalArgumentException("moleFractions must not be null");
        }
        if (moleFractions.length != numElements) {
            throw new IllegalArgumentException(
                "moleFractions length (" + moleFractions.length + 
                ") must equal numElements (" + numElements + ")");
        }
        if (cfBasisIndices == null) {
            throw new IllegalStateException("cfBasisIndices is null; cannot compute random CFs");
        }

        // Step 1: Get basis vectors for this K-component system
        double[] basis = RMatrixCalculator.buildBasis(numElements);

        // Step 2: Compute Kâˆ’1 point CFs (site operators) from composition
        // For binary: âŸ¨ÏƒÂ¹âŸ© = 2x_B - 1
        // For ternary: âŸ¨ÏƒÂ¹âŸ© = Î£_i x_i Â· t_iÂ¹,  âŸ¨ÏƒÂ²âŸ© = Î£_i x_i Â· t_iÂ²
        int nxcf = numElements - 1;
        double[] pointCFs = new double[nxcf];
        
        for (int k = 0; k < nxcf; k++) {
            int power = k + 1;  // ÏƒÂ¹ is power 1, ÏƒÂ² is power 2, etc.
            for (int i = 0; i < numElements; i++) {
                pointCFs[k] += moleFractions[i] * Math.pow(basis[i], power);
            }
        }

        // Step 3: For each non-point CF, compute product of site operators
        // cfBasisIndices[icf] tells us which site operators (1..K-1) to multiply
        int ncf = cfBasisIndices.length - nxcf;  // number of non-point CFs
        double[] uRandom = new double[ncf];
        
        for (int icf = 0; icf < ncf; icf++) {
            int[] indices = cfBasisIndices[icf];
            double val = 1.0;
            // Each index is 1-based: index 1 means ÏƒÂ¹ = pointCFs[0]
            for (int b : indices) {
                val *= pointCFs[b - 1];
            }
            uRandom[icf] = val;
        }

        return uRandom;
    }
}

