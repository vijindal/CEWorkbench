package org.ce.domain.cluster;

/**
 * Computes Kikuchi-Baker entropy coefficients for the CVM free energy functional.
 *
 * <h2>Physical meaning</h2>
 * <p>In the Cluster Variation Method, the configurational entropy is approximated
 * as a weighted sum over cluster probability distributions:</p>
 * <pre>
 *   S_CVM = -k_B Î£_t  kb[t] * m[t] * Î£_Ïƒ  Ï_t(Ïƒ) ln Ï_t(Ïƒ)
 * </pre>
 * <p>The Kikuchi-Baker coefficients {@code kb[t]} ensure that all sub-cluster
 * entropies that are already counted in the larger clusters are subtracted out,
 * avoiding double-counting.  For the maximal cluster {@code kb[0] = 1.0} by
 * definition.  For smaller clusters the coefficient is determined by the
 * inclusion-exclusion recurrence below.</p>
 *
 * <h2>Algorithm â€” the recurrence</h2>
 * <p>Clusters are ordered from largest ({@code t=0}) to smallest
 * ({@code t=tcdis-1}).  The KB coefficient for cluster {@code j} is:</p>
 * <pre>
 *   kb[j] = 1  -  Î£_{i < j}  nij[i][j] * kb[i]
 *
 *         = (m[j] - Î£_{i<j} m[i] * nij[i][j] * kb[i]) / m[j]
 * </pre>
 * <p>where {@code nij[i][j]} is the number of times cluster {@code j} appears
 * as a sub-cluster inside cluster {@code i} (from the Nij table), and
 * {@code m[i]} are the (possibly symbolic) multiplicities.
 * After computing the symbolic expression, multiplicities are substituted with
 * their numerical values, which is why two lists are supplied.</p>
 *
 * <h2>Mathematica equivalent</h2>
 * <pre>
 * kbdis = generateKikuchiBakerCoefficients[msdis, nijTable]
 * kbdis = kbdis /. Table[msdis[[i]] â†’ mhdis[[i]], {i, 1, Length[mhdis]-1}]
 *
 * generateKikuchiBakerCoefficients[mList_, nijTable_] := Module[...
 *   For[j=1, j<=Length[mList], j++,
 *     tempSum = Î£_{i<j}  (mList[[i]] * nijTable[[i]][[j]] * kbList[[i]]);
 *     kbList[[j]] = (mList[[j]] - tempSum) / mList[[j]];
 *   ]
 * ]
 * </pre>
 *
 * <p>In the Mathematica code, {@code msdis} are symbolic multiplicity variables
 * that are later substituted with numerical values {@code mhdis}.  Since in
 * Java we work directly with numerical values, a single {@code multiplicities}
 * array suffices.</p>
 *
 * <h2>Verification for A2-T (BCC tetrahedron)</h2>
 * <p>For the A2 binary system with a tetrahedron CVM approximation, the
 * expected clusters (in descending size order) and their KB coefficients are:
 * <pre>
 *   t=0: tetrahedron  (4 sites)   kb = +1
 *   t=1: triangle     (3 sites)   kb = -3
 *   t=2: 1NN pair     (2 sites)   kb = +3
 *   t=3: 2NN pair     (2 sites)   kb = +1
 *   t=4: point        (1 site)    kb = -1  (before empty cluster subtraction)
 * </pre>
 * These are the well-known Kikuchi-Baker coefficients for the BCC tetrahedron
 * approximation.</p>
 *
 * @see NijTableCalculator
 * @see ClusterIdentifier
 */
public class KikuchiBakerCalculator {

    private KikuchiBakerCalculator() {}

    /**
     * Computes the Kikuchi-Baker entropy coefficients for all HSP cluster types.
     *
     * <p>The input list {@code multiplicities} contains the normalised orbit-size
     * values {@code mhdis[t]} for each cluster type {@code t}.  Only the relative
     * ratios matter for the recurrence, so the normalization to point-cluster
     * multiplicity does not affect correctness.</p>
     *
     * @param multiplicities normalised multiplicity of each HSP cluster type,
     *                       ordered from largest to smallest (index 0 = maximal
     *                       cluster); must not be {@code null} and must have the
     *                       same length as one dimension of {@code nijTable}
     * @param nijTable       containment count matrix from {@link NijTableCalculator};
     *                       {@code nijTable[i][j]} counts how many times cluster
     *                       {@code j} appears as a sub-cluster of cluster {@code i}
     * @return array of KB coefficients of length {@code multiplicities.size()};
     *         {@code result[0] = 1.0} always (maximal cluster)
     * @throws IllegalArgumentException if sizes are inconsistent
     */
    public static double[] compute(double[] multiplicities, int[][] nijTable) {

        int tcdis = multiplicities.length;

        if (nijTable.length != tcdis || nijTable[0].length != tcdis) {
            throw new IllegalArgumentException(
                "Nij table dimensions " + nijTable.length + "Ã—" + nijTable[0].length +
                " do not match multiplicity count " + tcdis);
        }

        double[] kb = new double[tcdis];

        // Mathematica: For[j=1, j<=Length[mList], j++, ...]
        for (int j = 0; j < tcdis; j++) {

            double sumTerm = 0.0;

            // Î£_{i < j}  m[i] * nij[i][j] * kb[i]
            for (int i = 0; i < j; i++) {
                sumTerm += multiplicities[i] * nijTable[i][j] * kb[i];
            }

            // kb[j] = (m[j] - Î£) / m[j]
            kb[j] = (multiplicities[j] - sumTerm) / multiplicities[j];
        }

        return kb;
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a table of cluster multiplicities and their KB coefficients.
     *
     * @param multiplicities the input multiplicity array
     * @param kb             the computed KB coefficients
     */
    public static void printDebug(double[] multiplicities, double[] kb) {
        System.out.println("[KikuchiBakerCalculator]");
        System.out.printf("  %-6s %-12s %-14s%n", "Type", "mhdis", "KB coeff");
        for (int t = 0; t < kb.length; t++) {
            System.out.printf("  t=%-4d %-12.4f %-14.8f%n",
                    t, multiplicities[t], kb[t]);
        }

        // Sanity check: Î£ kb[t] * m[t] should equal 0 for sub-maximal clusters
        // (standard KB sum rule)
        double kbSum = 0.0;
        for (int t = 0; t < kb.length; t++) {
            kbSum += kb[t] * multiplicities[t];
        }
        System.out.printf("  Î£ kb[t]*m[t] = %.6f  (should be 0 for closed cluster set)%n",
                kbSum);
    }
}


