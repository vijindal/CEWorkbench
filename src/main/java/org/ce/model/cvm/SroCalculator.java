package org.ce.model.cvm;

import java.util.ArrayList;
import java.util.List;

/**
 * Cowley-Warren short-range order (SRO) parameters from CVM cluster variables.
 *
 * <p>Implements Jindal &amp; Lele, <i>Calphad</i> 89 (2025) 102825, Eqs. 40-41:</p>
 *
 * <pre>
 *   pair:        α_ij^PR    = 1 − p_ij^PR    / (x_P x_R)             (Eq. 40)
 *   multisite:   α_1234^MPRT = 1 − ρ_1234^MPRT / (x_M x_P x_R x_T)   (Eq. 41)
 * </pre>
 *
 * <p>As that paper notes, the CVM already produces the cluster probabilities during
 * Gibbs-energy minimization, so no extra machinery is needed — the cluster variables
 * {@code cv[t][j][v]} returned by
 * {@link CVMGibbsModel#clusterVariablesAt(double[], double[])} <em>are</em> the
 * probabilities p and ρ.</p>
 *
 * <h2>Sign convention</h2>
 * <ul>
 *   <li>{@code α < 0} — unlike pairs enriched relative to random: <b>ordering</b></li>
 *   <li>{@code α = 0} — random (ideal) solution</li>
 *   <li>{@code α > 0} — unlike pairs depleted: <b>clustering / phase separation</b></li>
 * </ul>
 *
 * <h2>Cluster-variable layout (verified empirically, BCC_A2 tetrahedron)</h2>
 * <p>Cluster types are ordered {@code t=0} tetrahedron, {@code t=1} triangle,
 * {@code t=2} 2nd-NN pair, {@code t=3} 1st-NN pair, {@code t=4} point. Note the
 * naming trap: the CVCF term {@code v21} (1st-NN) corresponds to {@code t=3}, while
 * {@code v22} (2nd-NN) corresponds to {@code t=2}.</p>
 *
 * <p>Pair CVs are laid out in upper-triangular species order — for K=3 that is
 * {@code [AA, AB, AC, BB, BC, CC]} — confirmed by evaluating the random state at a
 * skewed composition, where each entry reproduces {@code x_P·x_R} exactly.</p>
 */
public final class SroCalculator {

    /** Cluster-type index of the 1st-nearest-neighbour pair (CVCF name {@code v21…}). */
    public static final int T_PAIR_1NN = 3;
    /** Cluster-type index of the 2nd-nearest-neighbour pair (CVCF name {@code v22…}). */
    public static final int T_PAIR_2NN = 2;

    private SroCalculator() {}

    /** One Cowley-Warren pair parameter. */
    public static final class PairSro {
        /** Index of the first species in the system's canonical element order. */
        public final int i;
        /** Index of the second species. */
        public final int j;
        /** Cluster probability p_ij^PR from the CVM. */
        public final double probability;
        /** Random-state reference x_P·x_R. */
        public final double reference;
        /** α = 1 − p/(x_P x_R). */
        public final double alpha;

        PairSro(int i, int j, double probability, double reference, double alpha) {
            this.i = i;
            this.j = j;
            this.probability = probability;
            this.reference = reference;
            this.alpha = alpha;
        }
    }

    /**
     * Computes Cowley-Warren pair SRO parameters (Eq. 40) for one pair shell.
     *
     * <p>Returns entries for all distinct species pairs {@code i ≤ j} in the
     * upper-triangular order the cluster variables use. Like-pairs ({@code i == j})
     * are included; they carry the complementary information to the unlike pairs.</p>
     *
     * @param cv            cluster variables from
     *                      {@link CVMGibbsModel#clusterVariablesAt(double[], double[])}
     * @param moleFractions composition, length K, in canonical element order
     * @param clusterType   {@link #T_PAIR_1NN} or {@link #T_PAIR_2NN}
     * @throws IllegalArgumentException if the CV block does not have the K(K+1)/2
     *         entries a pair cluster must have — which would mean the layout
     *         assumption above no longer holds and the mapping cannot be trusted
     */
    public static List<PairSro> pairSro(double[][][] cv, double[] moleFractions, int clusterType) {
        int K = moleFractions.length;
        int expected = K * (K + 1) / 2;

        if (cv == null || clusterType < 0 || clusterType >= cv.length || cv[clusterType].length == 0)
            throw new IllegalArgumentException(
                    "No cluster variables for cluster type t=" + clusterType);

        double[] p = cv[clusterType][0];
        if (p.length != expected)
            throw new IllegalArgumentException(
                    "Pair cluster t=" + clusterType + " has " + p.length + " cluster variables, "
                    + "expected " + expected + " for K=" + K + " (upper-triangular species pairs). "
                    + "The CV layout assumed by SroCalculator does not hold for this system.");

        List<PairSro> out = new ArrayList<>(expected);
        int v = 0;
        for (int i = 0; i < K; i++) {
            for (int j = i; j < K; j++, v++) {
                double ref = moleFractions[i] * moleFractions[j];
                double alpha = (ref > 0.0) ? 1.0 - p[v] / ref : Double.NaN;
                out.add(new PairSro(i, j, p[v], ref, alpha));
            }
        }
        return out;
    }

    /**
     * Convenience: the unlike-pair ({@code i != j}) parameters only, which are the
     * ones normally reported and plotted.
     */
    public static List<PairSro> unlikePairSro(double[][][] cv, double[] moleFractions, int clusterType) {
        List<PairSro> all = pairSro(cv, moleFractions, clusterType);
        List<PairSro> out = new ArrayList<>();
        for (PairSro s : all) if (s.i != s.j) out.add(s);
        return out;
    }
}
