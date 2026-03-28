package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.ClusterVariableEvaluator;

import java.util.List;
import java.util.logging.Logger;

/**
 * Computes the CVM free-energy functional and its derivatives.
 *
 * <p>The Gibbs energy of mixing is:</p>
 * <pre>
 *   G = H âˆ’ T Â· S
 * </pre>
 *
 * <h2>Enthalpy (linear in non-point CVCF variables)</h2>
 * <pre>
 *   H            = sum_{l=0}^{ncf-1} ECI[l] * v[l]
 *   dH/dv[l]     = ECI[l]
 *   d2H/dv2      = 0
 * </pre>
 *
 * <h2>Entropy (CVM Kikuchiâ€“Barker formula)</h2>
 * <pre>
 *   S = âˆ’R Â· Î£_t kb[t] Â· ms[t] Â· Î£_j mh[t][j] Â· Î£_v wcv[t][j][incv] Â· CV Â· ln(CV)
 *
 *   Scu[l] = âˆ’R Â· Î£_t kb[t] Â· ms[t] Â· Î£_j mh[t][j]
 *            Â· Î£_v wcv[t][j][incv] Â· cmat[t][j][v][l] Â· ln(CV[t][j][v])
 *
 *   Scuu[l1][l2] = âˆ’R Â· Î£_t kb[t] Â· ms[t] Â· Î£_j mh[t][j]
 *                  Â· Î£_v wcv[t][j][incv] Â· cmat[t][j][v][l1] Â· cmat[t][j][v][l2] / CV[t][j][v]
 * </pre>
 *
 * <p>Note: in the Mathematica code, {@code ms[t] = mhdis[t]} (HSP multiplicities
 * are used as the "system parameter" multiplicities in the entropy).</p>
 *
 * <h2>Gibbs energy derivatives</h2>
 * <pre>
 *   Gcu  = Hcu âˆ’ T Â· Scu
 *   Gcuu = âˆ’T Â· Scuu          (since Hcuu = 0)
 * </pre>
 *
 * @see ClusterVariableEvaluator
 */
public final class CVMFreeEnergy {

    private static final Logger LOG = Logger.getLogger(CVMFreeEnergy.class.getName());

    /** Gas constant in J/(molÂ·K). */
    public static final double R_GAS = 8.3144598;

    private CVMFreeEnergy() { /* utility class */ }

    // =========================================================================
    // Result container
    // =========================================================================

    /**
     * Holds the evaluated free-energy functional and its first/second derivatives.
     */
    public static final class EvalResult {
        /** Gibbs energy of mixing G = H - T*S. */
        public final double G;
        /** Enthalpy of mixing. */
        public final double H;
        /** Entropy of mixing. */
        public final double S;
        /** Enthalpy gradient dH/du (length ncf). */
        public final double[] Hcu;
        /** Enthalpy Hessian d²H/du² (always 0, since H is linear in u). */
        public final double[][] Hcuu;
        /** Entropy gradient dS/du (length ncf). */
        public final double[] Scu;
        /** Entropy Hessian d²S/du² (ncf x ncf). */
        public final double[][] Scuu;
        /** Gibbs gradient dG/du (length ncf). */
        public final double[] Gcu;
        /** Gibbs Hessian d²G/du² (ncf x ncf). */
        public final double[][] Gcuu;

        public EvalResult(double G, double H, double S,
                         double[] Hcu, double[][] Hcuu,
                         double[] Scu, double[][] Scuu,
                         double[] Gcu, double[][] Gcuu) {
            this.G = G;
            this.H = H;
            this.S = S;
            this.Hcu = Hcu;
            this.Hcuu = Hcuu;
            this.Scu = Scu;
            this.Scuu = Scuu;
            this.Gcu = Gcu;
            this.Gcuu = Gcuu;
        }
    }

    // =========================================================================
    // Main evaluation
    // =========================================================================

    /**
     * Evaluates the CVM free-energy functional, gradient, and Hessian at the
     * given CVCF values.
     *
     * <p>In the CVCF basis ECIs already include cluster multiplicities, so the
     * enthalpy is simply {@code H = Σ eci[l] · v[l]}. The entropy formula is
     * unchanged — it runs through the transformed C-matrix exactly as before.</p>
     *
     * @param v             non-point CVCF values (length ncf) — the optimisation variables
     * @param moleFractions mole fractions of all K components (length K, Σ = 1)
     * @param temperature   temperature in Kelvin
     * @param eci           effective cluster interactions in CVCF basis (length ncf,
     *                      multiplicities already absorbed)
     * @param mhdis         HSP cluster multiplicities (size tcdis) — used in entropy
     * @param kb            Kikuchi-Baker entropy coefficients (length tcdis)
     * @param mh            normalised multiplicities: mh[t][j]
     * @param lc            ordered clusters per HSP type: lc[t]
     * @param cmat          CVCF C-matrix: cmat.get(t).get(j)[v][k], last col = constant
     * @param lcv           CV counts: lcv[t][j]
     * @param wcv           CV weights: wcv.get(t).get(j)[v]
     * @param tcdis         number of HSP cluster types
     * @param tcf           total number of CFs
     * @param ncf           number of non-point CFs
     * @return evaluation result containing G, H, S, gradient, and Hessian
     */
    public static EvalResult evaluate(
            double[] v,
            double[] moleFractions,
            double temperature,
            double[] eci,
            List<Double> mhdis,
            double[] kb,
            double[][] mh,
            int[] lc,
            List<List<double[][]>> cmat,
            int[][] lcv,
            List<List<int[]>> wcv,
            int tcdis,
            int tcf,
            int ncf) {

        LOG.finer("CVMFreeEnergy.evaluate — ENTER: T=" + temperature + ", ncf=" + ncf + ", tcf=" + tcf);

        // Build full CVCF vector: vFull = [v | moleFractions]
        double[] vFull = ClusterVariableEvaluator.buildFullCVCFVector(v, moleFractions, ncf, tcf);
        double[][][] cv = ClusterVariableEvaluator.evaluate(vFull, cmat, lcv, tcdis, lc);

        // --- Enthalpy ---
        // H = Σ_{l=0}^{ncf-1} eci[l] · v[l]
        // ECIs in CVCF basis already include cluster multiplicities.
        double Hval = 0.0;
        double[] Hcu = new double[ncf];
        for (int l = 0; l < ncf; l++) {
            Hcu[l] = eci[l];
            Hval += eci[l] * v[l];
        }

        // --- Entropy ---
        // S = âˆ’R Â· Î£_t kb[t] Â· ms[t] Â· Î£_j mh[t][j] Â· Î£_v w Â· cv Â· ln(cv)
        // Scu[l] = âˆ’R Â· Î£_t kb[t] Â· ms[t] Â· Î£_j mh[t][j] Â· Î£_v w Â· cmat[l] Â· ln(cv)
        // Scuu[l1][l2] = âˆ’R Â· Î£_t kb[t] Â· ms[t] Â· Î£_j mh[t][j] Â· Î£_v w Â· cmat[l1] Â· cmat[l2] / cv
        double Sval = 0.0;
        double[] Scu = new double[ncf];
        double[][] Scuu = new double[ncf][ncf];

        // Use the physical gas constant R_GAS = 8.3144598 J/(mol·K).
        // ECI must be provided in J/mol and temperature in Kelvin.
        double R = R_GAS;

        for (int t = 0; t < tcdis; t++) {
            double coeff_t = kb[t] * mhdis.get(t); // kb[t] Â· ms[t]
            for (int j = 0; j < lc[t]; j++) {
                double mh_tj = mh[t][j];
                double[][] cm = cmat.get(t).get(j);
                int[] w = wcv.get(t).get(j);
                int nv = lcv[t][j];

                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = cv[t][j][incv];
                    int wv = w[incv];

                    // Smooth entropy extension for CV â‰¤ EPS.
                    // For CV > EPS: use exact cvÂ·ln(cv), ln(cv), 1/cv.
                    // For CV â‰¤ EPS: use a CÂ² quadratic extension that creates
                    // a soft barrier pushing the solver toward positive CVs.
                    // This is critical for Kâ‰¥3 where the all-zero initial
                    // guess produces negative CVs.
                    final double EPS = 1.0e-6;
                    double sContrib; // cvÂ·ln(cv) or smooth extension
                    double logEff;   // effective ln(cv) for gradient
                    double invEff;   // effective 1/cv for Hessian

                    if (cvVal > EPS) {
                        double logCv = Math.log(cvVal);
                        sContrib = cvVal * logCv;
                        logEff = logCv;
                        invEff = 1.0 / cvVal;
                    } else {
                        double logEps = Math.log(EPS);
                        double d = cvVal - EPS;
                        sContrib = EPS * logEps + (1.0 + logEps) * d + 0.5 / EPS * d * d;
                        logEff = logEps + d / EPS;
                        invEff = 1.0 / EPS;
                    }

                    double prefix = coeff_t * mh_tj * wv;

                    // Entropy value
                    Sval -= R * prefix * sContrib;

                    // Gradient: only first ncf columns of cmat contribute
                    // (point-CF columns and constant column are not optimisation variables)
                    for (int l = 0; l < ncf; l++) {
                        double cml = cm[incv][l];
                        if (cml == 0.0) continue;
                        // d(cvÂ·ln(cv))/du_l = cmat[v][l] Â· (1 + ln(cv))
                        // But the "+1" term vanishes due to normalization constraint:
                        // Î£_v wcv Â· cmat[v][l] = 0 for l > 0 (non-trivial CFs)
                        // So: Scu[l] = âˆ’R Â· Î£ prefix Â· cmat[v][l] Â· ln(cv)
                        Scu[l] -= R * prefix * cml * logEff;
                    }

                    // Hessian
                    for (int l1 = 0; l1 < ncf; l1++) {
                        double cml1 = cm[incv][l1];
                        if (cml1 == 0.0) continue;
                        for (int l2 = l1; l2 < ncf; l2++) {
                            double cml2 = cm[incv][l2];
                            if (cml2 == 0.0) continue;
                            double val = -R * prefix * cml1 * cml2 * invEff;
                            Scuu[l1][l2] += val;
                            if (l1 != l2) {
                                Scuu[l2][l1] += val;
                            }
                        }
                    }
                }
            }
        }

        // --- Gibbs energy ---
        double Gval = Hval - temperature * Sval;
        double[] Gcu = new double[ncf];
        double[][] Gcuu = new double[ncf][ncf];

        for (int l = 0; l < ncf; l++) {
            Gcu[l] = Hcu[l] - temperature * Scu[l];
        }

        for (int l1 = 0; l1 < ncf; l1++) {
            for (int l2 = 0; l2 < ncf; l2++) {
                // Hcuu = 0, so Gcuu = âˆ’T Â· Scuu
                Gcuu[l1][l2] = -temperature * Scuu[l1][l2];
            }
        }

        LOG.finer("CVMFreeEnergy.evaluate — EXIT: G=" + String.format("%.8e", Gval)
                + ", H=" + String.format("%.8e", Hval) + ", S=" + String.format("%.8e", Sval));
        return new EvalResult(Gval, Hval, Sval, Hcu, new double[ncf][ncf], Scu, Scuu, Gcu, Gcuu);
    }

    /**
     * Convenience overload for binary systems.
     *
     * @param v           non-point CVCF values (length ncf)
     * @param composition mole fraction of component B (binary shorthand)
     * @param temperature temperature in Kelvin
     * @param eci         effective cluster interactions in CVCF basis (length ncf)
     * @param mhdis       HSP cluster multiplicities
     * @param kb          Kikuchi-Baker entropy coefficients
     * @param mh          normalised multiplicities: mh[t][j]
     * @param lc          ordered clusters per HSP type: lc[t]
     * @param cmat        CVCF C-matrix
     * @param lcv         CV counts
     * @param wcv         CV weights
     * @param tcdis       number of HSP cluster types
     * @param tcf         total number of CFs
     * @param ncf         number of non-point CFs
     * @return evaluation result
     */
    public static EvalResult evaluate(
            double[] v,
            double composition,
            double temperature,
            double[] eci,
            List<Double> mhdis,
            double[] kb,
            double[][] mh,
            int[] lc,
            List<List<double[][]>> cmat,
            int[][] lcv,
            List<List<int[]>> wcv,
            int tcdis,
            int tcf,
            int ncf) {

        return evaluate(v, new double[]{1.0 - composition, composition},
                temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv,
                tcdis, tcf, ncf);
    }
}

