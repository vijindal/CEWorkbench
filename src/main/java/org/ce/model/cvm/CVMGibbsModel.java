package org.ce.model.cvm;

import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.ClusterVariableEvaluator;
import org.ce.model.cluster.cvcf.CvCfBasis;
import org.ce.model.cluster.CFIdentificationResult;
import org.ce.model.cluster.CMatrix;
import org.ce.model.cluster.ClusterIdentificationResult;
import org.ce.model.cluster.ClusterMath;

import java.util.List;

/**
 * Physical model for the Cluster Variation Method (CVM).
 *
 * <p>
 * This class encapsulates the structural geometry, multiplicities, and
 * entropy coefficients. It provides the Gibbs free energy functional and its
 * derivatives (Gradient, Hessian) for a given thermodynamic state.
 * </p>
 */
public class CVMGibbsModel {


    public static final double R_GAS = 8.3144598;
    private static final double ENTROPY_SMOOTH_EPS = 1.0e-6;

    private final int tcdis, tcf, ncf;
    private final List<Double> mhdis;
    private final double[] kb;
    private final double[][] mh;
    private final int[] lc;
    private final List<List<double[][]>> cmat;
    private final int[][] lcv;
    private final List<List<int[]>> wcv;
    private final int[][] orthCfBasisIndices;
    private final CvCfBasis basis;

    /**
     * Calculated free energy and derivatives.
     */
    public static class ModelResult {
        public final double G, H, S;
        public final double[] Gu;
        public final double[][] Guu;
        public final double[] Hu;
        public final double[] Su;
        public final double[][] Suu;

        public ModelResult(double G, double H, double S,
                double[] Gu, double[][] Guu,
                double[] Hu, double[] Su, double[][] Suu) {
            this.G = G;
            this.H = H;
            this.S = S;
            this.Gu = Gu;
            this.Guu = Guu;
            this.Hu = Hu;
            this.Su = Su;
            this.Suu = Suu;
        }
    }

    /**
     * Main constructor that initializes from cluster identification results and basis.
     */
    public CVMGibbsModel(
            ClusterIdentificationResult stage1,
            CFIdentificationResult stage2,
            CMatrix.Result stage3,
            CvCfBasis basis) {
        this.tcdis = stage1.getTcdis();
        this.mhdis = stage1.getDisClusterData().getMultiplicities();
        this.kb = stage1.getKbCoefficients();
        this.mh = stage1.getMh();
        this.lc = stage1.getLc();

        this.cmat = stage3.getCmat();
        this.lcv = stage3.getLcv();
        this.wcv = stage3.getWcv();
        this.orthCfBasisIndices = stage3.getCfBasisIndices();

        this.basis = basis;
        this.ncf = basis.numNonPointCfs;
        this.tcf = basis.totalCfs();
    }

    /**
     * Internal physics evaluator shared by instance and static methods.
     */
    private static ModelResult evaluateInternal(
            double[] u, double[] moleFractions, double temperature, double[] eci,
            int tcdis, int tcf, int ncf, List<Double> mhdis, double[] kb, double[][] mh, int[] lc,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv) {

        double[] uFull = ClusterVariableEvaluator.buildFullCVCFVector(u, moleFractions, ncf, tcf);
        double[][][] cv = ClusterVariableEvaluator.evaluate(uFull, cmat, lcv, tcdis, lc);

        double Hval = 0.0;
        double Sval = 0.0;
        double[] Su = new double[ncf];
        double[][] Suu = new double[ncf][ncf];

        // Enthalpy is linear in non-point CVCF variables
        for (int l = 0; l < ncf; l++) {
            Hval += eci[l] * u[l];
        }

        for (int t = 0; t < tcdis; t++) {
            double coeff_t = kb[t] * mhdis.get(t);
            for (int j = 0; j < lc[t]; j++) {
                double mh_tj = mh[t][j];
                double[][] cm = cmat.get(t).get(j);
                int[] w = wcv.get(t).get(j);
                int nv = lcv[t][j];

                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = cv[t][j][incv];
                    int wv = w[incv];

                    double sContrib, logEff, invEff;
                    if (cvVal > ENTROPY_SMOOTH_EPS) {
                        double logCv = Math.log(cvVal);
                        sContrib = cvVal * logCv;
                        logEff = logCv;
                        invEff = 1.0 / cvVal;
                    } else {
                        double logEps = Math.log(ENTROPY_SMOOTH_EPS);
                        double d = cvVal - ENTROPY_SMOOTH_EPS;
                        sContrib = ENTROPY_SMOOTH_EPS * logEps + (1.0 + logEps) * d + 0.5 / ENTROPY_SMOOTH_EPS * d * d;
                        logEff = logEps + d / ENTROPY_SMOOTH_EPS;
                        invEff = 1.0 / ENTROPY_SMOOTH_EPS;
                    }

                    double prefix = coeff_t * mh_tj * wv;
                    Sval -= R_GAS * prefix * sContrib;

                    for (int l = 0; l < ncf; l++) {
                        double cml = cm[incv][l];
                        if (cml != 0.0)
                            Su[l] -= R_GAS * prefix * cml * logEff;
                    }

                    for (int l1 = 0; l1 < ncf; l1++) {
                        double cml1 = cm[incv][l1];
                        if (cml1 == 0.0)
                            continue;
                        for (int l2 = l1; l2 < ncf; l2++) {
                            double cml2 = cm[incv][l2];
                            if (cml2 == 0.0)
                                continue;
                            double val = -R_GAS * prefix * cml1 * cml2 * invEff;
                            Suu[l1][l2] += val;
                            if (l1 != l2)
                                Suu[l2][l1] += val;
                        }
                    }
                }
            }
        }

        double Gval = Hval - temperature * Sval;
        double[] Gu = new double[ncf];
        double[][] Guu = new double[ncf][ncf];
        for (int l = 0; l < ncf; l++) {
            Gu[l] = eci[l] - temperature * Su[l];
            for (int l2 = 0; l2 < ncf; l2++) {
                Guu[l][l2] = -temperature * Suu[l][l2];
            }
        }

        return new ModelResult(Gval, Hval, Sval, Gu, Guu, eci, Su, Suu);
    }

    /**
     * Evaluates the CVM Gibbs free energy at the given state.
     */
    public ModelResult evaluate(double[] u, double[] moleFractions, double temperature, double[] eci) {
        return evaluateInternal(u, moleFractions, temperature, eci,
                tcdis, tcf, ncf, mhdis, kb, mh, lc, cmat, lcv, wcv);
    }

    /**
     * Static evaluator for one-off calculations (e.g. from tests or CLI).
     */
    public static ModelResult evaluate(
            double[] u, double[] moleFractions, double temperature, double[] eci,
            int tcdis, int tcf, int ncf, List<Double> mhdis, double[] kb, double[][] mh, int[] lc,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv) {
        return evaluateInternal(u, moleFractions, temperature, eci,
                tcdis, tcf, ncf, mhdis, kb, mh, lc, cmat, lcv, wcv);
    }

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
     *       the T matrix.</li>
     *   <li><b>Transform to CVCF basis</b>: apply {@code Tinv} (inverse of the orthogonal-to-CVCF
     *       transformation matrix) to obtain the CVCF coordinates.</li>
     * </ol>
     *
     * @param moleFractions   mole fractions (length K = 2 or more, Σ = 1.0)
     * @return full CVCF vector of length {@link #getTcf()} = numNonPointCfs + K
     * @throws IllegalArgumentException if {@code orthCfBasisIndices} is null or empty
     */
    public double[] computeRandomCFs(double[] moleFractions) {
        int K    = moleFractions.length;
        int nxcf = K - 1;                         // # orthogonal point CFs
        int orthTcf = orthCfBasisIndices.length;   // total orthogonal CFs (non-point + point)
        int orthNcf = orthTcf - nxcf;             // # orthogonal non-point CFs

        // ── Step 1: K-1 orthogonal point CFs from composition ───────────────
        double[] basisVec = ClusterMath.buildBasis(K);
        double[] pointCF  = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            for (int i = 0; i < K; i++) {
                pointCF[k] += moleFractions[i] * Math.pow(basisVec[i], k + 1);
            }
        }

        // ── Step 2: orthogonal non-point random CFs ──────────────────────────
        double[] uNonPoint = new double[orthNcf];
        for (int col = 0; col < orthNcf; col++) {
            double val = 1.0;
            for (int b : orthCfBasisIndices[col]) val *= pointCF[b - 1];
            uNonPoint[col] = val;
        }

        // ── Step 3: full orthogonal vector [u_non-point | u_point | 1.0] ────
        double[] uPoint = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            int col = orthNcf + k;
            int power = orthCfBasisIndices[col][0];
            uPoint[k] = pointCF[power - 1];
        }

        double[][] tInv = basis.Tinv;
        int tRows = tInv[0].length;

        double[] uOrthFull = new double[tRows];
        System.arraycopy(uNonPoint, 0, uOrthFull, 0, orthNcf);
        System.arraycopy(uPoint,    0, uOrthFull, orthNcf, nxcf);
        if (tRows == orthTcf + 1) {
            uOrthFull[tRows - 1] = 1.0;
        } else if (tRows != orthTcf) {
            throw new IllegalStateException(
                "T row count mismatch: T.rows=" + tRows
                + ", orthTcf=" + orthTcf
                + "  (expected T.rows == orthTcf or orthTcf+1)");
        }

        // ── Step 4: transform to CVCF ────────────────────────────────────────
        double[] vFull = new double[tcf];
        for (int j = 0; j < tcf; j++) {
            double sum = 0.0;
            double[] row = tInv[j];
            for (int i = 0; i < row.length; i++) sum += row[i] * uOrthFull[i];
            vFull[j] = sum;
        }

        for (int i = 0; i < K; i++) vFull[ncf + i] = moleFractions[i];

        return vFull;
    }

    /**
     * Finds the maximum step size $\alpha \in (0, 1]$ such that $u + \alpha p$
     * maintains all cluster variables within $[0, 1]$.
     */
    public double calculateStepLimit(double[] u, double[] p, double[] moleFractions) {
        double fmin = 1.0;
        double[] uTrial = new double[ncf];
        for (int i = 0; i < ncf; i++)
            uTrial[i] = u[i] + p[i];

        double[][][] cv_old = updateCVInternal(u, moleFractions);
        double[][][] cv_new = updateCVInternal(uTrial, moleFractions);

        for (int i = 0; i < tcdis - 1; i++) {
            for (int j = 0; j < lc[i]; j++) {
                for (int v = 0; v < lcv[i][j]; v++) {
                    double vO = cv_old[i][j][v], vN = cv_new[i][j][v];
                    if (vN <= 0)
                        fmin = Math.min(fmin, Math.abs(vO / (vN - vO)));
                    if (vN >= 1)
                        fmin = Math.min(fmin, Math.abs((1.0 - vO) / (vN - vO)));
                }
            }
        }
        return (fmin >= 1.0) ? 1.0 : (0.1 * fmin);
    }

    private double[][][] updateCVInternal(double[] u, double[] moleFractions) {
        double[] uFull = ClusterVariableEvaluator.buildFullCVCFVector(u, moleFractions, ncf, tcf);
        return ClusterVariableEvaluator.evaluate(uFull, cmat, lcv, tcdis, lc);
    }

    // Accessors
    public CvCfBasis getBasis() {
        return basis;
    }

    public int getNcf() {
        return ncf;
    }

    public int getTcf() {
        return tcf;
    }
}
