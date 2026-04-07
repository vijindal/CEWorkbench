package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.ClusterVariableEvaluator;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.CMatrix;
import org.ce.domain.cluster.ClusterIdentificationResult;

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

    /**
     * Input contract for the CVM model.
     */
    public static final class CVMInput {
        private final ClusterIdentificationResult stage1;
        private final CFIdentificationResult stage2;
        private final CMatrix.Result stage3;
        private final String systemId;
        private final String systemName;
        private final int numComponents;
        private final CvCfBasis basis;

        public CVMInput(
                ClusterIdentificationResult stage1,
                CFIdentificationResult stage2,
                CMatrix.Result stage3,
                String systemId,
                String systemName,
                int numComponents,
                CvCfBasis basis) {
            this.stage1 = stage1;
            this.stage2 = stage2;
            this.stage3 = stage3;
            this.systemId = systemId;
            this.systemName = systemName;
            this.numComponents = numComponents;
            this.basis = basis;
        }

        public ClusterIdentificationResult getStage1()  { return stage1; }
        public CFIdentificationResult      getStage2()  { return stage2; }
        public CMatrix.Result               getStage3()  { return stage3; }
        public String                      getSystemId()      { return systemId; }
        public String                      getSystemName()    { return systemName; }
        public int                         getNumComponents() { return numComponents; }
        public CvCfBasis                   getBasis()         { return basis; }
    }

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

    public CVMGibbsModel(CVMInput input) {
        ClusterIdentificationResult stage1 = input.getStage1();
        this.tcdis = stage1.getTcdis();
        this.mhdis = stage1.getDisClusterData().getMultiplicities();
        this.kb = stage1.getKbCoefficients();
        this.mh = stage1.getMh();
        this.lc = stage1.getLc();

        CMatrix.Result stage3 = input.getStage3();
        this.cmat = stage3.getCmat();
        this.lcv = stage3.getLcv();
        this.wcv = stage3.getWcv();
        this.orthCfBasisIndices = stage3.getCfBasisIndices();

        this.basis = input.getBasis();
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
     * Gets the initial random state for the minimization.
     */
    public double[] getInitialGuess(double[] moleFractions) {
        return basis.computeRandomState(moleFractions, orthCfBasisIndices);
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
        return (fmin >= 1.0) ? 1.0 : (0.1 * fmin); // Margin factor
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
