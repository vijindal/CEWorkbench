package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.ClusterVariableEvaluator;
import org.ce.domain.cluster.LinearAlgebra;
import org.ce.domain.cluster.cvcf.CvCfBasis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.logging.Logger;

/**
 * Newton-Raphson solver for CVM free-energy minimization.
 *
 * <p>Restored to the 2012 working version with absolute error metrics,
 * boundary-hit detection, and 0.1 damping factor.</p>
 */
public final class NewtonRaphsonSolverSimple {

    private static final Logger LOG = Logger.getLogger(NewtonRaphsonSolverSimple.class.getName());

    /** Convergence thresholds from VJ logic. */
    private static final double TOLF = 1.0e-8;
    private static final double TOLX = 1.0e-12;

    /** Maximum iterations. */
    private static final int MAX_ITER = 400;

    /** Safety margin (0.1) as per VJ logic for stability. */
    private static final double STEP_MARGIN = 0.1;

    private NewtonRaphsonSolverSimple() { /* utility class */ }

    // =========================================================================
    // Result container
    // =========================================================================

    public static final class CVMSolverResult {

        public static final class IterationSnapshot {
            private final int iteration;
            private final double gibbsEnergy;
            private final double enthalpy;
            private final double entropy;
            private final double gradientNorm; // L1-norm as per VJ logic
            private final double[] cf;
            private final double[] dGdu;

            public IterationSnapshot(int iteration, double gibbsEnergy, double enthalpy,
                    double entropy, double gradientNorm, double[] cf, double[] dGdu) {
                this.iteration = iteration;
                this.gibbsEnergy = gibbsEnergy;
                this.enthalpy = enthalpy;
                this.entropy = entropy;
                this.gradientNorm = gradientNorm;
                this.cf = cf;
                this.dGdu = dGdu;
            }

            public int getIteration()       { return iteration; }
            public double getGibbsEnergy()  { return gibbsEnergy; }
            public double getEnthalpy()     { return enthalpy; }
            public double getEntropy()      { return entropy; }
            public double getGradientNorm() { return gradientNorm; }
            public double[] getCf()         { return cf; }
            public double[] getDGdu()       { return dGdu; }
        }

        private final double[] equilibriumCFs;
        private final double   gibbsEnergy;
        private final double   enthalpy;
        private final double   entropy;
        private final int      iterations;
        private final double   gradientNorm;
        private final boolean  converged;
        private final List<IterationSnapshot> iterationTrace;

        public CVMSolverResult(double[] equilibriumCFs, double gibbsEnergy, double enthalpy,
                double entropy, int iterations, double gradientNorm, boolean converged,
                List<IterationSnapshot> iterationTrace) {
            this.equilibriumCFs = equilibriumCFs;
            this.gibbsEnergy = gibbsEnergy;
            this.enthalpy = enthalpy;
            this.entropy = entropy;
            this.iterations = iterations;
            this.gradientNorm = gradientNorm;
            this.converged = converged;
            this.iterationTrace = new ArrayList<>(iterationTrace);
        }

        public double[] getEquilibriumCFs()            { return equilibriumCFs; }
        public double   getGibbsEnergy()               { return gibbsEnergy; }
        public double   getEnthalpy()                  { return enthalpy; }
        public double   getEntropy()                   { return entropy; }
        public int      getIterations()                { return iterations; }
        public double   getGradientNorm()              { return gradientNorm; }
        public boolean  isConverged()                  { return converged; }
        public List<IterationSnapshot> getIterationTrace() { return new ArrayList<>(iterationTrace); }
    }

    // =========================================================================
    // Data holder
    // =========================================================================

    public static final class CVMData {
        public final int tcdis, tcf, ncf;
        public final int[] lc;
        public final double[] kb;
        public final List<Double> mhdis;
        public final double[][] mh;
        public final int[][] lcv;
        public final List<List<int[]>> wcv;
        public final List<List<double[][]>> cmat;
        public final double[] eci;
        public final double temperature;
        public final double[] moleFractions;
        public final CvCfBasis basis;
        public final int[][] orthCfBasisIndices;

        public CVMData(int tcdis, int tcf, int ncf, int[] lc,
                double[] kb, List<Double> mhdis, double[][] mh,
                int[][] lcv, List<List<int[]>> wcv, List<List<double[][]>> cmat,
                double[] eci, double temperature,
                double[] moleFractions, CvCfBasis basis, int[][] orthCfBasisIndices) {
            this.tcdis = tcdis; this.tcf = tcf; this.ncf = ncf; this.lc = lc;
            this.kb = kb; this.mhdis = mhdis; this.mh = mh;
            this.lcv = lcv; this.wcv = wcv; this.cmat = cmat;
            this.eci = eci; this.temperature = temperature;
            this.moleFractions = moleFractions.clone();
            this.basis = basis;
            this.orthCfBasisIndices = orthCfBasisIndices;
        }
    }

    // =========================================================================
    // Public API
    // =========================================================================

    public static CVMSolverResult solve(
            double[] moleFractions, double temperature, double[] eci,
            List<Double> mhdis, double[] kb, double[][] mh, int[] lc,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv,
            int tcdis, int tcf, int ncf, CvCfBasis basis, int[][] orthCfBasisIndices,
            double tolerance) {

        return solve(moleFractions, temperature, eci, mhdis, kb, mh, lc,
                cmat, lcv, wcv, tcdis, tcf, ncf, basis, orthCfBasisIndices,
                MAX_ITER, tolerance);
    }

    public static CVMSolverResult solve(
            double[] moleFractions, double temperature, double[] eci,
            List<Double> mhdis, double[] kb, double[][] mh, int[] lc,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv,
            int tcdis, int tcf, int ncf, CvCfBasis basis, int[][] orthCfBasisIndices,
            int maxIter, double tolerance) {

        CVMData data = new CVMData(tcdis, tcf, ncf, lc, kb, mhdis, mh,
                lcv, wcv, cmat, eci, temperature, moleFractions, basis, orthCfBasisIndices);
        return minimize(data, maxIter, tolerance);
    }

    public static CVMSolverResult minimize(CVMData data) {
        return minimize(data, MAX_ITER, TOLF);
    }

    // =========================================================================
    // Implementation Loop (2012 Version)
    // =========================================================================

    private static CVMSolverResult minimize(CVMData data, int maxIter, double tolerance) {
        int n = data.ncf;
        int tcf = data.tcf;

        double xold[] = getURand(data);//Initializing with random values (full-sized tcf)
        double[][][] cvInit = updateCV(data, xold);

        // Diagnostic: Print Random State Cluster Variables (CVs)
        System.out.println("DEBUG: Initial Cluster Variables (CVs) from random CFs:");
        for (int t = 0; t < data.tcdis; t++) {
            for (int j = 0; j < data.lc[t]; j++) {
                System.out.print("    Type " + t + " Cluster " + j + ": [");
                for (int v = 0; v < data.lcv[t][j]; v++) {
                    System.out.printf("%.6f", cvInit[t][j][v]);
                    if (v < data.lcv[t][j] - 1) System.out.print(", ");
                }
                System.out.println("]");
            }
        }

        double xTrial[] = new double[tcf];
        double[] fvec = new double[n];//Gradient vector
        double[][] fjac = new double[n][n];

        System.arraycopy(xold, 0, xTrial, 0, tcf);

        // Diagnostic: Print Random State Thermodynamics
        CVMFreeEnergy.EvalResult result = CVMFreeEnergy.evaluate(
                xold, data.moleFractions, data.temperature, data.eci,
                data.mhdis, data.kb, data.mh, data.lc,
                data.cmat, data.lcv, data.wcv,
                data.tcdis, data.tcf, data.ncf);
        
        System.out.println("DEBUG: Random State Thermodynamics (x=" + Arrays.toString(data.moleFractions) + ", T=" + data.temperature + "K)");
        System.out.println("  H: " + result.H + " J/mol");
        System.out.println("  S: " + result.S + " J/mol-K");
        System.out.println("  G: " + result.G + " J/mol");
        
        System.out.println("  Per-Cluster Contributions:");
        for (int t = 0; t < data.tcdis; t++) {
            double hSub = result.hContribution[t];
            double sSub = result.sContribution[t];
            System.out.printf("    Type %d: H_sub=%.4f, S_sub=%.4f (gamma_kb=%.2f, m_s=%.2f)%n", 
                t, hSub, sSub, data.kb[t], data.mhdis.get(t));
        }

        double errf = 0, errx = 0;

        List<CVMSolverResult.IterationSnapshot> trace = new ArrayList<>();

        for (int its = 0; its < maxIter; its++) {
            if (Thread.currentThread().isInterrupted())
                throw new CancellationException("CVM cancelled at iter " + its);

            // 1. In-bounds check
            double[][][] cv_old = updateCV(data, xold);
            double cvMin = findMin(cv_old);
            if (cvMin <= 0) {
                LOG.warning("Boundary hit (cvMin=" + cvMin + ") at iter " + its);
                double[] vals = usrfun(data, xold, fvec, fjac);
                return new CVMSolverResult(xold, vals[0], vals[1], vals[2], its, errf, true, trace);
            }

            // 2. Evaluate
            double[] vals = usrfun(data, xold, fvec, fjac);
            double G = vals[0], H = vals[1], S = vals[2];

            // 3. errf convergence
            errf = 0;
            for (int i = 0; i < n; i++) errf += Math.abs(fvec[i]);
            
            trace.add(new CVMSolverResult.IterationSnapshot(
                    its, G, H, S, errf, xold.clone(), fvec.clone()));

            if (errf <= tolerance) {
                return new CVMSolverResult(xold, G, H, S, its, errf, true, trace);
            }

            // 4. Solve
            double[] p = new double[n];
            try {
                double[] negFvec = new double[n];
                for (int i = 0; i < n; i++) negFvec[i] = -fvec[i];
                p = LinearAlgebra.solve(fjac, negFvec);
            } catch (Exception e) {
                LOG.warning("Singular Hessian at iter " + its);
                return new CVMSolverResult(xold, G, H, S, its, errf, false, trace);
            }

            // 5. Step scaling & Update
            errx = 0;
            for (int i = 0; i < n; i++) {
                errx += Math.abs(p[i]);
                xTrial[i] = xold[i] + p[i];
            }
            
            double scalFactor = stpmx(data, xold, xTrial);

            for (int i = 0; i < n; i++) {
                xold[i] = xold[i] + scalFactor * p[i];
            }

            // 6. errx convergence
            if (errx <= TOLX) {
                vals = usrfun(data, xold, fvec, fjac);
                return new CVMSolverResult(xold, vals[0], vals[1], vals[2], its, errf, true, trace);
            }
        }

        double[] last = usrfun(data, xold, fvec, fjac);
        return new CVMSolverResult(xold, last[0], last[1], last[2], maxIter, errf, false, trace);
    }

    // =========================================================================
    // Utilities
    // =========================================================================

    private static double findMin(double[][][] cv) {
        double min = Double.MAX_VALUE;
        for (double[][] m : cv) for (double[] r : m) for (double v : r) if (v < min) min = v;
        return min;
    }

    private static double stpmx(CVMData data, double[] xold, double[] xTrial) {
        double fmin = 1.0;
        double[][][] cv_old = updateCV(data, xold);
        double[][][] cv_new = updateCV(data, xTrial);

        for (int i = 0; i < data.tcdis - 1; i++) {
            for (int j = 0; j < data.lc[i]; j++) {
                for (int v = 0; v < data.lcv[i][j]; v++) {
                    double vO = cv_old[i][j][v], vN = cv_new[i][j][v];
                    if (vN <= 0) fmin = Math.min(fmin, Math.abs(vO / (vN - vO)));
                    if (vN >= 1) fmin = Math.min(fmin, Math.abs((1.0 - vO) / (vN - vO)));
                }
            }
        }
        return (fmin >= 1.0) ? 1.0 : (STEP_MARGIN * fmin);
    }

    private static double[] usrfun(CVMData data, double[] u, double[] fvec, double[][] fjac) {
        CVMFreeEnergy.EvalResult r = CVMFreeEnergy.evaluate(
                u, data.moleFractions, data.temperature, data.eci, data.mhdis, data.kb, data.mh,
                data.lc, data.cmat, data.lcv, data.wcv, data.tcdis, data.tcf, data.ncf);
        System.arraycopy(r.Gcu, 0, fvec, 0, data.ncf);
        for (int i = 0; i < data.ncf; i++) System.arraycopy(r.Gcuu[i], 0, fjac[i], 0, data.ncf);
        return new double[]{r.G, r.H, r.S};
    }

    private static double[][][] updateCV(CVMData data, double[] u) {
        double[] uFull = ClusterVariableEvaluator.buildFullCVCFVector(u, data.moleFractions, data.ncf, data.tcf);
        return ClusterVariableEvaluator.evaluate(uFull, data.cmat, data.lcv, data.tcdis, data.lc);
    }

    private static double[] getURand(CVMData data) {
        return ClusterVariableEvaluator.generateFullRandomCVCFs(data.moleFractions, data.basis, data.orthCfBasisIndices);
    }
}
