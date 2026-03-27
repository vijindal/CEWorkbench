package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.ClusterVariableEvaluator;
import org.ce.domain.cluster.LinearAlgebra;
import org.ce.domain.cluster.cvcf.CvCfBasis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Newton-Raphson solver for CVM free-energy minimization.
 *
 * <p>Based on proven working implementation with simple 4-loop structure.
 * Minimizes G(u) = H(u) - T·S(u) over non-point correlation functions.</p>
 *
 * <p>Key features:</p>
 * <ul>
 *   <li>LU decomposition for solving linear system</li>
 *   <li>Step limiting to keep all CVs positive</li>
 *   <li>Random-state initial guess: u[icf] = (2xB - 1)^rank</li>
 * </ul>
 */
public final class NewtonRaphsonSolverSimple {

    private static final Logger LOG = Logger.getLogger(NewtonRaphsonSolverSimple.class.getName());

    /** Gas constant R = R_GAS (physical units, J/(mol·K)). Used only by legacy dead-code methods. */
    private static final double R = CVMFreeEnergy.R_GAS;

    /** Function tolerance for convergence. */
    private static final double TOLF = 1.0e-8;

    /** Step tolerance for convergence. */
    private static final double TOLX = 1.0e-12;

    /** Minimum CV value threshold. */
    private static final double CV_MIN = 1.0e-30;

    /** Maximum iterations. */
    private static final int MAX_ITER = 400;

    private NewtonRaphsonSolverSimple() { /* utility class */ }

    // =========================================================================
    // Result container
    // =========================================================================

    /** Immutable result of a CVM Newton-Raphson free-energy minimisation. */
    public static final class CVMSolverResult {

        /** Per-iteration Newton-Raphson diagnostics. */
        public static final class IterationSnapshot {
            private final int iteration;
            private final double gibbsEnergy;
            private final double enthalpy;
            private final double entropy;
            private final double gradientNorm;
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
                double entropy, int iterations, double gradientNorm, boolean converged) {
            this(equilibriumCFs, gibbsEnergy, enthalpy, entropy, iterations, gradientNorm,
                    converged, Collections.emptyList());
        }

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
    // Data holder for all CVM parameters
    // =========================================================================

    /**
     * Holds all CVM data needed for the N-R solver.
     */
    public static final class CVMData {
        public final int tcdis;           // number of HSP cluster types
        public final int tcf;             // total number of CFs
        public final int ncf;             // non-point CFs (optimization variables)
        public final int[] lc;            // ordered clusters per HSP type
        public final double[] kb;         // Kikuchi-Baker coefficients
        public final List<Double> msdis;  // HSP multiplicities
        public final double[][] m;        // normalized multiplicities mh[t][j]
        public final int[][] lcv;         // CV counts per (type, group)
        public final List<List<int[]>> wcv;     // CV weights
        public final List<List<double[][]>> cmat; // C-matrix
        public final double[] eci;        // effective cluster interactions
        public final double temperature;
        public final double xB;           // composition (mole fraction of B)
        public final double[] moleFractions;
        public final CvCfBasis basis;     // CVCF basis (holds Tinv for random init)

        public CVMData(
                int tcdis, int tcf, int ncf, int[] lc,
                double[] kb, List<Double> msdis, double[][] m,
                int[][] lcv, List<List<int[]>> wcv, List<List<double[][]>> cmat,
                double[] eci, double temperature, double xB,
                double[] moleFractions, CvCfBasis basis) {
            this.tcdis = tcdis;
            this.tcf = tcf;
            this.ncf = ncf;
            this.lc = lc;
            this.kb = kb;
            this.msdis = msdis;
            this.m = m;
            this.lcv = lcv;
            this.wcv = wcv;
            this.cmat = cmat;
            this.eci = eci;
            this.temperature = temperature;
            this.xB = xB;
            this.moleFractions = moleFractions.clone();
            this.basis = basis;
        }
    }

    // =========================================================================
    // Main solver
    // =========================================================================

    /**
     * Solves the CVM equilibrium problem.
     */
    public static CVMSolverResult solve(
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
            int ncf,
            CvCfBasis basis,
            int maxIter,
            double tolerance,
            double step) {

        // For binary: xB = moleFractions[1]
        double xB = moleFractions.length > 1 ? moleFractions[1] : moleFractions[0];

        CVMData data = new CVMData(tcdis, tcf, ncf, lc, kb, mhdis, mh,
                lcv, wcv, cmat, eci, temperature, xB, moleFractions, basis);

        return minimize(data, maxIter, tolerance);
    }

    /**
     * Convenience overload using default parameters.
     */
    public static CVMSolverResult solve(
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
            int ncf,
            CvCfBasis basis,
            double tolerance) {

        return solve(moleFractions, temperature, eci, mhdis, kb, mh, lc,
                cmat, lcv, wcv, tcdis, tcf, ncf, basis,
                MAX_ITER, tolerance, 1.0);
    }

    // =========================================================================
    // N-R minimization (from proven working code)
    // =========================================================================

    /**
     * Newton-Raphson minimization loop.
     */
    private static CVMSolverResult minimize(CVMData data, int maxIter, double tolerance) {
        int ncf = data.ncf;

        // Initialize CFs at random state
        double[] u = getURand(data);

        // Compute initial CVs
        double[][][] cv = updateCV(data, u);

        LOG.fine("NewtonRaphsonSolverSimple.minimize — ENTER: ncf=" + ncf
                + ", tcf=" + data.tcf + ", T=" + data.temperature
                + ", xB=" + data.xB + ", tolerance=" + tolerance);

        double[] Gu = new double[ncf];
        double[][] Guu = new double[ncf][ncf];

        double G = 0, H = 0, S = 0;
        double gradNorm = Double.MAX_VALUE;
        int iter = 0;
        List<CVMSolverResult.IterationSnapshot> trace = new ArrayList<>();

        for (iter = 1; iter <= maxIter; iter++) {
            if (Thread.currentThread().isInterrupted())
                throw new CancellationException("CVM cancelled at iteration " + iter);
            // Evaluate function and derivatives
            double[] vals = usrfun(data, u, Gu, Guu);
            G = vals[0];
            H = vals[1];
            S = vals[2];

            // Check convergence on gradient norm
            gradNorm = 0.0;
            for (int i = 0; i < ncf; i++) {
                gradNorm += Gu[i] * Gu[i];
            }
            gradNorm = Math.sqrt(gradNorm);

            trace.add(new CVMSolverResult.IterationSnapshot(
                    iter, G, H, S, gradNorm, u.clone(), Gu.clone()));

            if (LOG.isLoggable(Level.FINEST) && (iter == 1 || iter % 20 == 0 || gradNorm < 1e-6)) {
                LOG.finest(String.format("NewtonRaphsonSolverSimple — iter %3d: G=%.8e H=%.8e S=%.8e ||Gu||=%.8e",
                        iter, G, H, S, gradNorm));
            }

            if (gradNorm < tolerance) {
                LOG.fine(String.format("NewtonRaphsonSolverSimple — CONVERGED (gradient norm < tolerance): iter=%d, ||Gu||=%.4e",
                        iter, gradNorm));
                return new CVMSolverResult(u, G, H, S, iter, gradNorm, true, trace);
            }

            // Solve Guu · du = -Gu using Gaussian elimination
            double[] du;
            try {
                double[] negGu = new double[ncf];
                for (int i = 0; i < ncf; i++) negGu[i] = -Gu[i];
                du = LinearAlgebra.solve(Guu, negGu);
            } catch (IllegalArgumentException e) {
                LOG.warning("NewtonRaphsonSolverSimple — SINGULAR HESSIAN at iter=" + iter + ": " + e.getMessage());
                return new CVMSolverResult(u, G, H, S, iter, gradNorm, false, trace);
            }

            // Step limiting to keep CVs positive
            double stpmax = stpmx(data, u, du, cv);

            // Update CFs
            for (int i = 0; i < ncf; i++) {
                u[i] += stpmax * du[i];
            }

            // Update CVs
            cv = updateCV(data, u);

            // Check for small step (secondary convergence)
            double stepNorm = 0.0;
            for (int i = 0; i < ncf; i++) {
                stepNorm += du[i] * du[i];
            }
            stepNorm = Math.sqrt(stepNorm) * stpmax;

            if (stepNorm < TOLX) {
                // Recalculate final values and check true stationarity
                vals = usrfun(data, u, Gu, Guu);
                double finalGradNorm = 0.0;
                for (int i = 0; i < ncf; i++) {
                    finalGradNorm += Gu[i] * Gu[i];
                }
                finalGradNorm = Math.sqrt(finalGradNorm);

                if (finalGradNorm < tolerance) {
                    LOG.fine(String.format("NewtonRaphsonSolverSimple — CONVERGED (step < TOLX and gradient < tolerance): iter=%d, ||Gu||=%.4e",
                            iter, finalGradNorm));
                    return new CVMSolverResult(u, vals[0], vals[1], vals[2], iter, finalGradNorm, true, trace);
                }

                LOG.warning(String.format("NewtonRaphsonSolverSimple — STALLED (step < TOLX but ||Gu||=%.4e > tolerance=%.4e): iter=%d",
                        finalGradNorm, tolerance, iter));
                return new CVMSolverResult(u, vals[0], vals[1], vals[2], iter, finalGradNorm, false, trace);
            }
        }

        LOG.warning(String.format("NewtonRaphsonSolverSimple — NOT CONVERGED after %d iterations (||Gu||=%.4e)",
                maxIter, gradNorm));
        return new CVMSolverResult(u, G, H, S, maxIter, gradNorm, false, trace);
    }

    // =========================================================================
    // usrfun: evaluate G, Gu, Guu
    // =========================================================================

    /**
     * Evaluates Gibbs energy, gradient, and Hessian.
     *
     * @return [G, H, S]
     */
    private static double[] usrfun(CVMData data, double[] u,
                                    double[] Gu, double[][] Guu) {
        CVMFreeEnergy.EvalResult eval = CVMFreeEnergy.evaluate(
                u,
                data.moleFractions,
                data.temperature,
                data.eci,
                data.msdis,
                data.kb,
                data.m,
                data.lc,
                data.cmat,
                data.lcv,
                data.wcv,
                data.tcdis,
                data.tcf,
                data.ncf
        );

        System.arraycopy(eval.Gcu, 0, Gu, 0, data.ncf);
        for (int i = 0; i < data.ncf; i++) {
            System.arraycopy(eval.Gcuu[i], 0, Guu[i], 0, data.ncf);
        }

        return new double[]{eval.G, eval.H, eval.S};
    }

    // =========================================================================
    // CV calculation
    // =========================================================================

    /**
     * CV[t][j][v] = Σ_icf cmat[v][icf] × u[icf] + cmat[v][tcf] (constant term)
     */
    private static double[][][] updateCV(CVMData data, double[] u) {
        double[][][] cv = new double[data.tcdis][][];

        // Build full CVCF vector: [v | moleFractions]
        double[] uFull = ClusterVariableEvaluator.buildFullCVCFVector(
                u, data.moleFractions, data.ncf, data.tcf);

        for (int itc = 0; itc < data.tcdis; itc++) {
            cv[itc] = new double[data.lc[itc]][];
            for (int inc = 0; inc < data.lc[itc]; inc++) {
                int nv = data.lcv[itc][inc];
                cv[itc][inc] = new double[nv];
                double[][] cm = data.cmat.get(itc).get(inc);

                for (int incv = 0; incv < nv; incv++) {
                    double sum = cm[incv][data.tcf]; // constant term
                    for (int icf = 0; icf < data.tcf; icf++) {
                        sum += cm[incv][icf] * uFull[icf];
                    }
                    cv[itc][inc][incv] = sum;
                }
            }
        }
        return cv;
    }

    // =========================================================================
    // Initial guess: random state
    // =========================================================================

    /**
     * Returns random-state (disordered) CF values for all non-point CFs.
     */
    private static double[] getURand(CVMData data) {
        return ClusterVariableEvaluator.computeRandomCVCFs(data.moleFractions, data.basis);
    }

    // =========================================================================
    // Step limiting (stpmx)
    // =========================================================================

    /**
     * Find maximum step size that keeps all CVs positive.
     *
     * <p>For each CV: cv_new = cv_old + stpmax × Δcv
     * where Δcv = Σ_icf cmat[v][icf] × du[icf]</p>
     *
     * <p>If Δcv < 0, max step = -cv_old / Δcv (but slightly smaller)</p>
     */
    private static double stpmx(CVMData data, double[] u, double[] du, double[][][] cv) {
        double stpmax = 1.0;
        final double MARGIN = 0.99; // safety margin

        for (int itc = 0; itc < data.tcdis; itc++) {
            for (int inc = 0; inc < data.lc[itc]; inc++) {
                double[][] cm = data.cmat.get(itc).get(inc);
                int nv = data.lcv[itc][inc];

                for (int incv = 0; incv < nv; incv++) {
                    // Compute change in CV per unit step
                    double dcv = 0.0;
                    for (int icf = 0; icf < data.ncf; icf++) {
                        dcv += cm[incv][icf] * du[icf];
                    }

                    // If CV would decrease, limit step
                    if (dcv < -CV_MIN) {
                        double cvOld = cv[itc][inc][incv];
                        if (cvOld > CV_MIN) {
                            double maxStep = -MARGIN * cvOld / dcv;
                            if (maxStep < stpmax && maxStep > 0) {
                                stpmax = maxStep;
                            }
                        }
                    }
                }
            }
        }

        // Minimum step floor
        return Math.max(stpmax, 1.0e-6);
    }

    // =========================================================================
    // Utilities
    // =========================================================================

    /**
     * Count number of CFs for a given cluster type.
     * For binary: each type has 1 CF.
     */
    private static int countCFsForType(CVMData data, int typeIdx) {
        return data.lc[typeIdx];
    }
}
