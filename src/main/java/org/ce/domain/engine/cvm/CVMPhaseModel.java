package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.CMatrix;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.ClusterVariableEvaluator;
import org.ce.domain.cluster.LinearAlgebra;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.result.EquilibriumState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Central thermodynamic model for CVM free-energy calculations.
 */
public class CVMPhaseModel {

    // =========================================================================
    // Solver Result Containers
    // =========================================================================

    public static final class SolverResult {

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

        public SolverResult(double[] equilibriumCFs, double gibbsEnergy, double enthalpy,
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
    // Input contract
    // =========================================================================

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

    // =========================================================================

    private static final Logger LOG = Logger.getLogger(CVMPhaseModel.class.getName());
    private static final double DEFAULT_TOLERANCE = 1.0e-5;
    private static final double MAX_TOLERANCE = 1.0e-3;

    private final int tcdis, tcf, ncf;
    private final List<Double> mhdis;
    private final double[] kb;
    private final double[][] mh;
    private final int[] lc;
    private final List<List<double[][]>> cmat;
    private final int[][] lcv;
    private final List<List<int[]>> wcv;
    private final int[][] orthCfBasisIndices;
    private final String systemId;
    private final String systemName;
    private final int numComponents;
    private final CvCfBasis basis;

    private double[] eci;
    private double tolerance;
    private double temperature;
    private double[] moleFractions;

    private boolean isMinimized = false;
    private long lastMinimizationTimeNanos;
    private double[] equilibriumCFs;
    private CVMFreeEnergy.EvalResult equilibrium;

    private int lastIterations;
    private double lastGradientNorm;
    private String lastConvergenceStatus;
    private List<SolverResult.IterationSnapshot> lastIterationTrace = new ArrayList<>();

    public static CVMPhaseModel create(
            CVMInput input,
            double[] eci,
            double temperature,
            double[] moleFractions) throws Exception {
        CVMPhaseModel model = new CVMPhaseModel(input);
        model.setECI(eci);
        model.setTemperature(temperature);
        model.setMoleFractions(moleFractions);
        model.ensureMinimized();
        return model;
    }

    private CVMPhaseModel(CVMInput input) {
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

        this.systemId = input.getSystemId();
        this.systemName = input.getSystemName();
        this.numComponents = input.getNumComponents();
        this.basis = input.getBasis();
        this.ncf = basis.numNonPointCfs;
        this.tcf = basis.totalCfs();

        this.tolerance = DEFAULT_TOLERANCE;
    }

    public void setECI(double[] newECI) {
        this.eci = Arrays.copyOf(newECI, ncf);
        invalidateMinimization();
    }

    public void setTemperature(double T_K) {
        this.temperature = T_K;
        invalidateMinimization();
    }

    public void setMoleFractions(double[] fractions) {
        this.moleFractions = fractions.clone();
        invalidateMinimization();
    }

    public void setTolerance(double tol) {
        this.tolerance = tol;
        invalidateMinimization();
    }

    private void invalidateMinimization() {
        this.isMinimized = false;
        this.equilibriumCFs = null;
        this.equilibrium = null;
        this.lastIterationTrace = new ArrayList<>();
    }

    public synchronized void ensureMinimized() throws Exception {
        if (isMinimized) return;
        minimize();
        if (!isMinimized) throw new Exception("Minimization failed: " + lastConvergenceStatus);
    }

    private void minimize() {
        long startTime = System.nanoTime();
        try {
            SolverResult result = runNewtonRaphson();
            this.lastIterationTrace = result.getIterationTrace();
            if (result.isConverged()) {
                this.equilibriumCFs = result.getEquilibriumCFs();
                this.equilibrium = CVMFreeEnergy.evaluate(
                    equilibriumCFs, moleFractions, temperature, eci,
                    mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf
                );
                this.isMinimized = true;
                this.lastIterations = result.getIterations();
                this.lastGradientNorm = result.getGradientNorm();
                this.lastConvergenceStatus = "OK";
                this.lastMinimizationTimeNanos = System.nanoTime() - startTime;
            } else {
                this.isMinimized = false;
                this.lastConvergenceStatus = "Solver did not converge";
                this.lastIterations = result.getIterations();
                this.lastGradientNorm = result.getGradientNorm();
            }
        } catch (Exception e) {
            this.isMinimized = false;
            this.lastConvergenceStatus = e.getMessage();
        }
    }

    // =========================================================================
    // Solver Implementation
    // =========================================================================

    private static final int MAX_ITER = 400;
    private static final double TOLX = 1.0e-12;
    private static final double STEP_MARGIN = 0.1;

    private SolverResult runNewtonRaphson() {
        int n = ncf;
        double xold[] = ClusterVariableEvaluator.generateFullRandomCVCFs(moleFractions, basis, orthCfBasisIndices);
        double xTrial[] = new double[tcf];
        double[] fvec = new double[n];
        double[][] fjac = new double[n][n];
        System.arraycopy(xold, 0, xTrial, 0, tcf);

        double errf = 0, errx = 0;
        List<SolverResult.IterationSnapshot> trace = new ArrayList<>();

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted()) throw new CancellationException();

            double[][][] cv_old = updateCVInternal(xold);
            double cvMin = findMinInternal(cv_old);
            if (cvMin <= 0) {
                double[] vals = usrfunInternal(xold, fvec, fjac);
                return new SolverResult(xold, vals[0], vals[1], vals[2], its, errf, true, trace);
            }

            double[] vals = usrfunInternal(xold, fvec, fjac);
            errf = 0;
            for (int i = 0; i < n; i++) errf += Math.abs(fvec[i]);
            trace.add(new SolverResult.IterationSnapshot(its, vals[0], vals[1], vals[2], errf, xold.clone(), fvec.clone()));

            if (errf <= tolerance) return new SolverResult(xold, vals[0], vals[1], vals[2], its, errf, true, trace);

            try {
                double[] negFvec = new double[n];
                for (int i = 0; i < n; i++) negFvec[i] = -fvec[i];
                double[] p = LinearAlgebra.solve(fjac, negFvec);
                errx = 0;
                for (int i = 0; i < n; i++) {
                    errx += Math.abs(p[i]);
                    xTrial[i] = xold[i] + p[i];
                }
                double scalFactor = stpmxInternal(xold, xTrial);
                for (int i = 0; i < n; i++) xold[i] = xold[i] + scalFactor * p[i];
                if (errx <= TOLX) {
                    vals = usrfunInternal(xold, fvec, fjac);
                    return new SolverResult(xold, vals[0], vals[1], vals[2], its, errf, true, trace);
                }
            } catch (Exception e) {
                return new SolverResult(xold, vals[0], vals[1], vals[2], its, errf, false, trace);
            }
        }
        double[] last = usrfunInternal(xold, fvec, fjac);
        return new SolverResult(xold, last[0], last[1], last[2], MAX_ITER, errf, false, trace);
    }

    private double findMinInternal(double[][][] cv) {
        double min = Double.MAX_VALUE;
        for (double[][] m : cv) for (double[] r : m) for (double v : r) if (v < min) min = v;
        return min;
    }

    private double stpmxInternal(double[] xold, double[] xTrial) {
        double fmin = 1.0;
        double[][][] cv_old = updateCVInternal(xold);
        double[][][] cv_new = updateCVInternal(xTrial);
        for (int i = 0; i < tcdis - 1; i++) {
            for (int j = 0; j < lc[i]; j++) {
                for (int v = 0; v < lcv[i][j]; v++) {
                    double vO = cv_old[i][j][v], vN = cv_new[i][j][v];
                    if (vN <= 0) fmin = Math.min(fmin, Math.abs(vO / (vN - vO)));
                    if (vN >= 1) fmin = Math.min(fmin, Math.abs((1.0 - vO) / (vN - vO)));
                }
            }
        }
        return (fmin >= 1.0) ? 1.0 : (STEP_MARGIN * fmin);
    }

    private double[] usrfunInternal(double[] u, double[] fvec, double[][] fjac) {
        CVMFreeEnergy.EvalResult r = CVMFreeEnergy.evaluate(
                u, moleFractions, temperature, eci, mhdis, kb, mh,
                lc, cmat, lcv, wcv, tcdis, tcf, ncf);
        System.arraycopy(r.Gcu, 0, fvec, 0, ncf);
        for (int i = 0; i < ncf; i++) System.arraycopy(r.Gcuu[i], 0, fjac[i], 0, ncf);
        return new double[]{r.G, r.H, r.S};
    }

    private double[][][] updateCVInternal(double[] u) {
        double[] uFull = ClusterVariableEvaluator.buildFullCVCFVector(u, moleFractions, ncf, tcf);
        return ClusterVariableEvaluator.evaluate(uFull, cmat, lcv, tcdis, lc);
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    public double getEquilibriumG() throws Exception { ensureMinimized(); return equilibrium.G; }
    public double getEquilibriumH() throws Exception { ensureMinimized(); return equilibrium.H; }
    public double getEquilibriumS() throws Exception { ensureMinimized(); return equilibrium.S; }
    public double[] getEquilibriumCFs() throws Exception { ensureMinimized(); return equilibriumCFs.clone(); }
    public int getLastIterations() { return lastIterations; }
    public double getLastGradientNorm() { return lastGradientNorm; }
    public List<SolverResult.IterationSnapshot> getLastIterationTrace() { return new ArrayList<>(lastIterationTrace); }

    public double[] computeGHS(double[] nonPointCFs) {
        CVMFreeEnergy.EvalResult r = CVMFreeEnergy.evaluate(
            nonPointCFs, moleFractions, temperature, eci,
            mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);
        return new double[]{r.G, r.H, r.S};
    }
}
