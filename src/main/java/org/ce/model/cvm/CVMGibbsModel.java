package org.ce.model.cvm;

import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cluster.SpaceGroup;
import org.ce.model.storage.InputLoader;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.PhysicsConstants;
import org.ce.model.ProgressEvent;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.function.Consumer;

/**
 * Physical model for the Cluster Variation Method (CVM).
 *
 * <p>Encapsulates the structural geometry, multiplicities, and entropy coefficients.
 * Provides the Gibbs free energy functional, its derivatives, and owns the full
 * Newton-Raphson equilibrium minimisation loop.</p>
 */
public class CVMGibbsModel {

    private static final double ENTROPY_SMOOTH_EPS = 1.0e-6;
    private static final int    MAX_ITER           = 20;
    private static final double TOLX               = 1.0e-12;

    private String elements;
    private int numComponents;
    private int tcdis, tcf, ncf;
    private double[] mhdis;
    private double[] kb;
    private double[][] mh;
    private int[] lc;
    private List<List<double[][]>> cmat;
    private int[][] lcv;
    private List<List<int[]>> wcv;
    private int[][] orthCfBasisIndices;
    private CvCfBasis basis;
    private CECEntry cecEntry;

    // Cached minimisation result — invalidated when T or composition changes
    private boolean isMinimized = false;
    private double currentTemperature = -1.0;
    private double[] currentComposition = null;
    private EquilibriumResult lastResult = null;

    // =========================================================================
    // Inner result type
    // =========================================================================

    /** Result returned by {@link #getEquilibriumState}. */
    public static final class EquilibriumResult {
        /** Physics values at the equilibrium point. */
        public final ModelResult modelResult;
        /** Equilibrium non-point CVCF correlation functions (length = ncf). */
        public final double[] u;
        /** Convergence flag. Check before using modelResult. */
        public final boolean converged;
        /** Iteration count at convergence or failure. */
        public final int iterations;
        /** Final gradient norm ||∇G|| at exit. */
        public final double finalGradientNorm;

        public EquilibriumResult(ModelResult modelResult, double[] u, boolean converged,
                int iterations, double finalGradientNorm) {
            this.modelResult       = modelResult;
            this.u                 = u;
            this.converged         = converged;
            this.iterations        = iterations;
            this.finalGradientNorm = finalGradientNorm;
        }
    }

    // =========================================================================
    // Inner physics result type
    // =========================================================================

    /** Calculated free energy and derivatives at a given (u, T, x) point. */
    public static class ModelResult {
        public final double G, H, S;
        public final double[] Gu;
        public final double[][] Guu;
        public final double[] Hu;
        public final double[] Su;
        public final double[][] Suu;
        public final double[] cfs;

        public ModelResult(double G, double H, double S,
                double[] Gu, double[][] Guu,
                double[] Hu, double[] Su, double[][] Suu,
                double[] cfs) {
            this.G   = G;   this.H   = H;   this.S   = S;
            this.Gu  = Gu;  this.Guu = Guu;
            this.Hu  = Hu;  this.Su  = Su;  this.Suu = Suu;
            this.cfs = cfs;
        }
    }

    /** Default constructor for lazy initialization. */
    public CVMGibbsModel() {}

    // =========================================================================
    // Initialization
    // =========================================================================

    /**
     * Primary entry point for model formation. Orchestrates the 4-stage pipeline
     * to identify clusters, generate C-matrix, and resolve CVCF basis using
     * the system identity.
     */
    public void initialize(
            String elements,
            String structure,
            String model,
            CECEntry cecEntry,
            Consumer<String> progressSink) {

        this.elements = elements;
        this.numComponents = elements.split("-").length;
        this.cecEntry = cecEntry;

        org.ce.model.cluster.ClusterIdentificationRequest request =
                org.ce.model.cluster.ClusterIdentificationRequest.fromSystem(elements, structure, model);

        if (progressSink != null) progressSink.accept("\n[STAGE 0]: Loading Inputs...");
        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(request.getDisorderedClusterFile());
        disorderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup disorderedSG = InputLoader.parseSpaceGroup(request.getDisorderedSymmetryGroup());

        List<Cluster> orderedClusters = InputLoader.parseClusterFile(request.getOrderedClusterFile());
        orderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup orderedSG = InputLoader.parseSpaceGroup(request.getOrderedSymmetryGroup());

        org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult pr =
                org.ce.model.cluster.ClusterCFIdentificationPipeline.run(
                        disorderedClusters,
                        disorderedSG.getOperations(),
                        orderedClusters,
                        orderedSG.getOperations(),
                        request.getTransformationMatrix(),
                        new double[] {
                            request.getTranslationVector().getX(),
                            request.getTranslationVector().getY(),
                            request.getTranslationVector().getZ()
                        },
                        request.getNumComponents(),
                        progressSink);

        org.ce.model.cluster.CMatrixPipeline.CMatrixData mdOrth =
                org.ce.model.cluster.CMatrixPipeline.run(
                    pr.toClusterIdentificationResult(),
                    pr.toCFIdentificationResult(),
                    disorderedClusters,
                    pr.getNumComponents(),
                    progressSink);

        org.ce.model.cvm.CvCfBasis basisRef =
                org.ce.model.cvm.CvCfBasis.generate(structure, pr, mdOrth, model, progressSink);

        org.ce.model.cluster.CMatrixPipeline.CMatrixData mdCvcf = basisRef.cvcfCMatrixData;

        this.tcdis = pr.getTcdis();
        this.mhdis = pr.getMhdis();
        this.kb    = pr.getKbdis();
        this.mh    = pr.getMh();
        this.lc    = pr.getLc();

        this.cmat  = mdCvcf.getCmat();
        this.lcv   = mdCvcf.getLcv();
        this.wcv   = mdCvcf.getWcv();
        this.orthCfBasisIndices = mdCvcf.getCfBasisIndices();

        this.basis = basisRef;
        this.ncf   = basisRef.numNonPointCfs;
        this.tcf   = basisRef.totalCfs();
    }

    // =========================================================================
    // Equilibrium resolution (Newton-Raphson loop)
    // =========================================================================

    /**
     * Returns the equilibrium state at the given (T, x). Result is cached and
     * reused on repeated calls at the same conditions.
     */
    public EquilibriumResult getEquilibriumState(
            double temperature, double[] composition, double tolerance,
            Consumer<String> progressSink, Consumer<ProgressEvent> eventSink) {

        boolean compositionChanged = currentComposition == null
                || !Arrays.equals(currentComposition, composition);
        boolean temperatureChanged = Math.abs(currentTemperature - temperature) > 1.0e-5;

        if (temperatureChanged || compositionChanged) {
            if (progressSink != null)
                progressSink.accept("  [Model] Parameters updated (T=" + temperature + " K) — resetting equilibrium state.");
            isMinimized = false;
            currentTemperature = temperature;
            currentComposition = Arrays.copyOf(composition, composition.length);
        } else {
            if (progressSink != null)
                progressSink.accept("  [Model] Using currently cached equilibrium state.");
        }

        if (!isMinimized || lastResult == null) {
            if (progressSink != null)
                progressSink.accept("  [Model] Initiating internal minimization (Newton-Raphson loop)...");
            lastResult = minimize(composition, temperature, tolerance, progressSink, eventSink);
            if (progressSink != null) {
                progressSink.accept(lastResult.converged
                        ? "  [Model] ✓ Minimization converged in " + lastResult.iterations + " iterations."
                        : "  [Model] ⚠ Minimization FAILED to converge.");
            }
            isMinimized = true;
        }

        return lastResult;
    }

    /**
     * Newton-Raphson minimization of the CVM Gibbs free energy.
     * Fully self-contained: uses only fields of this instance.
     */
    private EquilibriumResult minimize(
            double[] moleFractions, double temperature, double tolerance,
            Consumer<String> progressSink, Consumer<ProgressEvent> eventSink) {

        if (eventSink != null)
            eventSink.accept(new ProgressEvent.EngineStart("CVM", 0));

        double[] u = computeRandomCFs(moleFractions);

        double errf = 0;
        ModelResult current = null;

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted())
                throw new CancellationException();

            current = evaluate(u, moleFractions, temperature);

            errf = 0;
            for (double g : current.Gu)
                errf += Math.abs(g);

            if (progressSink != null) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("    iter %3d  |∇G| = %.3e  G = %11.4f  H = %11.4f  S = %9.6f",
                        its, errf, current.G, current.H, current.S));
                if (current.cfs != null && current.cfs.length > 0) {
                    sb.append("  CFs: [");
                    for (int i = 0; i < current.cfs.length; i++) {
                        if (i > 0) sb.append(", ");
                        sb.append(String.format("%.4f", current.cfs[i]));
                    }
                    sb.append("]");
                }
                progressSink.accept(sb.toString());
            }
            if (eventSink != null)
                eventSink.accept(new ProgressEvent.CvmIteration(its, current.G, errf, current.H, current.S, u));

            if (errf <= tolerance)
                return new EquilibriumResult(current, u.clone(), true, its, errf);

            try {
                double[] negGu = new double[ncf];
                for (int i = 0; i < ncf; i++)
                    negGu[i] = -current.Gu[i];

                double[] p = LinearAlgebra.solve(current.Guu, negGu);

                double alpha = calculateStepLimit(u, p, moleFractions);

                double errx = 0;
                for (int i = 0; i < ncf; i++) {
                    double delta = alpha * p[i];
                    u[i] += delta;
                    errx += Math.abs(delta);
                }

                if (errx <= TOLX) {
                    current = evaluate(u, moleFractions, temperature);
                    return new EquilibriumResult(current, u.clone(), true, its, errf);
                }

            } catch (Exception e) {
                return new EquilibriumResult(current, u.clone(), false, its, errf);
            }
        }

        return new EquilibriumResult(current, u.clone(), false, MAX_ITER, errf);
    }

    // =========================================================================
    // Physics evaluation
    // =========================================================================

    /** Evaluates the CVM Gibbs free energy and derivatives at the given state. */
    public ModelResult evaluate(double[] u, double[] moleFractions, double temperature) {
        return evaluateInternal(u, moleFractions, temperature, cecEntry, basis,
                tcdis, ncf, mhdis, kb, mh, lc, cmat, lcv, wcv);
    }

    /** Static overload for callers that supply all parameters explicitly. */
    public static ModelResult evaluate(
            double[] u, double[] moleFractions, double temperature, CECEntry cecEntry, CvCfBasis basis,
            int tcdis, int ncf, double[] mhdis, double[] kb, double[][] mh, int[] lc,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv) {
        return evaluateInternal(u, moleFractions, temperature, cecEntry, basis,
                tcdis, ncf, mhdis, kb, mh, lc, cmat, lcv, wcv);
    }

    private static ModelResult evaluateInternal(
            double[] u, double[] moleFractions, double temperature, CECEntry cecEntry, CvCfBasis basis,
            int tcdis, int ncf, double[] mhdis, double[] kb, double[][] mh, int[] lc,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv) {

        double[] eci = CECEvaluator.evaluate(cecEntry, temperature, basis, "CVM");

        double[] uFull = CMatrixPipeline.buildFullCVCFVector(u, moleFractions, ncf);
        double[][][] cv = CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);

        double Hval = 0.0;
        double Sval = 0.0;
        double[] Su  = new double[ncf];
        double[][] Suu = new double[ncf][ncf];

        for (int l = 0; l < ncf; l++)
            Hval += eci[l] * u[l];

        for (int t = 0; t < tcdis; t++) {
            double coeff_t = kb[t] * mhdis[t];
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
                        logEff   = logCv;
                        invEff   = 1.0 / cvVal;
                    } else {
                        double logEps = Math.log(ENTROPY_SMOOTH_EPS);
                        double d = cvVal - ENTROPY_SMOOTH_EPS;
                        sContrib = ENTROPY_SMOOTH_EPS * logEps + (1.0 + logEps) * d + 0.5 / ENTROPY_SMOOTH_EPS * d * d;
                        logEff   = logEps + d / ENTROPY_SMOOTH_EPS;
                        invEff   = 1.0 / ENTROPY_SMOOTH_EPS;
                    }

                    double prefix = coeff_t * mh_tj * wv;
                    Sval -= PhysicsConstants.R_GAS * prefix * sContrib;

                    for (int l = 0; l < ncf; l++) {
                        double cml = cm[incv][l];
                        if (cml != 0.0)
                            Su[l] -= PhysicsConstants.R_GAS * prefix * cml * logEff;
                    }

                    for (int l1 = 0; l1 < ncf; l1++) {
                        double cml1 = cm[incv][l1];
                        if (cml1 == 0.0) continue;
                        for (int l2 = l1; l2 < ncf; l2++) {
                            double cml2 = cm[incv][l2];
                            if (cml2 == 0.0) continue;
                            double val = -PhysicsConstants.R_GAS * prefix * cml1 * cml2 * invEff;
                            Suu[l1][l2] += val;
                            if (l1 != l2) Suu[l2][l1] += val;
                        }
                    }
                }
            }
        }

        double Gval = Hval - temperature * Sval;
        double[] Gu  = new double[ncf];
        double[][] Guu = new double[ncf][ncf];
        for (int l = 0; l < ncf; l++) {
            Gu[l] = eci[l] - temperature * Su[l];
            for (int l2 = 0; l2 < ncf; l2++)
                Guu[l][l2] = -temperature * Suu[l][l2];
        }

        return new ModelResult(Gval, Hval, Sval, Gu, Guu, eci, Su, Suu, uFull);
    }

    // =========================================================================
    // Helpers used by the N-R loop
    // =========================================================================

    /**
     * Computes the full disordered-state (random) CVCF vector for N-R initialisation.
     */
    public double[] computeRandomCFs(double[] moleFractions) {
        return basis.computeRandomState(moleFractions, orthCfBasisIndices);
    }

    /**
     * Finds the maximum step size α ∈ (0, 1] such that u + α·p keeps all
     * cluster variables within [0, 1].
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
                    if (vN <= 0) fmin = Math.min(fmin, Math.abs(vO / (vN - vO)));
                    if (vN >= 1) fmin = Math.min(fmin, Math.abs((1.0 - vO) / (vN - vO)));
                }
            }
        }
        return (fmin >= 1.0) ? 1.0 : (0.1 * fmin);
    }

    private double[][][] updateCVInternal(double[] u, double[] moleFractions) {
        double[] uFull = CMatrixPipeline.buildFullCVCFVector(u, moleFractions, ncf);
        return CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    public int      getNumComponents()      { return numComponents; }
    public String   getElements()           { return elements; }
    public CvCfBasis getBasis()             { return basis; }
    public int      getNcf()               { return ncf; }
    public int      getTcf()               { return tcf; }
    public int      getTcdis()             { return tcdis; }
    public int[]    getLc()                { return lc; }
    public int[][]  getLcv()               { return lcv; }
    public int[][]  getOrthCfBasisIndices() { return orthCfBasisIndices; }

    /** Computes cluster variables cv[t][j][v] from the given non-point CFs and composition. */
    public double[][][] evaluateClusterVariables(double[] u, double[] moleFractions) {
        double[] uFull = CMatrixPipeline.buildFullCVCFVector(u, moleFractions, ncf);
        return CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);
    }
}
