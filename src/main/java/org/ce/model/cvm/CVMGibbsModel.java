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
 * <p>
 * Encapsulates the structural geometry, multiplicities, and entropy
 * coefficients.
 * Provides the Gibbs free energy functional, its derivatives, and owns the full
 * Newton-Raphson equilibrium minimisation loop.
 * </p>
 */
public class CVMGibbsModel {

    private static final double ENTROPY_SMOOTH_EPS = 1.0e-6;
    private static final int MAX_ITER = 20;
    private static final double TOLX = 1.0e-12;

    private String elements;
    private String structure;
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
    private double[] eci;
    private org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult pipelineResult;

    // --- Current State for Standalone Methods ---
    private double[] u;
    private double[] x_mole;
    private double temp;
    private double[][][] currentCv;

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
            this.modelResult = modelResult;
            this.u = u;
            this.converged = converged;
            this.iterations = iterations;
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
            this.G = G;
            this.H = H;
            this.S = S;
            this.Gu = Gu;
            this.Guu = Guu;
            this.Hu = Hu;
            this.Su = Su;
            this.Suu = Suu;
            this.cfs = cfs;
        }
    }

    /** Default constructor for lazy initialization. */
    public CVMGibbsModel() {
    }

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
        this.structure = structure;
        this.numComponents = elements.split("-").length;
        this.cecEntry = cecEntry;

        // --- INTERNAL PATH RESOLUTION ---
        String parentStructure = resolveParentStructure(structure);
        String disorderedFile = resolveClusterFile(parentStructure, model);
        String orderedFile = resolveClusterFile(structure, model);
        String disorderedSGName = resolveSymmetryGroup(parentStructure);
        String orderedSGName = resolveSymmetryGroup(structure);

        if (progressSink != null) {
            progressSink.accept("\n[STAGE 0]: Loading Inputs...");
            progressSink.accept(String.format("  > Elements:          %s", elements));
            progressSink.accept(String.format("  > Structure (Child): %s", structure));
            progressSink.accept(String.format("  > Structure (Parent):%s", parentStructure));
            progressSink.accept(String.format("  > Model:             %s", model));
            progressSink.accept(String.format("  > Components:        %d", numComponents));
            progressSink.accept(String.format("  > CEC Entry:         %s (%s)",
                    cecEntry != null ? cecEntry.elements : "null",
                    cecEntry != null ? cecEntry.structurePhase : "null"));

            if (cecEntry != null && cecEntry.cecTerms != null) {
                for (CECEntry.CECTerm term : cecEntry.cecTerms) {
                    progressSink.accept(String.format("    - %-10s: a = %10.6f, b = %10.6f",
                            term.name, term.a, term.b));
                }
            }
            progressSink.accept(
                    String.format("  > Files (Disord):    [clus: %s, sym: %s]", disorderedFile, disorderedSGName));
            progressSink
                    .accept(String.format("  > Files (Order):     [clus: %s, sym: %s]", orderedFile, orderedSGName));
        }

        // --- STAGE 0: LOADING ---
        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(disorderedFile);
        disorderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup disorderedSG = InputLoader.parseSpaceGroup(disorderedSGName);

        List<Cluster> orderedClusters = InputLoader.parseClusterFile(orderedFile);
        orderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup orderedSG = InputLoader.parseSpaceGroup(orderedSGName);

        // Extract transformation from the ordered phase mapping
        double[][] transformationMatrix = orderedSG.getRotateMat();
        double[] translationVector = orderedSG.getTranslateMat();

        // --- STAGE 1 & 2: IDENTIFICATION ---
        if (progressSink != null)
            progressSink.accept("\n[STAGE 1/2]: Running Identification Pipeline...");
        org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult pr = org.ce.model.cluster.ClusterCFIdentificationPipeline
                .run(
                        disorderedClusters,
                        disorderedSG.getOperations(),
                        orderedClusters,
                        orderedSG.getOperations(),
                        transformationMatrix,
                        translationVector,
                        numComponents,
                        progressSink);

        // --- STAGE 3: C-MATRIX (ORTHOGONAL) ---
        if (progressSink != null)
            progressSink.accept("\n[STAGE 3]: Running C-Matrix Pipeline...");
        org.ce.model.cluster.CMatrixPipeline.CMatrixData cmatOrth = org.ce.model.cluster.CMatrixPipeline.run(
                pr.toClusterIdentificationResult(),
                pr.toCFIdentificationResult(),
                orderedClusters,
                numComponents,
                progressSink);
        // 3b. Build and Verify Random State
        double[] x = new double[numComponents];
        java.util.Arrays.fill(x, 1.0 / numComponents);
        // double[] x = { 0.33, 0.33, 0.34 };

        System.out.println("\n[Stage 3b] Testing CMatrixPipeline.verifyRandomCVs...");
        CMatrixPipeline.verifyRandomCVs(x, pr, cmatOrth, progressSink);

        // --- STAGE 4: CVCF BASIS & FINAL DATA ---
        if (progressSink != null)
            progressSink.accept("\n[STAGE 4]: Basis Transformation...");
        org.ce.model.cvm.CvCfBasis basisRef = org.ce.model.cvm.CvCfBasis.generate(structure, pr, cmatOrth, model,
                progressSink);

        org.ce.model.cluster.CMatrixPipeline.CMatrixData cmatCvcf = basisRef.cvcfCMatrixData;

        this.mh = pr.getMh();
        this.lc = pr.getLc();

        this.cmat = cmatCvcf.getCmat();
        this.lcv = cmatCvcf.getLcv();
        this.wcv = cmatCvcf.getWcv();

        this.pipelineResult = pr;
        this.basis = basisRef;
        this.ncf = basisRef.numNonPointCfs;
        this.tcf = basisRef.totalCfs();
        this.orthCfBasisIndices = cmatCvcf.getCfBasisIndices();

        this.tcdis = pr.getTcdis();
        this.mhdis = pr.getMhdis();
        this.kb = pr.getKbdis();

        // --- Fetch and Initialize Interactions (ECIs) ---
        this.eci = CECEvaluator.evaluate(cecEntry, 300.0, basisRef, "CVM-INIT");

        if (progressSink != null && this.eci != null) {
            progressSink.accept("  > Final Mapped ECIs (Evaluated at 300K):");
            for (int i = 0; i < ncf; i++) {
                progressSink.accept(String.format("    - %-10s: %12.6f", basisRef.cfNames.get(i), eci[i]));
            }
        }

        if (progressSink != null)
            progressSink.accept("\n Intialization stage completed. ");
    }

    private String resolveParentStructure(String structure) {
        if (structure == null)
            return null;
        // Strip _CVCF if present (internal convention matching original Request logic)
        String base = structure.replace("_CVCF", "");
        if (base.equals("BCC_B2"))
            return "BCC_A2";
        if (base.equals("FCC_L12"))
            return "FCC_A1";
        return base;
    }

    private String resolveClusterFile(String structure, String model) {
        String mod = model != null ? model.replace("_CVCF", "") : "";
        return structure + "-" + mod + ".txt";
    }

    private String resolveSymmetryGroup(String structure) {
        return structure + "-SG";
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
                progressSink.accept(
                        "  [Model] Parameters updated (T=" + temperature + " K) — resetting equilibrium state.");
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

        setT(temperature);
        setX(moleFractions);

        double[] u = computeRandomCFs(moleFractions);

        double errf = 0;

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted())
                throw new CancellationException();

            setU(u);
            double[] Gu = calculateGu();
            double G = calculateG();
            double H = calculateH();
            double S = calculateS();

            errf = 0;
            for (double g : Gu)
                errf += Math.abs(g);

            if (progressSink != null) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("    iter %3d  |del_G| = %.3e  G = %11.4f  H = %11.4f  S = %9.6f",
                        its, errf, G, H, S));

                double[] cfs = calculateCfs();
                if (cfs != null && cfs.length > 0) {
                    sb.append("  CFs: [");
                    for (int i = 0; i < cfs.length; i++) {
                        if (i > 0)
                            sb.append(", ");
                        sb.append(String.format("%.4f", cfs[i]));
                    }
                    sb.append("]");
                }
                progressSink.accept(sb.toString());
            }
            if (eventSink != null)
                eventSink.accept(new ProgressEvent.CvmIteration(its, G, errf, H, S, u));

            if (errf <= tolerance) {
                ModelResult finalModelRes = new ModelResult(G, H, S, Gu, calculateGuu(), calculateHu(), calculateSu(),
                        calculateSuu(), calculateCfs());
                return new EquilibriumResult(finalModelRes, u.clone(), true, its, errf);
            }

            try {
                double[] negGu = new double[ncf];
                for (int i = 0; i < ncf; i++)
                    negGu[i] = -Gu[i];

                double[][] Guu = calculateGuu();
                double[] p = LinearAlgebra.solve(Guu, negGu);

                double alpha = calculateStepLimit(u, p, moleFractions);

                double errx = 0;
                for (int i = 0; i < ncf; i++) {
                    double delta = alpha * p[i];
                    u[i] += delta;
                    errx += Math.abs(delta);
                }

                if (errx <= TOLX) {
                    setU(u);
                    ModelResult finalModelRes = new ModelResult(calculateG(), calculateH(), calculateS(), calculateGu(),
                            calculateGuu(), calculateHu(), calculateSu(), calculateSuu(), calculateCfs());
                    return new EquilibriumResult(finalModelRes, u.clone(), true, its, errf);
                }

            } catch (Exception e) {
                ModelResult finalModelRes = new ModelResult(G, H, S, Gu, calculateGuu(), calculateHu(), calculateSu(),
                        calculateSuu(), calculateCfs());
                return new EquilibriumResult(finalModelRes, u.clone(), false, its, errf);
            }
        }

        setU(u);
        ModelResult finalModelRes = new ModelResult(calculateG(), calculateH(), calculateS(), calculateGu(),
                calculateGuu(), calculateHu(), calculateSu(), calculateSuu(), calculateCfs());
        return new EquilibriumResult(finalModelRes, u.clone(), false, MAX_ITER, errf);
    }

    // =========================================================================
    // Physics evaluation - Standalone Methods
    // =========================================================================

    /**
     * Set the current correlation functions (non-point).
     * Triggers a re-calculation of internal cluster variables.
     */
    public void setU(double[] u) {
        this.u = u.clone();
        if (this.x_mole != null) {
            syncCv();
        }
        this.isMinimized = false;
    }

    /**
     * Set the current mole fractions (composition).
     * Triggers a re-calculation of internal cluster variables.
     */
    public void setX(double[] x) {
        this.x_mole = x.clone();
        this.currentComposition = x.clone();
        if (this.u != null) {
            syncCv();
        }
        this.isMinimized = false;
    }

    /**
     * Set the current temperature.
     * Triggers a re-calculation of internal interactions (ECIs).
     */
    public void setT(double temperature) {
        this.temp = temperature;
        this.currentTemperature = temperature;
        this.eci = CECEvaluator.evaluate(cecEntry, temperature, basis, "CVM");
        this.isMinimized = false;
    }

    private void syncCv() {
        double[] uFull = CMatrixPipeline.buildFullCVCFVector(u, x_mole, ncf);
        this.currentCv = CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);
    }

    private void checkMinimized() {
        if (!isMinimized) {
            throw new IllegalStateException("CVM Model is not minimized. Please call getEquilibriumState() first.");
        }
    }

    public double calH() {
        checkMinimized();
        return calculateH();
    }

    public double[] calHu() {
        checkMinimized();
        return calculateHu();
    }

    public double[][] calHuu() {
        checkMinimized();
        return calculateHuu();
    }

    public double calS() {
        checkMinimized();
        return calculateS();
    }

    public double[] calSu() {
        checkMinimized();
        return calculateSu();
    }

    public double[][] calSuu() {
        checkMinimized();
        return calculateSuu();
    }

    public double calG() {
        checkMinimized();
        return calculateG();
    }

    public double[] calGu() {
        checkMinimized();
        return calculateGu();
    }

    public double[][] calGuu() {
        checkMinimized();
        return calculateGuu();
    }

    public double[] calCfs() {
        checkMinimized();
        return calculateCfs();
    }

    private double calculateH() {
        double Hval = 0.0;
        for (int l = 0; l < ncf; l++)
            Hval += eci[l] * u[l];
        return Hval;
    }

    private double[] calculateHu() {
        return eci.clone();
    }

    private double[][] calculateHuu() {
        return new double[ncf][ncf];
    }

    private double calculateS() {
        double Sval = 0.0;
        for (int t = 0; t < tcdis; t++) {
            double coeff_t = kb[t] * mhdis[t];
            for (int j = 0; j < lc[t]; j++) {
                double mh_tj = mh[t][j];
                int[] w = wcv.get(t).get(j);
                int nv = lcv[t][j];
                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = currentCv[t][j][incv];
                    double sContrib;
                    if (cvVal > ENTROPY_SMOOTH_EPS) {
                        sContrib = cvVal * Math.log(cvVal);
                    } else {
                        double logEps = Math.log(ENTROPY_SMOOTH_EPS);
                        double d = cvVal - ENTROPY_SMOOTH_EPS;
                        sContrib = ENTROPY_SMOOTH_EPS * logEps + (1.0 + logEps) * d + 0.5 / ENTROPY_SMOOTH_EPS * d * d;
                    }
                    Sval -= PhysicsConstants.R_GAS * coeff_t * mh_tj * w[incv] * sContrib;
                }
            }
        }
        return Sval;
    }

    private double[] calculateSu() {
        double[] Su = new double[ncf];
        for (int t = 0; t < tcdis; t++) {
            double coeff_t = kb[t] * mhdis[t];
            for (int j = 0; j < lc[t]; j++) {
                double mh_tj = mh[t][j];
                double[][] cm = cmat.get(t).get(j);
                int[] w = wcv.get(t).get(j);
                int nv = lcv[t][j];
                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = currentCv[t][j][incv];
                    double logEff;
                    if (cvVal > ENTROPY_SMOOTH_EPS) {
                        logEff = Math.log(cvVal);
                    } else {
                        double d = cvVal - ENTROPY_SMOOTH_EPS;
                        logEff = Math.log(ENTROPY_SMOOTH_EPS) + d / ENTROPY_SMOOTH_EPS;
                    }
                    double prefix = coeff_t * mh_tj * w[incv];
                    for (int l = 0; l < ncf; l++) {
                        double cml = cm[incv][l];
                        if (cml != 0.0)
                            Su[l] -= PhysicsConstants.R_GAS * prefix * cml * logEff;
                    }
                }
            }
        }
        return Su;
    }

    private double[][] calculateSuu() {
        double[][] Suu = new double[ncf][ncf];
        for (int t = 0; t < tcdis; t++) {
            double coeff_t = kb[t] * mhdis[t];
            for (int j = 0; j < lc[t]; j++) {
                double mh_tj = mh[t][j];
                double[][] cm = cmat.get(t).get(j);
                int[] w = wcv.get(t).get(j);
                int nv = lcv[t][j];
                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = currentCv[t][j][incv];
                    double invEff;
                    if (cvVal > ENTROPY_SMOOTH_EPS) {
                        invEff = 1.0 / cvVal;
                    } else {
                        invEff = 1.0 / ENTROPY_SMOOTH_EPS;
                    }
                    double prefix = coeff_t * mh_tj * w[incv];
                    for (int l1 = 0; l1 < ncf; l1++) {
                        double cml1 = cm[incv][l1];
                        if (cml1 == 0.0)
                            continue;
                        for (int l2 = l1; l2 < ncf; l2++) {
                            double cml2 = cm[incv][l2];
                            if (cml2 == 0.0)
                                continue;
                            double val = -PhysicsConstants.R_GAS * prefix * cml1 * cml2 * invEff;
                            Suu[l1][l2] += val;
                            if (l1 != l2)
                                Suu[l2][l1] += val;
                        }
                    }
                }
            }
        }
        return Suu;
    }

    private double calculateG() {
        return calculateH() - temp * calculateS();
    }

    private double[] calculateGu() {
        double[] Hu = calculateHu();
        double[] Su = calculateSu();
        double[] Gu = new double[ncf];
        for (int i = 0; i < ncf; i++) {
            Gu[i] = Hu[i] - temp * Su[i];
        }
        return Gu;
    }

    private double[][] calculateGuu() {
        double[][] Suu = calculateSuu();
        double[][] Guu = new double[ncf][ncf];
        for (int i = 0; i < ncf; i++) {
            for (int j = 0; j < ncf; j++) {
                Guu[i][j] = -temp * Suu[i][j];
            }
        }
        return Guu;
    }

    private double[] calculateCfs() {
        return CMatrixPipeline.buildFullCVCFVector(u, x_mole, ncf);
    }

    // =========================================================================
    // Physics evaluation - Core Newton-Raphson
    // =========================================================================

    /** Evaluates the CVM Gibbs free energy and derivatives at the given state. */
    public ModelResult evaluate(double[] u, double[] moleFractions, double temperature) {
        setT(temperature);
        setX(moleFractions);
        setU(u);

        return new ModelResult(
                calculateG(),
                calculateH(),
                calculateS(),
                calculateGu(),
                calculateGuu(),
                calculateHu(),
                calculateSu(),
                calculateSuu(),
                calculateCfs());
    }

    /** Static overload for callers that supply all parameters explicitly. */
    public static ModelResult evaluate(
            double[] u, double[] moleFractions, double temperature, CECEntry cecEntry, CvCfBasis basis,
            int tcdis, int ncf, double[] mhdis, double[] kb, double[][] mh, int[] lc,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv) {
        // NOTE: Static method doesn't have an instance to store state.
        // We could create a temporary instance or keep the legacy logic.
        // For simplicity and to avoid duplication, we'll create a temp instance.
        CVMGibbsModel tempModel = new CVMGibbsModel();
        // Manual initialization of required fields for temp model
        tempModel.cecEntry = cecEntry;
        tempModel.basis = basis;
        tempModel.tcdis = tcdis;
        tempModel.ncf = ncf;
        tempModel.mhdis = mhdis;
        tempModel.kb = kb;
        tempModel.mh = mh;
        tempModel.lc = lc;
        tempModel.cmat = cmat;
        tempModel.lcv = lcv;
        tempModel.wcv = wcv;

        return tempModel.evaluate(u, moleFractions, temperature);
    }
    // =========================================================================
    // Helpers used by the N-R loop
    // =========================================================================

    /**
     * Computes the full disordered-state (random) CVCF vector for N-R
     * initialisation.
     */
    public double[] computeRandomCFs(double[] moleFractions) {
        return basis.computeRandomCvcfCFs(moleFractions, pipelineResult);
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
        double[] uFull = CMatrixPipeline.buildFullCVCFVector(u, moleFractions, ncf);
        return CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    public int getNumComponents() {
        return numComponents;
    }

    public String getElements() {
        return elements;
    }

    public CvCfBasis getBasis() {
        return basis;
    }

    public int getNcf() {
        return ncf;
    }

    public int getTcf() {
        return tcf;
    }

    public int getTcdis() {
        return tcdis;
    }

    public int[] getLc() {
        return lc;
    }

    public int[][] getLcv() {
        return lcv;
    }

    public int[][] getOrthCfBasisIndices() {
        return orthCfBasisIndices;
    }

    /**
     * Computes cluster variables cv[t][j][v] from the given non-point CFs and
     * composition.
     */
    public double[][][] evaluateClusterVariables(double[] u, double[] moleFractions) {
        double[] uFull = CMatrixPipeline.buildFullCVCFVector(u, moleFractions, ncf);
        return CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);
    }
}
