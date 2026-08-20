package org.ce.scratch;

import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cvm.CvCfBasis;
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
/**
 * Frozen copy of {@code CVMGibbsModel} as it stood at commit 18b965d, before
 * it became a facade over {@code CvmEvaluator}.
 *
 * <p>Kept solely so {@link CvmEvaluatorParity} retains an independent second
 * implementation of the free-energy expressions. Once the production class
 * delegates to the evaluator, comparing against it would compare the evaluator
 * to itself and pass vacuously -- the same reason
 * {@code HardcodedLatticeStability} was frozen when {@code LatticeStability}
 * became a facade over {@code SgteDatabase}.</p>
 *
 * <p><b>Not for production use, and not to be "fixed".</b> It is a historical
 * artifact carrying the pre-refactor arithmetic verbatim.</p>
 */
public class PreFacadeCVMGibbsModel {

    private static final double ENTROPY_SMOOTH_EPS = 1.0e-6;
    private static final int MAX_ITER = 100;
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
    public PreFacadeCVMGibbsModel() {
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
            progressSink.accept("\n  [CVM-Setup] Stage 0: Loading Inputs...");
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
            progressSink.accept("  [CVM-Setup] Stage 1: Identification...");
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
            progressSink.accept("  [CVM-Setup] Stage 2: C-Matrix Generation...");
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
            progressSink.accept("  [CVM-Setup] Stage 3: Basis Transformation...");
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


        if (progressSink != null)
            progressSink.accept("  [CVM-Setup] ✓ Initialization complete.");
    }

    private String resolveParentStructure(String structure) {
        return org.ce.model.cluster.StructurePhaseRegistry.parentOf(structure);
    }

    private String resolveClusterFile(String structure, String model) {
        String mod = model != null ? model.replace("_CVCF", "") : "";
        return "clus/" + structure + "-" + mod + ".txt";
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
            Consumer<String> progressSink, Consumer<ProgressEvent> eventSink,
            org.ce.calculation.CalculationDescriptor.Property required) {

        boolean compositionChanged = currentComposition == null
                || !Arrays.equals(currentComposition, composition);
        boolean temperatureChanged = Math.abs(currentTemperature - temperature) > 1.0e-5;

        if (temperatureChanged || compositionChanged) {
            if (progressSink != null) {
                progressSink.accept(String.format(
                        "\n  [Model] Parameters updated: T = %.1f K, x = %s",
                        temperature, Arrays.toString(composition)));
            }
            isMinimized = false;
            currentTemperature = temperature;
            currentComposition = Arrays.copyOf(composition, composition.length);
        } else if (isMinimized) {
            if (progressSink != null)
                progressSink.accept("  [Model] Reusing cached equilibrium state for these parameters.");
        }

        if (!isMinimized || lastResult == null) {
            if (progressSink != null)
                progressSink.accept("  [Model] Initiating internal minimization (Newton-Raphson loop)...");
            lastResult = minimize(composition, temperature, tolerance, progressSink, eventSink, required);
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
     *
     * <p>Loop structure follows the original (pre-refactor) reference
     * solver, {@code CVMBINCE.minimize()} (which solves in the *orthogonal*
     * CF basis): an early exit the moment any cluster variable is
     * non-positive at the *current* point (not just a step-limit check on
     * the trial point; see {@link #minClusterVariable}), a gradient-norm
     * convergence check, then a Newton step whose size is limited by
     * {@link #calculateStepLimit} (port of the reference's {@code stpmx})
     * so no cluster variable leaves {@code [0, 1]}, followed by a
     * raw-Newton-step-size convergence check ({@code errx}). The reference's
     * additional pre-clamp on the trial correlation functions themselves
     * ({@code Utils.normalU}, restricting them to a fixed {@code [-1, 1]})
     * was tried here too but deliberately dropped: that bound is meaningful
     * for the orthogonal basis's Chebyshev-like CFs, but this solver works
     * in the CVCF basis, whose {@code u} components have a different,
     * non-uniform natural range (observed roughly 1e-4 to 0.2 in practice) --
     * a blind &plusmn;1 clamp there is a no-op, not a safeguard, and did not
     * change behavior when tested. The near-edge convergence stall this
     * porting effort was investigating (a dilute-composition Newton
     * direction that repeatedly re-approaches, but never crosses, a
     * cluster-variable boundary) remains open; see CLAUDE.md's note on
     * near-edge ternary solver fragility.</p>
     */
    private EquilibriumResult minimize(
            double[] moleFractions, double temperature, double tolerance,
            Consumer<String> progressSink, Consumer<ProgressEvent> eventSink,
            org.ce.calculation.CalculationDescriptor.Property required) {

        if (eventSink != null)
            eventSink.accept(new ProgressEvent.EngineStart("CVM", 0));

        setT(temperature, progressSink);
        setX(moleFractions);

        double[] u = computeRandomCFs(moleFractions);

        double errf = 0;

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted())
                throw new CancellationException();

            setU(u);

            // Early exit if the *current* point already has a non-positive
            // cluster variable -- mirrors the reference solver's cvMin<=0
            // check. This is NOT a failure: a cluster variable sitting at or
            // below zero at the current point (e.g. a configuration that is
            // genuinely disallowed at this composition/order) means there is
            // nowhere further for the Newton step to usefully go without
            // dividing by a near-zero probability, so the reference solver
            // accepts the current point as converged rather than attempting
            // a step. Runs before the gradient/Newton-step machinery each
            // iteration, not just as a trial-step limiter.
            double minCv = minClusterVariable();
            if (minCv <= 0) {
                double G0 = calculateGm(), H0 = calculateHm(), S0 = calculateSm();
                double[] Gu0 = calculateGmu();
                double errf0 = 0;
                for (double g : Gu0) errf0 += Math.abs(g);
                ModelResult finalModelRes = new ModelResult(G0, H0, S0, Gu0, calculateGmuu(), calculateHmu(),
                        calculateSmu(), calculateSmuu(), calculateCfs());
                return new EquilibriumResult(finalModelRes, u.clone(), true, its, errf0);
            }

            double[] Gu = calculateGmu();
            double G = calculateGm();
            double H = calculateHm();
            double S = calculateSm();

            errf = 0;
            for (double g : Gu)
                errf += Math.abs(g);

            if (eventSink != null)
                eventSink.accept(new ProgressEvent.CvmIteration(its, G, errf, H, S, u));

            if (errf <= tolerance) {
                ModelResult finalModelRes = new ModelResult(G, H, S, Gu, calculateGmuu(), calculateHmu(), calculateSmu(),
                        calculateSmuu(), calculateCfs());
                return new EquilibriumResult(finalModelRes, u.clone(), true, its, errf);
            }

            try {
                double[] negGu = new double[ncf];
                for (int i = 0; i < ncf; i++)
                    negGu[i] = -Gu[i];

                double[][] Guu = calculateGmuu();
                double[] p = LinearAlgebra.solve(Guu, negGu);

                double errx = 0;
                for (double v : p) errx += Math.abs(v);

                // The reference solver's Utils.normalU pre-clamps the trial
                // point to a fixed [-1,1] range before the cluster-variable
                // check below -- but that bound is meaningful for the
                // *orthogonal* CF basis that reference solves in (Chebyshev
                // -like polynomials, genuinely bounded by construction). This
                // solver works in the CVCF basis, whose u components have a
                // different, non-uniform natural range (observed ~1e-4 to
                // ~0.2 in practice) -- a blind +-1 clamp there is a no-op, not
                // a safeguard. So this only applies the cluster-variable-
                // space clamp (stpmx), directly on the unclamped u+p trial
                // point, same as before the reference-port attempt.
                double[] uTrial = new double[ncf];
                for (int i = 0; i < ncf; i++) uTrial[i] = u[i] + p[i];
                double alpha = calculateStepLimit(u, uTrial, moleFractions);

                for (int i = 0; i < ncf; i++) {
                    double delta = alpha * p[i];
                    u[i] += delta;
                }

                // errx small means the raw (unclamped) Newton step p was
                // already tiny -- a genuine sign of convergence, unlike a
                // small step caused purely by the boundary clamp above.
                if (errx <= TOLX) {
                    setU(u);
                    ModelResult finalModelRes = new ModelResult(calculateGm(), calculateHm(), calculateSm(), calculateGmu(),
                            calculateGmuu(), calculateHmu(), calculateSmu(), calculateSmuu(), calculateCfs());
                    return new EquilibriumResult(finalModelRes, u.clone(), true, its, errf);
                }

            } catch (Exception e) {
                ModelResult finalModelRes = new ModelResult(G, H, S, Gu, calculateGmuu(), calculateHmu(), calculateSmu(),
                        calculateSmuu(), calculateCfs());
                return new EquilibriumResult(finalModelRes, u.clone(), false, its, errf);
            }
        }

        setU(u);
        ModelResult finalModelRes = new ModelResult(calculateGm(), calculateHm(), calculateSm(), calculateGmu(),
                calculateGmuu(), calculateHmu(), calculateSmu(), calculateSmuu(), calculateCfs());
        return new EquilibriumResult(finalModelRes, u.clone(), false, MAX_ITER, errf);
    }

    /**
     * Minimum cluster variable value across all non-point cluster types, at
     * the current (u, x) state -- port of the reference solver's
     * {@code findMin}, which explicitly excludes the last (point) cluster
     * type ({@code tcdis-1}); the point type holds the mole fractions
     * themselves, not derived cluster probabilities, and including it here
     * would make this check spuriously sensitive to composition rather than
     * to how far the *solve* is from a degenerate cluster configuration.
     *
     * <p>Deliberately kept separate from the public, stateless {@link
     * #isValidParams}: this method reads the instance's already-set {@code
     * currentCv} (only valid after {@link #setU}/{@link #setX}) and is used
     * exclusively by {@code minimize()}'s own per-iteration early-exit
     * check; {@code isValidParams} evaluates an arbitrary trial point
     * without touching instance state, for external callers. Both check the
     * same physical condition over the same cluster types, just via
     * different data sources for different callers.</p>
     */
    private double minClusterVariable() {
        double minCv = Double.POSITIVE_INFINITY;
        for (int t = 0; t < tcdis - 1; t++) {
            double[][] tt = currentCv[t];
            if (tt == null) continue;
            for (double[] jj : tt) {
                if (jj == null) continue;
                for (double v : jj) minCv = Math.min(minCv, v);
            }
        }
        return minCv;
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
        setT(temperature, null);
    }

    public void setT(double temperature, java.util.function.Consumer<String> sink) {
        this.temp = temperature;
        this.currentTemperature = temperature;
        this.eci = CECEvaluator.evaluate(cecEntry, temperature, basis, "CVM", sink);
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

    public double calHm() {
        checkMinimized();
        return calculateHm();
    }

    public double[] calHmu() {
        checkMinimized();
        return calculateHmu();
    }

    public double[][] calHuu() {
        checkMinimized();
        return calculateHuu();
    }

    public double calSm() {
        checkMinimized();
        return calculateSm();
    }

    public double[] calSmu() {
        checkMinimized();
        return calculateSmu();
    }

    public double[][] calSmuu() {
        checkMinimized();
        return calculateSmuu();
    }

    public double calGm() {
        checkMinimized();
        return calculateGm();
    }

    public double[] calGmu() {
        checkMinimized();
        return calculateGmu();
    }

    public double[][] calGmuu() {
        checkMinimized();
        return calculateGmuu();
    }

    // =========================================================================
    // Reference energy (G0m) and the absolute total: G = G0m + Gm
    //
    //   G0m  reference energy of the mechanical mixture of pure elements,
    //        Sum_i x_i * G0(element_i, phase, T). Pure energy: linear in
    //        composition, independent of the CVCF variables u, and carrying
    //        no configurational entropy of its own.
    //
    //   Gm   the CVM mixing contribution -- the ECI energy Hm together with
    //        the configurational entropy of mixing Sm that the cluster
    //        variation method computes. This is what the Newton-Raphson loop
    //        minimises and what CLAUDE.md's documented verification values
    //        (e.g. -3480.5209063901 for Nb-Ti) are anchored to.
    //
    //   G    the absolute Gibbs energy, G = G0m + Gm.
    //
    // Because G0m depends only on (x, T) and not on u, every u-derivative of
    // the absolute quantity equals the mixing one exactly:
    //   dG0m/du = 0, d2G0m/du2 = 0, so calGmu/calGmuu serve both.
    // Only a widened gradient over uFull = [u ; x] differs, and only in its
    // trailing composition block, where dG0m/dx_i = G0(element_i, phase, T).
    // =========================================================================

    /**
     * Reference energy of the mechanical mixture of pure elements,
     * {@code G0m = Sum_i x_i * G0(element_i, phase, T)}, via
     * {@link org.ce.model.equilibrium.LatticeStability#g0m}.
     *
     * <p>Carries no configurational entropy, so {@code H0m = G0m} and
     * {@code S0m = 0}.</p>
     */
    public double calG0m() {
        checkMinimized();
        return calculateG0m();
    }

    /** {@code H0m = G0m} -- the reference term is pure energy. */
    public double calH0m() {
        checkMinimized();
        return calculateG0m();
    }

    /** {@code S0m = 0} -- a mechanical mixture of pure elements has no entropy of mixing. */
    public double calS0m() {
        checkMinimized();
        return 0.0;
    }

    /** Absolute Gibbs energy {@code G = G0m + Gm}. */
    public double calG() {
        checkMinimized();
        return calculateG0m() + calculateGm();
    }

    /** Absolute enthalpy {@code H = H0m + Hm = G0m + Hm}. */
    public double calH() {
        checkMinimized();
        return calculateG0m() + calculateHm();
    }

    /** Absolute entropy {@code S = S0m + Sm = Sm}, since {@code S0m = 0}. */
    public double calS() {
        checkMinimized();
        return calculateSm();
    }

    private double calculateG0m() {
        return org.ce.model.equilibrium.LatticeStability.g0m(
                java.util.List.of(elements.split("-")), structure, x_mole, temp);
    }

    /**
     * Widened gradient of G with respect to the full internal-parameter
     * vector {@code uFull = [u (non-point CFs) ; x (mole fractions)]},
     * length {@code ncf + K} -- for the Hillert multi-phase equilibrium
     * solver only (HILLERT_SOLVER_PLAN.md).
     *
     * <p><b>Deliberately separate from {@link #calGmu}/{@link #calGmuu} and
     * the private {@code calculateGmu}/{@code calculateGmuu} the Newton-Raphson
     * loop in {@code minimize()} actually solves against.</b> Those must stay
     * exactly {@code ncf}-length: {@code minimize()} builds
     * {@code negGu}/{@code Guu} sized {@code ncf} and hands them to
     * {@code LinearAlgebra.solve}, so widening them in place would silently
     * change the single-phase solver to solve a bigger, different linear
     * system every iteration -- not just report a bigger result. This method
     * (and {@link #calGuuFull}) exist alongside the untouched originals
     * specifically so the single-phase solve path is unaffected.</p>
     *
     * <p>The trailing {@code K} entries are the composition-block gradient
     * (chemical-potential-like quantities, {@code Gx} in the ported
     * Mathematica reference {@code delxGCVM}); the leading {@code ncf}
     * entries equal {@link #calGmu} exactly, since both are computed the same
     * way over the same {@code cm} columns, just over a different column
     * range.</p>
     */
    public double[] calGuFull() {
        checkMinimized();
        return calculateGuFull();
    }

    /**
     * Widened Hessian of G over the full {@code uFull} space, {@code (ncf+K)
     * x (ncf+K)} -- see {@link #calGuFull} for why this is separate from
     * {@link #calGmuu}/the Newton-Raphson loop's own {@code calculateGmuu}.
     */
    public double[][] calGuuFull() {
        checkMinimized();
        return calculateGuuFull();
    }

    public double[] calCfs() {
        checkMinimized();
        return calculateCfs();
    }

    private double calculateHm() {
        double Hval = 0.0;
        for (int l = 0; l < ncf; l++)
            Hval += eci[l] * u[l];
        return Hval;
    }

    private double[] calculateHmu() {
        return eci.clone();
    }

    private double[][] calculateHuu() {
        return new double[ncf][ncf];
    }

    private double calculateSm() {
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

    private double[] calculateSmu() {
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

    private double[][] calculateSmuu() {
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

    /**
     * Widened entropy gradient over {@code uFull = [u ; x]}, length
     * {@code ncf + numComponents} -- for {@link #calGuFull} (Hillert solver
     * only). Independent copy of {@link #calculateSmu}'s loop body, widened
     * to walk {@code cm}'s trailing composition columns too, rather than a
     * shared parameterized helper: {@code calculateSmu}/{@code calculateSmuu}
     * are the Newton-Raphson loop's own hot path (the target of this
     * session's convergence-criterion and step-damping fixes), and keeping
     * this method fully separate means nothing here can accidentally change
     * their behavior.
     */
    private double[] calculateSuFull() {
        int width = ncf + numComponents;
        double[] Su = new double[width];
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
                    for (int l = 0; l < width; l++) {
                        double cml = cm[incv][l];
                        if (cml != 0.0)
                            Su[l] -= PhysicsConstants.R_GAS * prefix * cml * logEff;
                    }
                }
            }
        }
        return Su;
    }

    /** Widened entropy Hessian over {@code uFull}, {@code (ncf+K) x (ncf+K)} -- see {@link #calculateSuFull}. */
    private double[][] calculateSuuFull() {
        int width = ncf + numComponents;
        double[][] Suu = new double[width][width];
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
                    for (int l1 = 0; l1 < width; l1++) {
                        double cml1 = cm[incv][l1];
                        if (cml1 == 0.0)
                            continue;
                        for (int l2 = l1; l2 < width; l2++) {
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

    /**
     * Widened gradient of G (= H - T*S) over {@code uFull = [u ; x]}. The
     * leading {@code ncf} entries equal {@link #calculateGmu} exactly (H is
     * linear in {@code u} only, so {@code Hu}'s contribution is identical
     * either way); the trailing {@code numComponents} entries are zero from
     * H (H has no x-dependence at all -- {@code eci[l]*u[l]} sums only over
     * {@code u}) and whatever {@code calculateSuFull}'s trailing block gives
     * from S.
     */
    private double[] calculateGuFull() {
        int width = ncf + numComponents;
        double[] SuFull = calculateSuFull();
        double[] Gu = new double[width];
        for (int i = 0; i < width; i++) {
            double hu = (i < ncf) ? eci[i] : 0.0; // H has no x-dependence
            Gu[i] = hu - temp * SuFull[i];
        }
        return Gu;
    }

    /** Widened Hessian of G over {@code uFull} -- see {@link #calculateGuFull}. H contributes zero everywhere (Huu = 0). */
    private double[][] calculateGuuFull() {
        int width = ncf + numComponents;
        double[][] SuuFull = calculateSuuFull();
        double[][] Guu = new double[width][width];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < width; j++) {
                Guu[i][j] = -temp * SuuFull[i][j];
            }
        }
        return Guu;
    }

    private double calculateGm() {
        return calculateHm() - temp * calculateSm();
    }

    private double[] calculateGmu() {
        double[] Hu = calculateHmu();
        double[] Su = calculateSmu();
        double[] Gu = new double[ncf];
        for (int i = 0; i < ncf; i++) {
            Gu[i] = Hu[i] - temp * Su[i];
        }
        return Gu;
    }

    private double[][] calculateGmuu() {
        double[][] Suu = calculateSmuu();
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
                calculateGm(),
                calculateHm(),
                calculateSm(),
                calculateGmu(),
                calculateGmuu(),
                calculateHmu(),
                calculateSmu(),
                calculateSmuu(),
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
        PreFacadeCVMGibbsModel tempModel = new PreFacadeCVMGibbsModel();
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

    /**
     * Result of {@link #solvePerPhaseStep}: the joint Newton step
     * {@code deltaY(mu)}, expressed as an <b>affine function of the trial
     * chemical-potential vector {@code mu}</b> rather than a value at one
     * fixed {@code mu} -- see the class-level note on {@link
     * #solvePerPhaseStep} for why this shape, not a single numeric result,
     * is what the outer Hillert solver actually needs.
     *
     * <p>{@code deltaY(mu) = deltaY0 + Σ_k mu[k]*deltaYSensitivity[k]}, and
     * likewise for {@code deltaComposition}/{@code lambda}. Evaluate at a
     * specific {@code mu} via {@link #evaluateAt}.</p>
     */
    public record PerPhaseStepResult(
            double[] deltaY0, double[][] deltaYSensitivity,
            double[] deltaComposition0, double[][] deltaCompositionSensitivity,
            double lambda0, double[] lambdaSensitivity) {

        /** Evaluates this affine result at a specific numeric {@code mu} -- convenience for testing/inspection. */
        public double[] deltaCompositionAt(double[] mu) {
            double[] result = deltaComposition0.clone();
            for (int k = 0; k < mu.length; k++) {
                for (int i = 0; i < result.length; i++) {
                    result[i] += mu[k] * deltaCompositionSensitivity[k][i];
                }
            }
            return result;
        }

        /** Evaluates the full joint deltaY (length ncf+K) at a specific numeric {@code mu} -- used by HillertSolver's inner loop to update uFull. */
        public double[] deltaYAt(double[] mu) {
            double[] result = deltaY0.clone();
            for (int k = 0; k < mu.length; k++) {
                for (int i = 0; i < result.length; i++) {
                    result[i] += mu[k] * deltaYSensitivity[k][i];
                }
            }
            return result;
        }
    }

    /**
     * One Hillert multi-phase equilibrium Newton step for this phase,
     * expressed as an <b>affine function of the trial chemical-potential
     * vector {@code mu}</b> -- port of the reference Mathematica
     * implementation's {@code delxGCVM} (specifically its {@code phaseEq}/
     * {@code ca}/{@code sol} linear-system assembly and solve).
     *
     * <p><b>Why affine-in-mu, not a fixed-mu numeric result (revised from
     * an earlier version of this method):</b> tracing {@code phaseq}'s
     * actual outer loop shows {@code delxGCVM} is called with {@code mu}
     * still <em>symbolic</em> (Mathematica's {@code Solve}/{@code
     * CoefficientArrays} return {@code delyN}/{@code delnN} as literal
     * formulas in {@code mu} -- confirmed directly against a real
     * {@code delxGCVM} output for Nb-Ti-V, whose printed result was e.g.
     * {@code 0.0000430165*(988.474+ΔT+0.618187*mu1-0.30448*mu2-...)}).
     * {@code genEqMat}'s outer mass-balance equations then substitute that
     * symbolic {@code delnN} in directly, so {@code mu} and {@code deltaN}
     * are solved <em>simultaneously</em> in one combined system -- the
     * per-phase step is never evaluated at a numeric {@code mu} on its own.
     * Since {@code deltaY} is provably affine in {@code mu} (the linear
     * system's matrix {@code A} doesn't depend on {@code mu}; only the
     * right-hand side does, and only in the x-block rows), the numeric
     * equivalent of Mathematica's symbolic substitution is: solve the same
     * system {@code K+1} times against basis right-hand-sides (once for
     * {@code mu=0}, once per unit vector {@code e_k}) to get the affine
     * coefficients directly, then let {@code HillertSolver.EquilibriumMatrix} fold that
     * affine form into its own equations before solving for {@code mu} --
     * exactly mirroring {@code genEqMat}'s substitution, entirely
     * numerically, with no symbolic algebra.</p>
     *
     * <p>Deliberately implemented on {@code PreFacadeCVMGibbsModel} rather than as a
     * separate class (HILLERT_SOLVER_PLAN.md): its only real inputs are
     * this class's own widened {@code Gu}/{@code Guu}, so a separate class
     * would mostly be a thin wrapper reaching into already-available state.
     * This mirrors {@code PreFacadeCVMGibbsModel} already being both an evaluator
     * and an optimizer (see {@code minimize()}). It is unrelated to
     * {@code getEquilibriumState}/{@code minimize()}'s single-phase
     * Newton-Raphson loop (fixed composition, stationary {@code G}) and
     * must never be confused with it: this solves for a stationary point of
     * {@code G} <em>relative to a trial {@code mu}</em>, with composition
     * itself among the unknowns.</p>
     *
     * <p><b>The linear system</b> (derived from {@code delxGCVM}'s {@code
     * phaseEq}, at fixed T/P so the {@code GxT*ΔT}/{@code GxP*ΔP} terms
     * vanish): unknowns are {@code deltaY[0..ncf+K-1]} (the joint step) and
     * {@code lambda} (Lagrange multiplier), {@code ncf+K+1} equations, for
     * a right-hand side {@code b(mu)}:</p>
     * <ul>
     *   <li>Rows {@code 0..ncf-1} (u-block): {@code Guu[i,:] . deltaY = -Gu[i]}
     *       -- the ordinary stationarity condition on the internal CFs,
     *       unconstrained by {@code mu}.</li>
     *   <li>Rows {@code ncf..ncf+K-1} (x-block): {@code Guu[i,:] . deltaY - lambda
     *       = mu[i-ncf] - Gu[i]} -- the only rows where {@code mu} appears,
     *       always with coefficient exactly {@code +1} on its own row --
     *       which is exactly why {@code deltaY} is affine in {@code mu} via
     *       one basis-vector solve per component.</li>
     *   <li>Row {@code ncf+K} (constraint): {@code sum(deltaY[ncf..ncf+K-1]) = 0}
     *       -- composition change stays on the simplex, independent of
     *       {@code mu} (the constraint row itself has no {@code mu} term).</li>
     * </ul>
     *
     * <p>Built entirely from analytic {@link #calculateGuFull}/{@link
     * #calculateGuuFull} -- no finite-differencing anywhere, per
     * HILLERT_SOLVER_PLAN.md's standing directive.</p>
     *
     * @param uFull current joint state {@code [u ; x]}, length {@code ncf+K}
     * @param temperature current temperature, K
     */
    public PerPhaseStepResult solvePerPhaseStep(double[] uFull, double temperature) {
        int width = ncf + numComponents;
        if (uFull.length != width) {
            throw new IllegalArgumentException(
                    "uFull.length=" + uFull.length + " != ncf+K=" + width);
        }

        double[] uBlock = new double[ncf];
        double[] xBlock = new double[numComponents];
        System.arraycopy(uFull, 0, uBlock, 0, ncf);
        System.arraycopy(uFull, ncf, xBlock, 0, numComponents);

        setT(temperature);
        setX(xBlock);
        setU(uBlock);

        double[] Gu = calculateGuFull();
        double[][] Guu = calculateGuuFull();

        int n = width + 1; // + lambda

        // Matrix A is the same for every right-hand side (mu doesn't appear
        // in it -- see class doc) -- build once, reuse for all K+1 solves.
        double[][] A = new double[n][n];
        for (int i = 0; i < width; i++) {
            System.arraycopy(Guu[i], 0, A[i], 0, width);
        }
        for (int i = ncf; i < width; i++) {
            A[i][width] = -1.0; // -lambda
        }
        for (int i = ncf; i < width; i++) {
            A[width][i] = 1.0; // sum(deltaX) = 0
        }

        // b0: the mu=0 right-hand side.
        double[] b0 = new double[n];
        for (int i = 0; i < ncf; i++) b0[i] = -Gu[i];
        for (int i = ncf; i < width; i++) b0[i] = -Gu[i];
        double[] sol0 = LinearAlgebra.solve(A, b0);

        // b_k: the k-th unit-mu right-hand side (mu[k]=1, all others 0) minus
        // b0 gives the pure sensitivity -- equivalently, solve with only the
        // mu[k] term present (all Gu terms zeroed) since the system is linear
        // and A is shared: deltaY(mu) - deltaY(0) is linear in mu, so solving
        // A*z = e_{ncf+k} directly gives d(deltaY)/d(mu_k).
        double[][] deltaYSens = new double[numComponents][];
        double[] lambdaSens = new double[numComponents];
        double[][] deltaCompSens = new double[numComponents][];
        for (int k = 0; k < numComponents; k++) {
            double[] ek = new double[n];
            ek[ncf + k] = 1.0;
            double[] solK = LinearAlgebra.solve(A, ek);
            double[] deltaYk = new double[width];
            System.arraycopy(solK, 0, deltaYk, 0, width);
            deltaYSens[k] = deltaYk;
            lambdaSens[k] = solK[width];
            double[] deltaCompK = new double[numComponents];
            System.arraycopy(deltaYk, ncf, deltaCompK, 0, numComponents);
            deltaCompSens[k] = deltaCompK;
        }

        double[] deltaY0 = new double[width];
        System.arraycopy(sol0, 0, deltaY0, 0, width);
        double[] deltaComposition0 = new double[numComponents];
        System.arraycopy(deltaY0, ncf, deltaComposition0, 0, numComponents);
        double lambda0 = sol0[width];

        return new PerPhaseStepResult(deltaY0, deltaYSens, deltaComposition0, deltaCompSens, lambda0, lambdaSens);
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
     * Finds the largest fraction α &isin; (0, 1] of the step from
     * {@code uOld} to {@code uTrial} that keeps every cluster variable
     * within [0, 1] -- port of the reference solver's {@code stpmx(uold,
     * unew)}. Unlike the old (pre-port) version of this method, it does not
     * compute the trial point itself: the caller is expected to hand in an
     * already-normalized trial point (see {@link #normalizeTrialU}), so this
     * method's only job is the cluster-variable-space clamp, matching the
     * reference's two-stage structure.
     */
    public double calculateStepLimit(double[] uOld, double[] uTrial, double[] moleFractions) {
        double fmin = 1.0;

        double[][][] cv_old = updateCVInternal(uOld, moleFractions);
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
        // Reference solver's stpmx: full step (alpha=1) if no cluster
        // variable would leave [0,1]; otherwise 0.1*fmin, backing well off
        // the boundary rather than landing exactly on it.
        return (fmin >= 1.0) ? 1.0 : (0.1 * fmin);
    }

    private double[][][] updateCVInternal(double[] u, double[] moleFractions) {
        double[] uFull = CMatrixPipeline.buildFullCVCFVector(u, moleFractions, ncf);
        return CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);
    }

    /**
     * True if every cluster variable at the given {@code (u, x)} point --
     * across <em>all</em> cluster types, including the point (composition)
     * block -- lies strictly inside {@code (0, 1)}. Public and side-effect-
     * free (does not call {@link #setU}/{@link #setX}, does not touch
     * {@code currentCv}) so external callers -- notably the Hillert
     * multi-phase solver (HILLERT_SOLVER_PLAN.md), which must validate a
     * trial joint step across several {@code PreFacadeCVMGibbsModel} instances
     * before committing to it -- can check a candidate point without
     * perturbing this instance's single-phase minimization state.
     *
     * <p>Port of the reference Hillert implementation's {@code
     * isValidParams}, confirmed against its actual definition: it checks
     * {@code cvt}, which the underlying {@code calGmcecvm} builds by
     * looping cluster type {@code i} over the <em>full</em> {@code 1..tcdis}
     * range (Mathematica 1-indexed) -- unlike {@link #minClusterVariable}'s
     * {@code findMin}, which explicitly excludes the last (point) type via
     * {@code tcdis-1}. So this method is deliberately <em>broader</em> than
     * {@code minClusterVariable}, not a stateless twin of it: it also
     * catches a trial point whose composition itself has drifted to a
     * pure-element boundary ({@code x_i = 0} or {@code 1}), which
     * {@code minClusterVariable}'s exclusion would miss. For a real mixture
     * with more than one component present this is normally trivially
     * satisfied, but it is not guaranteed -- e.g. exactly the kind of
     * boundary the Hillert outer/inner Newton step's backtracking is meant
     * to catch when a trial phase composition wanders to a corner. Do not
     * "simplify" this back to the {@code tcdis-1} exclusion; that would
     * silently narrow what this check catches.</p>
     *
     * <p>Evaluated the same way {@link #calculateStepLimit} already
     * evaluates a trial point (via {@link #updateCVInternal}), just
     * returning a boolean instead of a step-limit fraction. Strict
     * inequality matches the reference's {@code 0.0 < # < 1.0} predicate
     * exactly.</p>
     */
    public boolean isValidParams(double[] u, double[] moleFractions) {
        double[][][] cv = updateCVInternal(u, moleFractions);
        for (int t = 0; t < tcdis; t++) {
            double[][] tt = cv[t];
            if (tt == null) continue;
            for (double[] jj : tt) {
                if (jj == null) continue;
                for (double v : jj) {
                    if (v <= 0.0 || v >= 1.0) return false;
                }
            }
        }
        return true;
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
