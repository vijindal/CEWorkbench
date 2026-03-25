package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.CMatrixResult;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.ClusterVariableEvaluator;
import org.ce.domain.cluster.cvcf.CvCfIntegration;
import org.ce.domain.engine.cvm.NewtonRaphsonSolverSimple.CVMSolverResult;
import org.ce.domain.result.EquilibriumState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Central thermodynamic model for CVM free-energy calculations.
 *
 * <p><b>User's Mental Model:</b>
 * Create a CVM phase model for a given system, provide system parameters (CECs)
 * and macro parameters (T, composition), and query for equilibrium thermodynamic
 * quantities. When parameters change, the model automatically re-minimizes to find
 * equilibrium values.
 *
 * <p><b>Data Ownership:</b>
 * <ul>
 *   <li><b>Immutable:</b> Cluster data (Stages 1-3) from AllClusterData â€" fixed at creation</li>
 *   <li><b>Mutable:</b> System parameters (ECI) and macro parameters (T, x) â€" can change anytime</li>
 *   <li><b>Cached:</b> Equilibrium results â€" invalidated when parameters change</li>
 * </ul>
 *
 * <p><b>Typical Usage:</b>
 * <pre>{@code
 * // Create model with initial parameters
 * CVMModelInput input = ...;
 * CVMPhaseModel model = CVMPhaseModel.create(input, eci, 1000.0, 0.5);
 *
 * // Query equilibrium properties (auto-minimizes on first call)
 * double G = model.getEquilibriumG();
 *
 * // Change temperature
 * model.setTemperature(1100.0);
 * double newG = model.getEquilibriumG();  // Auto-minimizes at new T
 *
 * // Scan over composition
 * for (double x = 0; x <= 1.0; x += 0.1) {
 *     model.setComposition(x);
 *     System.out.println("x=" + x + " G=" + model.getEquilibriumG());
 * }
 * }</pre>
 *
 * @see CVMFreeEnergy
 * @see NewtonRaphsonSolverSimple
 */
public class CVMPhaseModel {

    // =========================================================================
    // Input contract (moved here from standalone CVMInput.java to avoid an
    // extra file — CVMPhaseModel is the sole consumer of this bundle)
    // =========================================================================

    /**
     * Immutable input bundle wrapping Stage 1-3 topology data plus system metadata.
     *
     * <p>Created by {@link org.ce.domain.engine.cvm.CVMEngine} and consumed by
     * {@link CVMPhaseModel#create} to decouple the engine from workbench containers.</p>
     */
    public static final class CVMInput {

        private final ClusterIdentificationResult stage1;
        private final CFIdentificationResult stage2;
        private final CMatrixResult stage3;
        private final String systemId;
        private final String systemName;
        private final int numComponents;

        public CVMInput(
                ClusterIdentificationResult stage1,
                CFIdentificationResult stage2,
                CMatrixResult stage3,
                String systemId,
                String systemName,
                int numComponents) {

            if (systemId == null || systemId.isBlank())
                throw new IllegalArgumentException("systemId must not be blank");
            if (systemName == null || systemName.isBlank())
                throw new IllegalArgumentException("systemName must not be blank");
            if (numComponents < 2)
                throw new IllegalArgumentException("numComponents must be >= 2");
            if (stage1 == null || stage2 == null || stage3 == null)
                throw new IllegalArgumentException("Stage 1/2/3 data must be non-null");

            this.stage1 = stage1;
            this.stage2 = stage2;
            this.stage3 = stage3;
            this.systemId = systemId;
            this.systemName = systemName;
            this.numComponents = numComponents;
        }

        public ClusterIdentificationResult getStage1()  { return stage1; }
        public CFIdentificationResult      getStage2()  { return stage2; }
        public CMatrixResult               getStage3()  { return stage3; }
        public String                      getSystemId()      { return systemId; }
        public String                      getSystemName()    { return systemName; }
        public int                         getNumComponents() { return numComponents; }
    }

    // =========================================================================

    private static final Logger LOG = Logger.getLogger(CVMPhaseModel.class.getName());

    // =========================================================================
    // IMMUTABLE: Cluster Data (from AllClusterData) — never changes
    // =========================================================================

    private final int tcdis;              // Number of cluster types
    private final int tcf;                // Total CF count (including point)
    private final int ncf;                // Non-point CF count (optimization vars)
    private final List<Double> mhdis;     // HSP cluster multiplicities
    private final double[] kb;            // Kikuchi-Baker entropy coefficients
    private final double[][] mh;          // Normalized multiplicities: mh[t][j]
    private final int[] lc;               // Cluster count per type: lc[t]
    private final int[][] lcf;            // CF count per (type, group): lcf[t][j]
    private final int[][] cfBasisIndices; // Per-CF basis decoration (Stage 3)
    private final List<List<double[][]>> cmat_old;  // C-matrix (old orthogonal basis): cmat[t][j][v][k]
    private final List<List<double[][]>> cmat_new;  // C-matrix (CVCF basis, optional): transformed via CvCfIntegration
    private final int[][] lcv;                      // CV counts: lcv[t][j]
    private final List<List<int[]>> wcv;            // CV weights: wcv.get(t).get(j)[v]
    private final String systemId;
    private final String systemName;
    private final int numComponents;
    private final CvCfIntegration cvcfIntegration;  // Optional CVCF basis adapter (null = use old basis)

    // =========================================================================
    // MUTABLE: System Parameters â€" can be changed at any time
    // =========================================================================

    private double[] eci;            // Effective Cluster Interactions (CECs)
    private double tolerance;        // Convergence criterion

    // =========================================================================
    // MUTABLE: Macro Parameters â€" can be changed at any time
    // =========================================================================

    private double temperature;      // Kelvin
    private double[] moleFractions;  // Composition (length numComponents)

    // =========================================================================
    // CACHED: Equilibrium results â€" updated when parameters change
    // =========================================================================

    private boolean isMinimized = false;
    private long lastMinimizationTimeNanos;
    private double[] equilibriumCFs;              // [ncf] non-point CFs
    private CVMFreeEnergy.EvalResult equilibrium; // G, H, S, âˆ‡G, âˆ‡Â²G

    // =========================================================================
    // DIAGNOSTICS: Solver convergence info
    // =========================================================================

    private int lastIterations;
    private double lastGradientNorm;
    private String lastConvergenceStatus;
    private List<CVMSolverResult.IterationSnapshot> lastIterationTrace = new ArrayList<>();

    // =========================================================================
    // FACTORY METHOD
    // =========================================================================

    /**
    * Factory creates a CVM phase model from stage-model input.
     *
     * <p>Loads all cluster data (Stages 1-3) from the input contract.
     * User provides initial CEC, temperature, and composition values.
     * Performs first minimization on creation.
     *
     * @param input CVM stage model input
     * @param eci initial effective cluster interactions
     * @param temperature initial temperature in Kelvin
     * @param composition initial composition (binary: mole fraction of component 2)
     * @return CVMPhaseModel ready for queries
     * @throws IllegalArgumentException if input is invalid or parameters invalid
     * @throws Exception if first minimization fails
     */
    /**
     * Creates a CVM phase model for any K-component system (old orthogonal basis).
     *
     * <p>This is the canonical factory method for old basis calculations.
     * Pass the full mole-fraction array {@code x[0..K-1]} — works for binary, ternary, and any K.</p>
     *
     * @param input         CVMInput carrying cluster topology data (Stages 1-3)
     * @param eci           ncf-length effective cluster interactions (J/mol)
     * @param temperature   initial temperature in Kelvin
     * @param moleFractions mole fractions x[0..K-1]; must sum to 1.0
     * @return CVMPhaseModel with first N-R minimization complete (old basis)
     * @throws Exception if minimization fails or input is invalid
     */
    public static CVMPhaseModel create(
            CVMInput input,
            double[] eci,
            double temperature,
            double[] moleFractions) throws Exception {

        return create(input, eci, temperature, moleFractions, null);
    }

    /**
     * Creates a CVM phase model with optional CVCF basis transformation.
     *
     * <p>Pass {@code cvcfIntegration = null} for old orthogonal basis (default).
     * Pass a {@link CvCfIntegration} instance to use CVCF basis internally.</p>
     *
     * @param input         CVMInput carrying cluster topology data (Stages 1-3)
     * @param eci           ncf-length effective cluster interactions (J/mol)
     * @param temperature   initial temperature in Kelvin
     * @param moleFractions mole fractions x[0..K-1]; must sum to 1.0
     * @param cvcfIntegration optional CVCF adapter (null = old basis)
     * @return CVMPhaseModel with first N-R minimization complete
     * @throws Exception if minimization fails or input is invalid
     */
    public static CVMPhaseModel create(
            CVMInput input,
            double[] eci,
            double temperature,
            double[] moleFractions,
            CvCfIntegration cvcfIntegration) throws Exception {

        if (input == null) {
            throw new IllegalArgumentException("CVMInput must not be null");
        }

        CVMPhaseModel model = new CVMPhaseModel(input, cvcfIntegration);
        model.setECI(eci);
        model.setTemperature(temperature);
        model.setMoleFractions(moleFractions);   // validates length and sum-to-1
        model.ensureMinimized();
        return model;
    }

    /**
     * Creates a CVM phase model using the binary B-fraction shorthand (K=2 only).
     *
     * <p>Convenience for callers that have a scalar x_B. Delegates to
     * {@link #create(CVMInput, double[], double, double[])} via
     * {@code [1-xB, xB]}. For K &ge; 3 use the array overload directly.</p>
     *
     * @deprecated Use {@link #create(CVMInput, double[], double, double[])} for all K
     */
    @Deprecated
    public static CVMPhaseModel create(
            CVMInput input,
            double[] eci,
            double temperature,
            double xB) throws Exception {

        return create(input, eci, temperature, new double[]{1.0 - xB, xB});
    }

    // =========================================================================
    // CONSTRUCTOR (Private â€" only called by factory)
    // =========================================================================

    /**
     * Private constructor. Extracts and caches all stage data from CVMInput.
     * Only called by the factory method.
     *
     * @param input CVMInput with cluster topology
     * @param cvcfIntegration optional CVCF basis adapter (null = use old basis)
     */
    private CVMPhaseModel(CVMInput input, CvCfIntegration cvcfIntegration) {
        // Extract Stage 1: Cluster Identification
        ClusterIdentificationResult stage1 = input.getStage1();
        this.tcdis = stage1.getTcdis();
        this.mhdis = stage1.getDisClusterData().getMultiplicities();
        this.kb = stage1.getKbCoefficients();
        this.mh = stage1.getMh();
        this.lc = stage1.getLc();

        // Extract Stage 2: CF Identification
        CFIdentificationResult stage2 = input.getStage2();
        this.tcf = stage2.getTcf();
        this.ncf = stage2.getNcf();
        this.lcf = stage2.getLcf();

        // Extract Stage 3: C-Matrix Structure
        CMatrixResult stage3 = input.getStage3();
        this.cfBasisIndices = stage3.getCfBasisIndices();
        this.cmat_old = stage3.getCmat();
        this.cmat_new = (cvcfIntegration != null) ? cvcfIntegration.getTransformedCmat().getCmat() : null;
        this.lcv = stage3.getLcv();
        this.wcv = stage3.getWcv();

        // System info
        this.systemId = input.getSystemId();
        this.systemName = input.getSystemName();
        this.numComponents = input.getNumComponents();

        // Optional CVCF integration
        this.cvcfIntegration = cvcfIntegration;

        // Initialize parameter storage (not yet set)
        this.eci = null;
        this.temperature = Double.NaN;
        this.moleFractions = null;
        this.tolerance = 1.0e-6;
        this.isMinimized = false;
    }

    // =========================================================================
    // PARAMETER SETTERS (Trigger Re-minimization)
    // =========================================================================

    /**
     * Sets cluster interaction energies (CECs).
     * Invalidates cached results â€" next query will re-minimize.
     *
     * @param newECI effective cluster interactions (must have length ncf)
     * @throws IllegalArgumentException if length mismatch
     */
    public void setECI(double[] newECI) throws IllegalArgumentException {
        if (newECI == null || newECI.length < this.ncf) {
            throw new IllegalArgumentException(
                "ECI too short: got " + (newECI == null ? 0 : newECI.length) +
                ", expected >= " + this.ncf);
        }
        // Accept ncf+1 or ncf+2 arrays (raw CEC with trailing point/empty cluster terms)
        // by trimming to ncf. This avoids an infrastructure dependency in the job layer.
        this.eci = Arrays.copyOf(newECI, this.ncf);
        invalidateMinimization();
    }

    /**
     * Sets temperature.
     * Invalidates cached results â€" next query will re-minimize.
     *
     * @param T_K temperature in Kelvin (must be positive)
     * @throws IllegalArgumentException if temperature invalid
     */
    public void setTemperature(double T_K) throws IllegalArgumentException {
        if (T_K <= 0) {
            throw new IllegalArgumentException("Temperature must be positive: " + T_K);
        }
        this.temperature = T_K;
        invalidateMinimization();
    }

    /**
     * Sets composition using the binary shorthand (only valid for K=2).
     *
     * <p>For K=2: builds {@code moleFractions = {1 − x_B, x_B}}.</p>
     *
     * <p>For K≥3 use {@link #setMoleFractions(double[])} instead. Calling this
     * method on a multi-component model throws {@link IllegalArgumentException}
     * because a single scalar cannot specify the full composition vector, and
     * silently zeroing the extra components gives physically wrong results.</p>
     *
     * @param x_B mole fraction of component B in [0,1]
     * @throws IllegalArgumentException if composition invalid or numComponents > 2
     */
    /**
     * @deprecated Use {@link #setMoleFractions(double[])} for all K
     */
    @Deprecated
    public void setComposition(double x_B) throws IllegalArgumentException {
        if (x_B < 0 || x_B > 1) {
            throw new IllegalArgumentException("Composition must be in [0,1]: " + x_B);
        }
        if (numComponents != 2) {
            throw new IllegalArgumentException(
                "setComposition(double) is only valid for binary (K=2) systems. "
                + "This model has K=" + numComponents + " components. "
                + "Use setMoleFractions(double[]) to specify the full composition vector.");
        }

        this.moleFractions = new double[]{1.0 - x_B, x_B};
        invalidateMinimization();
    }

    /**
     * Sets full mole fraction vector (for K â‰¥ 3).
     *
     * @param fractions mole fractions (length numComponents, summing to 1.0)
     * @throws IllegalArgumentException if length mismatch or sum != 1
     */
    public void setMoleFractions(double[] fractions) throws IllegalArgumentException {
        if (fractions == null || fractions.length != numComponents) {
            throw new IllegalArgumentException(
                "Mole fractions length mismatch: got " + (fractions == null ? 0 : fractions.length) +
                ", expected " + numComponents);
        }

        double sum = 0;
        for (double x : fractions) {
            if (x < 0 || x > 1) {
                throw new IllegalArgumentException("Mole fractions must be in [0,1]");
            }
            sum += x;
        }
        if (Math.abs(sum - 1.0) > 1.0e-9) {
            throw new IllegalArgumentException("Mole fractions don't sum to 1: " + sum);
        }

        this.moleFractions = fractions.clone();
        invalidateMinimization();
    }

    /**
     * Sets convergence tolerance.
     *
     * @param tol convergence tolerance (recommended: 1e-6 to 1e-3)
     * @throws IllegalArgumentException if tolerance out of reasonable range
     */
    public void setTolerance(double tol) throws IllegalArgumentException {
        if (tol <= 0 || tol > 1.0e-3) {
            throw new IllegalArgumentException("Tolerance out of range: " + tol);
        }
        this.tolerance = tol;
        invalidateMinimization();
    }

    /**
     * Invalidates cached minimization results.
     * Next call to getEquilibrium* will trigger re-minimization.
     */
    private void invalidateMinimization() {
        this.isMinimized = false;
        this.equilibriumCFs = null;
        this.equilibrium = null;
        this.lastIterationTrace = new ArrayList<>();
    }

    // =========================================================================
    // CMAT SELECTION (Old vs CVCF basis)
    // =========================================================================

    /**
     * Returns the active C-matrix (either old orthogonal or CVCF basis).
     * @return cmat_new if CVCF is active, otherwise cmat_old
     */
    private List<List<double[][]>> getActiveCmat() {
        return (cvcfIntegration != null) ? cmat_new : cmat_old;
    }

    // =========================================================================
    // AUTOMATIC LAZY MINIMIZATION
    // =========================================================================

    /**
     * Ensures model is minimized.
     * If parameters changed since last minimization, re-minimizes.
     * Safe to call multiple times â€" only computes if needed.
     *
     * @throws Exception if minimization fails
     * @throws IllegalStateException if required parameters not set
     */
    public synchronized void ensureMinimized() throws Exception {
        if (isMinimized && equilibriumCFs != null && equilibrium != null) {
            return;  // Already minimized, all results cached
        }

        // Validate parameters are set
        if (eci == null || Double.isNaN(temperature) || moleFractions == null) {
            throw new IllegalStateException(
                "Cannot minimize: missing parameters\n" +
                "  ECI set: " + (eci != null) + "\n" +
                "  T set: " + !Double.isNaN(temperature) + "\n" +
                "  x set: " + (moleFractions != null));
        }

        // Perform minimization
        LOG.fine("CVMPhaseModel.ensureMinimized — starting minimization: T=" + temperature
                + " K, x=" + Arrays.toString(moleFractions));
        minimize();

        if (!isMinimized) {
            throw new Exception("Minimization failed: " + lastConvergenceStatus);
        }
    }

    /**
     * Internal minimization routine (private).
     * Uses either old orthogonal basis (default) or CVCF basis depending on cvcfIntegration.
     */
    private void minimize() {
        long startTime = System.nanoTime();

        try {
            // Select active C-matrix (old or CVCF basis)
            List<List<double[][]>> activeCmat = getActiveCmat();

            // Run Newton-Raphson solver
            // (Solver generates initial guess internally)
            CVMSolverResult result = NewtonRaphsonSolverSimple.solve(
                moleFractions,
                numComponents,
                temperature,
                eci,
                mhdis,
                kb,
                mh,
                lc,
                activeCmat,
                lcv,
                wcv,
                tcdis,
                tcf,
                ncf,
                lcf,
                cfBasisIndices,
                tolerance
            );

            // 3. Cache results
            this.lastIterationTrace = result.getIterationTrace();

            if (result.isConverged()) {
                this.equilibriumCFs = result.getEquilibriumCFs();

                // Evaluate thermodynamics at equilibrium non-point CFs
                this.equilibrium = CVMFreeEnergy.evaluate(
                    equilibriumCFs,
                    moleFractions,
                    numComponents,
                    temperature,
                    eci,
                    mhdis,
                    kb,
                    mh,
                    lc,
                    activeCmat,
                    lcv,
                    wcv,
                    tcdis,
                    tcf,
                    ncf,
                    lcf,
                    cfBasisIndices
                );

                this.isMinimized = true;
                this.lastIterations = result.getIterations();
                this.lastGradientNorm = result.getGradientNorm();
                this.lastConvergenceStatus = "OK";
                this.lastMinimizationTimeNanos = System.nanoTime() - startTime;

                logMinimizationSuccess();
            } else {
                this.isMinimized = false;
                this.lastConvergenceStatus = "Solver did not converge";
                this.lastIterations = result.getIterations();
                this.lastGradientNorm = result.getGradientNorm();
                logMinimizationFailure();
            }
        } catch (Exception e) {
            this.isMinimized = false;
            this.lastConvergenceStatus = "Exception: " + e.getMessage();
            LOG.log(Level.WARNING, "CVMPhaseModel.minimize — EXCEPTION: " + e.getMessage(), e);
        }
    }



    /**
     * Build full CF vector [tcf] from non-point CFs [ncf].
     * Adds point CFs based on composition.
     */
    private double[] buildFullCFVector(double[] u_nonpoint) {
        return ClusterVariableEvaluator.buildFullCFVector(
            u_nonpoint, moleFractions, numComponents, cfBasisIndices, ncf, tcf);
    }

    private void logMinimizationSuccess() {
        LOG.fine("CVMPhaseModel.minimize — SUCCESS: T=" + temperature + " K"
                + ", x=" + Arrays.toString(moleFractions)
                + ", iterations=" + lastIterations
                + ", ||dG||=" + String.format("%.4e", lastGradientNorm)
                + ", G=" + String.format("%.6e", equilibrium != null ? equilibrium.G : Double.NaN) + " J/mol"
                + ", elapsed=" + (lastMinimizationTimeNanos / 1_000_000) + " ms");
    }

    private void logMinimizationFailure() {
        LOG.warning("CVMPhaseModel.minimize — FAILED: T=" + temperature + " K"
                + ", status=" + lastConvergenceStatus
                + ", iterations=" + lastIterations);
    }

    // =========================================================================
    // QUERY INTERFACE: Thermodynamic Quantities
    // =========================================================================

    /**
     * Gets Gibbs energy of mixing at current equilibrium.
     * Automatically minimizes if parameters changed.
     *
     * @return G in J/mol
     * @throws Exception if minimization fails
     */
    public double getEquilibriumG() throws Exception {
        ensureMinimized();
        return equilibrium.G;
    }

    /**
     * Gets enthalpy of mixing at current equilibrium.
     */
    public double getEquilibriumH() throws Exception {
        ensureMinimized();
        return equilibrium.H;
    }

    /**
     * Gets entropy of mixing at current equilibrium.
     */
    public double getEquilibriumS() throws Exception {
        ensureMinimized();
        return equilibrium.S;
    }

    /**
     * Gets Helmholtz free energy F = H - TÂ·S (alternative to G).
     */
    public double getHelmholtzF() throws Exception {
        ensureMinimized();
        return equilibrium.H - temperature * equilibrium.S;
    }

    /**
     * Gets equilibrium non-point correlation functions.
     *
     * @return array of length ncf
     */
    public double[] getEquilibriumCFs() throws Exception {
        ensureMinimized();
        return equilibriumCFs.clone();
    }

    /**
     * Gets full equilibrium CF vector including point CFs.
     *
     * @return array of length tcf
     */
    public double[] getEquilibriumCFsFull() throws Exception {
        ensureMinimized();
        return buildFullCFVector(equilibriumCFs);
    }

    /**
     * Gets cluster variables at equilibrium.
     *
     * @return cv[tcdis][lc[t]][lcv[t][j]]
     */
    public double[][][] getEquilibriumCVs() throws Exception {
        ensureMinimized();
        double[] uFull = getEquilibriumCFsFull();
        return ClusterVariableEvaluator.evaluate(uFull, cmat, lcv, tcdis, lc);
    }

    /**
     * Gets short-range order parameters.
     * SRO_t measures ordering in cluster type t.
     *
     * @return array of length tcdis-1 (excluding point type)
     */
    public double[] getSROs() throws Exception {
        ensureMinimized();
        double[][][] cvs = getEquilibriumCVs();
        double[] sro = new double[tcdis];

        for (int t = 0; t < tcdis - 1; t++) {  // Exclude point type
            for (int j = 0; j < lc[t]; j++) {
                for (int v = 0; v < lcv[t][j]; v++) {
                    sro[t] += wcv.get(t).get(j)[v] * cvs[t][j][v];
                }
            }
        }
        return sro;
    }

    /**
     * Gets gradient of Gibbs energy (should be ~0 at minimum).
     */
    public double[] getGradient() throws Exception {
        ensureMinimized();
        return equilibrium.Gcu.clone();
    }

    /**
     * Gets gradient norm ||âˆ‡G|| (convergence measure).
     */
    public double getGradientNorm() throws Exception {
        ensureMinimized();
        return getGradientNormValue(equilibrium.Gcu);
    }

    /**
     * Gets Gibbs energy Hessian (for stability analysis).
     */
    public double[][] getStabilityMatrix() throws Exception {
        ensureMinimized();
        return copyMatrix(equilibrium.Gcuu);
    }

    /**
     * Checks if phase is thermodynamically stable (local minimum of G).
     *
     * <p>A phase is stable if and only if the Hessian d²G/du² is positive
     * definite — i.e., all eigenvalues are strictly positive. This is verified
     * via Cholesky decomposition: the decomposition succeeds (all pivots > 0)
     * if and only if the matrix is positive definite.</p>
     *
     * <p>The previous diagonal-only check ({@code H[i][i] > 0}) was insufficient:
     * a matrix can have all positive diagonal entries yet still be indefinite
     * (saddle point) when off-diagonal terms are large.</p>
     *
     * @return true if stable (Hessian is positive definite), false otherwise
     */
    public boolean isStable() throws Exception {
        ensureMinimized();
        double[][] H = equilibrium.Gcuu;
        return isPositiveDefinite(H);
    }

    /**
     * Tests positive definiteness via Cholesky decomposition.
     *
     * <p>Attempts the decomposition H = L · Lᵀ in-place. Succeeds iff every
     * pivot is strictly positive ({@code > STAB_EPS}). No allocation of L is
     * needed — only the diagonal pivots are required for the test.</p>
     *
     * @param H symmetric matrix to test
     * @return true if H is positive definite
     */
    private static boolean isPositiveDefinite(double[][] H) {
        final double STAB_EPS = 1.0e-10;
        int n = H.length;
        // Work on a copy to avoid mutating the Hessian
        double[][] A = new double[n][n];
        for (int i = 0; i < n; i++) {
            A[i] = H[i].clone();
        }
        for (int j = 0; j < n; j++) {
            // Subtract contributions from previous columns
            for (int k = 0; k < j; k++) {
                A[j][j] -= A[j][k] * A[j][k];
            }
            if (A[j][j] <= STAB_EPS) {
                return false;  // Non-positive pivot → not positive definite
            }
            double diag = Math.sqrt(A[j][j]);
            for (int i = j + 1; i < n; i++) {
                for (int k = 0; k < j; k++) {
                    A[i][j] -= A[i][k] * A[j][k];
                }
                A[i][j] /= diag;
            }
        }
        return true;
    }

    /**
     * Returns the unified domain {@link EquilibriumState} at current equilibrium.
     *
     * <p>CVM-specific diagnostics (gradient ∇G and Hessian ∇²G) are stored in
     * {@link EngineMetrics.CvmMetrics} and accessible via pattern matching.</p>
     */
    public EquilibriumState getEquilibriumState() throws Exception {
        ensureMinimized();
        return new EquilibriumState(
            temperature,
            moleFractions.clone(),
            equilibrium.H,
            equilibrium.G
        );
    }

    // =========================================================================
    // QUERY INTERFACE: Model Introspection
    // =========================================================================

    public double getTemperature() { return temperature; }
    public double[] getMoleFractions() { return moleFractions == null ? null : moleFractions.clone(); }
    public double[] getECI() { return eci == null ? null : eci.clone(); }
    public int getNumComponents() { return numComponents; }
    public int getNumCFs() { return ncf; }
    public int getTotalCFs() { return tcf; }
    public int getNumClusterTypes() { return tcdis; }
    public double getTolerance() { return tolerance; }
    public String getSystemId() { return systemId; }
    public String getSystemName() { return systemName; }

    /**
     * Minimization status check.
     */
    public boolean isMinimized() { return isMinimized; }

    /**
     * Gets diagnostic information about last minimization.
     */
    public int getLastIterations() throws IllegalStateException {
        if (!isMinimized) throw new IllegalStateException("Not yet minimized");
        return lastIterations;
    }

    public double getLastGradientNorm() throws IllegalStateException {
        if (!isMinimized) throw new IllegalStateException("Not yet minimized");
        return lastGradientNorm;
    }

    public long getLastMinimizationTimeMs() {
        return lastMinimizationTimeNanos / 1_000_000;
    }

    /** Returns per-iteration N-R diagnostics (CFs and dG/du). */
    public List<CVMSolverResult.IterationSnapshot> getLastIterationTrace() {
        return new ArrayList<>(lastIterationTrace);
    }

    /**
     * Computes G, H, S for a given non-point CF vector using the current model state.
     * Used by CVMEngine to compute thermodynamic properties at each N-R iteration.
     *
     * @param nonPointCFs the non-point cluster fractions (length ncf)
     * @return {G, H, S} in J/mol and J/(mol·K)
     */
    double[] computeGHS(double[] nonPointCFs) {
        CVMFreeEnergy.EvalResult r = CVMFreeEnergy.evaluate(
            nonPointCFs, moleFractions, numComponents, temperature, eci,
            mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf, lcf, cfBasisIndices);
        return new double[]{r.G, r.H, r.S};
    }

    /**
     * Summary report.
     */
    public String getSummary() {
        StringBuilder sb = new StringBuilder();
        sb.append("=== CVMPhaseModel ===\n");
        sb.append("System: ").append(systemName).append(" (").append(systemId).append(")\n");
        sb.append("Components: ").append(numComponents).append("\n");
        sb.append("Cluster types: ").append(tcdis).append("\n");
        sb.append("CFs: ").append(ncf).append(" (total ").append(tcf).append(")\n");
        sb.append("\nCurrent Parameters:\n");
        sb.append("  T: ").append(temperature).append(" K\n");
        if (moleFractions != null) {
            sb.append("  x: ").append(Arrays.toString(moleFractions)).append("\n");
        }
        sb.append("  ECI length: ").append(eci == null ? "not set" : eci.length).append("\n");
        sb.append("\nMinimization Status:\n");
        sb.append("  Minimized: ").append(isMinimized).append("\n");
        if (isMinimized && equilibrium != null) {
            try {
                sb.append("  G_eq: ").append(String.format("%.6e", getEquilibriumG())).append(" J/mol\n");
                sb.append("  H_eq: ").append(String.format("%.6e", getEquilibriumH())).append(" J/mol\n");
                sb.append("  S_eq: ").append(String.format("%.6e", getEquilibriumS())).append(" J/(molÂ·K)\n");
                sb.append("  ||âˆ‡G||: ").append(String.format("%.2e", getLastGradientNorm())).append("\n");
                sb.append("  Iterations: ").append(getLastIterations()).append("\n");
                sb.append("  Time: ").append(getLastMinimizationTimeMs()).append(" ms\n");
            } catch (Exception e) {
                sb.append("  [Query failed: ").append(e.getMessage()).append("]\n");
            }
        }
        return sb.toString();
    }

    // =========================================================================
    // UTILITY METHODS
    // =========================================================================

    private double getGradientNormValue(double[] grad) {
        double norm = 0.0;
        for (double g : grad) {
            norm += g * g;
        }
        return Math.sqrt(norm);
    }

    private double[][] copyMatrix(double[][] mat) {
        double[][] copy = new double[mat.length][];
        for (int i = 0; i < mat.length; i++) {
            copy[i] = mat[i].clone();
        }
        return copy;
    }

}