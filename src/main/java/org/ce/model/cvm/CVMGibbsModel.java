package org.ce.model.cvm;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.ProgressEvent;

import java.util.Arrays;
import java.util.function.Consumer;

/**
 * Physical model for the Cluster Variation Method (CVM): the Newton-Raphson
 * equilibrium minimisation at fixed composition, plus the Hillert per-phase
 * Newton step.
 *
 * <h2>Now a façade over {@link CvmEvaluator}</h2>
 *
 * <p>This class formerly carried the pipeline setup, seventeen fields of
 * geometry, a mutable {@code (T, x, u)} point, and its own copy of every free
 * energy expression. Those are gone. Setup lives in {@link CvmGeometry},
 * evaluation in {@link CvmEvaluator}/{@link CvmState}, and what remains here is
 * the part that is genuinely a solver: the minimisation loop, its step
 * limiting, and the per-phase linear system.</p>
 *
 * <p>Two things motivated the split:</p>
 *
 * <ul>
 *   <li><b>Duplicated expressions.</b> The entropy sum existed here in four
 *       near-identical loops (gradient and Hessian, each in an {@code ncf} and
 *       an {@code ncf+K} variant) and again in {@code CvmState}. An audit of a
 *       formula that exists five times checks nothing.</li>
 *   <li><b>Temporal coupling.</b> Every {@code calculateXxx} read fields staged
 *       by a prior {@code setT}/{@code setX}/{@code setU} sequence, a contract
 *       the compiler cannot express and which {@code checkMinimized()} could
 *       only partly guard. A {@link CvmState} carries its own point.</li>
 * </ul>
 *
 * <p>The public API is unchanged and every caller still works. Results are
 * bit-identical: {@code org.ce.scratch.CvmEvaluatorParity} gates the evaluator
 * against the expressions this class used to hold.</p>
 *
 * <p><b>Still stateful, deliberately.</b> {@code setU}/{@code setX}/{@code setT}
 * and the {@code cal*} accessors remain, because {@code ThermodynamicWorkflow}
 * and the GUI depend on the "minimise, then read results off the model" shape.
 * The mutable point is now a single {@link CvmState} reference rather than five
 * loose fields, so it cannot fall out of sync with itself. Callers that want no
 * state at all should use {@link CvmEvaluator} directly, as the Hillert solver
 * will once it is migrated.</p>
 */
public class CVMGibbsModel {

    /** Immutable lattice combinatorics; replaces the seventeen setup fields. */
    private CvmGeometry geo;
    /** Geometry bound to a Hamiltonian; the source of every {@link CvmState}. */
    private CvmEvaluator evaluator;
    private CECEntry cecEntry;

    /**
     * The current thermodynamic point, or null before one is set. Replaces the
     * former {@code u}/{@code x_mole}/{@code temp}/{@code currentCv}/{@code eci}
     * fields: those could disagree with each other (notably {@code setX} synced
     * cluster variables only when {@code u} was already non-null, and vice
     * versa), whereas a state is computed whole or not at all.
     */
    private CvmState current;

    /** Pending point, staged by the setters until all three parts are present. */
    private double[] pendingU;
    private double[] pendingX;
    private double pendingT = Double.NaN;

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

    /**
     * Calculated free energy and derivatives at a given (u, T, x) point.
     *
     * <p><b>These are the mixing quantities</b> (Gm/Hm/Sm and their
     * derivatives), not the pure-element-anchored absolutes. The field names
     * keep their historical short spelling, but every producer fills them from
     * the mixing accessors. {@code HillertSolver} relies on this: it reads
     * {@code evaluate(...).G} as Gm and adds {@code LatticeStability.g0m}
     * itself, so filling these with absolute values would double-count G0m.</p>
     */
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
     * Primary entry point for model formation. Delegates the Stage 1-4 pipeline
     * to {@link CvmGeometry#build} and binds the result to {@code cecEntry}.
     *
     * <p>The geometry is Hamiltonian-independent -- the pipeline never reads
     * {@code cecEntry} -- so it is cacheable on (elements, structure, model)
     * alone.</p>
     */
    public void initialize(
            String elements,
            String structure,
            String model,
            CECEntry cecEntry,
            Consumer<String> progressSink) {

        this.cecEntry = cecEntry;

        if (progressSink != null) {
            progressSink.accept(String.format("  > CEC Entry:         %s (%s)",
                    cecEntry != null ? cecEntry.elements : "null",
                    cecEntry != null ? cecEntry.structurePhase : "null"));
            if (cecEntry != null && cecEntry.cecTerms != null) {
                for (CECEntry.CECTerm term : cecEntry.cecTerms) {
                    progressSink.accept(String.format("    - %-10s: a = %10.6f, b = %10.6f",
                            term.name, term.a, term.b));
                }
            }
        }

        this.geo = CvmGeometry.build(elements, structure, model, progressSink);
        this.geo.validate();
        this.evaluator = new CvmEvaluator(geo, cecEntry);
    }

    /** The immutable geometry this model was initialized against. */
    public CvmGeometry getGeometry() {
        return geo;
    }

    /** The pure evaluator underlying this model. */
    public CvmEvaluator getEvaluator() {
        return evaluator;
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
     * Runs the fixed-composition Newton-Raphson minimisation.
     *
     * <p>Delegates to {@link CvmNewtonSolver}, which owns the loop, its
     * convergence criteria and its step limiting. This method only adapts
     * between that solver's {@link CvmNewtonSolver.Result} and the
     * {@link EquilibriumResult}/{@link ModelResult} shape callers expect, and
     * leaves this model's current point at the converged state so the
     * {@code cal*} accessors read from it.</p>
     */
    private EquilibriumResult minimize(
            double[] moleFractions, double temperature, double tolerance,
            Consumer<String> progressSink, Consumer<ProgressEvent> eventSink,
            org.ce.calculation.CalculationDescriptor.Property required) {

        CvmNewtonSolver.Result result = new CvmNewtonSolver(evaluator)
                .solve(temperature, moleFractions, tolerance, progressSink, eventSink);

        // Leave the model sitting at the converged point: ThermodynamicWorkflow
        // and the GUI minimise first, then read G/H/S off the model.
        this.pendingT = temperature;
        this.pendingX = moleFractions.clone();
        this.pendingU = result.u().clone();
        this.current = result.state();

        return new EquilibriumResult(
                modelResult(), result.u(), result.converged(),
                result.iterations(), result.finalGradientNorm());
    }

    /** Bundles the current state into the result shape callers expect. */
    private ModelResult modelResult() {
        return new ModelResult(
                current.gm(), current.hm(), current.sm(),
                current.gmu(), current.gmuu(),
                current.hmu(), current.smu(), current.smuu(),
                current.cfs());
    }


    // =========================================================================
    // Current-point staging
    //
    // The setters accumulate the parts of a point and rebuild the CvmState as
    // soon as all three are present. The former fields could disagree with one
    // another -- setX synced cluster variables only when u was already set, and
    // setU only when x was -- whereas a state is computed whole or not at all.
    // =========================================================================

    /** Sets the current correlation functions (non-point). */
    public void setU(double[] u) {
        this.pendingU = u.clone();
        refreshState();
        this.isMinimized = false;
    }

    /** Sets the current mole fractions (composition). */
    public void setX(double[] x) {
        this.pendingX = x.clone();
        this.currentComposition = x.clone();
        refreshState();
        this.isMinimized = false;
    }

    /** Sets the current temperature. */
    public void setT(double temperature) {
        setT(temperature, null);
    }

    public void setT(double temperature, Consumer<String> sink) {
        this.pendingT = temperature;
        this.currentTemperature = temperature;
        refreshState();
        this.isMinimized = false;
    }

    private void refreshState() {
        if (pendingU != null && pendingX != null && !Double.isNaN(pendingT)) {
            this.current = evaluator.stateAt(pendingT, pendingX, pendingU);
        }
    }

    private void checkMinimized() {
        if (!isMinimized) {
            throw new IllegalStateException("CVM Model is not minimized. Please call getEquilibriumState() first.");
        }
    }

    // =========================================================================
    // Mixing quantities (Gm = Hm - T*Sm)
    // =========================================================================

    public double calHm() {
        checkMinimized();
        return current.hm();
    }

    public double[] calHmu() {
        checkMinimized();
        return current.hmu();
    }

    public double[][] calHuu() {
        checkMinimized();
        return current.hmuu();
    }

    public double calSm() {
        checkMinimized();
        return current.sm();
    }

    public double[] calSmu() {
        checkMinimized();
        return current.smu();
    }

    public double[][] calSmuu() {
        checkMinimized();
        return current.smuu();
    }

    public double calGm() {
        checkMinimized();
        return current.gm();
    }

    public double[] calGmu() {
        checkMinimized();
        return current.gmu();
    }

    public double[][] calGmuu() {
        checkMinimized();
        return current.gmuu();
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
    //        the configurational entropy of mixing Sm. This is what the
    //        Newton-Raphson loop minimises and what CLAUDE.md's documented
    //        verification values (e.g. -3480.5209063901 for Nb-Ti) are
    //        anchored to.
    //
    //   G    the absolute Gibbs energy, G = G0m + Gm.
    //
    // Because G0m depends only on (x, T) and not on u, every u-derivative of
    // the absolute quantity equals the mixing one exactly, so calGmu/calGmuu
    // serve both. Only a widened gradient over uFull = [u ; x] differs, and
    // only in its trailing composition block.
    // =========================================================================

    /** Reference energy of the mechanical mixture of pure elements. */
    public double calG0m() {
        checkMinimized();
        return current.g0m();
    }

    /** {@code H0m = G0m} -- the reference term is pure energy. */
    public double calH0m() {
        checkMinimized();
        return current.h0m();
    }

    /** {@code S0m = 0} -- a mechanical mixture of pure elements has no entropy of mixing. */
    public double calS0m() {
        checkMinimized();
        return current.s0m();
    }

    /** Absolute Gibbs energy {@code G = G0m + Gm}. */
    public double calG() {
        checkMinimized();
        return current.g();
    }

    /** Absolute enthalpy {@code H = H0m + Hm}. */
    public double calH() {
        checkMinimized();
        return current.h();
    }

    /** Absolute entropy {@code S = S0m + Sm = Sm}. */
    public double calS() {
        checkMinimized();
        return current.s();
    }

    /**
     * Widened gradient of Gm with respect to {@code uFull = [u ; x]}, length
     * {@code ncf + K} -- for the Hillert multi-phase equilibrium solver.
     *
     * <p>Distinct from {@link #calGmu}, which must stay exactly {@code ncf}
     * long: {@code minimize()} sizes its linear system from it, so widening
     * that in place would silently change the single-phase solver to solve a
     * different system. The leading {@code ncf} entries of this vector equal
     * {@link #calGmu} exactly.</p>
     *
     * <p>Note this returns the <em>mixing</em> gradient, matching the previous
     * behaviour of this method: {@code HillertSolver} adds the pure-element
     * reference itself. {@link CvmState#guFull()} is the absolute version.</p>
     */
    public double[] calGuFull() {
        checkMinimized();
        return current.gmuFull();
    }

    /** Widened Hessian over {@code uFull} -- see {@link #calGuFull}. */
    public double[][] calGuuFull() {
        checkMinimized();
        return current.gmuuFull();
    }

    public double[] calCfs() {
        checkMinimized();
        return current.cfs();
    }

    // =========================================================================
    // Physics evaluation
    //
    // There is deliberately no evaluate(u, x, T) here any more. It duplicated
    // CvmEvaluator.stateAt: a stateless-looking method on a stateful class that
    // silently moved this model's current point as a side effect of a read, and
    // computed all nine quantities whether the caller wanted one or all. It
    // also left isMinimized false, so its result could not be read back through
    // the cal* accessors -- a method at odds with its own class.
    //
    // For evaluation at an arbitrary point use getEvaluator().stateAt(T, x, u),
    // which builds an immutable CvmState and touches nothing here.
    // =========================================================================

    /**
     * Result of {@link #solvePerPhaseStep}: the joint Newton step
     * {@code deltaY(mu)}, expressed as an <b>affine function of the trial
     * chemical-potential vector {@code mu}</b> rather than a value at one
     * fixed {@code mu} -- see the note on {@link #solvePerPhaseStep} for why
     * this shape, not a single numeric result, is what the outer Hillert
     * solver actually needs.
     *
     * <p>{@code deltaY(mu) = deltaY0 + Σ_k mu[k]*deltaYSensitivity[k]}, and
     * likewise for {@code deltaComposition}/{@code lambda}.</p>
     */
    public record PerPhaseStepResult(
            double[] deltaY0, double[][] deltaYSensitivity,
            double[] deltaComposition0, double[][] deltaCompositionSensitivity,
            double lambda0, double[] lambdaSensitivity) {

        /** Evaluates this affine result at a specific numeric {@code mu}. */
        public double[] deltaCompositionAt(double[] mu) {
            double[] result = deltaComposition0.clone();
            for (int k = 0; k < mu.length; k++) {
                for (int i = 0; i < result.length; i++) {
                    result[i] += mu[k] * deltaCompositionSensitivity[k][i];
                }
            }
            return result;
        }

        /** Evaluates the full joint deltaY (length ncf+K) at a specific numeric {@code mu}. */
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
     * implementation's {@code delxGCVM}.
     *
     * <p><b>Why affine-in-mu, not a fixed-mu numeric result:</b> tracing
     * {@code phaseq}'s outer loop shows {@code delxGCVM} is called with
     * {@code mu} still <em>symbolic</em>, and {@code genEqMat}'s outer
     * mass-balance equations substitute that symbolic result in directly, so
     * {@code mu} and {@code deltaN} are solved <em>simultaneously</em> in one
     * combined system -- the per-phase step is never evaluated at a numeric
     * {@code mu} on its own. Since {@code deltaY} is provably affine in
     * {@code mu} (the matrix {@code A} does not depend on it; only the
     * right-hand side does, and only in the x-block rows), the numeric
     * equivalent is to solve the same system {@code K+1} times against basis
     * right-hand sides -- once for {@code mu=0}, once per unit vector -- and
     * let {@link org.ce.model.equilibrium.EquilibriumMatrix} fold the affine
     * form into its own equations, mirroring {@code genEqMat}'s substitution
     * entirely numerically.</p>
     *
     * <p>Unrelated to {@code minimize()}'s single-phase Newton-Raphson loop
     * (fixed composition, stationary {@code G}) and must not be confused with
     * it: this solves for a stationary point of {@code G} <em>relative to a
     * trial {@code mu}</em>, with composition itself among the unknowns.</p>
     *
     * <p><b>The linear system</b> (at fixed T/P, so the {@code GxT*ΔT} and
     * {@code GxP*ΔP} terms vanish): unknowns are {@code deltaY[0..ncf+K-1]}
     * and {@code lambda}, over {@code ncf+K+1} equations:</p>
     * <ul>
     *   <li>Rows {@code 0..ncf-1} (u-block): {@code Guu[i,:] . deltaY = -Gu[i]}
     *       -- ordinary stationarity on the internal CFs, unconstrained by
     *       {@code mu}.</li>
     *   <li>Rows {@code ncf..ncf+K-1} (x-block): {@code Guu[i,:] . deltaY -
     *       lambda = mu[i-ncf] - Gu[i]} -- the only rows where {@code mu}
     *       appears, always with coefficient exactly {@code +1} on its own
     *       row, which is why one basis solve per component suffices.</li>
     *   <li>Row {@code ncf+K}: {@code sum(deltaY[ncf..]) = 0} -- the
     *       composition change stays on the simplex.</li>
     * </ul>
     *
     * <p>Built entirely from analytic widened derivatives -- no
     * finite-differencing anywhere.</p>
     *
     * @param uFull current joint state {@code [u ; x]}, length {@code ncf+K}
     * @param temperature current temperature, K
     */
    public PerPhaseStepResult solvePerPhaseStep(double[] uFull, double temperature) {
        int ncf = geo.ncf;
        int numComponents = geo.numComponents;
        int width = ncf + numComponents;
        if (uFull.length != width) {
            throw new IllegalArgumentException(
                    "uFull.length=" + uFull.length + " != ncf+K=" + width);
        }

        CvmState state = evaluator.stateAtFull(temperature, uFull);
        double[] Gu = state.gmuFull();
        double[][] Guu = state.gmuuFull();

        int n = width + 1; // + lambda

        // Matrix A is the same for every right-hand side (mu does not appear
        // in it) -- build once, reuse for all K+1 solves.
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
        for (int i = 0; i < width; i++) b0[i] = -Gu[i];
        double[] sol0 = LinearAlgebra.solve(A, b0);

        // Solving A*z = e_{ncf+k} directly gives d(deltaY)/d(mu_k), since the
        // system is linear and A is shared across right-hand sides.
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
     * Correlation functions of the fully disordered state at this composition
     * -- the Newton loop's starting iterate.
     *
     * <p>Delegates to {@link CvmNewtonSolver}, which owns the minimisation, so
     * the starting point used here and the one the solver actually starts from
     * cannot drift apart. Returns the full CVCF vector (length {@code tcf}) for
     * backwards compatibility; the solver itself takes only the leading
     * {@code ncf} block.</p>
     */
    public double[] computeRandomCFs(double[] moleFractions) {
        return geo.basis.computeRandomCvcfCFs(moleFractions, geo.pipelineResult);
    }

    /**
     * Largest fraction of the step from {@code uOld} to {@code uTrial} that
     * keeps every cluster variable within {@code [0, 1]} -- the reference
     * solver's {@code stpmx}. Delegates to {@link CvmNewtonSolver#stepLimit},
     * which the minimisation itself uses.
     */
    public double calculateStepLimit(double[] uOld, double[] uTrial, double[] moleFractions) {
        return new CvmNewtonSolver(evaluator).stepLimit(uOld, uTrial, moleFractions);
    }

    /**
     * True if every cluster variable at the given {@code (u, x)} point --
     * across <em>all</em> cluster types, including the point (composition)
     * block -- lies strictly inside {@code (0, 1)}. Side-effect-free, so
     * external callers -- notably the Hillert multi-phase solver, which must
     * validate a trial joint step across several models before committing to
     * it -- can check a candidate point without perturbing this instance's
     * minimization state.
     *
     * <p>Port of the reference implementation's {@code isValidParams}, which
     * checks {@code cvt} over the <em>full</em> {@code 1..tcdis} range --
     * unlike {@link #minClusterVariable}'s {@code findMin}, which explicitly
     * excludes the last (point) type. So this is deliberately <em>broader</em>
     * than {@code minClusterVariable}, not a stateless twin of it: it also
     * catches a trial point whose composition itself has drifted to a
     * pure-element boundary ({@code x_i = 0} or {@code 1}), exactly the kind
     * of boundary the Hillert backtracking is meant to catch. Do not
     * "simplify" this back to the {@code tcdis-1} exclusion; that would
     * silently narrow what this check catches.</p>
     */
    public boolean isValidParams(double[] u, double[] moleFractions) {
        double[][][] cv = evaluateClusterVariables(u, moleFractions);
        for (int t = 0; t < geo.tcdis; t++) {
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
        return geo.numComponents;
    }

    public String getElements() {
        return geo.elements;
    }

    public CvCfBasis getBasis() {
        return geo.basis;
    }

    public int getNcf() {
        return geo.ncf;
    }

    public int getTcf() {
        return geo.tcf;
    }

    public int getTcdis() {
        return geo.tcdis;
    }

    public int[] getLc() {
        return geo.lc;
    }

    public int[][] getLcv() {
        return geo.lcv;
    }

    public int[][] getOrthCfBasisIndices() {
        return geo.orthCfBasisIndices;
    }

    /**
     * Computes cluster variables cv[t][j][v] from the given non-point CFs and
     * composition.
     */
    public double[][][] evaluateClusterVariables(double[] u, double[] moleFractions) {
        return geo.evaluateCVs(u, moleFractions);
    }
}
