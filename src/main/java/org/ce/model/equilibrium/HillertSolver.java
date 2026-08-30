package org.ce.model.equilibrium;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.ModelSession;
import org.ce.model.cvm.CVMGibbsModel;

/**
 * N-phase equilibrium solver (Hillert's method), ported from the reference
 * {@code phaseq} Mathematica implementation -- see HILLERT_SOLVER_PLAN.md.
 *
 * <p>{@code mu} throughout is the <b>absolute (SER-referenced) component
 * chemical potential</b>: {@link PhaseStep} solves the stationarity of the
 * absolute {@code G = G0m + Gm} (its x-block uses {@code dG/dx_i}, including the
 * pure-element reference {@code G0_i(T)}), and {@link EquilibriumMatrix} uses the
 * same absolute {@code G} in its Gibbs-Duhem rows. The two halves therefore
 * share one reference; this was verified against the {@code phaseq} reference
 * (see {@code org.ce.scratch.HillertAbsoluteMuValidation}).</p>
 *
 * <p>Outer loop, once per iteration:
 * <ol>
 *   <li>{@link PhaseStep} for every phase -- each
 *       phase's Newton step {@code deltaY(mu)} as an affine function of the
 *       shared trial chemical potential {@code mu}.</li>
 *   <li>Assemble {@link EquilibriumMatrix.PhaseContribution} for every
 *       currently-<em>active</em> phase ({@link Phase#active}) using
 *       {@code G = G0m + Gm} (absolute reference energy -- see
 *       {@link LatticeStability}) and solve for numeric {@code mu}/{@code
 *       deltaN}. If that solve returns a non-finite result, the run stops
 *       with {@link ConvergenceReason#NUMERICAL_BREAKDOWN} before any line
 *       search.</li>
 *   <li>Backtracking line search: starting at {@code lambda=1}, evaluate the
 *       trial step for every active phase and halve {@code lambda} until all
 *       active phases' trial state is valid ({@code State.isValidIncludingPoints})
 *       and every trial amount is finite and {@code >= 0}, up to
 *       {@code innerBacktrackTries} tries.</li>
 *   <li>Accept the trial state, then apply phase removal (STEP 8) and phase
 *       addition (STEP 12); if neither changed the phase set, check
 *       convergence via the per-phase step-norm.</li>
 * </ol>
 *
 * <p><b>Convergence rule and its one deliberate deviation from {@code phaseq}.</b>
 * The reference tests {@code Min[normList] <= 1e-10} where {@code normList[p]}
 * is the L2 norm of phase {@code p}'s full joint Newton step {@code deltaY_p(mu)}.
 * {@code Min} is sound <em>there</em> because {@code phaseq}'s phase-set
 * management keeps every phase in {@code normList} genuinely participating at a
 * true fixed point, so {@code Min}, {@code Max} and "all below tol" coincide.
 * This solver has phase removal (STEP 8) and phase addition (STEP 12), but its
 * candidate search is a coarse composition grid (see below), not the reference's
 * exact driving-force test -- so a stable phase whose favourable composition
 * falls between grid points could be missed, and a non-participating or
 * near-dead phase could still drive {@code Min} below {@code tol} while a
 * genuinely stable phase is far from equilibrium (a false positive). So
 * convergence here uses {@code max} over the currently <em>active</em> phases:
 *
 * <pre>
 *   converged  iff  a step was accepted this iteration
 *              AND   no phase was added or removed this iteration
 *              AND   max_{p active} ||deltaY_p(mu)||_2  <= tol
 *              AND   outer-solve relative residual      <= max(1e-8, tol)
 *              AND   mass-balance residual              <= massRelGate   (STEP 10)
 * </pre>
 *
 * <p>This accepts the identical set of true-converged states {@code phaseq}
 * accepts (all norms ~0 at equilibrium, so their max is ~0 too) and additionally
 * rejects only the false positives {@code Min} would have let through. It is not
 * "Max because Max seems better" -- it is the minimum change that makes the
 * criterion sound without the exact phase-set management {@code Min} depends
 * on.</p>
 *
 * <p>Non-converged runs report a {@link ConvergenceReason}. Besides
 * {@link ConvergenceReason#CONVERGED} the currently implemented reasons are:
 * {@link ConvergenceReason#MAX_ITERATIONS iteration cap},
 * {@link ConvergenceReason#LINE_SEARCH_FAILED line-search failure} (no feasible
 * step from the current state),
 * {@link ConvergenceReason#STALLED numerical stall} (two consecutive accepted
 * steps below the machine/tolerance floor),
 * {@link ConvergenceReason#SINGULAR_OUTER_SYSTEM a singular outer system}
 * (exactly rank-deficient -- a redundant phase set),
 * {@link ConvergenceReason#NUMERICAL_BREAKDOWN a non-finite outer solve}
 * (NaN/Inf {@code mu} or residual),
 * {@link ConvergenceReason#INITIAL_MASS_BALANCE an invalid problem definition}
 * (bad {@code overallAmounts}, or a seed that does not represent it), and
 * {@link ConvergenceReason#MASS_BALANCE_DRIFT a converged step-norm whose
 * represented inventory cannot be reconciled with the target}. A failed line
 * search can no longer be mistaken for convergence: the rule requires a step to
 * have been accepted.</p>
 *
 * <h2>Phase-set management (STEP 8)</h2>
 *
 * <p>Each {@link Phase} carries an {@link Phase#active} flag; only active phases
 * contribute a Gibbs-Duhem row, a {@code deltaN} unknown and mass-balance
 * terms, and only active phases are stepped by the line search. A phase seeded
 * with {@code amount > 0} starts active.</p>
 *
 * <p><b>Phase removal</b> is implemented. When an active phase's amount reaches
 * numerical zero -- below {@code 1e-9} of the total represented amount while
 * still shrinking, after at least {@code 3} accepted iterations of residence --
 * it is removed: its amount is set to exactly zero, its state frozen, its
 * {@link Phase#active} flag cleared, and a {@link PhaseSetEvent} of type
 * {@link PhaseSetEventType#PHASE_REMOVED} recorded. The next iteration rebuilds
 * the outer system on the reduced active set and the run continues (this is
 * <em>not</em> a {@link ConvergenceReason#STALLED} outcome). Removing a
 * zero-amount phase changes the represented total {@code sum_p N_p x^p_i} by
 * exactly {@code N*x = 0}, so removal is exactly mass-conserving under this
 * solver's self-referential (Cauchy) mass balance and needs no
 * overall-composition target.</p>
 *
 * <p><b>Phase addition (STEP 12).</b> After an accepted active-set step and the
 * removal block, if there is an inactive candidate phase and the active set is
 * near its own fixed point, each inactive candidate is scored by the absolute
 * tangent-plane driving force
 * {@code dGf_beta = sum_i mu_i x^beta_i - G^beta} at a <em>relaxed</em>
 * candidate state -- its composition searched on a coarse grid (1-D for a
 * binary, a small barycentric grid for a ternary; higher K is not searched and
 * the candidate is reported unevaluated) and its internal CVCF variables
 * minimised at each trial composition by {@link CvmNewtonSolver}. The candidate
 * {@link Phase} is <b>not mutated</b> during this scan. At most one phase is
 * added per outer iteration: the candidate with the largest {@code dGf_beta}
 * above {@code addThreshold = max(1.0, 1e-6 * |G^beta|)} J/mol is activated with
 * amount {@code epsilon = 1e-6 * sum(overallAmounts)} and its {@code uFull} set
 * from the relaxed candidate state; a {@link PhaseSetEventType#PHASE_ADDED}
 * event is recorded. The inserted mass is <b>not</b> hand-redistributed -- the
 * next iteration's STEP-10 target-aware mass-balance rows carry
 * {@code r_i = b_i - sum_p N_p x^p_i = -epsilon*x^beta_i} and drive the
 * redistribution (verified on an analytic two-phase system in STEP 11 --
 * {@code org.ce.scratch.HillertAbsoluteMuValidation} Part 9). An iteration that
 * adds a phase is <b>not</b> a converged iteration: the stall/convergence check
 * is skipped for it, exactly as for a removal.</p>
 *
 * <p><b>Anti-oscillation.</b> A just-added or just-removed phase carries a small
 * {@link Phase#phaseSetCooldown} (set to {@code MIN_RESIDENCE}); a phase may not
 * be added while its cooldown is nonzero, and the cooldown decrements once per
 * outer iteration. This is the minimum guard against an add/remove cycle; it is
 * not a general phase-set state machine.</p>
 *
 * <h2>V1 scope</h2>
 *
 * <p><b>Implemented and validated:</b></p>
 * <ul>
 *   <li>N-phase equilibrium for <b>single-site, one-atom-per-site disordered
 *       phases</b> (BCC_A2 / FCC_A1 / HCP_A3), where a phase's mole fraction
 *       {@code x^p_i} equals its moles of component {@code i} per formula unit,
 *       so {@code target_i = sum_p N_p x^p_i}.</li>
 *   <li><b>Fixed T and P</b> ({@code GxT}/{@code GxP} terms vanish).</li>
 *   <li>An <b>explicit {@code overallAmounts} target</b> (not normalised;
 *       {@code sum_i} is the total system amount). The seed must represent it
 *       within the entry tolerance or the run is rejected without iterating.</li>
 *   <li><b>Phase removal</b> (STEP 8) and <b>phase addition</b> (STEP 12), one
 *       phase-set change per iteration, with a {@code MIN_RESIDENCE} cooldown.</li>
 *   <li>Phase-addition <b>candidate search for K = 2 (binary) and K = 3
 *       (ternary)</b>: a coarse composition grid, {@link CvmNewtonSolver} for
 *       the internal relaxation at each trial point, absolute driving force
 *       {@code sum_i mu_i x_i - G}.</li>
 * </ul>
 *
 * <p><b>Deliberate V1 limitations</b> (reported as a non-converged
 * {@link ConvergenceReason}, never as a silently wrong answer):</p>
 * <ul>
 *   <li>Phase-addition candidate search for <b>K &gt;= 4</b> is not implemented;
 *       such a candidate is reported unevaluated and never activated.</li>
 *   <li>The coarse candidate grid gives <b>no global guarantee</b>: a stable
 *       phase whose favourable composition lies between grid points behind a
 *       concave barrier can be missed. A missed phase surfaces as
 *       {@link ConvergenceReason#MASS_BALANCE_DRIFT}, not a wrong equilibrium.</li>
 *   <li>The known near-edge (dilute-composition) Newton fragility of the inner
 *       CVM solve is inherited here; near-edge compositions can fail to
 *       converge (reported), and are out of V1 scope.</li>
 *   <li>{@code HillertSolver} does not itself poll {@link Thread#interrupted()};
 *       cooperative cancellation is observed only where the run calls
 *       {@link CvmNewtonSolver} (the candidate scan), which does.</li>
 * </ul>
 *
 * <p><b>Future enhancements (not in V1):</b> multi-sublattice / ordered phases
 * (an {@code M_A} generalisation of mass balance); {@code GxT}/{@code GxP} for a
 * scanning T or P; a global (non-grid) phase-stability test; a dedicated
 * {@code PHASE_SET_OSCILLATION} reason.</p>
 */
public final class HillertSolver {

    private HillertSolver() {}

    // =========================================================================
    // Inputs and outputs
    //
    // Nested here for the same reason CvmNewtonSolver.Result is nested in its
    // solver: these are this solver's working state and its output, used
    // nowhere else. Keeping them alongside the loop that drives them means the
    // whole multi-phase contract reads in one place.
    // =========================================================================

    /**
     * One candidate phase's mutable state during the solve.
     *
     * <p>The {@link CVMGibbsModel} it carries is a pure evaluator and holds no
     * per-point state, so it may safely be shared between phases of the same
     * system; what is mutable here is {@link #amount} and {@link #uFull}, which
     * the outer loop updates. (An earlier version required a separate model
     * instance per phase, because the model then carried a current
     * {@code (T, x, u)} internally.)</p>
     *
     * <p>Distinct from a single-phase {@link ModelSession}-driven calculation:
     * a Hillert phase's composition is itself an unknown, solved for jointly
     * with its internal CVM parameters by {@link PhaseStep} -- not
     * a fixed input the way {@link org.ce.calculation.Conditions} treats it for
     * {@code CalculationService.calculate}.</p>
     */
    public static final class Phase {

        /**
         * Amount below which a phase is treated as unstable. Matches the
         * reference's amount-sign-only check ({@code amount > 0}), not a
         * rigorous Gibbs phase rule.
         */
        private static final double STABILITY_THRESHOLD = 0.0;

        public final String label;
        public final ModelSession session;
        public final CVMGibbsModel model;
        public final int ncf;
        public final int numComponents;

        /** Current phase amount (moles of formula units, "N" in the reference). */
        public double amount;

        /** Current joint internal-parameter vector {@code uFull = [u ; x]}, length {@code ncf+K}. */
        public double[] uFull;

        /**
         * Whether this phase currently participates in the outer equilibrium
         * system (contributes a Gibbs-Duhem row, a {@code deltaN} unknown,
         * mass-balance terms, and is stepped by the line search).
         *
         * <p>An <b>inactive candidate</b> ({@code active == false}) keeps its
         * model, label and {@code uFull} (hence {@code composition()}) but
         * contributes nothing to a Newton step -- it is frozen until phase-set
         * management activates it. Seeded {@code true} from
         * {@code initialAmount > 0}. Cleared to {@code false} by phase removal
         * (STEP 8) when a phase's amount reaches numerical zero. Set back to
         * {@code true} by phase addition (STEP 12) when an inactive candidate's
         * relaxed tangent-plane driving force exceeds {@code addThreshold} --
         * its {@code uFull} is then overwritten with the relaxed candidate state
         * and its {@code amount} set to {@code 1e-6 * sum(overallAmounts)}. The
         * outer loop filters exclusively on this flag (not on
         * {@code amount > 0}); the two can differ for one iteration around a
         * phase-set change.</p>
         */
        public boolean active;

        /**
         * Number of accepted outer iterations this phase has been continuously
         * active for -- gates removal ({@code MIN_RESIDENCE}) so a phase that
         * dips low transiently and recovers is not removed prematurely.
         */
        int activeResidence;

        /**
         * Outer iterations remaining before this phase may take part in a
         * phase-set event again -- set to {@code MIN_RESIDENCE} when the phase
         * is added or removed, decremented once per outer iteration (STEP 12
         * anti-oscillation). A phase with {@code phaseSetCooldown > 0} is not a
         * candidate for addition. Removal ignores it -- a phase that genuinely
         * hits zero amount must still be removable immediately.
         */
        int phaseSetCooldown;

        public Phase(String label, ModelSession session, CVMGibbsModel model,
                double initialAmount, double[] initialUFull) {
            this.label = label;
            this.session = session;
            this.model = model;
            this.ncf = model.ncf();
            this.numComponents = initialUFull.length - ncf;
            this.amount = initialAmount;
            this.uFull = initialUFull.clone();
            this.active = initialAmount > STABILITY_THRESHOLD;
            this.activeResidence = 0;
            this.phaseSetCooldown = 0;
        }

        /**
         * Current composition -- the trailing {@code K} entries of
         * {@link #uFull}. Port of the reference's {@code updateComp}:
         * composition is always exactly this slice, never a separate inversion.
         */
        public double[] composition() {
            double[] x = new double[numComponents];
            System.arraycopy(uFull, ncf, x, 0, numComponents);
            return x;
        }

        /**
         * True if this phase's amount is strictly positive -- the reference's
         * amount-sign-only physical check. Distinct from {@link #active}: a
         * just-removed phase has {@code amount == 0} and {@code active == false};
         * the two agree for every phase the current version produces, but
         * {@code active} is what the outer loop filters on.
         */
        public boolean isStable() {
            return amount > STABILITY_THRESHOLD;
        }
    }

    /** Why {@link #solve} stopped iterating. Exactly one applies to a {@link Result}. */
    public enum ConvergenceReason {
        /** {@code max} per-stable-phase Newton-step norm fell to {@code tol} on an
         *  iteration whose step was accepted and whose outer solve was trustworthy. */
        CONVERGED,
        /** The outer iteration cap was hit without meeting the convergence rule. */
        MAX_ITERATIONS,
        /** The final executed iteration's backtracking line search found no valid
         *  step at any {@code lambda}; the state could not move and convergence
         *  was not independently established. */
        LINE_SEARCH_FAILED,
        /** Two consecutive accepted iterations moved the joint state by a
         *  relative amount below the machine/tolerance floor while still above
         *  {@code tol} -- no meaningful progress is being made. */
        STALLED,
        /** {@link EquilibriumMatrix#solve} could not solve the outer system
         *  (two stable phases with the same composition -> identical Gibbs-Duhem
         *  rows -> singular matrix). A physically redundant phase set; resolving
         *  it is phase-set management, not this loop's job. */
        SINGULAR_OUTER_SYSTEM,
        /** The outer solve returned a non-finite result -- a {@code NaN} or
         *  {@code Inf} chemical potential or backward-error residual (e.g. a
         *  phase driven to a degenerate state whose widened derivatives are not
         *  finite fed a non-finite coefficient into the outer system, which
         *  {@link LinearAlgebra#solveChecked}'s exact-singularity pivot guard
         *  does not catch). Distinct from {@link #SINGULAR_OUTER_SYSTEM} (an
         *  exactly rank-deficient matrix, which <em>throws</em>): here the solve
         *  completes but its output is numerically meaningless. The loop stops
         *  immediately -- no {@link PhaseStep}, no line search -- and the
         *  {@link Result}'s {@code mu} is the last finite value (all-zero if the
         *  first iteration broke down); it is never a {@code NaN}/{@code Inf}
         *  vector presented as a plausible state. */
        NUMERICAL_BREAKDOWN,
        /** The problem definition is invalid: the {@code overallAmounts} target
         *  is malformed (wrong length, a negative / non-finite entry, or a
         *  non-positive sum), or the seeded phase states do not represent it
         *  ({@code sum_p N_p x^p_i != target_i} beyond the mass-balance
         *  tolerance). The solver does not iterate -- the caller must fix the
         *  problem. No initial amount or composition is silently changed. */
        INITIAL_MASS_BALANCE,
        /** The step-norm criterion was met but the represented overall
         *  composition {@code sum_p N_p x^p_i} drifted off {@code overallAmounts}
         *  by more than the mass-balance tolerance -- converged in the Newton
         *  sense, but not to the prescribed inventory. The iteration equations
         *  are unchanged; this is a reported downgrade of what would otherwise
         *  be {@link #CONVERGED}. */
        MASS_BALANCE_DRIFT
    }

    /** A change to the active phase set during a {@link #solve} run. */
    public enum PhaseSetEventType {
        /** An active phase's amount reached numerical zero and it was removed. */
        PHASE_REMOVED,
        /** An inactive candidate's absolute tangent-plane driving force
         *  {@code sum_i mu_i x^beta_i - G^beta} exceeded {@code addThreshold}
         *  and it was activated with amount {@code 1e-6 * sum(overallAmounts)}
         *  from its relaxed candidate state (STEP 12). The following iteration's
         *  STEP-10 mass-balance rows redistribute the inserted mass. */
        PHASE_ADDED
    }

    /**
     * One active-phase-set change. Reported through {@link ConvergenceReport}
     * and the progress sink -- not a framework, just a record so a caller can
     * assert on what happened.
     *
     * @param drivingForce {@code NaN} for a removal (there is no driving force
     *                     to report -- the phase simply reached zero amount);
     *                     for an addition, the selected candidate's absolute
     *                     tangent-plane driving force {@code sum_i mu_i x_i - G}
     *                     at its relaxed state (in J/mol, always
     *                     {@code > addThreshold > 0}).
     */
    public record PhaseSetEvent(
            PhaseSetEventType type, String label, int iteration,
            double oldAmount, double newAmount, double drivingForce) {
    }

    /**
     * Structured convergence diagnostics for one {@link #solve} run -- the
     * independent quantities the outer loop can measure, so a caller can see
     * <em>why</em> a run did or did not converge rather than just a single
     * norm. None of these is a raw combined L2 norm over everything; each is on
     * its own scale.
     *
     * @param reason              why the loop stopped
     * @param maxPhaseStepNorm    {@code max_{p active} ||deltaY_p(mu)||_2} at the
     *                            final executed iteration -- the reference's
     *                            {@code normList} entries, but taken as the
     *                            <em>maximum</em> over participating phases, not
     *                            the minimum (see {@link #solve}'s class note on
     *                            the deliberate deviation from {@code phaseq}'s
     *                            {@code Min}). Compared to {@code tol}.
     * @param maxPhaseAmountStep  {@code max_p |lambda_accepted * deltaN_p| / sum_p N_p}
     *                            at the final iteration -- the largest relative
     *                            phase-amount move. {@code NaN} if no step was
     *                            accepted.
     * @param linearSolveResidual {@code eqStep.relativeResidual()} of the outer
     *                            {@code (K+np)} solve at the final iteration --
     *                            backward error, dimensionless.
     * @param lastLambda          the accepted {@code lambda} at the final
     *                            iteration, or {@code NaN} if the line search
     *                            failed there.
     * @param lastStepAccepted    whether a step was actually applied on the
     *                            final iteration.
     * @param allStatesValid      whether every phase's final {@code uFull}
     *                            passes {@link CVMGibbsModel.State#isValidIncludingPoints()}.
     * @param iterationsRun       number of outer iterations actually executed
     *                            (no off-by-one compensation).
     * @param phaseSetEvents      active-set changes during the run, in order.
     *                            Empty for the common no-event case.
     * @param massBalance         the prescribed vs achieved overall component
     *                            inventory and its residuals -- always populated.
     */
    public record ConvergenceReport(
            ConvergenceReason reason,
            double maxPhaseStepNorm,
            double maxPhaseAmountStep,
            double linearSolveResidual,
            double lastLambda,
            boolean lastStepAccepted,
            boolean allStatesValid,
            int iterationsRun,
            List<PhaseSetEvent> phaseSetEvents,
            MassBalanceReport massBalance) {
    }

    /**
     * The system's conserved component inventory: what was prescribed
     * ({@code overallAmounts}) versus what the accepted phase states represent
     * ({@code sum_{p active} N_p x^p_i}), with residuals reconstructed
     * <em>independently</em> from the phase states -- never from the Newton
     * right-hand side.
     *
     * <p><b>Scope.</b> This contract holds for single-site, one-atom-per-site
     * disordered phases (BCC_A2 / FCC_A1 / HCP_A3), where a phase's mole
     * fraction {@code x^p_i} is exactly its moles of component {@code i} per
     * formula unit, so {@code target_i = sum_p N_p x^p_i}. Ordered /
     * multi-sublattice phases would need an {@code M_A} generalization that is
     * out of scope.</p>
     *
     * @param targetOverall     prescribed overall composition, normalized to
     *                          {@code sum_i = 1} ({@code overallAmounts[i] /
     *                          sum_j overallAmounts[j]}), for readability.
     * @param calculatedOverall achieved overall composition, normalized:
     *                          {@code (sum_p N_p x^p_i) / sum_p N_p}.
     * @param maxAbsResidual    {@code max_i |overallAmounts[i] - sum_p N_p x^p_i|}
     *                          at the FINAL state, in moles -- the primary
     *                          quantity, not normalized, so a real inventory
     *                          error is not hidden.
     * @param maxRelResidual    {@code maxAbsResidual / max(sum_i overallAmounts[i], 1)}
     *                          -- dimensionless, for a scale-independent gate.
     * @param residualBeforeLastStep {@code max_i |r_i|} at the START of the last
     *                          outer iteration (before its Newton step), in
     *                          moles. Together with {@link #maxAbsResidual}
     *                          (after the step) this shows whether the STEP-10
     *                          mass-balance Newton correction is reducing the
     *                          residual per iteration (quadratically near a
     *                          solution). {@code NaN} if no iteration ran.
     */
    public record MassBalanceReport(
            double[] targetOverall,
            double[] calculatedOverall,
            double maxAbsResidual,
            double maxRelResidual,
            double residualBeforeLastStep) {
    }

    /**
     * Immutable output of {@link #solve} -- the multi-phase counterpart to
     * {@code ThermodynamicResult} for the single-phase path.
     *
     * <p><b>Check {@link #overallConverged} (or
     * {@link #convergenceReport}{@code .reason()}) before using any value.</b> A
     * non-converged run still returns plausible-looking numbers.</p>
     *
     * <p>{@link #overallConverged()}, {@link #outerIterations()} and
     * {@link #finalResidualNorm()} are kept as convenience accessors; they now
     * read straight off {@link #convergenceReport} so the two can never
     * disagree.</p>
     */
    public record Result(
            List<PhaseResult> phases,
            double[] mu,
            ConvergenceReport convergenceReport) {

        /** True iff {@link #convergenceReport}{@code .reason() == CONVERGED}. */
        public boolean overallConverged() {
            return convergenceReport.reason() == ConvergenceReason.CONVERGED;
        }

        /** Outer iterations actually executed -- {@link ConvergenceReport#iterationsRun()}. */
        public int outerIterations() {
            return convergenceReport.iterationsRun();
        }

        /**
         * The convergence metric at the final iteration --
         * {@link ConvergenceReport#maxPhaseStepNorm()} (the <em>max</em>
         * per-stable-phase Newton-step norm, not the old {@code min}).
         */
        public double finalResidualNorm() {
            return convergenceReport.maxPhaseStepNorm();
        }
    }

    /**
     * One phase's outcome: amount, composition, and energetics at the final
     * iterate.
     *
     * <p>{@link #state} is that phase's model evaluated at its converged joint
     * point {@code uFull = [u ; x]} -- the same object {@code model.atFull(T,
     * uFull)} would produce, retained rather than discarded so every other
     * property is reachable without re-solving or re-evaluating:</p>
     *
     * <pre>
     *   for (PhaseResult p : eq.phases()) {
     *       double s    = p.state().sm();              // entropy of this phase
     *       double[] dg = p.state().gmuFull();         // its widened gradient
     *       var sro     = p.state().pairSroByShell();  // its short-range order
     *   }
     * </pre>
     *
     * <p>{@link #g} is kept as its own field because it is the quantity the
     * outer equilibrium assembly actually solved with -- the absolute
     * {@code G = G0m + Gm}, which must share one zero across phases for
     * chemical potentials to be comparable. It equals {@code state().g()}; the
     * field records what the solve used, the state offers everything else.</p>
     */
    public record PhaseResult(
            String label,
            double amount,
            double[] composition,
            double g,
            CVMGibbsModel.State state,
            boolean phaseConverged) {

        /** The model this phase was evaluated against. */
        public CVMGibbsModel model() {
            return state.model();
        }
    }


    /**
     * One phase's inner Newton step -- the per-phase half of the Hillert
     * iteration, solved once per phase per outer step.
     *
     * <p>Holds a {@link CVMGibbsModel} and drives it from the outside, exactly
     * as {@link CvmNewtonSolver} does for the fixed-composition minimisation.
     * The two differ only in which unknowns they solve for, and therefore in
     * which block of the same evaluated state they read:</p>
     *
     * <table border="1">
     *   <caption>The two solvers over one model</caption>
     *   <tr><th></th><th>CvmNewtonSolver</th><th>this</th></tr>
     *   <tr><td>composition</td><td>fixed constraint</td><td>an unknown</td></tr>
     *   <tr><td>reads</td><td>{@code gmu} / {@code gmuu} ({@code ncf})</td>
     *       <td>{@code guFull} (absolute {@code dG/dY}) / {@code gmuuFull}
     *       ({@code ncf+K}; {@code == guuFull} since {@code G0m} is linear in x)</td></tr>
     *   <tr><td>solves for</td><td>stationary mixing {@code Gm} at fixed x</td>
     *       <td>stationary absolute {@code G = G0m + Gm} relative to a trial mu;
     *       mu is the absolute (SER-referenced) chemical potential</td></tr>
     * </table>
     *
     * <p>Nested rather than standing alone because it is used by exactly one
     * caller -- {@link HillertSolver#solve} -- and is meaningless outside it:
     * a per-phase step expressed as an affine function of a chemical potential
     * has no use unless something downstream is solving for that potential.
     * Its {@link Step} output exists solely to be folded into
     * {@link EquilibriumMatrix}'s equations.</p>
     *
     * <p>It began life on {@code CVMGibbsModel} itself, on the reasoning that
     * its only inputs were that class's own widened derivatives. True, but it
     * made the model both an evaluator and two different solvers. The model
     * now evaluates; the solvers solve.</p>
     */
    public static final class PhaseStep {


        private final CVMGibbsModel model;

        public PhaseStep(CVMGibbsModel model) {
            if (model == null) {
                throw new IllegalArgumentException("model must not be null");
            }
            this.model = model;
        }

        /** The model this solver drives. */
        public CVMGibbsModel model() {
            return model;
        }

        /**
         * Result of {@link #step}: the joint Newton step
         * {@code deltaY(mu)}, expressed as an <b>affine function of the trial
         * chemical-potential vector {@code mu}</b> rather than a value at one
         * fixed {@code mu} -- see the note on {@link #step} for why
         * this shape, not a single numeric result, is what the outer Hillert
         * solver actually needs.
         *
         * <p>{@code deltaY(mu) = deltaY0 + Σ_k mu[k]*deltaYSensitivity[k]}, and
         * likewise for {@code deltaComposition}/{@code lambda}.</p>
         *
         * @param maxRelativeResidual the worst {@link LinearAlgebra.Solution#relativeResidual()}
         *        across the {@code K+1} linear solves this step required (one
         *        for {@code mu=0}, one per component's sensitivity column) --
         *        all share the same coefficient matrix, so this is a single
         *        ill-conditioning signal for the whole affine step, not just
         *        one of its solves. Near machine epsilon for a well-behaved
         *        state; orders of magnitude larger signals the widened
         *        Hessian's diagonal rescaling was not enough to protect this
         *        step's accuracy from round-off (see {@link LinearAlgebra#solveChecked}).
         */
        public record Step(
                double[] deltaY0, double[][] deltaYSensitivity,
                double[] deltaComposition0, double[][] deltaCompositionSensitivity,
                double lambda0, double[] lambdaSensitivity,
                double maxRelativeResidual) {

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
         * let {@link EquilibriumMatrix} fold the affine
         * form into its own equations, mirroring {@code genEqMat}'s substitution
         * entirely numerically.</p>
         *
         * <p>Unrelated to {@code CvmNewtonSolver}'s single-phase Newton-Raphson loop
         * (fixed composition, stationary mixing {@code Gm}) and must not be
         * confused with it: this solves for a stationary point of the
         * <b>absolute Gibbs energy {@code G = G0m + Gm}</b> <em>relative to a
         * trial {@code mu}</em>, with composition itself among the unknowns.</p>
         *
         * <p><b>Absolute stationarity, not mixing.</b> {@code Gu} is the
         * <em>absolute</em> widened gradient {@code state.guFull() = dG/dY}, so
         * the x-block residual is {@code dG/dx_i - mu_i - lambda}, with
         * {@code dG/dx_i = dGm/dx_i + G0_i(T)} carrying the pure-element
         * reference. Consequently the {@code mu} that {@link EquilibriumMatrix}
         * solves for -- which also uses the absolute {@code G = g()} in its
         * Gibbs-Duhem rows -- is the <b>absolute (SER-referenced) component
         * chemical potential</b>, directly comparable across phases and to a
         * CALPHAD {@code mu_i}. Using {@code gmuFull()} (mixing-only
         * {@code dGm/dY}) here while {@code EquilibriumMatrix} used absolute
         * {@code G} was an inconsistent hybrid; it was corrected to
         * {@code guFull()} after the Mo-Nb-Ta / 1000 K {@code phaseq} reference
         * chemical potentials were shown to be reproduced by {@code guFull()} to
         * ~0.03 J/mol and not by {@code gmuFull()} (see
         * {@code org.ce.scratch.HillertAbsoluteMuValidation}).</p>
         *
         * <p><b>The Hessian is unchanged by this.</b> {@code G0m} is linear in
         * composition, so {@code d2G0m/dY2 = 0} and
         * {@code state.guuFull() == state.gmuuFull()} exactly -- the reference
         * term contributes to the gradient's x-block only, never to {@code Guu}.</p>
         *
         * <p><b>The linear system</b> (at fixed T/P, so the {@code GxT*ΔT} and
         * {@code GxP*ΔP} terms vanish): unknowns are {@code deltaY[0..ncf+K-1]}
         * and {@code lambda}, over {@code ncf+K+1} equations:</p>
         * <ul>
         *   <li>Rows {@code 0..ncf-1} (u-block): {@code Guu[i,:] . deltaY = -Gu[i]}
         *       -- ordinary stationarity on the internal CFs, unconstrained by
         *       {@code mu}. (Here {@code Gu[i] = dGm/du_i} since {@code G0m} does
         *       not depend on {@code u}, so this block is the same whether the
         *       gradient is absolute or mixing.)</li>
         *   <li>Rows {@code ncf..ncf+K-1} (x-block): {@code Guu[i,:] . deltaY -
         *       lambda = mu[i-ncf] - Gu[i]} with {@code Gu[i] = dG/dx_{i-ncf}}
         *       (absolute) -- the only rows where {@code mu} appears, always with
         *       coefficient exactly {@code +1} on its own row, which is why one
         *       basis solve per component suffices.</li>
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
        public Step step(double[] uFull, double temperature) {
            int ncf = model.ncf();
            int numComponents = model.numComponents();
            int width = ncf + numComponents;
            if (uFull.length != width) {
                throw new IllegalArgumentException(
                        "uFull.length=" + uFull.length + " != ncf+K=" + width);
            }

            CVMGibbsModel.State state = model.atFull(temperature, uFull);
            // Absolute gradient dG/dY = d(G0m + Gm)/dY -- the x-block carries the
            // pure-element reference term dG0m/dx_i = G0_i(T), so mu solved
            // against it is the absolute (SER-referenced) component chemical
            // potential, consistent with EquilibriumMatrix's absolute G = g().
            // Using gmuFull() (mixing only) here was an inconsistent hybrid; see
            // the step() Javadoc.
            double[] Gu = state.guFull();
            // Hessian is unchanged: d2G0m/dY2 = 0 (G0m is linear in x), so
            // guuFull() == gmuuFull() exactly and the reference term contributes
            // nothing to the KKT matrix.
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
            LinearAlgebra.Solution sol0Checked = LinearAlgebra.solveChecked(A, b0);
            double[] sol0 = sol0Checked.x();
            double maxRelativeResidual = sol0Checked.relativeResidual();

            // Solving A*z = e_{ncf+k} directly gives d(deltaY)/d(mu_k), since the
            // system is linear and A is shared across right-hand sides.
            double[][] deltaYSens = new double[numComponents][];
            double[] lambdaSens = new double[numComponents];
            double[][] deltaCompSens = new double[numComponents][];
            for (int k = 0; k < numComponents; k++) {
                double[] ek = new double[n];
                ek[ncf + k] = 1.0;
                LinearAlgebra.Solution solKChecked = LinearAlgebra.solveChecked(A, ek);
                double[] solK = solKChecked.x();
                maxRelativeResidual = Math.max(maxRelativeResidual, solKChecked.relativeResidual());
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

            return new Step(deltaY0, deltaYSens, deltaComposition0, deltaCompSens, lambda0, lambdaSens,
                    maxRelativeResidual);
        }


    }

    /**
     * Outer-loop mass-balance + chemical-potential-equality assembly for
     * Hillert's method — port of the reference Mathematica implementation's
     * {@code genEqMat}, as actually called from {@code phaseq}'s outer
     * iteration (not a reconstruction from first principles; the full
     * {@code genEqMat} source was obtained and ported directly).
 *
     * <p><b>Where this sits relative to {@link PhaseStep}
     * (HILLERT_SOLVER_PLAN.md):</b> {@code PhaseStep.step} solves one
     * phase's local Newton step as an <em>affine function of the trial
     * chemical potential {@code mu}</em>, not at one fixed numeric {@code mu}
     * — tracing {@code phaseq}'s actual outer loop showed {@code delxGCVM} is
     * called with {@code mu} still symbolic, and {@code genEqMat}'s mass-
     * balance equations substitute that symbolic {@code delnN} in directly
     * before solving, so {@code mu} and {@code deltaN} are solved
     * <em>simultaneously</em> in one combined system. This class is the
     * numeric equivalent of that symbolic substitution: it takes each phase's
     * affine {@code deltaComposition(mu) = deltaComposition0 +
     * Σ_k mu[k]*deltaCompositionSensitivity[k]} and folds the {@code mu}
     * coefficients directly into the mass-balance rows' left-hand side,
     * producing one {@code (K+np)×(K+np)} system solved for {@code mu} and
     * {@code deltaN} together — exactly mirroring {@code genEqMat}'s
     * substitution, with no symbolic algebra.</p>
 *
     * <p><b>The linear system</b>, from {@code genEqMat}'s {@code gExpr}/
     * {@code nExpr} (fixed T, P, so the {@code GT*ΔT}/{@code GP*ΔP} terms
     * vanish): unknowns are {@code mu[0..K-1]} (chemical potentials) and
     * {@code deltaN[0..np-1]} (phase amount changes), {@code K+np} equations:</p>
     * <ul>
     *   <li><b>Gibbs-Duhem rows</b> (one per phase, {@code np} rows):
     *       {@code sum_i mu[i]     * composition[phase][i] = G[phase]} — a phase's
     *       Gibbs energy must equal the composition-weighted sum of chemical
     *       potentials, the standard partial-molar-quantity relation. {@code G}
     *       doesn't depend on {@code mu}, so these rows are unaffected by the
     *       affine substitution below.</li>
     *   <li><b>Mass-balance rows</b> (one per component, {@code K} rows) --
     *       the Newton linearization of the nonlinear inventory constraint
     *       {@code F_i(N,Y) = sum_p N_p x^p_i - b_i = 0} (STEP 10). Linearizing
     *       {@code F_i} about the current iterate:
     *       {@code F_i + sum_p x^p_i deltaN_p + sum_p N_p deltaComposition^p_i = 0},
     *       i.e. {@code sum_p x^p_i deltaN_p + sum_p N_p deltaComposition^p_i(mu)
     *       = r_i}, where {@code r_i = b_i - sum_p N_p x^p_i} is the
     *       <b>current mass-balance residual</b> (zero only when the iterate
     *       already satisfies the target -- see the {@code target} parameter of
     *       {@link #solve(List, double[])}). Substituting the affine form and
     *       collecting {@code mu}'s coefficients onto the left-hand side:
     *       {@code sum_phase deltaN[phase]*x[phase][i]
     *       + Σ_k mu[k]*(Σ_phase amount[phase]*deltaCompositionSensitivity[phase][k][i])
     *       = r_i - Σ_phase amount[phase]*deltaComposition0[phase][i]}.
     *       <b>When {@code r_i = 0} this is bit-identical to the pre-STEP-10 RHS
     *       {@code -Σ N c0}</b>; the {@code +r_i} term only bites when the
     *       represented inventory has drifted off (or been perturbed from) the
     *       target -- e.g. after a phase removal. Only the RHS changes; every
     *       matrix coefficient is exactly as before, and {@link PhaseStep} is
     *       untouched.</li>
     * </ul>
 *
     * <p><b>Why this is not on a model or on {@link PhaseStep}.</b> Its inputs are
     * aggregated quantities from <em>several</em> phases at once, not one model's
     * own state, so no per-phase object is a natural owner. It belongs to the outer
     * loop, which is the only thing that sees every phase at once -- hence its
     * place here rather than beside the model.</p>
 */
    public static final class EquilibriumMatrix {


        private EquilibriumMatrix() {}

        /**
         * One phase's contribution to the outer system: current amount,
         * composition, {@code G}, and its affine {@code deltaComposition(mu)}
         * (from {@link PhaseStep#step}) --
         * {@code deltaComposition0} is the {@code mu=0} value,
         * {@code deltaCompositionSensitivity[k]} is {@code d(deltaComposition)/d(mu_k)}.
         */
        public record PhaseContribution(
                double amount, double[] composition, double g,
                double[] deltaComposition0, double[][] deltaCompositionSensitivity) {}

        /**
         * Solved outer step: updated chemical potentials and each phase's amount
         * change, plus the backward error of the outer linear solve.
         *
         * @param relativeResidual {@code ||A z - b|| / ||b||} of the outer
         *        {@code (K+np)} system, from {@link LinearAlgebra#solveChecked}.
         *        Near machine epsilon for any physically distinct phase set (the
         *        Jacobi scaling inside {@code solveChecked} absorbs the O(1e4)
         *        magnitude spread between the O(1) composition rows and the
         *        O(1e-5) {@code N*sensitivity} rows). Orders of magnitude larger
         *        means the Gibbs-Duhem rows are near-parallel -- two stable
         *        phases with nearly identical composition, i.e. an almost
         *        rank-deficient (physically near-redundant) phase set. An
         *        exactly redundant phase set makes the outer matrix singular and
         *        {@code solveChecked} throws instead; {@link HillertSolver#solve}
         *        catches that and reports a non-converged {@link Result}.
         */
        public record EquilibriumStepResult(double[] mu, double[] deltaN, double relativeResidual) {}

        /**
         * Solves the combined outer Gibbs-Duhem + mass-balance system for one
         * outer iteration, with each phase's {@code mu}-dependence folded in
         * algebraically (see class doc).
         *
         * @param phases K-component contributions from every phase currently
         *               treated as stable (the reference's {@code
         *               unStablePhaseRules} excludes non-positive-amount phases
         *               from this assembly entirely — callers must filter
         *               before calling, not pass all phases and expect
         *               filtering here)
         * @param numComponents K, the number of system components
         * @param target the prescribed component inventory {@code b_i} (moles),
         *               length K. The mass-balance rows are the Newton
         *               linearization of {@code sum_p N_p x^p_i = b_i}, so their
         *               RHS carries the current residual
         *               {@code r_i = b_i - sum_p N_p x^p_i} (STEP 10). Pass
         *               {@code null} to keep the pre-STEP-10 behaviour
         *               ({@code r_i} treated as 0) -- used only by callers that
         *               do not prescribe an inventory.
         */
        public static EquilibriumStepResult solve(
                List<PhaseContribution> phases, int numComponents, double[] target) {
            int np = phases.size();
            int n = numComponents + np; // unknowns: mu[0..K-1], deltaN[0..np-1]

            double[][] A = new double[n][n];
            double[] b = new double[n];

            // Gibbs-Duhem rows (0..np-1): sum_i mu[i]*x[phase][i] = G[phase]
            // (G has no mu-dependence, so this block is unaffected by the affine substitution below)
            for (int p = 0; p < np; p++) {
                PhaseContribution phase = phases.get(p);
                for (int i = 0; i < numComponents; i++) {
                    A[p][i] = phase.composition()[i];
                }
                b[p] = phase.g();
            }

            // Mass-balance rows (np..np+K-1): Newton linearization of
            // F_i = sum_p N_p x^p_i - b_i = 0. RHS = r_i - sum_p N_p c0_p[i],
            // r_i = b_i - sum_p N_p x^p_i (STEP 10). When the current iterate
            // already meets the target (r_i = 0) this equals the pre-STEP-10 RHS
            // -sum_p N_p c0_p[i] exactly. The matrix coefficients are unchanged:
            // mu's coefficients on the LHS (columns 0..K-1), deltaN's on
            // columns K..K+np-1.
            for (int i = 0; i < numComponents; i++) {
                int row = np + i;
                double represented = 0.0;
                double rhs = 0.0;
                for (int p = 0; p < np; p++) {
                    PhaseContribution phase = phases.get(p);
                    A[row][numComponents + p] = phase.composition()[i];
                    for (int k = 0; k < numComponents; k++) {
                        A[row][k] += phase.amount() * phase.deltaCompositionSensitivity()[k][i];
                    }
                    rhs -= phase.amount() * phase.deltaComposition0()[i];
                    represented += phase.amount() * phase.composition()[i];
                }
                double ri = (target != null) ? target[i] - represented : 0.0;
                b[row] = ri + rhs;
            }

            // solveChecked, not solve: same algorithm and same exact answer
            // (solve() is literally solveChecked().x()), but it also returns the
            // backward error so the outer loop can flag a near-degenerate phase
            // set the same way PhaseStep already flags an ill-conditioned inner
            // step. An exactly singular matrix (two stable phases with identical
            // composition and G) still throws IllegalArgumentException here;
            // HillertSolver.solve catches it.
            LinearAlgebra.Solution checked = LinearAlgebra.solveChecked(A, b);
            double[] sol = checked.x();

            double[] mu = new double[numComponents];
            System.arraycopy(sol, 0, mu, 0, numComponents);
            double[] deltaN = new double[np];
            System.arraycopy(sol, numComponents, deltaN, 0, np);

            return new EquilibriumStepResult(mu, deltaN, checked.relativeResidual());
        }

    }

    /**
     * Runs the outer/inner Hillert iteration to equilibrium.
     *
     * @param phases per-phase mutable state (updated in place on the accepted-step path)
     * @param overallAmounts the system's conserved component inventory:
     *            {@code overallAmounts[i]} = total moles of component {@code i},
     *            length {@code K}, every entry {@code >= 0} and finite, with
     *            {@code sum_i overallAmounts[i] > 0}. <b>Not normalized</b> --
     *            pass a vector summing to 1 for a unit-total system; the total
     *            system amount is {@code sum_i overallAmounts[i]} and is not
     *            independently fixed. Same units as {@link Phase#amount} (moles
     *            of formula units; for the supported single-site
     *            one-atom-per-site disordered phases, moles of atoms). The seeded
     *            phase states must already represent this target
     *            ({@code sum_p N_p x^p_i == overallAmounts[i]} within the
     *            mass-balance tolerance) or the run is rejected with
     *            {@link ConvergenceReason#INITIAL_MASS_BALANCE} -- no initial
     *            amount or composition is silently changed. The Newton
     *            mass-balance equations are <b>unchanged</b>; the target enters
     *            only validation and the independent residual report.
     * @param temperature fixed temperature (K) -- GxT/GxP not yet supported, v1 is fixed-T,P
     * @param maxOuterIterations outer Newton iteration cap
     * @param innerBacktrackTries max lambda-halving tries per outer iteration
     * @param tol convergence tolerance on {@code max_{p active} ||deltaY_p(mu)||}
     *            (see the class note on why {@code max}, not {@code phaseq}'s {@code min})
     */
    public static Result solve(
            List<Phase> phases,
            double[] overallAmounts,
            double temperature,
            int maxOuterIterations,
            int innerBacktrackTries,
            double tol,
            Consumer<String> progressSink) {

        int numComponents = phases.get(0).numComponents;
        double[] mu = new double[numComponents];

        // ---- Mass-balance target and its two tolerances (STEP 9) ----
        // Mass balance is in moles; the phase-step tol is dimensionless, so the
        // two are not the same quantity and are not reused for each other.
        //
        // ENTRY tolerance (massTolEntry = massRelEntry * massScale, massRelEntry
        // = 1e-9): the seed is exact by construction, so this can be tight --
        // seven orders below the O(0.1-1) phase-amount scale (a genuine caller
        // error, seeding phases that represent a different overall composition,
        // is off by O(0.01)+ and always caught) and seven orders above the
        // ~1e-16 round-off of the sum_p N_p x^p_i dot product (a correctly
        // specified problem is never falsely rejected).
        //
        // EXIT / CONVERGENCE-GATE tolerance (massRelGate = max(1e-7, 10*tol)):
        // STEP 10 made the mass-balance constraint part of the Newton system --
        // its rows now carry the current residual r_i = b_i - sum_p N_p x^p_i,
        // so each step actively drives r_i -> 0 (quadratically near a solution,
        // ~(1-lambda) per damped step) instead of merely "holding" the drifted
        // value. The residual is therefore a genuine component of what the
        // iteration solves, and a HARD part of the CONVERGED gate: a converged
        // Result must have maxRelResidual <= massRelGate. This is an order
        // tighter than STEP-9's post-hoc massRelExit (the correction earns it),
        // still loose enough for the linearization's per-step quadratic residual
        // and the damped-step tail. If the target cannot be represented by the
        // active phase set, the iteration cannot reduce r_i -> it ends
        // non-converged with a large, clearly reported residual (MASS_BALANCE_DRIFT).
        final double massRelEntry = 1.0e-9;
        final double massRelGate = Math.max(1.0e-7, 10.0 * tol);
        double[] target = overallAmounts == null ? null : overallAmounts.clone();
        double targetSum = 0.0;
        boolean targetValid = target != null && target.length == numComponents;
        if (targetValid) {
            for (double v : target) {
                if (!Double.isFinite(v) || v < 0.0) { targetValid = false; break; }
                targetSum += v;
            }
            if (targetSum <= 0.0) targetValid = false;
        }
        final double massScale = Math.max(targetSum, 1.0);
        final double massTolEntry = massRelEntry * massScale;

        if (!targetValid) {
            if (progressSink != null) {
                progressSink.accept("Hillert: invalid overallAmounts target "
                        + (target == null ? "(null)"
                           : "(length=" + target.length + ", expected " + numComponents
                             + "; must be finite, >= 0, sum > 0)")
                        + " -- rejecting problem, not iterating.");
            }
            return rejectedResult(phases, temperature, mu, ConvergenceReason.INITIAL_MASS_BALANCE,
                    massBalance(phases, target, targetSum, massScale));
        }

        // Initial-state mass balance: reconstruct sum_p N_p x^p_i from the seed
        // phase states (independent of any Newton quantity) and compare to target.
        MassBalanceReport initialMb = massBalance(phases, target, targetSum, massScale);
        if (initialMb.maxAbsResidual() > massTolEntry) {
            if (progressSink != null) {
                progressSink.accept(String.format(
                        "Hillert: seeded phases do not represent overallAmounts "
                        + "(max |target_i - sum_p N_p x^p_i| = %.3e > massTolEntry %.3e) "
                        + "-- rejecting problem; fix the seed, no silent correction.",
                        initialMb.maxAbsResidual(), massTolEntry));
            }
            return rejectedResult(phases, temperature, mu,
                    ConvergenceReason.INITIAL_MASS_BALANCE, initialMb);
        }

        // A step below this fraction of the joint-state scale is not progress:
        // 16 ulps is "no change" in double precision; 1e-3*tol is a step three
        // orders below the convergence tolerance, which cannot plausibly carry
        // the solve across that tolerance in a sane number of iterations.
        final double stallRel = Math.max(16.0 * Math.ulp(1.0), 1.0e-3 * tol);

        // Trustworthy-outer-solve bound: reuse the STEP-5 backward-error signal,
        // loose so it never blocks a genuinely converged run.
        final double outerResidualBound = Math.max(1.0e-8, tol);

        // ---- Phase-set management (STEP 8) ----
        // Remove an active phase once its amount is numerically zero relative to
        // the total represented amount. REMOVE_REL = 1e-9: seven orders above
        // the outer solve's deltaN round-off floor (~1e-16 backward error, so it
        // can never be tripped by noise in a converging phase) and seven orders
        // below the O(0.1) lever-rule fraction a genuinely stable phase holds;
        // also the order of the tightest tol any caller uses. A phase reaches
        // this bound geometrically (STEP 7's line search leaves a shrinking
        // phase at >= 50% of the remaining distance to zero each step), so the
        // crossing is a bounded, deterministic ~33 iterations from N ~ N_total/2
        // -- not chatter. MIN_RESIDENCE = 3 accepted iterations before a phase
        // may be removed, so a transient dip that recovers is not acted on.
        // Removing a zero-amount phase changes the represented total by exactly
        // N*x = 0, so removal is exactly mass-conserving with respect to the
        // explicit overallAmounts target -- removal alone never perturbs it, so
        // it needs no help from the STEP-10 target-aware rows.
        //
        // Phase ADDITION (STEP 12): after the removal block, an inactive
        // candidate whose relaxed absolute driving force
        // dGf_beta = sum_i mu_i x^beta_i - G^beta exceeds addThreshold is
        // activated with amount epsilon = 1e-6 * sum(overallAmounts). The
        // inserted mass is NOT hand-redistributed -- STEP 10's target-aware
        // mass-balance rows carry r_i = -epsilon*x^beta_i on the next iteration
        // and drive the redistribution (STEP 11 proved this on an analytic
        // two-phase system). At most one phase per outer iteration; the
        // iteration that adds is skipped for convergence, as for a removal.
        final double removeRel = 1.0e-9;
        final int minResidence = 3;
        final double epsilonAdd = 1.0e-6 * targetSum;   // initial amount for an inserted phase
        List<PhaseSetEvent> phaseSetEvents = new ArrayList<>();

        ConvergenceReason reason = ConvergenceReason.MAX_ITERATIONS;
        int iterationsRun = 0;
        double lastMaxStepNorm = Double.POSITIVE_INFINITY;
        double lastMaxAmountStep = Double.NaN;
        double lastLinearResidual = Double.NaN;
        double lastLambda = Double.NaN;
        boolean lastAccepted = false;
        int consecutiveNegligible = 0;
        double lastMassResidualBefore = Double.NaN;   // r_i at the start of the last iteration (STEP 13)

        for (int outerIter = 1; outerIter <= maxOuterIterations; outerIter++) {
            iterationsRun = outerIter;

            // STEP 12 anti-oscillation: age every phase's phase-set cooldown by
            // one outer iteration. A phase added or removed within the last
            // MIN_RESIDENCE iterations is not a candidate for addition.
            for (Phase phase : phases) {
                if (phase.phaseSetCooldown > 0) {
                    phase.phaseSetCooldown--;
                }
            }

            List<PhaseStep.Step> steps = new ArrayList<>(phases.size());
            for (int p = 0; p < phases.size(); p++) {
                Phase phase = phases.get(p);
                PhaseStep.Step step = new PhaseStep(phase.model).step(phase.uFull, temperature);
                steps.add(step);
                // See LinearAlgebra.solveChecked: a widened Hessian near a
                // dilute/near-boundary composition can span 10+ orders of
                // magnitude on its diagonal and still clear the elimination's
                // singularity guard while losing real accuracy to round-off.
                // Surface that here rather than let it silently degrade the
                // Newton step -- this is exactly the failure mode that made
                // this solver's trace diverge from the reference phaseq port
                // on such compositions.
                if (step.maxRelativeResidual() > 1e-6 && progressSink != null) {
                    progressSink.accept(String.format(
                            "Hillert outer iter %d, phase '%s': ill-conditioned Newton step "
                            + "(relative residual %.3e) -- this state's widened Hessian may be "
                            + "too badly scaled for a reliable step here.",
                            outerIter, phase.label, step.maxRelativeResidual()));
                }
            }

            // Active phase set for this iteration, with an explicit
            // phase-index <-> active-slot mapping (STEP 10): rebuilt from
            // phase.active every iteration so a removal cannot corrupt it, and
            // deterministic (ascending phase index). activeIndices.get(slot) is
            // the phase index of active slot `slot`; activeSlotOf[p] is p's slot
            // or -1 if inactive.
            List<EquilibriumMatrix.PhaseContribution> contributions = new ArrayList<>();
            List<Integer> activeIndices = new ArrayList<>();
            int[] activeSlotOf = new int[phases.size()];
            java.util.Arrays.fill(activeSlotOf, -1);
            for (int p = 0; p < phases.size(); p++) {
                Phase phase = phases.get(p);
                if (!phase.active) {
                    continue;
                }
                PhaseStep.Step step = steps.get(p);
                double g = currentG(phase, temperature);
                contributions.add(new EquilibriumMatrix.PhaseContribution(
                        phase.amount, phase.composition(), g,
                        step.deltaComposition0(), step.deltaCompositionSensitivity()));
                activeSlotOf[p] = activeIndices.size();
                activeIndices.add(p);
            }

            // Mass-balance residual r_i = b_i - sum_{p active} N_p x^p_i at the
            // CURRENT iterate -- reconstructed independently from the phase
            // states, the "before-step" residual (STEP 13 diagnostic). Fed to
            // EquilibriumMatrix.solve as the Newton RHS so the step drives it
            // toward zero (STEP 10). When it is already ~0 the outer system is
            // bit-identical to the pre-STEP-10 formulation.
            MassBalanceReport beforeStepMb = massBalance(phases, target, targetSum, massScale);
            lastMassResidualBefore = beforeStepMb.maxAbsResidual();

            EquilibriumMatrix.EquilibriumStepResult eqStep;
            try {
                eqStep = EquilibriumMatrix.solve(contributions, numComponents, target);
            } catch (RuntimeException ex) {
                // Singular outer matrix: two stable phases share one composition
                // (and G), so their Gibbs-Duhem rows are identical -- a genuinely
                // rank-deficient, physically redundant phase set. There is
                // nothing to iterate towards (the phase set will not change on
                // its own; that is phase-set management, not this loop's job), so
                // stop and return a non-converged Result from the last good
                // state rather than letting an unchecked exception abort the
                // whole solve with no Result at all.
                if (progressSink != null) {
                    progressSink.accept("Hillert outer iter " + outerIter
                            + ": outer equilibrium matrix is singular (" + ex.getMessage()
                            + ") -- likely two stable phases with the same composition; "
                            + "stopping, result is non-converged.");
                }
                reason = ConvergenceReason.SINGULAR_OUTER_SYSTEM;
                break;
            }

            // N1: solveChecked's pivot guard only catches an EXACTLY singular
            // matrix; a NaN/Inf coefficient (from a phase whose widened
            // derivatives went non-finite) sails through elimination and yields
            // a NaN mu vector with a NaN relativeResidual. Catch that here
            // BEFORE assigning `mu`, so the returned Result keeps the last
            // finite mu instead of presenting NaN/Inf as a solver state, and
            // stop the loop without attempting a PhaseStep or line search.
            if (!isFinite(eqStep.mu()) || !Double.isFinite(eqStep.relativeResidual())) {
                if (progressSink != null) {
                    progressSink.accept("Hillert outer iter " + outerIter
                            + ": outer equilibrium solve returned a non-finite result "
                            + "(mu or backward-error residual is NaN/Inf) -- numerical breakdown; "
                            + "stopping, result is non-converged, mu left at its last finite value.");
                }
                reason = ConvergenceReason.NUMERICAL_BREAKDOWN;
                lastAccepted = false;
                break;
            }

            mu = eqStep.mu();
            double[] deltaN = eqStep.deltaN();
            lastLinearResidual = eqStep.relativeResidual();
            // Same ill-conditioning signal PhaseStep already surfaces for the
            // inner step: a large backward error means the Gibbs-Duhem rows are
            // near-parallel (two stable phases with nearly identical
            // composition), so mu/deltaN from this solve may be unreliable.
            if (eqStep.relativeResidual() > 1e-6 && progressSink != null) {
                progressSink.accept(String.format(
                        "Hillert outer iter %d: ill-conditioned outer equilibrium solve "
                        + "(relative residual %.3e) -- two stable phases may be nearly "
                        + "identical in composition.",
                        outerIter, eqStep.relativeResidual()));
            }

            // ---- Backtracking line search (STEP 7 feasibility hardening) ----
            // Only ACTIVE phases (activeIndices) are stepped and
            // validity-checked. An inactive phase was excluded from the outer
            // solve, so stepping its uFull against the shared mu is meaningless
            // and its trial state must not veto the active phases' step.
            // Inactive phases are frozen (uFull, amount untouched).
            //
            // Feasibility of an accepted trial now requires, in addition to the
            // reference's CV validity check:
            //   (i)  every stepped phase's trial amount N + lambda*deltaN >= 0
            //        and finite -- the model needs N >= 0 and the old code never
            //        checked it (committed negatives);
            //   (ii) the trial actually moves the joint state by more than the
            //        machine/tolerance floor -- a "valid" step of relative size
            //        ~1e-13 is not progress (see the class note on the STEP-6
            //        stall threshold, reused here per-attempt).
            //
            // The lambda=1,1/2,1/4,... cadence is unchanged; only its starting
            // value is capped so the first probe is already amount-feasible:
            //   lambda_amount = min over stepped p with deltaN<0 of (-N/deltaN)
            //   lambda_start   = min(1, 0.5 * lambda_amount)
            // This is a single linear ratio per phase -- it has none of the
            // dilute-composition stalling of an analytic cluster-variable bound
            // (that approach was tried and rejected; see HILLERT_SOLVER_PLAN.md).
            double lambdaAmount = Double.POSITIVE_INFINITY;
            for (int slot = 0; slot < activeIndices.size(); slot++) {
                double dN = deltaN[slot];
                if (dN < 0.0) {
                    double n = phases.get(activeIndices.get(slot)).amount;
                    lambdaAmount = Math.min(lambdaAmount, -n / dN);
                }
            }
            double lambda = Math.min(1.0, 0.5 * lambdaAmount);

            double maxUFullScale = 0.0;
            double preStepMaxNorm = 0.0;
            for (int slot = 0; slot < activeIndices.size(); slot++) {
                int p = activeIndices.get(slot);
                maxUFullScale = Math.max(maxUFullScale, l2Norm(phases.get(p).uFull));
                preStepMaxNorm = Math.max(preStepMaxNorm, l2Norm(steps.get(p).deltaYAt(mu)));
            }
            double stateScale = Math.max(maxUFullScale, 1.0);
            // Near the fixed point the Newton step itself is tiny, so lambda*deltaY
            // is legitimately below stallRel -- that is convergence, not a stall.
            // Only enforce the "step must be meaningful" gate when we are still
            // short of the convergence tolerance in step norm.
            boolean enforceMeaningfulStep = preStepMaxNorm > tol;

            double[][] trialUFull = new double[phases.size()][];
            double[] trialAmount = new double[phases.size()];
            boolean accepted = false;
            boolean feasibleButNegligibleSeen = false;
            double maxRelStateChange = 0.0;

            for (int tries = 0; tries < innerBacktrackTries; tries++) {
                boolean allFeasible = true;
                double relChange = 0.0;
                for (int slot = 0; slot < activeIndices.size(); slot++) {
                    int p = activeIndices.get(slot);
                    Phase phase = phases.get(p);
                    double[] deltaY = steps.get(p).deltaYAt(mu);
                    double[] u = new double[phase.uFull.length];
                    for (int i = 0; i < u.length; i++) {
                        u[i] = phase.uFull[i] + lambda * deltaY[i];
                    }
                    trialUFull[p] = u;
                    trialAmount[p] = phase.amount + lambda * deltaN[slot];

                    relChange = Math.max(relChange, l2Norm(scaled(deltaY, lambda)));

                    if (!(trialAmount[p] >= 0.0) || !Double.isFinite(trialAmount[p])) {
                        allFeasible = false; // negative or non-finite phase amount
                    }
                    double[] uOnly = java.util.Arrays.copyOfRange(u, 0, phase.ncf);
                    double[] xOnly = java.util.Arrays.copyOfRange(u, phase.ncf, u.length);
                    if (!phase.model.at(temperature, xOnly, uOnly).isValidIncludingPoints()) {
                        allFeasible = false;
                    }
                }
                relChange /= stateScale;

                if (allFeasible) {
                    if (!enforceMeaningfulStep || relChange >= stallRel) {
                        accepted = true;
                        maxRelStateChange = relChange;
                        break;
                    }
                    // A feasible step exists but, while still short of
                    // convergence, it is physically negligible -- do not accept
                    // it; keep halving in case a larger valid lambda is found,
                    // but remember this so a total failure here is classified as
                    // a stall, not an infeasibility.
                    feasibleButNegligibleSeen = true;
                }
                lambda *= 0.5;
            }

            // Reference (phaseq) has no special handling for exhausting all
            // itr backtracking tries without success: the state is left
            // unchanged and the outer loop proceeds. Port that literally -- do
            // not abort here. STEP 6 classifies the terminal outcome; STEP 7
            // adds the "every feasible step is negligible" -> stall distinction.
            if (accepted) {
                for (int slot = 0; slot < activeIndices.size(); slot++) {
                    int p = activeIndices.get(slot);
                    phases.get(p).uFull = trialUFull[p];
                    phases.get(p).amount = trialAmount[p];
                    phases.get(p).activeResidence++;
                }
            } else {
                if (feasibleButNegligibleSeen) {
                    // Force STEP-6's stall path immediately: a feasible step was
                    // available every attempt, it was just too small to matter.
                    consecutiveNegligible++;
                }
                if (progressSink != null) {
                    progressSink.accept("Hillert outer iter " + outerIter
                            + (feasibleButNegligibleSeen
                                    ? ": every feasible line-search step is numerically negligible; "
                                      + "state unchanged."
                                    : ": backtracking exhausted (no feasible step), state unchanged, "
                                      + "continuing."));
                }
            }

            // ---- Phase removal (STEP 8) ----
            // After an accepted step, an active phase whose amount has reached
            // numerical zero (relative to the total represented amount) and that
            // has been active long enough (MIN_RESIDENCE) is removed: amount set
            // to exact 0, frozen, marked inactive. This is a PHASE_REMOVED event
            // -- distinct from a solver stall. The next iteration rebuilds the
            // outer system on the reduced active set and the iteration
            // continues. Removing a zero-amount phase changes the represented
            // total by N*x = 0, so it is exactly mass-conserving.
            boolean removedThisIter = false;
            if (accepted) {
                double activeTotal = 0.0;
                for (int slot = 0; slot < activeIndices.size(); slot++) {
                    activeTotal += phases.get(activeIndices.get(slot)).amount;
                }
                double removeThreshold = removeRel * Math.max(activeTotal, 1.0);
                for (int slot = 0; slot < activeIndices.size(); slot++) {
                    int p = activeIndices.get(slot);
                    Phase phase = phases.get(p);
                    boolean shrinking = deltaN[slot] < 0.0;
                    boolean atZero = phase.amount <= 0.0
                            || (phase.amount < removeThreshold && shrinking);
                    if (atZero && phase.activeResidence >= minResidence
                            && activeIndices.size() >= 2) {
                        double old = phase.amount;
                        phase.amount = 0.0;
                        phase.active = false;
                        phase.activeResidence = 0;
                        phase.phaseSetCooldown = minResidence;   // STEP 12: no immediate re-add
                        removedThisIter = true;
                        phaseSetEvents.add(new PhaseSetEvent(
                                PhaseSetEventType.PHASE_REMOVED, phase.label, outerIter,
                                old, 0.0, Double.NaN));
                        if (progressSink != null) {
                            progressSink.accept(String.format(
                                    "Hillert outer iter %d: phase '%s' removed (amount %.3e -> 0, "
                                    + "threshold %.3e); continuing with %d active phase(s).",
                                    outerIter, phase.label, old, removeThreshold,
                                    activeIndices.size() - 1));
                        }
                    }
                }
                if (removedThisIter) {
                    // Removal IS progress -- do not let this iteration trip the
                    // stall detector or (with a still-nonzero step) declare
                    // convergence on the now-stale metric. Reset and re-solve
                    // the reduced system next iteration.
                    consecutiveNegligible = 0;
                    lastMaxStepNorm = Double.POSITIVE_INFINITY;
                    lastAccepted = true;
                    continue;
                }
            }

            // ---- Phase addition (STEP 12) ----
            // After an accepted step and the removal block, if the active set is
            // near its own fixed point and an inactive candidate exists, score
            // each inactive candidate by its relaxed absolute driving force
            // dGf = sum_i mu_i x^cand_i - G^cand (mu from the accepted step) and
            // add the single most-favourable one above addThreshold. The
            // candidate Phase is NOT mutated by the scan (relaxCandidate works
            // on copies); only the selected phase is committed. The inserted
            // mass is left for STEP 10's target-aware Newton rows to
            // redistribute next iteration -- no hand redistribution here.
            boolean addedThisIter = false;
            if (accepted && target != null) {
                double activeSetNorm = 0.0;
                boolean anyInactiveCandidate = false;
                for (int p = 0; p < phases.size(); p++) {
                    Phase phase = phases.get(p);
                    if (phase.active) {
                        activeSetNorm = Math.max(activeSetNorm,
                                l2Norm(steps.get(p).deltaYAt(mu)));
                    } else if (phase.phaseSetCooldown == 0) {
                        anyInactiveCandidate = true;
                    }
                }
                // Only chase candidates once the active set has (nearly) stopped
                // moving -- a good mu estimate is a precondition for a
                // meaningful driving force (STEP 11 PART 10 / PART K step 4).
                boolean activeSetSettled = activeSetNorm <= Math.max(1.0e-4, 100.0 * tol);

                if (anyInactiveCandidate && activeSetSettled) {
                    int bestIdx = -1;
                    double bestDrivingForce = 0.0;
                    CandidatePhaseState bestState = null;
                    int evaluated = 0;
                    for (int p = 0; p < phases.size(); p++) {
                        Phase phase = phases.get(p);
                        if (phase.active || phase.phaseSetCooldown > 0) {
                            continue;
                        }
                        CandidatePhaseState cand = relaxCandidate(phase, mu, temperature);
                        if (!cand.valid()) {
                            if (progressSink != null) {
                                progressSink.accept(String.format(
                                        "Hillert outer iter %d: inactive candidate '%s' could not be "
                                        + "evaluated (no converged relaxed state on the search grid) "
                                        + "-- not added.",
                                        outerIter, phase.label));
                            }
                            continue;
                        }
                        evaluated++;
                        double addThreshold = Math.max(1.0, 1.0e-6 * Math.abs(cand.gibbsEnergy()));
                        if (cand.drivingForce() > addThreshold) {
                            if (cand.drivingForce() > bestDrivingForce) {
                                bestIdx = p;
                                bestDrivingForce = cand.drivingForce();
                                bestState = cand;
                            }
                        } else if (progressSink != null) {
                            progressSink.accept(String.format(
                                    "Hillert outer iter %d: inactive candidate '%s' evaluated, "
                                    + "dGf=%.3f J/mol <= addThreshold %.3f -- not favourable, not added.",
                                    outerIter, phase.label, cand.drivingForce(), addThreshold));
                        }
                    }

                    if (bestIdx >= 0 && epsilonAdd > 0.0 && Double.isFinite(epsilonAdd)
                            && isFinite(bestState.uFull())) {
                        Phase beta = phases.get(bestIdx);
                        double old = beta.amount;   // 0 for an inactive phase
                        beta.uFull = bestState.uFull().clone();
                        beta.amount = epsilonAdd;
                        beta.active = true;
                        beta.activeResidence = 0;
                        beta.phaseSetCooldown = minResidence;
                        addedThisIter = true;
                        phaseSetEvents.add(new PhaseSetEvent(
                                PhaseSetEventType.PHASE_ADDED, beta.label, outerIter,
                                old, epsilonAdd, bestDrivingForce));
                        if (progressSink != null) {
                            progressSink.accept(String.format(
                                    "Hillert outer iter %d: phase '%s' added (amount 0 -> %.3e, "
                                    + "dGf=%.3f J/mol > threshold); now %d active phase(s). "
                                    + "STEP-10 mass balance will redistribute next iteration.",
                                    outerIter, beta.label, epsilonAdd, bestDrivingForce,
                                    activeIndices.size() + 1));
                        }
                    } else if (evaluated == 0 && progressSink != null) {
                        progressSink.accept(String.format(
                                "Hillert outer iter %d: no inactive candidate could be evaluated; "
                                + "convergence (if reached) is on the current active set only.",
                                outerIter));
                    }
                }
            }
            if (addedThisIter) {
                // Addition IS progress -- same handling as a removal: do not
                // let this iteration trip the stall detector or declare
                // convergence. Re-solve the widened system next iteration.
                consecutiveNegligible = 0;
                lastMaxStepNorm = Double.POSITIVE_INFINITY;
                lastAccepted = true;
                continue;
            }

            // Convergence metric: MAX (not phaseq's MIN) over the currently
            // stable phases' full joint Newton step -- see the class note.
            // activeIndices is in ascending, unique order, so its position is
            // exactly the phase's deltaN slot.
            double maxNorm = 0.0;
            double maxAmountStep = 0.0;
            double totalAmount = 0.0;
            for (int slot = 0; slot < activeIndices.size(); slot++) {
                int p = activeIndices.get(slot);
                maxNorm = Math.max(maxNorm, l2Norm(steps.get(p).deltaYAt(mu)));
                maxAmountStep = Math.max(maxAmountStep,
                        Math.abs((accepted ? lambda : 1.0) * deltaN[slot]));
                totalAmount += phases.get(p).amount;
            }
            double relAmountStep = totalAmount > 0.0 ? maxAmountStep / totalAmount : maxAmountStep;

            lastMaxStepNorm = maxNorm;
            lastMaxAmountStep = accepted ? relAmountStep : Double.NaN;
            lastLambda = accepted ? lambda : Double.NaN;
            lastAccepted = accepted;

            if (progressSink != null) {
                double afterStepMass = accepted
                        ? massBalance(phases, target, targetSum, massScale).maxAbsResidual()
                        : lastMassResidualBefore;
                progressSink.accept(String.format(
                        "Hillert outer iter %d: lambda=%s accepted=%s maxNorm=%.3e outerResid=%.3e "
                        + "massResid %.3e -> %.3e",
                        outerIter, accepted ? Double.toString(lambda) : "n/a",
                        accepted, maxNorm, lastLinearResidual,
                        lastMassResidualBefore, afterStepMass));
            }

            // Numerical stall. Two sources feed consecutiveNegligible:
            //   - STEP 7's line search already incremented it above when every
            //     feasible trial this iteration was below stallRel (a feasible
            //     step exists but is physically meaningless);
            //   - an ACCEPTED step whose relative state change is below stallRel
            //     while still short of convergence. (With STEP 7's per-attempt
            //     stallRel gate this second case can no longer occur -- an
            //     accepted step is >= stallRel by construction -- but the check
            //     is kept as a defensive backstop.)
            // Either way: two consecutive such iterations -> STALLED.
            boolean negligibleThisIter =
                    (accepted && maxNorm > tol && maxRelStateChange < stallRel)
                    || (!accepted && feasibleButNegligibleSeen);
            if (negligibleThisIter) {
                if (accepted && maxNorm > tol && maxRelStateChange < stallRel) {
                    consecutiveNegligible++; // STEP 7's failed-search path already counted
                }
                if (consecutiveNegligible >= 2 && maxNorm > tol) {
                    reason = ConvergenceReason.STALLED;
                    if (progressSink != null) {
                        progressSink.accept("Hillert outer iter " + outerIter
                                + ": stalled (no meaningful progress possible for two iterations "
                                + "while maxNorm=" + maxNorm + " > tol); stopping, non-converged.");
                    }
                    break;
                }
            } else {
                consecutiveNegligible = 0;
            }

            // Convergence requires: a step was actually applied this iteration
            // (a failed line search is not convergence), the max stable-phase
            // step norm is below tol, and the outer solve was trustworthy.
            if (accepted && maxNorm <= tol && lastLinearResidual <= outerResidualBound) {
                reason = ConvergenceReason.CONVERGED;
                break;
            }
        }

        // If the loop ran to the cap (or ended) with the last executed iteration
        // failing its line search, that is a distinct non-converged condition,
        // not "max iterations".
        if (reason == ConvergenceReason.MAX_ITERATIONS && !lastAccepted) {
            reason = ConvergenceReason.LINE_SEARCH_FAILED;
        }

        // Final mass-balance check (STEP 10). The mass-balance constraint is now
        // part of the Newton system (its rows carry r_i, driven to zero), so
        // maxRelResidual <= massRelGate is a HARD part of the CONVERGED gate,
        // not an after-the-fact downgrade. If the step-norm criterion was met
        // but r_i is still above massRelGate, the represented inventory cannot
        // be reconciled with the target by the current active phase set (e.g. a
        // phase was removed and no single-phase equilibrium at the target
        // exists) -- report MASS_BALANCE_DRIFT, an explicit non-converged
        // condition. The residual is reconstructed independently from the
        // accepted phase states, never from a Newton quantity.
        MassBalanceReport rawMb = massBalance(phases, target, targetSum, massScale);
        MassBalanceReport finalMb = new MassBalanceReport(
                rawMb.targetOverall(), rawMb.calculatedOverall(),
                rawMb.maxAbsResidual(), rawMb.maxRelResidual(), lastMassResidualBefore);
        if (reason == ConvergenceReason.CONVERGED && finalMb.maxRelResidual() > massRelGate) {
            reason = ConvergenceReason.MASS_BALANCE_DRIFT;
            if (progressSink != null) {
                progressSink.accept(String.format(
                        "Hillert: step-norm converged but the mass-balance residual r_i is still "
                        + "%.3e relative (> massRelGate %.3e) -- the target is not representable by "
                        + "the current active phase set; reporting MASS_BALANCE_DRIFT.",
                        finalMb.maxRelResidual(), massRelGate));
            }
        }

        boolean converged = reason == ConvergenceReason.CONVERGED;

        // Evaluate each phase once at its final point and keep the state, so a
        // caller can read any further property (entropy, gradients, SRO) on
        // demand rather than re-evaluating. G comes from that same state, so
        // the reported energy and anything derived from it cannot disagree.
        boolean allValid = true;
        List<PhaseResult> entries = new ArrayList<>();
        for (Phase phase : phases) {
            CVMGibbsModel.State state = phase.model.atFull(temperature, phase.uFull);
            allValid &= state.isValidIncludingPoints();
            entries.add(new PhaseResult(
                    phase.label, phase.amount, phase.composition(), state.g(), state, converged));
        }

        ConvergenceReport report = new ConvergenceReport(
                reason, lastMaxStepNorm, lastMaxAmountStep, lastLinearResidual,
                lastLambda, lastAccepted, allValid, iterationsRun,
                List.copyOf(phaseSetEvents), finalMb);
        return new Result(entries, mu, report);
    }

    // =========================================================================
    // Phase-addition candidate scan (STEP 12)
    // =========================================================================

    /**
     * A relaxed inactive-candidate state produced by {@link #relaxCandidate} --
     * the composition and internal CVCF vector at which the candidate phase's
     * tangent-plane driving force was best, plus that driving force. Immutable;
     * the candidate {@link Phase} is never touched to build one.
     *
     * @param x            candidate composition, length {@code K}
     * @param uFull        joint {@code [u ; x]} at that composition with the
     *                     internal block minimised by {@link CvmNewtonSolver}
     * @param gibbsEnergy  absolute {@code G = G0m + Gm} there
     * @param drivingForce {@code sum_i mu_i x_i - gibbsEnergy} (absolute tangent
     *                     driving force); {@code > 0} favours addition
     * @param valid        whether at least one grid point produced a converged
     *                     relaxed state -- {@code false} means the candidate
     *                     could not be evaluated and must not be added
     */
    private record CandidatePhaseState(
            double[] x, double[] uFull, double gibbsEnergy,
            double drivingForce, boolean valid) {

        static CandidatePhaseState invalid() {
            return new CandidatePhaseState(null, null, Double.NaN, Double.NaN, false);
        }
    }

    /**
     * Scores one inactive candidate phase by its absolute tangent-plane driving
     * force at a <em>relaxed</em> state, WITHOUT mutating the candidate
     * {@link Phase} or the caller's arrays (STEP 12 PART 3).
     *
     * <p>Search (STEP 12 PART 2 -- deliberately the smallest thing that works
     * for the supported CVM scope, not a global optimiser):</p>
     * <ul>
     *   <li><b>K = 2</b> (binary): a coarse 1-D composition grid,
     *       {@code x0 in {1/(G+1), 2/(G+1), ..., G/(G+1)}}.</li>
     *   <li><b>K = 3</b> (ternary): a small barycentric grid, all
     *       {@code (i,j,k)/G} with {@code i+j+k = G}, interior points only.</li>
     *   <li><b>K >= 4</b>: not searched -- returns
     *       {@link CandidatePhaseState#invalid()} (the caller logs it). Phase
     *       addition for quaternary+ is a documented v1 limitation.</li>
     * </ul>
     *
     * <p>At each grid composition the internal CVCF variables are minimised by
     * {@link CvmNewtonSolver#solve} (the project's hardened fixed-composition
     * minimiser; it returns an immutable result and mutates nothing external).
     * Non-converged grid points are skipped. Among the converged ones the one
     * with the largest {@code sum_i mu_i x_i - G} is returned. If no grid point
     * converges, the result is {@link CandidatePhaseState#invalid()}.</p>
     *
     * @param candidate the inactive phase to score -- read only
     * @param mu        the active set's absolute chemical potentials
     * @param temperature K
     */
    private static CandidatePhaseState relaxCandidate(
            Phase candidate, double[] mu, double temperature) {
        int k = candidate.numComponents;
        int ncf = candidate.ncf;
        CVMGibbsModel model = candidate.model;

        List<double[]> grid = candidateCompositionGrid(k);
        if (grid.isEmpty()) {
            return CandidatePhaseState.invalid();   // K >= 4: unsupported search
        }

        boolean any = false;
        double bestDrivingForce = Double.NEGATIVE_INFINITY;
        double[] bestX = null;
        double[] bestUFull = null;
        double bestG = Double.NaN;

        for (double[] x : grid) {
            // CvmNewtonSolver.solve handles its OWN numerical failures (a
            // singular inner Hessian, a degenerate cluster-variable point)
            // internally -- it returns a non-converged Result, it does not
            // throw for them. The one exception it does throw is
            // CancellationException, on Thread interruption; that is
            // cooperative cancellation and MUST propagate out of the whole
            // solve, not be turned into "candidate could not be evaluated".
            // Anything else escaping here is an unexpected programming error
            // (NPE, AIOOBE, ...) and must NOT be silently swallowed -- let it
            // propagate so it is seen, rather than masked as an unfavourable
            // candidate.
            CvmNewtonSolver.Result r =
                    new CvmNewtonSolver(model).solve(temperature, x, 1.0e-10, null, null);
            if (!r.converged() || !r.state().isValidIncludingPoints()) {
                continue;   // expected: this grid point did not yield a usable relaxed state
            }
            any = true;
            double g = r.state().g();                 // absolute G0m + Gm
            double dgf = 0.0;
            for (int i = 0; i < k; i++) {
                dgf += mu[i] * x[i];
            }
            dgf -= g;
            if (dgf > bestDrivingForce) {
                bestDrivingForce = dgf;
                bestX = x.clone();
                double[] uFull = new double[ncf + k];
                System.arraycopy(r.u(), 0, uFull, 0, ncf);
                System.arraycopy(x, 0, uFull, ncf, k);
                bestUFull = uFull;
                bestG = g;
            }
        }

        if (!any) {
            return CandidatePhaseState.invalid();
        }
        return new CandidatePhaseState(bestX, bestUFull, bestG, bestDrivingForce, true);
    }

    /**
     * Coarse composition grid for the candidate search -- interior points only
     * (a zero mole fraction is a lower-order subsystem, not this phase). Empty
     * for {@code K >= 4} (unsupported; see {@link #relaxCandidate}).
     */
    private static List<double[]> candidateCompositionGrid(int k) {
        List<double[]> grid = new ArrayList<>();
        if (k == 2) {
            int g = 12;
            for (int i = 1; i < g; i++) {
                double x0 = i / (double) g;
                grid.add(new double[] { x0, 1.0 - x0 });
            }
        } else if (k == 3) {
            int g = 8;
            for (int i = 1; i < g; i++) {
                for (int j = 1; j < g - i; j++) {
                    int m = g - i - j;
                    if (m < 1) {
                        continue;
                    }
                    grid.add(new double[] { i / (double) g, j / (double) g, m / (double) g });
                }
            }
        }
        return grid;   // k >= 4 -> empty
    }

    // =========================================================================
    // Mass-balance helpers (STEP 9)
    // =========================================================================

    /**
     * Independently reconstructs the represented component inventory
     * {@code sum_{p active} N_p x^p_i} from the phase states (NOT from any
     * Newton right-hand side) and forms the residual against {@code target}.
     *
     * @param target      the prescribed {@code overallAmounts}, or {@code null}
     *                    / wrong length -- then the report carries {@code NaN}
     *                    residuals and an empty target so a rejection still has
     *                    a populated {@link MassBalanceReport}.
     * @param targetSum   {@code sum_i target_i} (0 if target is unusable).
     * @param massScale   {@code max(targetSum, 1)} -- the relative-residual denominator.
     */
    private static MassBalanceReport massBalance(
            List<Phase> phases, double[] target, double targetSum, double massScale) {
        int k = phases.get(0).numComponents;

        double[] calc = new double[k];
        double calcSum = 0.0;
        for (Phase p : phases) {
            if (!p.active) {
                continue;
            }
            double[] x = p.composition();
            for (int i = 0; i < k; i++) {
                calc[i] += p.amount * x[i];
            }
            calcSum += p.amount;
        }

        double[] targetOverall = new double[k];
        double[] calcOverall = new double[k];
        if (targetSum > 0.0) {
            for (int i = 0; i < k; i++) targetOverall[i] = target[i] / targetSum;
        }
        if (calcSum > 0.0) {
            for (int i = 0; i < k; i++) calcOverall[i] = calc[i] / calcSum;
        }

        double maxAbs = Double.NaN;
        double maxRel = Double.NaN;
        if (target != null && target.length == k && targetSum > 0.0) {
            maxAbs = 0.0;
            for (int i = 0; i < k; i++) {
                maxAbs = Math.max(maxAbs, Math.abs(target[i] - calc[i]));
            }
            maxRel = maxAbs / massScale;
        }
        // residualBeforeLastStep is filled in by solve() at report time (it is
        // per-run state, not derivable here); NaN as the neutral default.
        return new MassBalanceReport(targetOverall, calcOverall, maxAbs, maxRel, Double.NaN);
    }

    /**
     * Builds a {@link Result} for a problem rejected before iterating -- state
     * echoed back exactly as passed (no silent change), reason and mass-balance
     * report supplied by the caller.
     */
    private static Result rejectedResult(List<Phase> phases, double temperature,
            double[] mu, ConvergenceReason reason, MassBalanceReport mb) {
        List<PhaseResult> entries = new ArrayList<>();
        for (Phase phase : phases) {
            CVMGibbsModel.State state = phase.model.atFull(temperature, phase.uFull);
            entries.add(new PhaseResult(
                    phase.label, phase.amount, phase.composition(), state.g(), state, false));
        }
        ConvergenceReport report = new ConvergenceReport(
                reason, Double.POSITIVE_INFINITY, Double.NaN, Double.NaN, Double.NaN,
                false, false, 0, List.of(), mb);
        return new Result(entries, mu, report);
    }

    /**
     * {@code G = G0m + Gm} at this phase's current joint state -- the absolute,
     * pure-element-anchored energy the outer equilibrium assembly needs, since
     * chemical potentials are equalised across phases and every phase's G must
     * therefore share one zero.
     *
     * <p>Evaluated through {@link org.ce.model.cvm.CVMGibbsModel}, which composes
     * the reference and mixing terms itself. This previously called
     * {@code model.evaluate(...).G} for Gm and added
     * {@code LatticeStability.g0m} separately -- a second place where
     * {@code G = G0m + Gm} was spelled out, and one that mutated the phase's
     * model as a side effect of a read.</p>
     */
    private static double currentG(Phase phase, double temperature) {
        return phase.model.atFull(temperature, phase.uFull).g();
    }

    private static boolean isFinite(double[] v) {
        for (double x : v) {
            if (!Double.isFinite(x)) {
                return false;
            }
        }
        return true;
    }

    private static double l2Norm(double[] v) {
        double sum = 0.0;
        for (double value : v) {
            sum += value * value;
        }
        return Math.sqrt(sum);
    }

    private static double[] scaled(double[] v, double factor) {
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            out[i] = v[i] * factor;
        }
        return out;
    }
}
