package org.ce.calculation.workflow.thermo;

import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.debug.MCSDebug;
import org.ce.model.ModelSession;
import org.ce.model.PhysicsConstants;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.equilibrium.CvmNewtonSolver;
import org.ce.model.mcs.MCSRunner;
import org.ce.model.mcs.MCSGeometry;
import org.ce.model.mcs.MetropolisMC.MCResult;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Unified engine for thermodynamic calculations (CVM and MCS).
 * Handles caching of expensive model state (minimization loops, geometry).
 */
public class ThermodynamicWorkflow {

    private static final Logger LOG = Logger.getLogger(ThermodynamicWorkflow.class.getName());

    // Caching records for value-based identity
    /**
     * The model for one session, together with the last minimisation result.
     *
     * <p>The (T, x) memoisation lives here rather than on the model: it is a
     * workflow concern -- avoiding a repeated solve for the same requested
     * point -- not physics. The model itself is a pure evaluator and holds no
     * current point.</p>
     */
    private static final class CvmCache {
        final ModelSession session;
        final CVMGibbsModel model;
        double lastT = Double.NaN;
        double[] lastX;
        CvmNewtonSolver.Result lastResult;

        CvmCache(ModelSession session, CVMGibbsModel model) {
            this.session = session;
            this.model = model;
        }

        boolean validFor(ModelSession s) {
            return session == s;
        }

        CVMGibbsModel model() {
            return model;
        }

        boolean hasResultFor(double t, double[] x) {
            return lastResult != null
                    && Math.abs(lastT - t) <= 1.0e-5
                    && lastX != null
                    && Arrays.equals(lastX, x);
        }

        void store(double t, double[] x, CvmNewtonSolver.Result result) {
            this.lastT = t;
            this.lastX = x.clone();
            this.lastResult = result;
        }
    }

    private record McsCache(ModelSession session, int L, MCSGeometry geo,
            double T, MCSRunner runner) {
        boolean geoValidFor(ModelSession s, int l) {
            return session == s && L == l;
        }

        boolean runnerValidFor(ModelSession s, int l, double t) {
            return geoValidFor(s, l) && Double.compare(T, t) == 0;
        }
    }

    private CvmCache cvmCache;
    private McsCache mcsCache;

    /**
     * Request describing a single calculation point.
     */
    public static record Request(
            double temperature,
            double[] composition,
            Property property,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink,
            int mcsL,
            int mcsNEquil,
            int mcsNAvg,
            double[] fixedCorrelations) {
        public Request(double T, double[] x) {
            this(T, x, Property.GIBBS_ENERGY, null, null, 4, 1000, 2000, null);
        }
    }

    public ThermodynamicResult runCalculation(ModelSession session, Request request) throws Exception {
        validateInputs(session, request.temperature, request.composition);

        emit(request.progressSink, String.format(
                "\nCALCULATION START [%s] — T=%.1f K, x=%s (%s)",
                session.engineConfig, request.temperature, Arrays.toString(request.composition), request.property));

        ThermodynamicResult result = switch (session.engineConfig) {
            case CVM -> runCvm(session, request);
            case MCS -> runMcs(session, request);
        };
        printResultSummary(request.progressSink, request.temperature, request.composition, result);
        return result;
    }

    private void validateInputs(ModelSession session, double T, double[] x) {
        if (T < 0)
            throw new IllegalArgumentException("Temperature cannot be negative: " + T);
        if (x == null || x.length == 0)
            throw new IllegalArgumentException("Composition array missing");
        if (x.length != session.numComponents())
            throw new IllegalArgumentException(
                    "Composition length " + x.length + " != numComponents " + session.numComponents());
        double sum = 0;
        for (double val : x)
            sum += val;
        if (Math.abs(sum - 1.0) > 1e-4) {
            throw new IllegalArgumentException("Composition does not sum to 1.0: " + Arrays.toString(x));
        }
    }

    // ── CVM Engine ────────────────────────────────────────────────────────────

    private ThermodynamicResult runCvm(ModelSession session, Request request) throws Exception {
        if (cvmCache == null || !cvmCache.validFor(session)) {
            CVMGibbsModel model = CVMGibbsModel.of(
                    session.systemId.elements(), session.systemId.structure(),
                    session.systemId.model(), session.cecEntry, request.progressSink);
            cvmCache = new CvmCache(session, model);
        }

        CVMGibbsModel model = cvmCache.model();
        double T = request.temperature;
        double[] x = request.composition;

        CvmNewtonSolver.Result eq;
        if (cvmCache.hasResultFor(T, x)) {
            emit(request.progressSink, "  [Model] Reusing cached equilibrium state for these parameters.");
            eq = cvmCache.lastResult;
        } else {
            emit(request.progressSink, String.format(
                    "\n  [Model] Parameters updated: T = %.1f K, x = %s", T, Arrays.toString(x)));
            emit(request.progressSink, "  [Model] Initiating internal minimization (Newton-Raphson loop)...");
            eq = new CvmNewtonSolver(model).solve(
                    T, x, 1e-5, request.progressSink(), request.eventSink());
            emit(request.progressSink, eq.converged()
                    ? "  [Model] \u2713 Minimization converged in " + eq.iterations() + " iterations."
                    : "  [Model] \u26a0 Minimization FAILED to converge.");
            cvmCache.store(T, x, eq);
        }

        CVMGibbsModel.State state = eq.state();

        if (!eq.converged()) {
            emit(request.progressSink, String.format(
                    "  [WARNING] CVM minimization did NOT converge at T=%.1f K, x=%s "
                    + "(%d iterations, final ||grad G||=%.3e). Results are unreliable.",
                    T, Arrays.toString(x), eq.iterations(), eq.finalGradientNorm()));
            LOG.warning(String.format("CVM non-convergence at T=%.1f K, x=%s: %d iters, gradNorm=%.3e",
                    T, Arrays.toString(x), eq.iterations(), eq.finalGradientNorm()));
        }

        // Mixing quantities (Gm/Hm/Sm), not the pure-element-anchored
        // absolutes (calG/calH/calS): these feed ThermodynamicResult, which
        // the CLI/GUI/JSON API report and which CLAUDE.md's documented
        // verification values (e.g. G = -3480.5209063901 for Nb-Ti) are
        // anchored to. Switching to absolute here would change every reported
        // energy without changing any physics.
        double g = Double.NaN, h = Double.NaN, s = Double.NaN;
        switch (request.property) {
            case GIBBS_ENERGY -> {
                g = state.gm();
                h = state.hm();
                s = state.sm();
            }
            case ENTHALPY -> h = state.hm();
            case ENTROPY -> s = state.sm();
        }

        ThermodynamicResult result = new ThermodynamicResult(
                T, x, g, h, s,
                Double.NaN, // stdEnthalpy
                Double.NaN, // heatCapacity
                state.cfs(),
                null, // avgCFs
                null, // stdCFs
                request.property());

        return result
                .withConvergence(eq.converged(), eq.iterations(), eq.finalGradientNorm())
                .withSro(computeSro(model, state, x, request.progressSink));
    }

    /**
     * Cowley-Warren SRO parameters at the converged point (Jindal &amp; Lele
     * 2025, Eq. 40), read off the state that already holds the cluster
     * probabilities. Returns null if the calculation fails outright — SRO is
     * supplementary, so an unsupported geometry must not fail the whole
     * calculation; individual unsupported shells are already skipped by
     * {@link CVMGibbsModel.State#pairSroByShell()}.
     */
    private Map<String, List<CVMGibbsModel.PairSro>> computeSro(
            CVMGibbsModel model, CVMGibbsModel.State state, double[] x, Consumer<String> sink) {
        try {
            return state.pairSroByShell();
        } catch (RuntimeException e) {
            emit(sink, "  [SRO] not computed: " + e.getMessage());
            return null;
        }
    }

    // ── MCS Engine ────────────────────────────────────────────────────────────

    private ThermodynamicResult runMcs(ModelSession session, Request request) throws Exception {
        // Route debug output to the GUI progress sink
        MCSDebug.setSink(request.progressSink());
        try {
            return runMcsInternal(session, request);
        } finally {
            MCSDebug.clearSink();
        }
    }

    private ThermodynamicResult runMcsInternal(ModelSession session, Request request) throws Exception {
        int L = request.mcsL;

        if (mcsCache == null || !mcsCache.geoValidFor(session, L)) {
            emit(request.progressSink, String.format("  [MCS Geometry] rebuilding for L=%d...", L));
            MCSGeometry geo = MCSGeometry.build(session, L, request.progressSink());
            mcsCache = new McsCache(session, L, geo, Double.NaN, null);
        }

        if (mcsCache.runner() == null || !mcsCache.runnerValidFor(session, L, request.temperature)) {
            emit(request.progressSink,
                    String.format("  [MCS Model] evaluating ECIs at T=%.1f K...", request.temperature));
            MCSRunner runner = MCSRunner.forTemperature(mcsCache.geo(), session, request.temperature, request.progressSink());
            mcsCache = new McsCache(session, L, mcsCache.geo(), request.temperature, runner);
        } else {
            emit(request.progressSink, String.format("  [MCS Model] reusing cached model (T=%.1f K, L=%d)",
                    request.temperature, L));
        }

        if (request.eventSink() != null) {
            request.eventSink()
                    .accept(new org.ce.model.ProgressEvent.EngineStart("MCS", request.mcsNEquil() + request.mcsNAvg()));
        }

        MCResult r = mcsCache.runner().run(
                request.composition(),
                request.mcsNEquil(),
                request.mcsNAvg(),
                System.currentTimeMillis(),
                request.progressSink(),
                mcUpdate -> {
                    if (request.eventSink() != null) {
                        request.eventSink().accept(new org.ce.model.ProgressEvent.McSweep(
                                mcUpdate.getStep(),
                                request.mcsNEquil() + request.mcsNAvg(),
                                mcUpdate.getE_total() / mcsCache.geo().nSites(),
                                mcUpdate.getAcceptanceRate(),
                                mcUpdate.getPhase() == org.ce.model.mcs.MetropolisMC.MCSUpdate.Phase.EQUILIBRATION,
                                mcUpdate.getCfs()));
                    }
                },
                () -> false // cancellationCheck
        ).result;

        // ── MCS-DBG: MCResult → ThermodynamicResult mapping ──
        if (MCSDebug.ENABLED) {
            MCSDebug.separator("MCResult → ThermodynamicResult MAPPING");
            MCSDebug.log("FLOW", "MCResult.hmixPerSite   = %.10f  → ThermodynamicResult.enthalpy", r.getHmixPerSite());
            MCSDebug.log("FLOW", "MCResult.energyPerSite = %.10f  (incremental running E/N, not used)", r.getEnergyPerSite());
            MCSDebug.log("FLOW", "MCResult.acceptRate    = %.4f", r.getAcceptRate());
            MCSDebug.vector("FLOW", "MCResult.avgCFs → ThermodynamicResult.avgCFs", r.getAvgCFs());
        }

        return new ThermodynamicResult(
                request.temperature(),
                request.composition(),
                Double.NaN, // G not directly available in single-point MCS
                r.getHmixPerSite(),  // ⟨Hmix⟩/site = Σ eci·⟨vCvcf⟩, matches CVM calH()
                Double.NaN, // S not directly available
                Double.NaN, // stdEnthalpy (σ) - not in MCResult
                Double.NaN, // heatCapacity - requires FSS
                null, // optimizedCFs
                r.getAvgCFs(),
                null, // stdCFs
                request.property());
    }

    private void printResultSummary(Consumer<String> sink, double T, double[] x, ThermodynamicResult r) {
        if (sink == null)
            return;
        emit(sink, "  RESULTS AT " + T + " K (" + Arrays.toString(x) + ")");
        emit(sink, "  " + "-".repeat(60));
        if (!Double.isNaN(r.gibbsEnergy))
            emit(sink, String.format("  Gibbs Energy (G): %15.6f J/mol", r.gibbsEnergy));
        if (!Double.isNaN(r.enthalpy))
            emit(sink, String.format("  Enthalpy (H):     %15.6f J/mol", r.enthalpy));
        if (!Double.isNaN(r.entropy))
            emit(sink, String.format("  Entropy (S):      %15.6f J/mol\u00B7K", r.entropy));
        emit(sink, "  " + "-".repeat(60));
    }

    private void emit(Consumer<String> sink, String msg) {
        if (sink != null)
            sink.accept(msg);
    }
}
