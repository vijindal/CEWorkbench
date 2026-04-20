package org.ce.calculation.workflow.thermo;

import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.model.ModelSession;
import org.ce.model.PhysicsConstants;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.mcs.MCSRunner;
import org.ce.model.mcs.MCSGeometry;
import org.ce.model.mcs.MetropolisMC.MCResult;

import java.util.Arrays;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Unified engine for thermodynamic calculations (CVM and MCS).
 * Handles caching of expensive model state (minimization loops, geometry).
 */
public class ThermodynamicWorkflow {

    private static final Logger LOG = Logger.getLogger(ThermodynamicWorkflow.class.getName());

    // Caching records for value-based identity
    private record CvmCache(ModelSession session, CVMGibbsModel model) {
        boolean validFor(ModelSession s) {
            return session == s;
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
        validateInputs(request.temperature, request.composition);

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

    private void validateInputs(double T, double[] x) {
        if (T < 0)
            throw new IllegalArgumentException("Temperature cannot be negative: " + T);
        if (x == null || x.length == 0)
            throw new IllegalArgumentException("Composition array missing");
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
            CVMGibbsModel model = new CVMGibbsModel();
            model.initialize(session.systemId.elements(), session.systemId.structure(), session.systemId.model(),
                    session.cecEntry, request.progressSink);
            cvmCache = new CvmCache(session, model);
        }

        CVMGibbsModel model = cvmCache.model();
        double T = request.temperature;
        double[] x = request.composition;

        // Extract result via procedural bridge logic (formerly ThermodynamicMethods)
        model.getEquilibriumState(T, x, 1e-5, request.progressSink(), request.eventSink(), request.property());

        double g = Double.NaN, h = Double.NaN, s = Double.NaN;
        switch (request.property) {
            case GIBBS_ENERGY -> g = model.calG();
            case ENTHALPY -> h = model.calH();
            case ENTROPY -> s = model.calS();
        }

        return new ThermodynamicResult(
                T, x, g, h, s,
                Double.NaN, // stdEnthalpy
                Double.NaN, // heatCapacity
                model.calCfs(),
                null, // avgCFs
                null, // stdCFs
                request.property());
    }

    // ── MCS Engine ────────────────────────────────────────────────────────────

    private ThermodynamicResult runMcs(ModelSession session, Request request) throws Exception {
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

        return new ThermodynamicResult(
                request.temperature(),
                request.composition(),
                Double.NaN, // G not directly available in single-point MCS
                r.getEnergyPerSite(),
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
