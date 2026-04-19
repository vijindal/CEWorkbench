package org.ce.calculation.workflow.thermo;

import org.ce.model.PhysicsConstants;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;
import org.ce.model.ProgressEvent;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CVMGibbsModel.ModelResult;
import org.ce.model.cvm.CVMGibbsModel.EquilibriumResult;
import org.ce.model.mcs.MCSRunner;
import org.ce.model.mcs.MetropolisMC.MCResult;
import org.ce.model.mcs.MetropolisMC.MCSUpdate;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;

import java.util.Arrays;
import java.util.logging.Logger;
import java.util.function.Consumer;

/**
 * Calculation-layer workflow for thermodynamic calculations.
 *
 * <p>Accepts a pre-built {@link ModelSession} (which already holds cluster data,
 * Hamiltonian, and resolved basis) and a {@link ThermodynamicRequest} (which
 * carries only calculation parameters: T, composition, MCS params, progress sinks).
 * No cluster identification or Hamiltonian loading occurs here.</p>
 *
 * <p>Dispatches directly to model-layer optimizers (CVMSolver for CVM, MCSRunner for MCS).</p>
 */
public class ThermodynamicWorkflow {

    private static final Logger LOG = Logger.getLogger(ThermodynamicWorkflow.class.getName());

    private CVMGibbsModel cachedGibbsModel;
    private ModelSession cachedSession;

    public ThermodynamicWorkflow() {
        // No injected engines; dispatch directly to model layer
    }

    /**
     * Runs a thermodynamic calculation using pre-built session state.
     *
     * <p>No cluster identification or Hamiltonian loading occurs here — both are
     * already available on {@code session}.</p>
     *
     * @param session pre-built model session holding cluster data, ECI, and basis
     * @param request calculation parameters (T, composition, MCS params, sinks)
     */
    public ThermodynamicResult runCalculation(
            ModelSession session,
            ThermodynamicRequest request) throws Exception {

        LOG.info("ThermodynamicWorkflow.runCalculation — ENTER");
        LOG.info("  session: " + session.label());
        LOG.info("  temperature: " + request.temperature + " K");
        LOG.info("  composition: " + Arrays.toString(request.composition));
        LOG.info("  engineType: " + session.engineConfig);

        emit(request.progressSink, "STAGE 3: Using pre-loaded Hamiltonian '"
                + session.resolvedHamiltonianId + "'");

        ThermodynamicResult result = switch (session.engineConfig) {
            case CVM -> runCvm(session, request);
            case MCS -> runMcs(session, request);
        };

        if (request.progressSink != null) {
            emit(request.progressSink, "");
            emit(request.progressSink, "  RESULTS AT " + request.temperature
                    + " K (composition: " + Arrays.toString(request.composition) + ")");
            emit(request.progressSink, "  " + "-".repeat(60));
            emit(request.progressSink, String.format("  Gibbs Energy (G): %15.6f J/mol", result.gibbsEnergy));
            emit(request.progressSink, String.format("  Enthalpy (H):     %15.6f J/mol", result.enthalpy));
            if (!Double.isNaN(result.entropy)) {
                emit(request.progressSink, String.format("  Entropy (S):      %15.6f J/(mol·K)", result.entropy));
            }
            if (result.optimizedCFs != null) {
                emit(request.progressSink, "  Equilibrium CFs (non-point):");
                for (int i = 0; i < result.optimizedCFs.length; i++) {
                    emit(request.progressSink, String.format("    CF[%d] = %15.10f", i, result.optimizedCFs[i]));
                }
            }
            emit(request.progressSink, "  " + "-".repeat(60));
            emit(request.progressSink, "  ✓ Calculation successful");
            emit(request.progressSink, "");
        }

        LOG.info("ThermodynamicWorkflow.runCalculation — EXIT: G="
                + String.format("%.4e", result.gibbsEnergy) + " J/mol");
        return result;
    }

    private ThermodynamicResult runCvm(ModelSession session, ThermodynamicRequest request)
            throws Exception {
        Consumer<String> progressSink = request.progressSink;
        printCvmHeader(progressSink);

        double temperature = request.temperature;
        double[] composition = request.composition;

        printInputParameters(progressSink, session.resolvedHamiltonianId, session.label(),
                            session.cecEntry.structurePhase, temperature, composition);

        validateInputs(temperature, composition);

        if (cachedGibbsModel == null || cachedSession != session) {
            cachedGibbsModel = new CVMGibbsModel();
            cachedGibbsModel.initialize(session.systemId.elements, session.systemId.structure, session.systemId.model, session.cecEntry, progressSink);
            cachedSession = session;
        }

        // 1. Resolve Equilibrium State (Minimization)
        EquilibriumResult solverResult = cachedGibbsModel.getEquilibriumState(
                temperature, composition, 1.0e-5,
                progressSink, request.eventSink);

        validateConvergence(solverResult, progressSink);

        // 2. Extract requested property via ThermodynamicMethods (Procedural Bridge)
        ThermodynamicMethods methods = new ThermodynamicMethods(cachedGibbsModel);
        double G = Double.NaN, H = Double.NaN, S = Double.NaN;
        
        // Prepare modDataIn array [output, T, x...]
        double[] modData = new double[2 + composition.length];
        modData[1] = temperature;
        System.arraycopy(composition, 0, modData, 2, composition.length);

        switch (request.property) {
            case ENTHALPY -> {
                methods.calHm(modData);
                H = modData[0];
            }
            case ENTROPY -> {
                methods.calSm(modData);
                S = modData[0];
            }
            case GIBBS_ENERGY -> {
                methods.calGm(modData);
                G = modData[0];
                H = cachedGibbsModel.calH();
                S = cachedGibbsModel.calS();
            }
            default -> {
                // For other properties, fall back to ModelResult if necessary
                ModelResult mr = solverResult.modelResult;
                G = mr.G; H = mr.H; S = mr.S;
            }
        }

        ThermodynamicResult result = new ThermodynamicResult(
                temperature, composition.clone(),
                G, H, S,
                Double.NaN, Double.NaN,
                solverResult.u,
                null, null,
                request.property
        );

        emit(progressSink, "================================================================================");
        return result;
    }

    private ThermodynamicResult runMcs(ModelSession session, ThermodynamicRequest request)
            throws Exception {

        Consumer<String> progressSink = request.progressSink;
        ClusCoordListData clusterData =
                session.clusterData.getDisorderedClusterResult().getDisClusterData();

        // Validate C-matrix dimensions match basis
        org.ce.model.cluster.CMatrixPipeline.CMatrixData matrixData = session.clusterData.getMatrixData();
        int cmatCols = (matrixData.getCmat().isEmpty() || matrixData.getCmat().get(0).isEmpty() 
                        || matrixData.getCmat().get(0).get(0).length == 0) 
                        ? 0 : matrixData.getCmat().get(0).get(0)[0].length;
        
        if (cmatCols != session.cvcfBasis.totalCfs()) {
            throw new IllegalStateException("C-matrix dimension mismatch (cmatCols=" + cmatCols 
                + ", basis.totalCfs=" + session.cvcfBasis.totalCfs() + ")");
        }
        LOG.fine("✓ C-matrix dimensions valid: " + cmatCols + " columns");

        double[] eci = CECEvaluator.evaluate(session.cecEntry, request.temperature, session.cvcfBasis, "MCS");

        int L           = request.mcsL;
        int nEquil      = request.mcsNEquil;
        int nAvg        = request.mcsNAvg;
        int totalSweeps = nEquil + nAvg;
        int N           = 2 * L * L * L;  // BCC sites

        MCSRunner.Builder builder = MCSRunner.builder()
                .clusterData(clusterData)
                .eci(eci)
                .numComp(request.composition.length)
                .T(request.temperature)
                .composition(request.composition)
                .nEquil(nEquil)
                .nAvg(nAvg)
                .L(L)
                .seed(System.currentTimeMillis())
                .R(PhysicsConstants.R_GAS)
                .basis(session.cvcfBasis)
                .matrixData(matrixData)
                .lcf(session.clusterData.getLcf())
                .cancellationCheck(Thread.currentThread()::isInterrupted);

        Consumer<String> strSink = progressSink;
        Consumer<ProgressEvent> evtSink = request.eventSink;

        if (strSink != null || evtSink != null) {
            if (strSink != null) {
                strSink.accept(String.format("  MCS started: L=%d (%d sites), %d equil + %d avg sweeps",
                        L, N, nEquil, nAvg));
            }
            // EngineStart emitted here — before build().run() — so the chart clears before sweeps arrive
            if (evtSink != null) {
                evtSink.accept(new ProgressEvent.EngineStart("MCS", totalSweeps));
            }
            final int sitesCount = N;
            final int sweepCount = totalSweeps;
            builder.updateListener(mcUpdate -> {
                if (strSink != null && (mcUpdate.getStep() % 100 == 0 || mcUpdate.getStep() == sweepCount)) {
                    strSink.accept(String.format("  [%-5s] sweep %4d/%-4d  <E>/site=%9.5f  accept=%5.1f%%",
                            mcUpdate.getPhase(),
                            mcUpdate.getStep(), sweepCount,
                            mcUpdate.getE_total() / sitesCount,
                            mcUpdate.getAcceptanceRate() * 100));
                }
                if (evtSink != null) {
                    evtSink.accept(new ProgressEvent.McSweep(
                            mcUpdate.getStep(), sweepCount,
                            mcUpdate.getE_total() / sitesCount,
                            mcUpdate.getAcceptanceRate(),
                            mcUpdate.getPhase() == MCSUpdate.Phase.EQUILIBRATION,
                            null));  // CFs no longer charted
                }
            });
        }

        MCSRunner.MCSRunResult runResult = builder.build().run();
        MCResult mcResult = runResult.result;

        // [DEBUG] Raw simulation output
        MCSRunner.Debug.printMcsSummary(mcResult);

        // Compute statistics from raw time series in the calculation layer
        MCSStatisticsProcessor statsProcessor = new MCSStatisticsProcessor(
                mcResult.getNSites(),
                PhysicsConstants.R_GAS,
                request.temperature,
                mcResult.getSeriesHmix(),
                mcResult.getSeriesE(),
                mcResult.getSeriesCF());
        statsProcessor.computeStatistics();

        // Post-run statistics summary with computed error bars
        if (strSink != null) {
            strSink.accept(
                "  ── MCS Statistics (Post-Processed) ────────────────────────");
            strSink.accept(String.format(
                "  τ_int = %.1f sweeps  s = %.3f  n_eff = %d  (block size = %d, %d blocks)",
                statsProcessor.getTauInt(), statsProcessor.getStatInefficiency(),
                statsProcessor.getNEff(), statsProcessor.getBlockSizeUsed(), statsProcessor.getNBlocks()));
            strSink.accept(String.format(
                "  ⟨H⟩/site = %.6f ± %.6f  J/mol",
                mcResult.getHmixPerSite(), nanZero(statsProcessor.getStdHmixPerSite())));
            strSink.accept(String.format(
                "  ⟨E⟩/site = %.6f ± %.6f  J/mol",
                mcResult.getEnergyPerSite(), nanZero(statsProcessor.getStdEnergyPerSite())));
            strSink.accept(String.format(
                "  Cv       = %.6f ± %.6f  J/(mol·K)  [jackknife]",
                statsProcessor.getCvJackknife(), nanZero(statsProcessor.getCvStdErr())));

            // Correlation functions with error bars
            double[] cfs    = mcResult.getAvgCFs();
            double[] stdCFs = statsProcessor.getStdCFs();
            if (cfs != null && cfs.length > 0) {
                strSink.accept(
                    "  ── Correlation Functions ────────────────────────────────");
                for (int i = 0; i < cfs.length; i++) {
                    double sigma = (stdCFs != null && i < stdCFs.length) ? stdCFs[i] : Double.NaN;
                    strSink.accept(String.format(
                        "  CF[%2d] = %+.6f ± %.6f",
                        i, cfs[i], nanZero(sigma)));
                }
            }
        }

        return new ThermodynamicResult(
                mcResult.getTemperature(),
                mcResult.getComposition(),
                Double.NaN,                           // gibbsEnergy — not available from canonical MCS
                mcResult.getHmixPerSite(),            // enthalpy
                Double.NaN,                           // entropy — not available from canonical MCS
                statsProcessor.getStdHmixPerSite(),   // stdEnthalpy (from post-processing)
                statsProcessor.getCvJackknife(),      // heat capacity from jackknife
                null,                                 // optimizedCFs (CVM only)
                mcResult.getAvgCFs(),                 // mean correlation functions
                statsProcessor.getStdCFs()            // CF standard errors from post-processing
        );
    }

    // ---- CVM helpers ----

    private void validateConvergence(EquilibriumResult result, Consumer<String> sink) {
        if (!result.converged) {
            emit(sink, "  [!] CVM minimization FAILED to converge!");
            throw new RuntimeException(
                    "CVM minimization failed to converge within " + result.iterations + " iterations.");
        }
        emit(sink, String.format("  - CONVERGED in %d iterations (||Gu||=%.4e)",
                result.iterations, result.finalGradientNorm));
    }

    private void printCvmHeader(Consumer<String> sink) {
        String border = "================================================================================";
        emit(sink, border);
        emit(sink, "                       CVM THERMODYNAMIC CALCULATION");
        emit(sink, border);
    }

    private void printInputParameters(Consumer<String> sink, String systemId, String systemName,
                                     String structurePhase, double temperature, double[] composition) {
        emit(sink, "\nINPUT PARAMETERS");
        emit(sink, "-----------------");
        emit(sink, "  - System ID:      " + systemId);
        emit(sink, "  - System Name:    " + systemName);
        emit(sink, "  - Structure:      " + structurePhase);
        emit(sink, "  - Temperature:    " + temperature + " K");
        emit(sink, "  - Composition:    [" + formatArray(composition) + "]");
    }

    private void validateInputs(double temperature, double[] composition) {
        if (temperature <= 0) {
            throw new IllegalArgumentException("Temperature must be positive: " + temperature);
        }
        if (composition == null || composition.length < 2) {
            throw new IllegalArgumentException("Invalid composition array");
        }
        double sum = 0.0;
        for (double x : composition) {
            if (x < 0 || x > 1) {
                throw new IllegalArgumentException("Composition values must be in [0,1]");
            }
            sum += x;
        }
        if (Math.abs(sum - 1.0) > 1e-9) {
            throw new IllegalArgumentException("Composition must sum to 1.0, got: " + sum);
        }
    }

    private static String formatArray(double[] arr) {
        if (arr == null || arr.length == 0)
            return "";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            if (i > 0)
                sb.append(", ");
            sb.append(String.format("%.6f", arr[i]));
        }
        return sb.toString();
    }

    private static double nanZero(double value) {
        return Double.isNaN(value) ? 0.0 : value;
    }

    private static void emit(Consumer<String> sink, String line) {
        if (sink != null) sink.accept(line);
    }

}
