package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.CMatrix;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.engine.cvm.CVMGibbsModel.CVMInput;
import org.ce.domain.result.EquilibriumState;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECEvaluator;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.workflow.ClusterIdentificationWorkflow;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;

/**
 * CVM thermodynamic engine.
 *
 * Implements the ThermodynamicEngine interface using the
 * Cluster Variation Method (CVM) for thermodynamic calculations.
 *
 * <p>
 * Responsibilities:
 * </p>
 * <ul>
 * <li>Create CVMInput (topology) from AllClusterData</li>
 * <li>Evaluate ECI at temperature</li>
 * <li>Validate input consistency</li>
 * <li>Run CVM minimization via CVMPhaseModel</li>
 * <li>Return EquilibriumState</li>
 * </ul>
 */
public class CVMEngine implements ThermodynamicEngine {

    private static final Logger LOG = Logger.getLogger(CVMEngine.class.getName());

    @Override
    public EquilibriumState compute(ThermodynamicInput input) throws Exception {

        // Resolve cluster data (on-the-fly or input-provided)
        AllClusterData clusterData = resolveClusterData(input);

        emit(input.progressSink, "================================================================================");
        emit(input.progressSink, "                       CVM THERMODYNAMIC CALCULATION");
        emit(input.progressSink, "================================================================================");

        CECEntry cec = input.cec;
        double temperature = input.temperature;
        double[] composition = input.composition;
        String systemId = input.systemId;
        String systemName = input.systemName;
        String structurePhase = cec.structurePhase;

        emit(input.progressSink, "\nINPUT PARAMETERS");
        emit(input.progressSink, "-----------------");
        emit(input.progressSink, "  - System ID:      " + systemId);
        emit(input.progressSink, "  - System Name:    " + systemName);
        emit(input.progressSink, "  - Structure:      " + structurePhase);
        emit(input.progressSink, "  - Temperature:    " + temperature + " K");
        emit(input.progressSink, "  - Composition:    [" + formatArray(composition) + "]");

        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get(structurePhase, cec.model, composition.length);

        // 1. Validate structural consistency
        validateCmatEciConsistency(clusterData.getCMatrixResult(), basis);

        // 2. Validate inputs
        validateInputs(temperature, composition);

        /*
         * 3. Create CVMInput (the static model topology class)
         */
        CVMInput cvmInput = new CVMInput(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                clusterData.getCMatrixResult(),
                systemId,
                systemName,
                composition.length,
                basis
        );

        /*
         * 4. Evaluate ECI at temperature
         */
        emit(input.progressSink, "\n  - Evaluating Hamiltonian at T=" + temperature + " K...");
        double[] eci = CECEvaluator.evaluate(cec, temperature, basis, "CVM");

        /*
         * 5. Run CVM minimization
         */
        emit(input.progressSink, "\n  - Running CVM N-R minimization...");
        CVMGibbsModel gibbsModel = new CVMGibbsModel(cvmInput);
        CVMSolver solver = new CVMSolver();

        CVMSolver.EquilibriumResult solverResult = solver.minimize(
                gibbsModel,
                composition,
                temperature,
                eci,
                1.0e-5, // default tolerance
                input.progressSink,
                input.eventSink
        );

        if (!solverResult.converged) {
            emit(input.progressSink, "  [!] CVM minimization FAILED to converge!");
            throw new RuntimeException(
                    "CVM minimization failed to converge within " + solverResult.iterations + " iterations.");
        }

        emit(input.progressSink, String.format("  - CONVERGED in %d iterations (||Gu||=%.4e)",
                solverResult.iterations, solverResult.finalGradientNorm));

        emit(input.progressSink, "\nEQUILIBRIUM RESULTS");
        emit(input.progressSink, "-------------------");
        emit(input.progressSink, String.format("  - G(eq): %.6e J/mol", solverResult.modelValues.G));
        emit(input.progressSink, String.format("  - H(eq): %.6e J/mol", solverResult.modelValues.H));
        emit(input.progressSink, String.format("  - S(eq): %.6e J/(mol·K)", solverResult.modelValues.S));

        /*
         * 6. Iteration trace (already printed live by solver)
         */

        EquilibriumState result = new EquilibriumState(
                temperature,
                composition.clone(),
                solverResult.modelValues.H,
                solverResult.modelValues.G,
                solverResult.modelValues.S,
                Double.NaN, // enthalpyStdErr (MCS only)
                Double.NaN, // heatCapacity (MCS only)
                solverResult.u,
                null, // avgCFs (MCS only)
                null // stdCFs (MCS only)
        );
        emit(input.progressSink, "================================================================================");
        return result;
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

    private void emit(java.util.function.Consumer<String> sink, String line) {
        if (sink != null)
            sink.accept(line);
        LOG.fine(line);
    }

    private static void validateCmatEciConsistency(CMatrix.Result cmat, CvCfBasis basis) {
        List<String> cmatNames = cmat.getCmatCfNames();
        if (cmatNames == null)
            return; // orthogonal path, skip

        int ncf = basis.numNonPointCfs;
        for (int i = 0; i < ncf; i++) {
            String cmatLabel = cmatNames.get(i);
            String basisLabel = basis.cfNames.get(i);
            if (!cmatLabel.equals(basisLabel)) {
                throw new IllegalStateException(
                        "C-matrix column " + i + " is labelled '" + cmatLabel
                                + "' but basis expects '" + basisLabel
                                + "'. The cluster_data.json and CVCF basis are out of sync. "
                                + "Re-run type1a to regenerate cluster_data.json.");
            }
        }
    }

    /**
     * Resolves the required structural cluster data.
     *
     * <p>Uses the provided input.clusterData if available, otherwise triggers
     * on-the-fly identification using standard naming conventions.</p>
     */
    private AllClusterData resolveClusterData(ThermodynamicInput input) {
        // Strict Always-Fresh Identification Workflow
        emit(input.progressSink, "  [NOTE] Starting Always Fresh Structural Identification...");
 
        ClusterIdentificationRequest request = new ClusterIdentificationRequest(input);
 
        // Run fresh identification workflow (No caching)
        return ClusterIdentificationWorkflow.identify(request, input.progressSink);
    }

    /**
     * Formats a double array for logging.
     */
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
}
