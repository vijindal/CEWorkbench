package org.ce.calculation.engine.cvm;

import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.cvcf.CvCfBasis;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.calculation.engine.ThermodynamicEngine;
import org.ce.model.ThermodynamicResult;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;

import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * CVM thermodynamic engine.
 *
 * Implements the ThermodynamicEngine interface using the
 * Cluster Variation Method (CVM) for thermodynamic calculations.
 *
 * <p>
 * Orchestrates the CVM calculation pipeline:
 * </p>
 * <ul>
 * <li>Identify cluster structures and basis functions</li>
 * <li>Evaluate ECI (effective cluster interactions) at temperature</li>
 * <li>Initialize CVMGibbsModel</li>
 * <li>Run CVM minimization via CVMSolver</li>
 * <li>Return EquilibriumState with thermodynamic properties</li>
 * </ul>
 */
public class CVMEngine implements ThermodynamicEngine {

    private static final Logger LOG = Logger.getLogger(CVMEngine.class.getName());

    @Override
    public ThermodynamicResult compute(ThermodynamicEngine.Input input) throws Exception {
        printHeader(input.progressSink);

        double temperature = input.temperature;
        double[] composition = input.composition;
        CECEntry cec = input.cec;

        printInputParameters(input.progressSink, input.systemId, input.systemName,
                            cec.structurePhase, temperature, composition);

        validateInputs(temperature, composition);
        CvCfBasis basis = getBasis(cec, composition.length);
        AllClusterData clusterData = input.clusterData;
        CVMGibbsModel gibbsModel = new CVMGibbsModel(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                clusterData.getCMatrixResult(),
                basis);

        double[] eci = evaluateEci(cec, temperature, basis, input.progressSink);
        CVMSolver.EquilibriumResult solverResult = runSolver(gibbsModel, composition, temperature, eci, input);

        validateConvergence(solverResult, input.progressSink);
        printEquilibriumResults(input.progressSink, solverResult);

        ThermodynamicResult result = buildEquilibriumState(temperature, composition, solverResult);
        emit(input.progressSink, "================================================================================");
        return result;
    }

    private CvCfBasis getBasis(CECEntry cec, int nComponents) {
        return CvCfBasis.Registry.INSTANCE.get(cec.structurePhase, cec.model, nComponents);
    }

    private double[] evaluateEci(CECEntry cec, double temperature, CvCfBasis basis, Consumer<String> sink) {
        emit(sink, "\n  - Evaluating Hamiltonian at T=" + temperature + " K...");
        return CECEvaluator.evaluate(cec, temperature, basis, "CVM");
    }

    private CVMSolver.EquilibriumResult runSolver(CVMGibbsModel gibbsModel, double[] composition,
                                                   double temperature, double[] eci, ThermodynamicEngine.Input input) {
        emit(input.progressSink, "\n  - Running CVM N-R minimization...");
        return new CVMSolver().minimize(
                gibbsModel,
                composition,
                temperature,
                eci,
                1.0e-5,
                input.progressSink,
                input.eventSink
        );
    }

    private ThermodynamicResult buildEquilibriumState(double temperature, double[] composition,
                                                    CVMSolver.EquilibriumResult result) {
        return new ThermodynamicResult(
                temperature,
                composition.clone(),
                result.modelValues.G,  // gibbsEnergy
                result.modelValues.H,  // enthalpy
                result.modelValues.S,  // entropy
                Double.NaN,            // stdEnthalpy
                Double.NaN,            // heatCapacity
                result.u,
                null,
                null
        );
    }

    private void printHeader(Consumer<String> sink) {
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


    private void validateConvergence(CVMSolver.EquilibriumResult result, Consumer<String> sink) {
        if (!result.converged) {
            emit(sink, "  [!] CVM minimization FAILED to converge!");
            throw new RuntimeException(
                    "CVM minimization failed to converge within " + result.iterations + " iterations.");
        }
        emit(sink, String.format("  - CONVERGED in %d iterations (||Gu||=%.4e)",
                result.iterations, result.finalGradientNorm));
    }

    private void printEquilibriumResults(Consumer<String> sink, CVMSolver.EquilibriumResult result) {
        emit(sink, "\nEQUILIBRIUM RESULTS");
        emit(sink, "-------------------");
        emit(sink, String.format("  - G(eq): %.6e J/mol", result.modelValues.G));
        emit(sink, String.format("  - H(eq): %.6e J/mol", result.modelValues.H));
        emit(sink, String.format("  - S(eq): %.6e J/(mol·K)", result.modelValues.S));
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

    private void emit(Consumer<String> sink, String line) {
        if (sink != null)
            sink.accept(line);
        LOG.fine(line);
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
