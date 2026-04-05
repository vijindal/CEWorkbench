package org.ce.domain.engine.cvm;

import static org.ce.domain.cluster.AllClusterData.ClusterData;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.CMatrix;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.engine.ProgressEvent;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.engine.cvm.CVMPhaseModel.CVMInput;
import org.ce.domain.engine.cvm.CVMPhaseModel.SolverResult;
import org.ce.domain.engine.cvm.CVMPhaseModel.SolverResult.IterationSnapshot;
import org.ce.domain.result.EquilibriumState;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECEvaluator;

import java.util.List;
import java.util.logging.Logger;

/**
 * CVM thermodynamic engine.
 *
 * Implements the ThermodynamicEngine interface using the
 * Cluster Variation Method (CVM) for thermodynamic calculations.
 *
 * <p>Responsibilities:</p>
 * <ul>
 *   <li>Create CVMInput (topology) from AllClusterData</li>
 *   <li>Evaluate ECI at temperature</li>
 *   <li>Validate input consistency</li>
 *   <li>Run CVM minimization via CVMPhaseModel</li>
 *   <li>Return EquilibriumState</li>
 * </ul>
 */
public class CVMEngine implements ThermodynamicEngine {

    private static final Logger LOG = Logger.getLogger(CVMEngine.class.getName());


    @Override
    public EquilibriumState compute(ThermodynamicInput input) throws Exception {

        LOG.info("=== CVMEngine.compute() START ===");
        AllClusterData clusterData = input.clusterData;
        CECEntry cec = input.cec;
        double temperature = input.temperature;
        double[] composition = input.composition;
        String systemId = input.systemId;
        String systemName = input.systemName;
        String structurePhase = cec.structurePhase;

        LOG.info("Input parameters:");
        LOG.info("  systemId: " + systemId);
        LOG.info("  systemName: " + systemName);
        LOG.info("  structurePhase: " + structurePhase);
        LOG.info("  temperature: " + temperature + " K");
        LOG.info("  composition: [" + formatArray(composition) + "]");
        LOG.info("  numComponents: " + composition.length);

        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get(structurePhase, composition.length);

        /*
         * 1. Validate CMatrix exists
         */
        LOG.fine("Validating C-Matrix...");
        if (clusterData.getCMatrixResult() == null) {
            LOG.severe("C-Matrix not found in AllClusterData!");
            throw new IllegalStateException(
                "CMatrix not found in AllClusterData. " +
                "Ensure ClusterIdentificationWorkflow generated Stage 3."
            );
        }
        LOG.fine("✓ C-Matrix found");

        // Validate C-matrix dimensions match basis
        CMatrix.Result cmatResult = clusterData.getCMatrixResult();
        cmatResult.validateCols(
                basis.totalCfs(),
                "C-matrix dimension mismatch (basis.numNonPointCfs=" + basis.numNonPointCfs
                + " + " + basis.numComponents + " point variables)"
        );
        LOG.fine("✓ C-matrix dimensions valid: " + basis.totalCfs() + " columns");
        
        // Fix 4: Structural consistency check (labels vs basis)
        validateCmatEciConsistency(cmatResult, basis);

        /*
         * 2. Validate composition
         */
        LOG.fine("Validating composition...");
        if (composition == null || composition.length < 2) {
            LOG.severe("Invalid composition array: " + (composition == null ? "null" : "length=" + composition.length));
            throw new IllegalArgumentException("Invalid composition array");
        }

        double sum = 0.0;
        for (double x : composition) {
            if (x < 0 || x > 1) {
                LOG.severe("Composition value out of range [0,1]: " + x);
                throw new IllegalArgumentException("Composition values must be in [0,1]");
            }
            sum += x;
        }

        if (Math.abs(sum - 1.0) > 1e-9) {
            LOG.severe("Composition sum != 1.0: " + sum);
            throw new IllegalArgumentException("Composition must sum to 1.0, got: " + sum);
        }
        LOG.fine("✓ Composition valid (sum=" + String.format("%.9f", sum) + ")");

        /*
         * 3. Validate temperature
         */
        LOG.fine("Validating temperature...");
        if (temperature <= 0) {
            LOG.severe("Invalid temperature: " + temperature);
            throw new IllegalArgumentException("Temperature must be positive: " + temperature);
        }
        LOG.fine("✓ Temperature valid");

        /*
         * 4. Extract Stage 1-3 data from AllClusterData
         */
        LOG.fine("STAGE 4a: Extract Stage 1 (Cluster Identification)...");
        ClusterIdentificationResult stage1 = clusterData.getDisorderedClusterResult();
        LOG.fine("  ✓ tcdis=" + stage1.getTcdis() + " (cluster types)");
        LOG.fine("  ✓ kb coefficients, mh multiplicities, lc counts loaded");

        LOG.fine("STAGE 4b: Extract Stage 2 (CF Identification)...");
        CFIdentificationResult stage2 = clusterData.getDisorderedCFResult();
        LOG.fine("  ✓ tcf=" + stage2.getTcf() + ", ncf=" + stage2.getNcf() + " (CF counts)");
        LOG.fine("  ✓ lcf array, CF basis indices loaded");

        LOG.fine("STAGE 4c: Extract Stage 3 (C-Matrix - CVCF basis)...");
        org.ce.domain.cluster.CMatrix.Result clusterCMatrix = clusterData.getCMatrixResult();
        LOG.fine("  ✓ cmat[t][j][v][k] transformation matrix loaded");
        LOG.fine("  ✓ lcv (CV counts), wcv (CV weights) loaded");
        LOG.fine("  ✓ cfBasisIndices: " + (clusterCMatrix.getCfBasisIndices() == null ? "null (CVCF)" : "present"));

        /*
         * 5. Create CVMInput with complete topology (Stages 1-3)
         */
        LOG.fine("STAGE 4d: Create CVMInput bundle...");
        // TODO: resolve basis from systemId + numComponents via registry (future)
        CVMInput cvmInput = new CVMInput(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                clusterCMatrix,
                systemId,
                systemName,
                composition.length,
                basis
        );
        LOG.fine("  ✓ CVMInput created with numComponents=" + composition.length);

        /*
         * 6. Evaluate ECI at temperature using CVCF names.
         *    ECIs are indexed by CF name (e4AB → v4AB at col 0, etc.).
         *    Missing terms default to 0 (direct inheritance property of CVCF basis).
        */
        LOG.fine("Evaluating ECI at T=" + temperature + " K...");
        double[] eci = CECEvaluator.evaluate(cec, temperature, basis, "CVM");
        LOG.fine("✓ ECI evaluated (" + eci.length + " non-point terms)");

        /*
         * 7. Create CVMPhaseModel and run N-R minimization
         */
        LOG.info("Running CVM N-R minimization...");
        CVMPhaseModel model = CVMPhaseModel.create(
                cvmInput,
                eci,
                temperature,
                composition
        );
        LOG.info("✓ CVM minimization converged in " + model.getLastIterations()
                + " iterations (||Gu||=" + String.format("%.4e", model.getLastGradientNorm()) + ")");
        LOG.info("  G(eq) = " + String.format("%.8e", model.getEquilibriumG()) + " J/mol");
        LOG.info("  H(eq) = " + String.format("%.8e", model.getEquilibriumH()) + " J/mol");
        LOG.info("  S(eq) = " + String.format("%.8e", model.getEquilibriumS()) + " J/(mol·K)");

        /*
         * 8. Emit post-hoc N-R iteration trace as structured chart events.
         *    CVMPhaseModel.create() throws on convergence failure, so this block
         *    only runs when the solver succeeded — the chart always shows a converged trace.
         */
        if (input.eventSink != null) {
            LOG.fine("Emitting N-R iteration trace to eventSink...");
            input.eventSink.accept(new ProgressEvent.EngineStart("CVM", 0));
            List<SolverResult.IterationSnapshot> trace = model.getLastIterationTrace();
            for (SolverResult.IterationSnapshot snap : trace) {
                double[] ghs = model.computeGHS(snap.getCf());
                input.eventSink.accept(new ProgressEvent.CvmIteration(
                        snap.getIteration(),
                        snap.getGibbsEnergy(),
                        snap.getGradientNorm(),
                        ghs[1],      // enthalpy (H)
                        ghs[2],      // entropy (S) in J/(mol·K)
                        snap.getCf())); // CFs for logging
            }
            LOG.fine("✓ Emitted " + trace.size() + " iteration snapshots");
        }

        EquilibriumState result = new EquilibriumState(
            temperature,
            composition.clone(),
            model.getEquilibriumH(),
            model.getEquilibriumG(),
            model.getEquilibriumS(),
            Double.NaN, // enthalpyStdErr (MCS only)
            Double.NaN, // heatCapacity (MCS only)
            model.getEquilibriumCFs(),
            null,       // avgCFs (MCS only)
            null        // stdCFs (MCS only)
        );
        LOG.info("=== CVMEngine.compute() SUCCESS ===");
        return result;
    }

    private static void validateCmatEciConsistency(CMatrix.Result cmat, CvCfBasis basis) {
        List<String> cmatNames = cmat.getCmatCfNames();
        if (cmatNames == null) return; // orthogonal path, skip

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
     * Formats a double array for logging.
     */
    private static String formatArray(double[] arr) {
        if (arr == null || arr.length == 0) return "";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            if (i > 0) sb.append(", ");
            sb.append(String.format("%.6f", arr[i]));
        }
        return sb.toString();
    }
}
