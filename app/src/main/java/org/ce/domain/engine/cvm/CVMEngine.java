package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.cvcf.BccA2CvCfTransformations;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.cluster.cvcf.CvCfCFRegistry;
import org.ce.domain.engine.ProgressEvent;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.engine.cvm.CVMPhaseModel.CVMInput;
import org.ce.domain.engine.cvm.NewtonRaphsonSolverSimple.CVMSolverResult;
import org.ce.domain.result.EquilibriumState;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;

import java.util.List;

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

    @Override
    public EquilibriumState compute(ThermodynamicInput input) throws Exception {

        AllClusterData clusterData = input.clusterData;
        CECEntry cec = input.cec;
        double temperature = input.temperature;
        double[] composition = input.composition;
        String systemId = input.systemId;
        String systemName = input.systemName;

        /*
         * 1. Validate CMatrix exists
         */
        if (clusterData.getCMatrixResult() == null) {
            throw new IllegalStateException(
                "CMatrix not found in AllClusterData. " +
                "Ensure ClusterIdentificationWorkflow generated Stage 3."
            );
        }

        /*
         * 2. Validate composition
         */
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

        /*
         * 3. Validate temperature
         */
        if (temperature <= 0) {
            throw new IllegalArgumentException("Temperature must be positive: " + temperature);
        }

        /*
         * 4. Use precomputed C-Matrix (from Type-1)
         */
        org.ce.domain.cluster.CMatrixResult clusterCMatrix = clusterData.getCMatrixResult();

        /*
         * 5. Create CVMInput with complete topology (Stages 1-3)
         */
        // TODO: resolve basis from systemId + numComponents via registry (future)
        CVMInput cvmInput = new CVMInput(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                clusterCMatrix,
                systemId,
                systemName,
                composition.length,
                BccA2CvCfTransformations.binaryBasis()
        );

        /*
         * 3. Evaluate ECI at temperature using CVCF names.
         *    ECIs are indexed by CF name (e4AB → v4AB at col 0, etc.).
         *    Missing terms default to 0 (direct inheritance property of CVCF basis).
         */
        CvCfBasis basis = BccA2CvCfTransformations.binaryBasis();
        double[] eci = evaluateECI(cec, temperature, basis);

        /*
         * 7. Create CVMPhaseModel and return equilibrium
         */
        CVMPhaseModel model = CVMPhaseModel.create(
                cvmInput,
                eci,
                temperature,
                composition
        );

        /*
         * 8. Emit post-hoc N-R iteration trace as structured chart events.
         *    CVMPhaseModel.create() throws on convergence failure, so this block
         *    only runs when the solver succeeded — the chart always shows a converged trace.
         */
        if (input.eventSink != null) {
            input.eventSink.accept(new ProgressEvent.EngineStart("CVM", 0));
            List<CVMSolverResult.IterationSnapshot> trace = model.getLastIterationTrace();
            for (CVMSolverResult.IterationSnapshot snap : trace) {
                double[] ghs = model.computeGHS(snap.getCf());
                input.eventSink.accept(new ProgressEvent.CvmIteration(
                        snap.getIteration(),
                        snap.getGibbsEnergy(),
                        snap.getGradientNorm(),
                        ghs[1],      // enthalpy (H)
                        ghs[2],      // entropy (S) in J/(mol·K)
                        snap.getCf())); // CFs for logging
            }
        }

        return model.getEquilibriumState();
    }

    /**
     * Evaluates ECI at given temperature using CVCF CF names.
     *
     * <p>Each term in the CEC entry has a name like {@code e4AB}, {@code e22AB}.
     * The leading {@code e} is replaced by {@code v} to get the CF name ({@code v4AB},
     * {@code v22AB}), which is looked up in the registry to find the column index.
     * Missing terms default to 0 — this is the direct-inheritance property of the
     * CVCF basis.</p>
     *
     * @param cec         CEC entry with named CVCF terms
     * @param temperature temperature in Kelvin
     * @param basis       CVCF basis defining ncf and CF name order
     * @return eci[ncf] array; entry l corresponds to CF basis.cfNames.get(l)
     */
    private double[] evaluateECI(CECEntry cec, double temperature, CvCfBasis basis) {

        int ncf = basis.numNonPointCfs;
        double[] eci = new double[ncf];  // defaults to 0 (missing ECIs = no interaction)

        if (cec.cecTerms == null || cec.cecTerms.length == 0) {
            return eci;
        }

        // TODO: resolve registry from basis (future — when ternary/quaternary supported)
        CvCfCFRegistry registry = CvCfCFRegistry.forBccA2Binary();

        for (CECTerm term : cec.cecTerms) {
            if (term.name == null || !term.name.startsWith("e")) continue;
            String cfName = "v" + term.name.substring(1);  // e4AB → v4AB
            if (!registry.contains(cfName)) continue;
            int idx = registry.indexOf(cfName);
            if (idx < ncf) {
                eci[idx] = term.a + term.b * temperature;
            }
        }

        return eci;
    }
}
