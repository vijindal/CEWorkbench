package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.engine.ProgressEvent;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.engine.cvm.CVMPhaseModel.CVMInput;
import org.ce.domain.engine.cvm.NewtonRaphsonSolverSimple.CVMSolverResult;
import org.ce.domain.result.EquilibriumState;
import org.ce.domain.hamiltonian.CECEntry;

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
        CVMInput cvmInput = new CVMInput(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                clusterCMatrix,
                systemId,
                systemName,
                composition.length
        );

        /*
         * 3. Evaluate ECI at temperature: eci[l] = a[l] + b[l] * T
         */
        double[] eci = evaluateECI(cec, temperature);

        /*
         * 5. Validate ECI length matches expected ncf exactly
         */
        int expectedNcf = clusterData.getDisorderedCFResult().getNcf();

        if (eci.length != expectedNcf) {
            throw new IllegalStateException(
                    "ECI length (" + eci.length + ") does not match ncf (" + expectedNcf + ")"
            );
        }

        /*
         * 6. Validate CEC ncf field matches expected ncf
         * (enforces: G = Σ (ECI_l * CF_l) → length(ECI) == ncf == number of CFs)
         */
        if (cec.ncf != expectedNcf) {
            throw new IllegalArgumentException(
                    "CEC ncf field (" + cec.ncf + ") does not match cluster data ncf (" + expectedNcf + ")"
            );
        }

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
            double T = input.temperature;
            for (CVMSolverResult.IterationSnapshot snap : trace) {
                double[] ghs = model.computeGHS(snap.getCf());
                input.eventSink.accept(new ProgressEvent.CvmIteration(
                        snap.getIteration(),
                        snap.getGibbsEnergy(),
                        snap.getGradientNorm(),
                        ghs[1],      // enthalpy (H)
                        -T * ghs[2], // −T·S
                        null));      // cfs unused
            }
        }

        return model.getEquilibriumState();
    }

    /**
     * Evaluates ECI at given temperature using the formula: eci[l] = a[l] + b[l] * T
     */
    private double[] evaluateECI(CECEntry cec, double temperature) {

        if (cec.cecTerms == null || cec.cecTerms.length == 0) {
            throw new IllegalArgumentException("CECEntry has no terms");
        }

        double[] eci = new double[cec.cecTerms.length];

        for (int i = 0; i < cec.cecTerms.length; i++) {
            double a = cec.cecTerms[i].a;
            double b = cec.cecTerms[i].b;
            eci[i] = a + b * temperature;
        }

        return eci;
    }
}
