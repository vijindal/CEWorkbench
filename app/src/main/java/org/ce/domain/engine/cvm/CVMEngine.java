package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CMatrixBuilder;
import org.ce.domain.cluster.CMatrixResult;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.result.EquilibriumState;
import org.ce.storage.cec.CECEntry;

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
         * 1. Build C-Matrix (Stage 3)
         * TODO: maxClusters - determine source for this parameter
         */
        CMatrixResult cmatrix = CMatrixBuilder.build(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                null,  // TODO: maxClusters
                composition.length
        );

        /*
         * 2. Create CVMInput with complete topology (Stages 1-3)
         */
        CVMInput cvmInput = new CVMInput(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                cmatrix,
                systemId,
                systemName,
                composition.length
        );

        /*
         * 3. Evaluate ECI at temperature: eci[l] = a[l] + b[l] * T
         */
        double[] eci = evaluateECI(cec, temperature);

        /*
         * 4. Validate ECI length matches ncf (prevents silent truncation)
         */
        int ncf = cvmInput.getStage2().getNcf();
        if (eci.length < ncf) {
            throw new IllegalArgumentException(
                    "ECI array length (" + eci.length + ") < ncf (" + ncf + "). " +
                    "CEC database is incomplete or mismatched with cluster data.");
        }

        /*
         * 5. Create CVMPhaseModel and run minimization
         */
        CVMPhaseModel model = CVMPhaseModel.create(
                cvmInput,
                eci,
                temperature,
                composition
        );

        /*
         * 6. Initialize model
         */
        boolean initOk = CVMPhaseModelExecutor.initializeModel(model);
        if (!initOk) {
            throw new RuntimeException("CVM initialization failed");
        }

        /*
         * 7. Query equilibrium and return state
         */
        EquilibriumState state = CVMPhaseModelExecutor.queryModel(model);
        if (state == null) {
            throw new RuntimeException("CVM query failed");
        }

        return state;
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
