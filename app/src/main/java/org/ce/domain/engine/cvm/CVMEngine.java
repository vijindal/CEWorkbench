package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicResult;
import org.ce.domain.engine.cvm.CVMPhaseModelExecutor;
import org.ce.storage.cec.CECEntry;

/**
 * CVM thermodynamic engine.
 *
 * This class adapts the legacy CVM implementation to the
 * new ThermodynamicEngine interface.
 */
public class CVMEngine implements ThermodynamicEngine {

    @Override
    public ThermodynamicResult compute(
            AllClusterData clusterData,
            CECEntry cec,
            double temperature,
            double[] composition) {

        double[] eci = computeECI(cec, temperature);

        /*
         * Create CVM solver
         */
        CVMPhaseModelExecutor executor = new CVMPhaseModelExecutor();

        /*
         * Run CVM calculation
         */
        double freeEnergy = executor.compute(
                clusterData,
                eci,
                temperature,
                composition
        );

        /*
         * For now we return only free energy.
         * Later we will also return equilibrium cluster probabilities.
         */
        return new ThermodynamicResult(
                freeEnergy,
                new double[0]
        );
    }

    private double[] computeECI(CECEntry cec, double temperature) {

        int n = cec.cecTerms.length;

        double[] eci = new double[n];

        for (int i = 0; i < n; i++) {

            double a = cec.cecTerms[i].a;
            double b = cec.cecTerms[i].b;

            eci[i] = a + b * temperature;
        }

        return eci;
    }
}
