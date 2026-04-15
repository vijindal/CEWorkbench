package org.ce.calculation.workflow.thermo;

import org.ce.model.cvm.CVMGibbsModel;
import java.io.IOException;

/**
 * Procedural bridge for specific thermodynamic algorithms.
 * Wraps the modern model-layer evaluators to provide property-specific
 * methods (calHm, calSm, calGm) as requested.
 */
public class ThermodynamicMethods {

    private final CVMGibbsModel model;

    public ThermodynamicMethods(CVMGibbsModel model) {
        this.model = model;
    }

    /**
     * Calculates Enthalpy (Hm) at given conditions.
     * modDataIn: [output, T, x]
     */
    public void calHm(double[] modDataIn) throws IOException {
        double T_local = modDataIn[1];
        double x_local = modDataIn[2];
        
        // Using binary-friendly setter if composition is passed as a single value
        double[] composition = (modDataIn.length > 3) 
            ? java.util.Arrays.copyOfRange(modDataIn, 2, modDataIn.length)
            : new double[]{1.0 - x_local, x_local};

        model.setT(T_local);
        model.setX(composition);
        
        // Ensure minimization occurs if values are stale
        model.getEquilibriumState(T_local, composition, 1e-5, null, null);
        
        modDataIn[0] = model.calH();
    }

    /**
     * Calculates Entropy (Sm) at given conditions.
     * modDataIn: [output, T, x]
     */
    public void calSm(double[] modDataIn) throws IOException {
        double T_local = modDataIn[1];
        double x_local = modDataIn[2];
        
        double[] composition = (modDataIn.length > 3) 
            ? java.util.Arrays.copyOfRange(modDataIn, 2, modDataIn.length)
            : new double[]{1.0 - x_local, x_local};

        model.setT(T_local);
        model.setX(composition);
        
        model.getEquilibriumState(T_local, composition, 1e-5, null, null);
        
        modDataIn[0] = model.calS();
    }

    /**
     * Calculates Gibbs Free Energy (Gm) at given conditions.
     * modDataIn: [output, T, x]
     */
    public void calGm(double[] modDataIn) throws IOException {
        double T_local = modDataIn[1];
        double x_local = modDataIn[2];
        
        double[] composition = (modDataIn.length > 3) 
            ? java.util.Arrays.copyOfRange(modDataIn, 2, modDataIn.length)
            : new double[]{1.0 - x_local, x_local};

        model.setT(T_local);
        model.setX(composition);
        
        model.getEquilibriumState(T_local, composition, 1e-5, null, null);
        
        modDataIn[0] = model.calG();
    }

    /**
     * (Stub) Iterative common tangent solver.
     */
    public void calMGX1(double[] modDataIn) throws IOException {
        // TODO: Implement iterative common tangent logic once 
        // chemical potentials are available in CVMGibbsModel.
        throw new UnsupportedOperationException("calMGX1 not yet implemented in this architecture.");
    }

    /**
     * (Stub) Consolute point calculation.
     */
    public void calCOPT(double[] modData) throws IOException {
        // TODO: Implement consolute point logic once 
        // higher-order derivatives (Gxx, Gxxx) are available in CVMGibbsModel.
        throw new UnsupportedOperationException("calCOPT not yet implemented in this architecture.");
    }
}
