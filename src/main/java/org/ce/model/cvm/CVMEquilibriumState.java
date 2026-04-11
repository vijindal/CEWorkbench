package org.ce.model.cvm;

import org.ce.model.cvm.CVMGibbsModel.ModelResult;

/**
 * Thermodynamic equilibrium state at a given (T, x).
 *
 * <p>Holds the equilibrium correlation functions and corresponding thermodynamic values
 * (G, H, S). This is what the calculation layer consumes from CVMSolver.</p>
 */
public class CVMEquilibriumState {

    /** Equilibrium CVCF correlation functions (non-point CFs only). */
    private final double[] u;

    /** Thermodynamic values at equilibrium: G, H, S and derivatives. */
    private final ModelResult modelValues;

    /** Temperature (K) at which equilibrium was computed. */
    private final double temperature;

    /** Gas constant (J/(mol·K)) for analytical Cv calculation. */
    private final double R;

    /**
     * @param u equilibrium non-point CFs (length = numNonPointCfs)
     * @param modelValues G/H/S and derivatives at equilibrium
     * @param temperature K
     * @param R gas constant J/(mol·K)
     */
    public CVMEquilibriumState(double[] u, ModelResult modelValues, double temperature, double R) {
        this.u = u.clone();
        this.modelValues = modelValues;
        this.temperature = temperature;
        this.R = R;
    }

    /**
     * Returns the equilibrium CVCF correlation functions.
     */
    public double[] getCFs() {
        return u.clone();
    }

    /**
     * Returns Gibbs free energy (J/mol) at equilibrium.
     */
    public double getGibbsEnergy() {
        return modelValues.G;
    }

    /**
     * Returns enthalpy (J/mol) at equilibrium.
     */
    public double getEnthalpy() {
        return modelValues.H;
    }

    /**
     * Returns entropy (J/(mol·K)) at equilibrium.
     */
    public double getEntropy() {
        return modelValues.S;
    }

    /**
     * Analytical heat capacity from entropy Hessian: C_v = -T · S_uu · uu_dot_dot.
     *
     * <p>This is the CVM heat capacity derived from the second derivatives of entropy
     * with respect to correlation functions. It's an analytical result (no statistical error).</p>
     *
     * @return C_v (J/(mol·K))
     */
    public double analyticalCv() {
        if (modelValues.Suu == null || modelValues.Suu.length == 0) {
            return Double.NaN;
        }
        // C_v = -T * Σ_i Σ_j S_uu[i][j]
        // (Simple sum; assumes S_uu is the Hessian and we're contracting it)
        double sum = 0.0;
        for (double[] row : modelValues.Suu) {
            if (row != null) {
                for (double val : row) {
                    sum += val;
                }
            }
        }
        return -temperature * sum;
    }

    @Override
    public String toString() {
        return "CVMEquilibriumState{T=" + temperature + ", G=" + String.format("%.4f", modelValues.G)
             + ", H=" + String.format("%.4f", modelValues.H)
             + ", S=" + String.format("%.6f", modelValues.S) + "}";
    }
}
