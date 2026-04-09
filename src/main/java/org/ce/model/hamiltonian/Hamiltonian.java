package org.ce.domain.hamiltonian;

/**
 * Represents a thermodynamic Hamiltonian used by engines (CVM / MCS).
 *
 * Implementations provide the energy evaluation for a given
 * set of correlation functions.
 */
public interface Hamiltonian {

    /**
     * Number of interaction terms (ECIs).
     */
    int size();

    /**
     * Returns a copy of the ECI vector.
     */
    double[] getECI();

    /**
     * Computes the configurational energy.
     *
     * E = Σ Jα Φα
     *
     * @param correlationFunctions correlation function vector
     * @return energy
     */
    double energy(double[] correlationFunctions);
}
