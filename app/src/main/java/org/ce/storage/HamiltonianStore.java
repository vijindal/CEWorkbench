package org.ce.storage;

/**
 * Storage service for Hamiltonian parameters (ECI).
 */
public class HamiltonianStore {

    private final Workspace workspace;

    public HamiltonianStore(Workspace workspace) {
        this.workspace = workspace;
    }

    public double[] loadECI(String systemId) {

        // Placeholder implementation.
        // Will later load ECI parameters from JSON.

        return new double[0];
    }
}
