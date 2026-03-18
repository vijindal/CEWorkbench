package org.ce.domain.hamiltonian;

/**
 * Represents a Cluster Expansion Coefficient (CEC) database entry.
 *
 * This corresponds directly to the contents of a cec.json file.
 */
public class CECEntry {

    /**
     * Elements string (example: "Nb-Ti").
     */
    public String elements;

    /**
     * Structure + phase identifier (example: "BCC_A2").
     */
    public String structurePhase;

    /**
     * Model identifier (example: "T").
     */
    public String model;

    /**
     * Array of cluster expansion terms.
     */
    public CECTerm[] cecTerms;

    /**
     * Units of the coefficients.
     */
    public String cecUnits;

    /**
     * Reference for the parameterization.
     */
    public String reference;

    /**
     * Additional notes.
     */
    public String notes;

    /**
     * Number of correlation functions (ncf).
     */
    public int ncf;

}
