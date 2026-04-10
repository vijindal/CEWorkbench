package org.ce.calculation.workflow.cec;

/**
 * Holds metadata about a correlation function extracted from cluster data.
 * Used during Hamiltonian scaffolding to enrich CECTerm with descriptive information.
 */
public class CFMetadata {

    public int numSites;
    public double multiplicity;
    public String description;

    public CFMetadata(int numSites, double multiplicity, String description) {
        this.numSites = numSites;
        this.multiplicity = multiplicity;
        this.description = description;
    }
}
