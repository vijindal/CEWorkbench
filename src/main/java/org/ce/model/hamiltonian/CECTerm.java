package org.ce.model.hamiltonian;

import static org.ce.model.cluster.ClusterPrimitives.*;

/**
 * Represents one Cluster Expansion coefficient term in the CEC database.
 *
 * Energy coefficient is temperature dependent:
 *
 *      J(T) = a + bT
 *
 * where
 *      a = constant coefficient
 *      b = temperature coefficient
 *
 * <p>Optional metadata fields (numSites, multiplicity, description) are
 * populated during scaffolding to help users identify which CF represents
 * which cluster type, and can be edited by users for documentation.</p>
 */
public class CECTerm {

    /**
     * Human-readable name of the correlation function
     * (e.g. "CF_0", "CF_1").
     */
    public String name;

    /**
     * Optional user-editable description identifying the cluster type
     * (e.g. "Point cluster (1 site, multiplicity 1)").
     */
    public String description;

    /**
     * Optional number of sites in this correlation function's cluster.
     * Populated during scaffolding to help users identify CF types.
     */
    public Integer numSites;

    /**
     * Optional multiplicity of this correlation function.
     * Populated during scaffolding to help users identify CF types.
     */
    public Double multiplicity;

    /**
     * Constant coefficient.
     */
    public double a;

    /**
     * Temperature coefficient.
     */
    public double b;

    /**
     * Validates that {@link #a} and {@link #b} are finite numbers.
     *
     * @throws IllegalStateException if either value is NaN or infinite
     */
    public void validate() {
        if (!Double.isFinite(a)) throw new IllegalStateException(
                "CECTerm '" + name + "': a=" + a + " is not finite");
        if (!Double.isFinite(b)) throw new IllegalStateException(
                "CECTerm '" + name + "': b=" + b + " is not finite");
    }

}
