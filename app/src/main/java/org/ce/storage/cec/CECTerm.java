package org.ce.storage.cec;

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
 */
public class CECTerm {

    /**
     * Human-readable name of the correlation function
     * (e.g. "Pair rc=[2] mult=4").
     */
    public String name;

    /**
     * Constant coefficient.
     */
    public double a;

    /**
     * Temperature coefficient.
     */
    public double b;

}
