package org.ce.domain.cluster;

/**
 * Builds p-operator expansion rules using the Inden R-matrix.
 */
public final class PRulesBuilder {

    private PRulesBuilder() {}

    /**
     * Builds a PRules instance for the given number of sites and elements.
     *
     * @param numSites number of sites in the maximal cluster
     * @param numElements number of components
     * @return PRules container for p-operator expansions
     */
    public static PRules build(int numSites, int numElements) {
        if (numSites <= 0) {
            throw new IllegalArgumentException("numSites must be > 0");
        }
        if (numElements < 2) {
            throw new IllegalArgumentException("numElements must be >= 2");
        }

        double[][] rMat = RMatrixCalculator.buildRMatrix(numElements);
        return new PRules(numSites, numElements, rMat);
    }
}

