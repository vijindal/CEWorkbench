package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.Arrays;

/**
 * Immutable container for p-operator expansion coefficients.
 */
public final class PRules {

    private final int numSites;
    private final int numElements;
    private final double[][] rMat;

    private PRules(int numSites, int numElements, double[][] rMat) {
        this.numSites = numSites;
        this.numElements = numElements;
        this.rMat = copyMatrix(rMat);
    }

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

    public int getNumSites() {
        return numSites;
    }

    public int getNumElements() {
        return numElements;
    }

    /**
     * Returns coefficients for p(site, element) in the basis
     * [1, s1(site), s2(site), ...]. Indices are 0-based.
     */
    public double[] coefficientsFor(int siteIndex, int elementIndex) {
        if (siteIndex < 0 || siteIndex >= numSites) {
            throw new IllegalArgumentException("siteIndex out of range");
        }
        if (elementIndex < 0 || elementIndex >= numElements) {
            throw new IllegalArgumentException("elementIndex out of range");
        }
        return Arrays.copyOf(rMat[elementIndex], numElements);
    }

    /**
     * Returns the full R-matrix used by these rules.
     */
    public double[][] getRMatrix() {
        return copyMatrix(rMat);
    }

    private static double[][] copyMatrix(double[][] src) {
        double[][] copy = new double[src.length][];
        for (int i = 0; i < src.length; i++) {
            copy[i] = Arrays.copyOf(src[i], src[i].length);
        }
        return copy;
    }
}

