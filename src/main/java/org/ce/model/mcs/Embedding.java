package org.ce.model.mcs;

import java.util.Arrays;

/** A single embedding of an abstract cluster type onto specific lattice sites. */
public class Embedding {

    private final int   clusterType;
    private final int   orbitMemberIndex;
    private final int[] siteIndices;
    private final int[] alphaIndices;

    public Embedding(int clusterType, int orbitMemberIndex, int[] siteIndices, int[] alphaIndices) {
        this.clusterType      = clusterType;
        this.orbitMemberIndex = orbitMemberIndex;
        this.siteIndices      = siteIndices;
        this.alphaIndices     = alphaIndices;
    }

    public int   getClusterType()      { return clusterType; }
    public int   getOrbitMemberIndex() { return orbitMemberIndex; }
    public int[] getSiteIndices()      { return siteIndices; }
    public int   size()                { return siteIndices.length; }
    public int[] getAlphaIndices()     { return alphaIndices; }

    @Override
    public String toString() {
        return "Embedding{type=" + clusterType + ", orbit=" + orbitMemberIndex
             + ", sites=" + Arrays.toString(siteIndices) + "}";
    }
}
