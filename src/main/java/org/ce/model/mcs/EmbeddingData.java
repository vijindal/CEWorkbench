package org.ce.model.mcs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Container for all Embedding instances generated for a supercell. */
public class EmbeddingData {

    private final List<Embedding>   allEmbeddings;
    private final List<Embedding>[] siteToEmbeddings;

    public EmbeddingData(List<Embedding> allEmbeddings, List<Embedding>[] siteToEmbeddings) {
        this.allEmbeddings    = allEmbeddings;
        this.siteToEmbeddings = siteToEmbeddings;
    }

    public List<Embedding>   getAllEmbeddings()    { return allEmbeddings; }
    public List<Embedding>[] getSiteToEmbeddings() { return siteToEmbeddings; }
    public int totalEmbeddingCount()               { return allEmbeddings.size(); }
    public int siteCount()                         { return siteToEmbeddings.length; }

    public double[] computeHmixCoeff(double[] eci, int tc) {
        int N = siteCount();
        int[] counts = new int[tc];
        int[] sizes  = new int[tc];
        for (Embedding e : allEmbeddings) {
            int t = e.getClusterType();
            if (t < tc) { counts[t]++; sizes[t] = e.size(); }
        }
        double[] coeff = new double[tc];
        for (int t = 0; t < tc && t < eci.length; t++) {
            if (sizes[t] > 1) coeff[t] = eci[t] * counts[t] / ((double) sizes[t] * N);
        }
        return coeff;
    }

    public int[] multiSiteEmbedCountsPerType(int tc) {
        int[] counts = new int[tc];
        for (Embedding e : allEmbeddings) {
            int t = e.getClusterType();
            if (t < tc && e.size() > 1) counts[t]++;
        }
        return counts;
    }

    public List<Embedding> getEmbeddingsForTypeAtSite(int clusterType, int siteIndex) {
        List<Embedding> result = new ArrayList<>();
        for (Embedding e : siteToEmbeddings[siteIndex]) {
            if (e.getClusterType() == clusterType) result.add(e);
        }
        return result;
    }

    /** A single embedding of an abstract cluster type onto specific lattice sites. */
    public static class Embedding {

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
}
