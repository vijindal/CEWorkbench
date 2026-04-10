package org.ce.model.mcs;

import java.util.ArrayList;
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
}
