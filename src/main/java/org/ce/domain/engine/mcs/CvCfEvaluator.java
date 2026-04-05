package org.ce.domain.engine.mcs;

import java.util.List;

/**
 * Utility for evaluating CVCF cluster variables directly from configuration.
 */
public final class CvCfEvaluator {

    private CvCfEvaluator() {}

    /**
     * Measures CVCF cluster variables directly from the current lattice configuration
     * by accumulating embedding products for each decorated CF pattern.
     *
     * @param config        current lattice configuration
     * @param cfEmbeddings  per-CF-column embedding lists (built once at startup)
     * @param basisMatrix   basis values: basisMatrix[occ][alpha-1]
     * @param ncf           number of non-point CFs to measure
     * @return v[0..ncf-1]: measured CVCF cluster variables
     */
    public static double[] measureCVsFromConfig(
            LatticeConfig config,
            List<List<Embedding>> cfEmbeddings,
            double[][] basisMatrix,
            int ncf) {

        double[] v = new double[ncf];
        for (int l = 0; l < ncf; l++) {
            List<Embedding> embs = cfEmbeddings.get(l);
            if (embs == null || embs.isEmpty()) {
                v[l] = 0.0;
                continue;
            }

            double sum = 0.0;
            for (Embedding e : embs) {
                int[] sites  = e.getSiteIndices();
                int[] alphas = e.getAlphaIndices();
                double prod  = 1.0;
                for (int k = 0; k < sites.length; k++) {
                    int occ = config.getOccupation(sites[k]);
                    // basisMatrix[occ][alpha-1]
                    prod *= basisMatrix[occ][alphas[k] - 1];
                }
                sum += prod;
            }
            v[l] = sum / embs.size();
        }
        return v;
    }

    /**
     * Helper to build the basis value lookup table needed for measureCVsFromConfig.
     */
    public static double[][] buildBasisValues(int numComp) {
        double[] sequence = org.ce.domain.cluster.ClusterMath.buildBasis(numComp);
        double[][] values = new double[numComp][numComp - 1];
        for (int occ = 0; occ < numComp; occ++) {
            for (int alpha = 1; alpha <= numComp - 1; alpha++) {
                values[occ][alpha - 1] = Math.pow(sequence[occ], alpha);
            }
        }
        return values;
    }
}
