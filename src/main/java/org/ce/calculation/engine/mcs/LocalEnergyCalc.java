package org.ce.calculation.engine.mcs;

import org.ce.model.mcs.EmbeddingData;
import org.ce.model.mcs.LatticeConfig;
import org.ce.model.mcs.SiteOperatorBasis;
import org.ce.model.cluster.Cluster;

import java.util.List;

/** Static utility for computing cluster-expansion energy contributions. */
public final class LocalEnergyCalc {

    private LocalEnergyCalc() {}

    public static double clusterProduct(EmbeddingData.Embedding e, LatticeConfig config) {
        SiteOperatorBasis basis  = config.getBasis();
        double prod  = 1.0;
        int[]  idx    = e.getSiteIndices();
        int[]  alphas = e.getAlphaIndices();
        for (int k = 0; k < idx.length; k++)
            prod *= basis.evaluate(alphas[k], config.getOccupation(idx[k]));
        return prod;
    }

    public static double totalEnergy(LatticeConfig config,
                                     EmbeddingData emb,
                                     double[] eci,
                                     List<List<Cluster>> orbits) {
        double sum = 0.0;
        for (EmbeddingData.Embedding e : emb.getAllEmbeddings()) {
            int size = e.size();
            if (size > 0)
                sum += eci[e.getClusterType()] * clusterProduct(e, config) / size;
            else
                sum += eci[e.getClusterType()] * clusterProduct(e, config);
        }
        return sum;
    }

    public static double deltaESingleSite(int i,
                                          int newOcc,
                                          LatticeConfig config,
                                          EmbeddingData emb,
                                          double[] eci,
                                          List<List<Cluster>> orbits) {
        int oldOcc = config.getOccupation(i);
        if (oldOcc == newOcc) return 0.0;

        double dE    = 0.0;
        SiteOperatorBasis basis = config.getBasis();

        for (EmbeddingData.Embedding e : emb.getSiteToEmbeddings()[i]) {
            int    t      = e.getClusterType();
            int[]  idx    = e.getSiteIndices();
            int[]  alphas = e.getAlphaIndices();
            double restProduct = 1.0;
            int    alphaI = -1;

            for (int k = 0; k < idx.length; k++) {
                if (idx[k] == i) alphaI = alphas[k];
                else restProduct *= basis.evaluate(alphas[k], config.getOccupation(idx[k]));
            }

            if (alphaI < 0) continue;
            double phiOld = basis.evaluate(alphaI, oldOcc);
            double phiNew = basis.evaluate(alphaI, newOcc);
            int size = e.size();
            double energyCont = eci[t] * (phiNew - phiOld) * restProduct;
            dE += (size > 0) ? (energyCont / size) : energyCont;
        }
        return dE;
    }

    public static double deltaEExchange(int i, int j,
                                        LatticeConfig config,
                                        EmbeddingData emb,
                                        double[] eci,
                                        List<List<Cluster>> orbits) {
        int occI = config.getOccupation(i);
        int occJ = config.getOccupation(j);
        if (occI == occJ) return 0.0;

        double dEi = deltaESingleSite(i, occJ, config, emb, eci, orbits);

        config.setOccupation(i, occJ);
        double dEj;
        try {
            dEj = deltaESingleSite(j, occI, config, emb, eci, orbits);
        } finally {
            config.setOccupation(i, occI);
        }
        return dEi + dEj;
    }
}
