package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterResults.*;
import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Unified utility for ordered-phase cluster operations.
 */
public class OrderedClusterOps {

    private OrderedClusterOps() {}

    // =========================================================
    // OrderedClusterClassifier logic
    // =========================================================
    public static ClassifiedClusterResult classify(
            ClusCoordListResult disClusData,
            ClusCoordListResult clusData,
            List<Cluster> transformedClusCoordList) {

        List<List<Cluster>> disClusOrbitList = disClusData.getOrbitList();
        List<Cluster> clusCoordList1 = clusData.getClusCoordList();
        List<Double> clusMList = clusData.getMultiplicities();
        List<List<Cluster>> clusOrbitList = clusData.getOrbitList();
        List<List<Integer>> clusRcList = clusData.getRcList();

        int tcdis = disClusOrbitList.size();
        int tc = transformedClusCoordList.size();

        List<List<Cluster>> classifiedClusCoordList = new ArrayList<>();
        List<List<Double>> classifiedClusMList = new ArrayList<>();
        List<List<List<Cluster>>> classifiedClusOrbitList = new ArrayList<>();
        List<List<List<Integer>>> classifiedClusRcList = new ArrayList<>();

        for (int j = 0; j < tcdis; j++) {
            classifiedClusCoordList.add(new ArrayList<>());
            classifiedClusMList.add(new ArrayList<>());
            classifiedClusOrbitList.add(new ArrayList<>());
            classifiedClusRcList.add(new ArrayList<>());

            for (int i = 0; i < tc; i++) {
                Cluster flattened = flattenCluster(transformedClusCoordList.get(i));

                if (ClusterUtils.isContained(disClusOrbitList.get(j), flattened)) {
                    classifiedClusCoordList.get(j).add(clusCoordList1.get(i));
                    classifiedClusMList.get(j).add(clusMList.get(i));
                    classifiedClusOrbitList.get(j).add(clusOrbitList.get(i));
                    classifiedClusRcList.get(j).add(clusRcList.get(i));
                }
            }
        }

        return new ClassifiedClusterResult(
                classifiedClusCoordList,
                classifiedClusMList,
                classifiedClusOrbitList,
                classifiedClusRcList
        );
    }

    private static Cluster flattenCluster(Cluster cluster) {
        // 1. Flatten
        List<Site> flatSites = new ArrayList<>(cluster.getAllSites());
        // 2. Sort (Mathematica sortClusCoord equivalent)
        flatSites.sort(Cluster::compareSites);
        // 3. Wrap into single sublattice
        List<Sublattice> singleSub = new ArrayList<>();
        singleSub.add(new Sublattice(flatSites));
        return new Cluster(singleSub);
    }

    // =========================================================
    // OrderedToDisorderedTransformer logic
    // =========================================================
    /**
     * Transforms a list of clusters using the given rotation matrix and translation vector.
     *
     * @param rotateMat 3x3 rotation matrix
     * @param translateMat 3-element translation vector
     * @param clusCoordCart clusters to transform
     * @return transformed clusters with new coordinates
     */
    public static List<Cluster> transformToDisordered(
            double[][] rotateMat,
            double[] translateMat,
            List<Cluster> clusCoordCart) {

        List<Cluster> result = new ArrayList<>();

        for (Cluster cluster : clusCoordCart) {
            List<Sublattice> newSublattices = new ArrayList<>();
            for (Sublattice sub : cluster.getSublattices()) {
                List<Site> newSites = new ArrayList<>();
                for (Site site : sub.getSites()) {
                    Position r = site.getPosition();

                    double x = rotateMat[0][0] * r.getX() +
                               rotateMat[0][1] * r.getY() +
                               rotateMat[0][2] * r.getZ() +
                               translateMat[0];

                    double y = rotateMat[1][0] * r.getX() +
                               rotateMat[1][1] * r.getY() +
                               rotateMat[1][2] * r.getZ() +
                               translateMat[1];

                    double z = rotateMat[2][0] * r.getX() +
                               rotateMat[2][1] * r.getY() +
                               rotateMat[2][2] * r.getZ() +
                               translateMat[2];

                    newSites.add(new Site(new Position(x, y, z), site.getSymbol()));
                }
                newSublattices.add(new Sublattice(newSites));
            }
            result.add(new Cluster(newSublattices));
        }

        return result;
    }
}
