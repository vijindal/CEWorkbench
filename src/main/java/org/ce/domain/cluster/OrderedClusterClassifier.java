package org.ce.domain.cluster;

import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.Site;
import org.ce.domain.cluster.Sublattice;
import org.ce.domain.cluster.ClassifiedClusterResult;
import org.ce.domain.cluster.ClusCoordListResult;
import org.ce.domain.cluster.OrbitUtils;

import java.util.ArrayList;
import java.util.List;

public class OrderedClusterClassifier {

    public static ClassifiedClusterResult classify(
            ClusCoordListResult disClusData,
            ClusCoordListResult clusData,
            List<Cluster> transformedClusCoordList) {

        List<List<Cluster>> disClusOrbitList =
                disClusData.getOrbitList();

        List<Cluster> clusCoordList1 =
                clusData.getClusCoordList();

        List<Double> clusMList =
                clusData.getMultiplicities();

        List<List<Cluster>> clusOrbitList =
                clusData.getOrbitList();

        List<List<Integer>> clusRcList =
                clusData.getRcList();

        int tcdis = disClusOrbitList.size();
        int tc = transformedClusCoordList.size();

        List<List<Cluster>> classifiedClusCoordList =
                new ArrayList<>();

        List<List<Double>> classifiedClusMList =
                new ArrayList<>();

        List<List<List<Cluster>>> classifiedClusOrbitList =
                new ArrayList<>();

        List<List<List<Integer>>> classifiedClusRcList =
                new ArrayList<>();

        for (int j = 0; j < tcdis; j++) {

            classifiedClusCoordList.add(new ArrayList<>());
            classifiedClusMList.add(new ArrayList<>());
            classifiedClusOrbitList.add(new ArrayList<>());
            classifiedClusRcList.add(new ArrayList<>());

            for (int i = 0; i < tc; i++) {

                // ðŸ”¹ Mathematica equivalent:
                // flattenClusCoord = sortClusCoord[Flatten[clusCoordList[[i]],1]]

                Cluster flattened =
                        flattenCluster(transformedClusCoordList.get(i));

                if (OrbitUtils.isContained(
                        disClusOrbitList.get(j),
                        flattened)) {

                    classifiedClusCoordList.get(j)
                            .add(clusCoordList1.get(i));

                    classifiedClusMList.get(j)
                            .add(clusMList.get(i));

                    classifiedClusOrbitList.get(j)
                            .add(clusOrbitList.get(i));

                    classifiedClusRcList.get(j)
                            .add(clusRcList.get(i));
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

    // =========================================================
    // Equivalent to:
    // sortClusCoord[Flatten[clusCoord,1]]
    // =========================================================
    private static Cluster flattenCluster(Cluster cluster) {

        // 1ï¸âƒ£ Flatten
        List<Site> flatSites =
                new ArrayList<>(cluster.getAllSites());

        // 2ï¸âƒ£ Sort (Mathematica sortClusCoord equivalent)
        flatSites.sort(Cluster::compareSites);

        // 3ï¸âƒ£ Wrap into single sublattice
        List<Sublattice> singleSub = new ArrayList<>();
        singleSub.add(new Sublattice(flatSites));

        return new Cluster(singleSub);
    }
}





