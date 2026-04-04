package org.ce.domain.cluster;

import java.util.ArrayList;
import java.util.List;

public class CFGroupGenerator {

    public static GroupedCFResult groupCFData(
            ClusCoordListResult disClusData,
            ClusCoordListResult disCFData,
            ClassifiedClusterResult ordCFData,
            List<String> basisSymbolList) {

        List<Cluster> disClusCoordData =
                disClusData.getClusCoordList();

        List<Cluster> disCFCoordData =
                disCFData.getClusCoordList();

        List<List<Cluster>> disClusOrbitList =
                disClusData.getOrbitList();

        List<List<Cluster>> ordCFCoordData =
                ordCFData.getCoordList();

        List<List<Double>> ordCFMData =
                ordCFData.getMultiplicityList();

        List<List<List<Cluster>>> ordCFOrbitData =
                ordCFData.getOrbitList();

        List<List<List<Integer>>> ordCFRData =
                ordCFData.getRcList();

        // ---------------------------------------------------
        // STEP 1: Remove decorations from disordered CF clusters
        // Equivalent to transDisCFCoordData in Mathematica
        // ---------------------------------------------------

        List<Cluster> transDisCFCoordData = new ArrayList<>();

        for (Cluster cf : disCFCoordData) {

            List<Sublattice> newSubs = new ArrayList<>();

            for (Sublattice sub : cf.getSublattices()) {

                List<Site> newSites = new ArrayList<>();

                for (Site s : sub.getSites()) {

                    newSites.add(
                        new Site(
                            s.getPosition(),
                            basisSymbolList.get(0) // canonical symbol
                        )
                    );
                }

                newSubs.add(new Sublattice(newSites));
            }

            transDisCFCoordData.add(new Cluster(newSubs));
        }

        // ---------------------------------------------------
        // STEP 2: Group ordered CFs into disordered clusters
        // ---------------------------------------------------

        List<List<List<Cluster>>> groupedCFCoordData =
                new ArrayList<>();

        List<List<List<Double>>> groupedCFMData =
                new ArrayList<>();

        List<List<List<List<Cluster>>>> groupedCFOrbitData =
                new ArrayList<>();

        List<List<List<List<Integer>>>> groupedCFRData =
                new ArrayList<>();

        int tcdis = disClusCoordData.size();

        for (int i = 0; i < tcdis; i++) {

            groupedCFCoordData.add(new ArrayList<>());
            groupedCFMData.add(new ArrayList<>());
            groupedCFOrbitData.add(new ArrayList<>());
            groupedCFRData.add(new ArrayList<>());

            List<Cluster> disOrbit =
                    disClusOrbitList.get(i);

            for (int j = 0; j < transDisCFCoordData.size(); j++) {

                Cluster undecoratedCF =
                        transDisCFCoordData.get(j);

                if (OrbitUtils.isContained(disOrbit, undecoratedCF)) {

                    groupedCFCoordData.get(i)
                            .add(ordCFCoordData.get(j));

                    groupedCFMData.get(i)
                            .add(ordCFMData.get(j));

                    groupedCFOrbitData.get(i)
                            .add(ordCFOrbitData.get(j));

                    groupedCFRData.get(i)
                            .add(ordCFRData.get(j));
                }
            }
        }

        return new GroupedCFResult(
                groupedCFCoordData,
                groupedCFMData,
                groupedCFOrbitData,
                groupedCFRData
        );
    }
}




