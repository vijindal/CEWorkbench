package org.ce.domain.cluster;

import org.ce.domain.cluster.Cluster;
import java.util.List;

public class ClassifiedClusterResult {

    private final List<List<Cluster>> coordList;
    private final List<List<Double>> multiplicityList;
    private final List<List<List<Cluster>>> orbitList;
    private final List<List<List<Integer>>> rcList;

    public ClassifiedClusterResult(
            List<List<Cluster>> coordList,
            List<List<Double>> multiplicityList,
            List<List<List<Cluster>>> orbitList,
            List<List<List<Integer>>> rcList) {

        this.coordList = coordList;
        this.multiplicityList = multiplicityList;
        this.orbitList = orbitList;
        this.rcList = rcList;
    }

    public List<List<Cluster>> getCoordList() { return coordList; }
    public List<List<Double>> getMultiplicityList() { return multiplicityList; }
    public List<List<List<Cluster>>> getOrbitList() { return orbitList; }
    public List<List<List<Integer>>> getRcList() { return rcList; }
}



