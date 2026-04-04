package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import org.ce.domain.cluster.Cluster;
import java.util.List;

public class ClassifiedClusterResult {

    private final List<List<Cluster>> coordList;
    private final List<List<Double>> multiplicityList;
    private final List<List<List<Cluster>>> orbitList;
    private final List<List<List<Integer>>> rcList;

    @JsonCreator
    public ClassifiedClusterResult(
            @JsonProperty("coordList")        List<List<Cluster>>        coordList,
            @JsonProperty("multiplicityList") List<List<Double>>          multiplicityList,
            @JsonProperty("orbitList")        List<List<List<Cluster>>>   orbitList,
            @JsonProperty("rcList")           List<List<List<Integer>>>   rcList) {

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



