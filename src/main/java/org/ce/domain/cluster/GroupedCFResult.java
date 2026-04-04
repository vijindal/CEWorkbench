package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;

public class GroupedCFResult {

    private final List<List<List<Cluster>>> coordData;
    private final List<List<List<Double>>> multiplicityData;
    private final List<List<List<List<Cluster>>>> orbitData;
    private final List<List<List<List<Integer>>>> rcData;

    @JsonCreator
    public GroupedCFResult(
            @JsonProperty("coordData")        List<List<List<Cluster>>>          coordData,
            @JsonProperty("multiplicityData") List<List<List<Double>>>            multiplicityData,
            @JsonProperty("orbitData")        List<List<List<List<Cluster>>>>     orbitData,
            @JsonProperty("rcData")           List<List<List<List<Integer>>>>     rcData) {

        this.coordData = coordData;
        this.multiplicityData = multiplicityData;
        this.orbitData = orbitData;
        this.rcData = rcData;
    }

    public List<List<List<Cluster>>> getCoordData() {
        return coordData;
    }

    public List<List<List<Double>>> getMultiplicityData() {
        return multiplicityData;
    }

    public List<List<List<List<Cluster>>>> getOrbitData() {
        return orbitData;
    }

    public List<List<List<List<Integer>>>> getRcData() {
        return rcData;
    }
}



