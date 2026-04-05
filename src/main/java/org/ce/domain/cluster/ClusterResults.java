package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;

/**
 * Unified container for cluster-related result objects.
 */
public final class ClusterResults {

    private ClusterResults() {
        // Utility container
    }

    /**
     * Immutable result object returned by ClusCoordListGenerator.
     */
    public static class ClusCoordListResult {
        private final List<Cluster> clusCoordList;
        private final List<Double> multiplicities;
        private final List<List<Cluster>> orbitList;
        private final List<List<Integer>> rcList;
        private final int tc;
        private final int numPointSubClusFound;

        @JsonCreator
        public ClusCoordListResult(
                @JsonProperty("clusCoordList") List<Cluster> clusCoordList,
                @JsonProperty("multiplicities") List<Double> multiplicities,
                @JsonProperty("orbitList") List<List<Cluster>> orbitList,
                @JsonProperty("rcList") List<List<Integer>> rcList,
                @JsonProperty("tc") int tc,
                @JsonProperty("numPointSubClusFound") int numPointSubClusFound) {
            this.clusCoordList = clusCoordList;
            this.multiplicities = multiplicities;
            this.orbitList = orbitList;
            this.rcList = rcList;
            this.tc = tc;
            this.numPointSubClusFound = numPointSubClusFound;
        }

        public List<Cluster> getClusCoordList() { return clusCoordList; }
        public List<Double> getMultiplicities() { return multiplicities; }
        public List<List<Cluster>> getOrbitList() { return orbitList; }
        public List<List<Integer>> getRcList() { return rcList; }
        public int getTc() { return tc; }
        public int getNumPointSubClusFound() { return numPointSubClusFound; }

        public void printDebug() {
            System.out.println("[ClusCoordListResult]");
            System.out.println("  total cluster types (tc)  : " + tc);
            System.out.println("  point sub-clusters found  : " + numPointSubClusFound);
            System.out.println("  number of entries         : " + clusCoordList.size());
            for (int i = 0; i < clusCoordList.size(); i++) {
                List<ClusterPrimitives.Site> sites = clusCoordList.get(i).getAllSites();
                System.out.printf("  type[%d] : sites=%d, mult=%.4f, rc=%s, orbitSize=%d%n",
                        i, sites.size(), multiplicities.get(i), rcList.get(i), orbitList.get(i).size());
                StringBuilder sb = new StringBuilder("    positions:");
                for (ClusterPrimitives.Site s : sites) {
                    sb.append(" ").append(s.getPosition());
                    if (s.getSymbol() != null) sb.append("(").append(s.getSymbol()).append(")");
                }
                System.out.println(sb);
            }
        }
    }

    /**
     * Immutable result object for classified clusters.
     */
    public static class ClassifiedClusterResult {
        private final List<List<Cluster>> coordList;
        private final List<List<Double>> multiplicityList;
        private final List<List<List<Cluster>>> orbitList;
        private final List<List<List<Integer>>> rcList;

        @JsonCreator
        public ClassifiedClusterResult(
                @JsonProperty("coordList") List<List<Cluster>> coordList,
                @JsonProperty("multiplicityList") List<List<Double>> multiplicityList,
                @JsonProperty("orbitList") List<List<List<Cluster>>> orbitList,
                @JsonProperty("rcList") List<List<List<Integer>>> rcList) {
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

    /**
     * Immutable result object for grouped CF data.
     */
    public static class GroupedCFResult {
        private final List<List<List<Cluster>>> coordData;
        private final List<List<List<Double>>> multiplicityData;
        private final List<List<List<List<Cluster>>>> orbitData;
        private final List<List<List<List<Integer>>>> rcData;

        @JsonCreator
        public GroupedCFResult(
                @JsonProperty("coordData") List<List<List<Cluster>>> coordData,
                @JsonProperty("multiplicityData") List<List<List<Double>>> multiplicityData,
                @JsonProperty("orbitData") List<List<List<List<Cluster>>>> orbitData,
                @JsonProperty("rcData") List<List<List<List<Integer>>>> rcData) {
            this.coordData = coordData;
            this.multiplicityData = multiplicityData;
            this.orbitData = orbitData;
            this.rcData = rcData;
        }

        public List<List<List<Cluster>>> getCoordData() { return coordData; }
        public List<List<List<Double>>> getMultiplicityData() { return multiplicityData; }
        public List<List<List<List<Cluster>>>> getOrbitData() { return orbitData; }
        public List<List<List<List<Integer>>>> getRcData() { return rcData; }
    }
}
