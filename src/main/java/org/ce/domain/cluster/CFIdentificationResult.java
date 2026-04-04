package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonAlias;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;

/**
 * Placeholder for correlation function identification result.
 * Will be fully implemented later.
 */
public class CFIdentificationResult {

    private final ClusCoordListResult disCFData;
    private final int tcfdis;
    private final ClusCoordListResult phaseCFData;
    private final ClassifiedClusterResult ordCFData;
    private final GroupedCFResult groupedCFData;
    private final int[][] lcf;
    private final int tcf;
    private final int nxcf;
    private final int ncf;
    private final List<String> uNames;
    private final List<String> eoNames;


    @JsonCreator
    public CFIdentificationResult(
            @JsonProperty("disCFData")     ClusCoordListResult     disCFData,
            @JsonProperty("tcfdis")        int                     tcfdis,
            @JsonProperty("phaseCFData")   ClusCoordListResult     phaseCFData,
            @JsonProperty("ordCFData")     ClassifiedClusterResult ordCFData,
            @JsonProperty("groupedCFData") GroupedCFResult         groupedCFData,
            @JsonProperty("lcf")           int[][]                 lcf,
            @JsonProperty("tcf")           int                     tcf,
            @JsonProperty("nxcf")          int                     nxcf,
            @JsonProperty("ncf")           int                     ncf,
            @JsonProperty("uNames") @JsonAlias("unames") List<String> uNames,
            @JsonProperty("eoNames") @JsonAlias("eonames") List<String> eoNames) {

        this.disCFData = disCFData;
        this.tcfdis = tcfdis;
        this.phaseCFData = phaseCFData;
        this.ordCFData = ordCFData;
        this.groupedCFData = groupedCFData;
        this.lcf = lcf;
        this.tcf = tcf;
        this.nxcf = nxcf;
        this.ncf = ncf;
        this.uNames = uNames;
        this.eoNames = eoNames;
    }

    @JsonProperty("disCFData")
    public ClusCoordListResult getDisCFData() { return disCFData; }
    @JsonProperty("tcfdis")
    public int getTcfdis() { return tcfdis; }
    @JsonProperty("phaseCFData")
    public ClusCoordListResult getPhaseCFData() { return phaseCFData; }
    @JsonProperty("ordCFData")
    public ClassifiedClusterResult getOrdCFData() { return ordCFData; }
    @JsonProperty("groupedCFData")
    public GroupedCFResult getGroupedCFData() { return groupedCFData; }
    @JsonProperty("lcf")
    public int[][] getLcf() { return lcf; }
    @JsonProperty("tcf")
    public int getTcf() { return tcf; }
    @JsonProperty("nxcf")
    public int getNxcf() { return nxcf; }
    @JsonProperty("ncf")
    public int getNcf() { return ncf; }

    @JsonProperty("uNames")
    public List<String> getUNames() { return uNames; }
    @JsonProperty("eoNames")
    public List<String> getEONames() { return eoNames; }
}
