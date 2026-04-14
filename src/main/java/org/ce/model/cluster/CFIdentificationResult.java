package org.ce.model.cluster;

import com.fasterxml.jackson.annotation.JsonAlias;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClassifiedData;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.GroupedCFData;

/**
 * Placeholder for correlation function identification result.
 * Will be fully implemented later.
 */
public class CFIdentificationResult {

    private final ClusCoordListData disCFData;
    private final int tcfdis;
    private final ClusCoordListData phaseCFData;
    private final ClassifiedData ordCFData;
    private final GroupedCFData groupedCFData;
    private final int[][] lcf;
    private final int tcf;
    private final int nxcf;
    private final int ncf;
    private final List<String> uNames;
    private final List<String> eoNames;


    @JsonCreator
    public CFIdentificationResult(
            @JsonProperty("disCFData")     ClusCoordListData     disCFData,
            @JsonProperty("tcfdis")        int                     tcfdis,
            @JsonProperty("phaseCFData")   ClusCoordListData     phaseCFData,
            @JsonProperty("ordCFData")     ClassifiedData        ordCFData,
            @JsonProperty("groupedCFData") GroupedCFData         groupedCFData,
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
    public ClusCoordListData getDisCFData() { return disCFData; }
    @JsonProperty("tcfdis")
    public int getTcfdis() { return tcfdis; }
    @JsonProperty("phaseCFData")
    public ClusCoordListData getPhaseCFData() { return phaseCFData; }
    @JsonProperty("ordCFData")
    public ClassifiedData getOrdCFData() { return ordCFData; }
    @JsonProperty("groupedCFData")
    public GroupedCFData getGroupedCFData() { return groupedCFData; }
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

    public void printSummary(java.util.function.Consumer<String> sink) {
        sink.accept("================================================================================");
        sink.accept("                       CF IDENTIFICATION RESULT");
        sink.accept("================================================================================");

        sink.accept("\nSTAGE 2a: HSP CFs (Disordered Phase)");
        sink.accept("------------------------------------");
        sink.accept(String.format("  - Total CF types (tcfdis): %d", tcfdis));
        disCFData.printSummary(sink);

        sink.accept("\nSTAGE 2b: Ordered Phase CFs (grouped)");
        sink.accept("--------------------------------------");
        sink.accept(String.format("  - Total CF types (tcf): %d", tcf));
        sink.accept(String.format("  - Point CF types (nxcf): %d", nxcf));
        sink.accept(String.format("  - Total CF variables (ncf): %d", ncf));
        groupedCFData.printSummary(sink);
        sink.accept("================================================================================");
    }
}
