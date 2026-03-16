package org.ce.domain.cluster;

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

    public CFIdentificationResult() {
        this(null, 0, null, null, null, null, 0, 0, 0);
    }

    public CFIdentificationResult(
            ClusCoordListResult disCFData,
            int tcfdis,
            ClusCoordListResult phaseCFData,
            ClassifiedClusterResult ordCFData,
            GroupedCFResult groupedCFData,
            int[][] lcf,
            int tcf,
            int nxcf,
            int ncf) {

        this.disCFData = disCFData;
        this.tcfdis = tcfdis;
        this.phaseCFData = phaseCFData;
        this.ordCFData = ordCFData;
        this.groupedCFData = groupedCFData;
        this.lcf = lcf;
        this.tcf = tcf;
        this.nxcf = nxcf;
        this.ncf = ncf;
    }

    public ClusCoordListResult getDisCFData() { return disCFData; }
    public int getTcfdis() { return tcfdis; }
    public ClusCoordListResult getPhaseCFData() { return phaseCFData; }
    public ClassifiedClusterResult getOrdCFData() { return ordCFData; }
    public GroupedCFResult getGroupedCFData() { return groupedCFData; }
    public int[][] getLcf() { return lcf; }
    public int getTcf() { return tcf; }
    public int getNxcf() { return nxcf; }
    public int getNcf() { return ncf; }
}
