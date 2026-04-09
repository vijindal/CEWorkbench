package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterResults.*;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;

/**
 * Immutable result of the two-stage Cluster Identification pipeline.
 *
 * <p>Annotated with @JsonIgnoreProperties(ignoreUnknown = true) to remain
 * robust against diagnostic fields (like nijTable) that may be added to the 
 * Mathematica-exported JSON but aren't strictly required for thermodynamic 
 * calculations.</p>
 *
 * <h2>Conceptual role</h2>
 * <p>Cluster identification determines <em>which clusters exist</em> for a
 * given CVM approximation and structure, independently of how many chemical
 * components decorate the sites.  It answers the question:
 * "How many distinct clusters are there, what are their multiplicities, and
 * what are their Kikuchi-Baker entropy coefficients?"</p>
 *
 * <h2>Two stages</h2>
 * <ol>
 *   <li><b>Stage 1a â€” HSP clusters</b> ({@code disClusterData})<br>
 *       Clusters of the Highest Symmetric Phase (e.g. A2 for a B2 ternary
 *       system), enumerated using A2 symmetry with a binary basis.  These
 *       define the canonical cluster types whose multiplicities and KB
 *       coefficients govern the CVM entropy functional.</li>
 *
 *   <li><b>Stage 1b â€” Phase clusters</b> ({@code ordClusterData})<br>
 *       Clusters of the actual ordered phase (e.g. B2), enumerated using B2
 *       symmetry with a binary basis, then classified back into the HSP
 *       cluster types.  This gives {@code lc[t]} â€” how many ordered-phase
 *       clusters map to HSP cluster type {@code t} â€” and the per-cluster
 *       normalized multiplicities {@code mh}.</li>
 * </ol>
 *
 * <h2>Mathematica equivalents</h2>
 * <pre>
 * disClusterData  â†â†’  disClusData  (from genClusCoordList[disMaxClusCoord, disSymOpList, basisBin])
 * phaseClusterData â†â†’ clusData     (from genClusCoordList[maxClusCoord, symOpList, basisBin])
 * ordClusterData  â†â†’  ordClusData  (from transClusCoordList[disClusData, clusData, ...])
 * nijTable        â†â†’  nijTable     (from getNijTable[...])
 * kbCoefficients  â†â†’  kbdis        (from generateKikuchiBakerCoefficients[...])
 * lc              â†â†’  lc           (= Map[Length, ordClusCoordList])
 * disClusterData  â† â†’  disClusData  (from genClusCoordList[disMaxClusCoord, disSymOpList, basisBin])
 * phaseClusterData â† â†’ clusData     (from genClusCoordList[maxClusCoord, symOpList, basisBin])
 * ordClusterData  â† â†’  ordClusData  (from transClusCoordList[disClusData, clusData, ...])
 * nijTable        â† â†’  nijTable     (from getNijTable[...])
 * kbCoefficients  â† â†’  kbdis        (from generateKikuchiBakerCoefficients[...])
 * lc              â† â†’  lc           (= Map[Length, ordClusCoordList])
 * mh              â† â†’  mh           (= ordClusMList / mhdis)
 * tcdis           â† â†’  tcdis        (= disClusData[[5]] - 1)
 * nxcdis          â† â†’  nxcdis       (always 1 for binary HSP)
 * tc              â† â†’  tc           (total ordered clusters)
 * nxc             â† â†’  nxc          (point clusters in ordered phase)
 * </pre>
 *
 * @see ClusterIdentifier
 */
public class ClusterIdentificationResult {

    // ---- Stage 1a outputs ----

    /** HSP cluster types (A2 binary): tc, multiplicities, orbits, rc. */
    private final ClusCoordListResult disClusterData;

    /**
     * Nij containment table.
     * {@code nijTable[i][j]} = number of times HSP cluster {@code j} appears
     * as a sub-cluster of HSP cluster {@code i}.
     */
    private final int[][]             nijTable;

    /**
     * Kikuchi-Baker entropy coefficients for each HSP cluster type.
     * {@code kbCoefficients[t]} gives the KB weight of cluster type {@code t}
     * in the CVM entropy functional.
     * <p>For the maximal cluster {@code kbCoefficients[0] = 1.0} always.</p>
     */
    private final double[]            kbCoefficients;

    // ---- Stage 1b outputs ----

    /** Phase cluster types (B2 binary), before classification. */
    private final ClusCoordListResult phaseClusterData;

    /**
     * Ordered-phase clusters classified into HSP cluster types.
     * {@code ordClusterData.getCoordList().get(t)} is the list of ordered-phase
     * cluster representatives that map to HSP cluster type {@code t}.
     */
    private final ClassifiedClusterResult ordClusterData;

    /**
     * Number of ordered-phase cluster representatives per HSP cluster type.
     * {@code lc[t] = ordClusterData.getCoordList().get(t).size()}.
     * Equivalent to Mathematica {@code lc = Map[Length, ordClusCoordList]}.
     */
    private final int[]               lc;

    /**
     * Normalized multiplicities of ordered-phase clusters within each HSP type.
     * {@code mh[t][j] = ordClusMList[t][j] / mhdis[t]}.
     * Equivalent to Mathematica {@code mh = ordClusMList / mhdis}.
     */
    private final double[][]          mh;

    // ---- Derived scalars ----

    /** Total number of HSP cluster types (excluding empty cluster). */
    private final int tcdis;

    /** Number of HSP point-cluster types (always 1 for disordered phase). */
    private final int nxcdis;

    /** Total distinct ordered-phase cluster types (excluding empty). */
    private final int tc;

    /** Number of ordered-phase point-cluster types. */
    private final int nxc;

    // -------------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------------

    @JsonCreator
    public ClusterIdentificationResult(
            @JsonProperty("disClusterData")   ClusCoordListResult     disClusterData,
            @JsonProperty("nijTable")         int[][]                 nijTable,
            @JsonProperty("kbCoefficients")   double[]                kbCoefficients,
            @JsonProperty("phaseClusterData") ClusCoordListResult     phaseClusterData,
            @JsonProperty("ordClusterData")   ClassifiedClusterResult ordClusterData,
            @JsonProperty("lc")               int[]                   lc,
            @JsonProperty("mh")               double[][]              mh,
            @JsonProperty("tcdis")            int                     tcdis,
            @JsonProperty("nxcdis")           int                     nxcdis,
            @JsonProperty("tc")               int                     tc,
            @JsonProperty("nxc")              int                     nxc) {

        this.disClusterData   = disClusterData;
        this.nijTable         = nijTable;
        this.kbCoefficients   = kbCoefficients;
        this.phaseClusterData = phaseClusterData;
        this.ordClusterData   = ordClusterData;
        this.lc               = lc;
        this.mh               = mh;
        this.tcdis            = tcdis;
        this.nxcdis           = nxcdis;
        this.tc               = tc;
        this.nxc              = nxc;
    }

    public int[][] getNijTable() { return nijTable; }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    @JsonProperty("disClusterData")
    public ClusCoordListResult getDisClusterData() { return disClusterData; }
    @JsonProperty("kbCoefficients")
    public double[] getKbCoefficients() { return kbCoefficients; }
    @JsonProperty("phaseClusterData")
    public ClusCoordListResult getPhaseClusterData() { return phaseClusterData; }
    @JsonProperty("ordClusterData")
    public ClassifiedClusterResult getOrdClusterData() { return ordClusterData; }
    @JsonProperty("lc")
    public int[] getLc() { return lc; }
    @JsonProperty("mh")
    public double[][] getMh() { return mh; }
    @JsonProperty("tcdis")
    public int getTcdis() { return tcdis; }
    @JsonProperty("nxcdis")
    public int getNxcdis() { return nxcdis; }
    @JsonProperty("tc")
    public int getTc() { return tc; }
    @JsonProperty("nxc")
    public int getNxc() { return nxc; }

    // -------------------------------------------------------------------------
    // Debug print
    // -------------------------------------------------------------------------

    public void printSummary(java.util.function.Consumer<String> sink) {
        sink.accept("================================================================================");
        sink.accept("                       CLUSTER IDENTIFICATION RESULT");
        sink.accept("================================================================================");

        sink.accept("\nSTAGE 1a: HSP Clusters (Disordered Phase)");
        sink.accept("------------------------------------------");
        disClusterData.printSummary(sink);

        sink.accept("\nSTAGE 1b: Ordered Phase Clusters (classified)");
        sink.accept("----------------------------------------------");
        sink.accept(String.format("  - Total ordered types (tc): %d", tc));
        sink.accept(String.format("  - Point clusters (nxc):     %d", nxc));
        ordClusterData.printSummary(sink);
        sink.accept("================================================================================");
    }

    public void printDebug() {
        System.out.println("================================================================================");
        System.out.println("                       CLUSTER IDENTIFICATION RESULT");
        System.out.println("================================================================================");

        // --- Stage 1a ---
        System.out.println("\n------ Stage 1a: HSP Clusters ------");
        System.out.println("  tcdis  = " + tcdis + "  (HSP cluster types, excl. empty)");
        System.out.println("  nxcdis = " + nxcdis + "  (HSP point-cluster types)");
        System.out.printf("  %-6s %-12s %-10s %-20s%n",
                "Type", "mhdis", "rcdis", "KB coeff");
        List<Double>       mhdis = disClusterData.getMultiplicities();
        List<List<Integer>> rcdis = disClusterData.getRcList();
        for (int t = 0; t < tcdis; t++) {
            System.out.printf("  t=%-4d %-12.4f %-10s %-20.8f%n",
                    t, mhdis.get(t), rcdis.get(t), kbCoefficients[t]);
        }

        // --- Stage 1b ---
        System.out.println("\n------ Stage 1b: Ordered Phase Clusters (classified) ------");
        System.out.println("  tc   = " + tc + "  (total ordered cluster types, excl. empty)");
        System.out.println("  nxc  = " + nxc + "  (ordered point-cluster types)");
        System.out.printf("  %-6s %-8s %s%n", "HSP-t", "lc[t]", "mh[t][j]");
        for (int t = 0; t < tcdis; t++) {
            System.out.printf("  t=%-4d lc=%-4d mh=", t, lc[t]);
            for (int j = 0; j < lc[t]; j++) {
                System.out.printf("%.4f ", mh[t][j]);
            }
            System.out.println();
        }

        System.out.println("================================================================================");
    }
}


