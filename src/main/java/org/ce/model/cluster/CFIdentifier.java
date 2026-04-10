package org.ce.model.cluster;

import static org.ce.model.cluster.ClusterResults.*;

import static org.ce.model.cluster.SpaceGroup.SymmetryOperation;

import static org.ce.model.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Orchestrates Stage 2 (Correlation Function Identification) of the CVM pipeline.
 *
 * <h2>What this stage does</h2>
 * <p>CF identification determines which correlation functions (CFs) exist for a
 * given CVM approximation, structure, and number of chemical components.
 * It is component-dependent: adding a third or fourth species creates new
 * CFs (new ways to decorate the same cluster geometry with different symbols).</p>
 *
 * <h2>Stage 2a â€" HSP CFs (n-component basis)</h2>
 * <p>Enumerate all distinct decorated clusters of the HSP using the full
 * n-component basis (e.g. {s1, s2} for ternary â€" two independent site
 * operators for 3 species).  This gives all correlation functions of the
 * disordered (highest-symmetry) reference phase.  Their multiplicities
 * ({@code mcfdis}) are used to normalize the CF basis.</p>
 *
 * <p>Mathematica equivalent:</p>
 * <pre>
 * basisN     = genBasisSymbolList[numComp, siteOpSymbol]   // {s[1], s[2]} for ternary
 * disCFData  = genClusCoordList[disMaxClusCoord, disSymOps, basisN]
 * tcfdis     = disCFData[[5]] - 1
 * </pre>
 *
 * <h2>Stage 2b â€" Phase CFs (n-component basis)</h2>
 * <p>Enumerate all distinct decorated clusters of the ordered phase using the
 * n-component basis.  Classify them into HSP cluster types.  Group the resulting
 * CFs using {@link ClusterBuilders#groupCFData} so that each CF knows which
 * HSP cluster type it belongs to, and within that type, which ordered-phase
 * cluster group it belongs to.</p>
 *
 * <p>Mathematica equivalent:</p>
 * <pre>
 * CFData       = genClusCoordList[maxClusCoord, symOps, basisN]
 * CFCoordList  = ordToDisordCoord[rotateMat, translateMat, CFData[[1]])
 * ordCFData    = transClusCoordList[disCFData, CFData, CFCoordList]
 * cfData       = groupCFData[disClusData, disCFData, ordCFData, basisBin]
 * lcf          = readLength[cfData[[1]]]
 * tcf          = Sum[Sum[lcf[i][j], j], i]
 * nxcf         = Sum[lcf[tcdis][j], j]
 * ncf          = tcf - nxcf
 * </pre>
 *
 * <h2>Dependency on Stage 1</h2>
 * <p>Stage 2 requires the {@link ClusterIdentificationResult} from Stage 1
 * because {@link CFGroupGenerator#groupCFData} needs:
 * <ul>
 *   <li>{@code disClusData} â€" HSP cluster orbits (from Stage 1a)</li>
 *   <li>{@code basisBin} â€" for stripping decorations in the grouping step</li>
 * </ul>
 * </p>
 *
 * @see CFIdentificationResult
 * @see ClusterBuilders
 * @see org.ce.domain.identification.cluster.ClusterIdentifier
 */
public class CFIdentifier {

    private static final Logger LOG = Logger.getLogger(CFIdentifier.class.getName());

    private CFIdentifier() {}

    /**
     * Runs the full two-stage CF identification.
     *
     * @param clusterResult   result from Stage 1 {@link org.ce.domain.identification.cluster.ClusterIdentifier}
     * @param disMaxClusCoord maximal clusters of the HSP
     * @param disSymOps       symmetry operations of the HSP
     * @param maxClusCoord    maximal clusters of the ordered phase
     * @param symOps          symmetry operations of the ordered phase
     * @param rotateMat       3Ã—3 rotation from ordered to HSP frame
     * @param translateMat    translation from ordered to HSP frame
     * @param numComp         number of chemical components (2 = binary, 3 = ternary, ...)
     * @return fully populated {@link CFIdentificationResult}
     */
    public static CFIdentificationResult identify(
            ClusterIdentificationResult clusterResult,
            List<Cluster>               disMaxClusCoord,
            List<SymmetryOperation>     disSymOps,
            List<Cluster>               maxClusCoord,
            List<SymmetryOperation>     symOps,
            double[][]                  rotateMat,
            double[]                    translateMat,
            int                         numComp) {

        int tcdis = clusterResult.getTcdis();

        // ================================================================
        // STAGE 2a: HSP CFs (n-component basis, HSP symmetry)
        // ================================================================
        LOG.fine("CFIdentifier.identify — ENTER: tcdis=" + tcdis + ", numComp=" + numComp);

        // Build n-component basis: {s1, s2, ..., s[numComp-1]}
        // Mathematica: basisN = genBasisSymbolList[numComp, siteOpSymbol]
        List<String> basisN = buildBasis(numComp);

        // Enumerate HSP CFs
        ClusCoordListResult disCFData =
                ClusCoordListGenerator.generate(disMaxClusCoord, disSymOps, basisN);

        int tcfdis = ClusterUtils.countNonEmpty(disCFData);

        // ================================================================
        // STAGE 2b: Phase CFs (n-component basis, phase symmetry)
        // ================================================================
        LOG.fine("CFIdentifier.identify — Stage 2a: tcfdis=" + tcfdis + " HSP CFs found");

        // 2b-1. Enumerate phase CFs
        // Mathematica: CFData = genClusCoordList[maxClusCoord, symOps, basisN]
        ClusCoordListResult phaseCFData =
                ClusCoordListGenerator.generate(maxClusCoord, symOps, basisN);

        int tcPhase = ClusterUtils.countNonEmpty(phaseCFData);

        // 2b-2. Transform phase CF coordinates to HSP frame
        // Mathematica: CFCoordList = ordToDisordCoord[rotateMat, translateMat, CFData[[1]]]
        List<Cluster> transformedCFList =
                OrderedClusterOps.transformToDisordered(
                        rotateMat,
                        translateMat,
                        phaseCFData.getClusCoordList());

        // 2b-3. Classify phase CFs into HSP CF types
        // Mathematica: ordCFData = transClusCoordList[disCFData, CFData, CFCoordList]
        ClusCoordListResult disCFNonEmpty   = ClusterUtils.trimToNonEmpty(disCFData,    tcfdis);
        ClusCoordListResult phaseCFNonEmpty = ClusterUtils.trimToNonEmpty(phaseCFData,  tcPhase);

        ClassifiedClusterResult ordCFData =
                OrderedClusterOps.classify(
                        disCFNonEmpty,
                        phaseCFNonEmpty,
                        transformedCFList.subList(0, tcPhase));

        // 2b-4. Group CFs by HSP cluster type
        // Mathematica: cfData = groupCFData[disClusData, disCFData, ordCFData, basisBin]
        // Note: groupCFData needs the *cluster* (not CF) disordered data for orbit lookup,
        //       and the binary basis for decoration stripping.
        List<String> basisBin = List.of("s1");

        ClusCoordListResult disClusDataForGrouping =
                ClusterUtils.trimToNonEmpty(clusterResult.getDisClusterData(), tcdis);

        GroupedCFResult groupedCFData =
                ClusterBuilders.groupCFData(
                        disClusDataForGrouping,
                        disCFNonEmpty,
                        ordCFData,
                        basisBin);

        // 2b-5. Compute lcf, tcf, nxcf, ncf
        // Mathematica: lcf = readLength[cfCoordList]
        List<List<List<Cluster>>> cfCoordData = groupedCFData.getCoordData();
        int[][] lcf = readLength(cfCoordData, tcdis);

        int tcf = 0;
        for (int t = 0; t < tcdis; t++) {
            for (int j = 0; j < lcf[t].length; j++) {
                tcf += lcf[t][j];
            }
        }

        // nxcf = Î£_j lcf[tcdis-1][j]  (point-cluster group is the last HSP type)
        int nxcf = 0;
        for (int j = 0; j < lcf[tcdis - 1].length; j++) {
            nxcf += lcf[tcdis - 1][j];
        }
        int ncf = tcf - nxcf;

        // Build lcf summary string for the log
        StringBuilder lcfStr = new StringBuilder();
        for (int t = 0; t < tcdis; t++) {
            lcfStr.append("t=").append(t).append("[");
            for (int j = 0; j < lcf[t].length; j++) {
                lcfStr.append(lcf[t][j]).append(j < lcf[t].length - 1 ? "," : "");
            }
            lcfStr.append("] ");
        }
        LOG.fine("CFIdentifier.identify — EXIT: tcf=" + tcf + ", nxcf=" + nxcf
                + ", ncf=" + ncf + ", lcf=" + lcfStr.toString().trim());

        // Generate symbolic names for uList and eoList
        List<String> uNames = new ArrayList<>();
        List<String> eoNames = new ArrayList<>();
        for (int t = 0; t < tcdis; t++) {
            for (int j = 0; j < lcf[t].length; j++) {
                int numCfs = lcf[t][j];
                for (int p = 0; p < numCfs; p++) {
                    // Orthogonal CF names: u[t+1][j+1][p+1]
                    String uName = String.format("u[%d][%d][%d]", t + 1, j + 1, p + 1);
                    uNames.add(uName);

                    // Orthogonal ECI names (exclude point clusters)
                    if (t < tcdis - 1) {
                        String eName = String.format("e[%d][%d][%d]", t + 1, j + 1, p + 1);
                        eoNames.add(eName);
                    }
                }
            }
        }

        // Add empty cluster symbol (constant term) as the last element: tcdis + 1
        // Following Mathematica convention, it is the next available cluster type
        String emptyClusterSymbol = String.format("u[%d][1][1]", tcdis + 1);
        uNames.add(emptyClusterSymbol);

        return new CFIdentificationResult(
                disCFData,
                tcfdis,
                phaseCFData,
                ordCFData,
                groupedCFData,
                lcf,
                tcf,
                nxcf,
                ncf,
                uNames,
                eoNames);
    }

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    /**
     * Builds the n-component basis symbol list.
     * For numComp=3: returns ["s1", "s2"] (numComp-1 independent operators).
     */
    private static List<String> buildBasis(int numComp) {
        List<String> basis = new ArrayList<>();
        for (int i = 1; i <= numComp - 1; i++) {
            basis.add("s" + i);
        }
        return basis;
    }

    /**
     * Reads the lcf array from grouped CF data.
     * Equivalent to Mathematica {@code readLength[cfCoordList]}.
     * {@code lcf[t][j] = cfCoordData[t][j].size()} = number of CFs in group (t,j).
     */
    private static int[][] readLength(
            List<List<List<Cluster>>> cfCoordData, int tcdis) {

        int[][] lcf = new int[tcdis][];
        for (int t = 0; t < tcdis; t++) {
            List<List<Cluster>> groupT = cfCoordData.get(t);
            lcf[t] = new int[groupT.size()];
            for (int j = 0; j < groupT.size(); j++) {
                lcf[t][j] = groupT.get(j).size();
            }
        }
        return lcf;
    }
}






