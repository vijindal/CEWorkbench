package org.ce.model.cluster;

import java.util.*;
import java.util.function.Consumer;

import static org.ce.model.cluster.ClusterPrimitives.*;
import static org.ce.model.cluster.SpaceGroup.SymmetryOperation;

/**
 * Self-contained Java translation of the Mathematica cluster/CF identification
 * pipeline (Stages 1a, 1b, 2a, 2b).
 *
 * <p>Each Mathematica function is translated 1-to-1 as a Java method.
 * No logic has been added or removed relative to the Mathematica source.</p>
 *
 * <h2>Mathematica functions → Java methods</h2>
 * <pre>
 *   genBasisSymbolList              → {@link #genBasisSymbolList}
 *   sortClusCoord                   → {@link #sortClusCoord}
 *   applySymOpPoint                 → {@link #applySymOpPoint}
 *   applySymOpClus                  → {@link #applySymOpClus}
 *   genOrbit                        → {@link #genOrbit}
 *   isTranslated                    → {@link #isTranslated}
 *   isContained                     → {@link #isContained}
 *   genSubClusCoord (undecorated)   → {@link #genSubClusCoord(Cluster)}
 *   genSubClusCoord (decorated)     → {@link #genSubClusCoord(Cluster, List)}
 *   genClusCoordList                → {@link #genClusCoordList}
 *   getNijTable                     → {@link #getNijTable}
 *   generateKikuchiBakerCoefficients→ {@link #generateKikuchiBakerCoefficients}
 *   ordToDisordCoord                → {@link #ordToDisordCoord}
 *   transClusCoordList              → {@link #transClusCoordList}
 *   groupCFData                     → {@link #groupCFData}
 *   readLength                      → {@link #readLength}
 *   genConfig                       → {@link #genConfig}
 *   calMultiplicity                 → {@link #calMultiplicity}
 *   groupSubClus                    → {@link #groupSubClus}
 * </pre>
 */
public final class ClusterCFIdentificationPipeline {

    private ClusterCFIdentificationPipeline() {}

    private static final double DELTA = 1e-6;

    // =====================================================================
    //  Result container: genClusCoordList output
    //  Mathematica: {clusCoordList, mList, orbitList, rcList, tc, numPoint}
    // =====================================================================

    public static final class ClusCoordListData {
        private final List<Cluster> clusCoordList;
        private final List<Double> multiplicities;
        private final List<List<Cluster>> orbitList;
        private final List<List<Integer>> rcList;
        private final int tc;
        private final int numPointSubClusFound;

        ClusCoordListData(List<Cluster> clusCoordList, List<Double> multiplicities,
                          List<List<Cluster>> orbitList, List<List<Integer>> rcList,
                          int tc, int numPointSubClusFound) {
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

        public void printSummary(Consumer<String> sink) {
            if (sink == null) return;
            sink.accept(String.format("  - Total cluster types (tc): %d (excl. empty)", tc));
            sink.accept(String.format("  - Point sub-clusters found: %d", numPointSubClusFound));
            for (int i = 0; i < tc; i++) {
                sink.accept(String.format("    t=%-3d nc=%-2d mult=%-8.4f rc=%s",
                    i, clusCoordList.get(i).getAllSites().size(), multiplicities.get(i), rcList.get(i)));
            }
        }
    }

    // =====================================================================
    //  Result container: transClusCoordList output
    //  Mathematica: {classifiedCoordList, classifiedMList,
    //                classifiedOrbitList, classifiedRcList}
    // =====================================================================

    public static final class ClassifiedData {
        private final List<List<Cluster>> coordList;
        private final List<List<Double>> multiplicityList;
        private final List<List<List<Cluster>>> orbitList;
        private final List<List<List<Integer>>> rcList;

        ClassifiedData(List<List<Cluster>> coordList, List<List<Double>> multiplicityList,
                       List<List<List<Cluster>>> orbitList, List<List<List<Integer>>> rcList) {
            this.coordList = coordList;
            this.multiplicityList = multiplicityList;
            this.orbitList = orbitList;
            this.rcList = rcList;
        }

        public List<List<Cluster>> getCoordList() { return coordList; }
        public List<List<Double>> getMultiplicityList() { return multiplicityList; }
        public List<List<List<Cluster>>> getOrbitList() { return orbitList; }
        public List<List<List<Integer>>> getRcList() { return rcList; }

        public void printSummary(Consumer<String> sink) {
            if (sink == null) return;
            for (int t = 0; t < coordList.size(); t++) {
                for (int j = 0; j < coordList.get(t).size(); j++) {
                    sink.accept(String.format("    t=%-3d j=%-2d nc=%-2d mult=%-8.4f rc=%s",
                            t, j, coordList.get(t).get(j).getAllSites().size(),
                            multiplicityList.get(t).get(j), rcList.get(t).get(j)));
                }
            }
        }
    }

    // =====================================================================
    //  Result container: groupCFData output
    //  Mathematica: {groupedCFCoordData, groupedCFMData,
    //                groupedCFOrbitData, groupedCFRData}
    // =====================================================================

    public static final class GroupedCFData {
        private final List<List<List<Cluster>>> coordData;
        private final List<List<List<Double>>> multiplicityData;
        private final List<List<List<List<Cluster>>>> orbitData;
        private final List<List<List<List<Integer>>>> rcData;

        GroupedCFData(List<List<List<Cluster>>> coordData,
                      List<List<List<Double>>> multiplicityData,
                      List<List<List<List<Cluster>>>> orbitData,
                      List<List<List<List<Integer>>>> rcData) {
            this.coordData = coordData;
            this.multiplicityData = multiplicityData;
            this.orbitData = orbitData;
            this.rcData = rcData;
        }

        public List<List<List<Cluster>>> getCoordData() { return coordData; }
        public List<List<List<Double>>> getMultiplicityData() { return multiplicityData; }
        public List<List<List<List<Cluster>>>> getOrbitData() { return orbitData; }
        public List<List<List<List<Integer>>>> getRcData() { return rcData; }

        public void printSummary(Consumer<String> sink) {
            if (sink == null) return;
            for (int t = 0; t < coordData.size(); t++) {
                for (int j = 0; j < coordData.get(t).size(); j++) {
                    for (int p = 0; p < coordData.get(t).get(j).size(); p++) {
                        sink.accept(String.format("    t=%-3d j=%-2d p=%-2d mult=%-8.4f rc=%s",
                                t, j, p, multiplicityData.get(t).get(j).get(p), rcData.get(t).get(j).get(p)));
                    }
                }
            }
        }
    }

    // =====================================================================
    //  Full pipeline result
    // =====================================================================

    public static final class PipelineResult {
        // Stage 1a
        private final ClusCoordListData disClusData;
        private final int tcdis;
        private final int nxcdis;
        private final int ncdis;
        private final double[] mhdis;
        private final int[][] nijTable;
        private final double[] kbdis;

        // Stage 1b
        private final ClusCoordListData phaseClusterData;
        private final ClassifiedData ordClusData;
        private final int[] lc;
        private final int tc;
        private final int nxc;
        private final int nc;
        private final double[][] mh;
        private final List<List<List<Integer>>> rc;

        // Stage 2a
        private final ClusCoordListData disCFData;
        private final int tcfdis;

        // Stage 2b
        private final ClusCoordListData phaseCFDataRaw;
        private final ClassifiedData ordCFData;
        private final GroupedCFData cfData;
        private final int[][] lcf;
        private final int tcf;
        private final int nxcf;
        private final int ncf;

        // CF basis indices: cfBasisIndices[col] = array of 1-based σ-powers for that CF column
        private final int[][] cfBasisIndices;

        // Number of components
        private final int numComponents;

        // Expanded fields for consolidated bundle
        private final CMatrixPipeline.CMatrixData matrixData;
        private final List<String> uList;
        private final List<String> vList;
        private final List<String> eOList;
        private final List<String> eList;
        private final double[] equiatomicCVCF;

        PipelineResult(
                ClusCoordListData disClusData, int tcdis, int nxcdis, int ncdis,
                double[] mhdis, int[][] nijTable, double[] kbdis,
                ClusCoordListData phaseClusterData, ClassifiedData ordClusData,
                int[] lc, int tc, int nxc, int nc, double[][] mh, List<List<List<Integer>>> rc,
                ClusCoordListData disCFData, int tcfdis,
                ClusCoordListData phaseCFDataRaw, ClassifiedData ordCFData,
                GroupedCFData cfData, int[][] lcf, int tcf, int nxcf, int ncf,
                int[][] cfBasisIndices, int numComponents,
                CMatrixPipeline.CMatrixData matrixData,
                List<String> uList, List<String> vList,
                List<String> eOList, List<String> eList,
                double[] equiatomicCVCF) {
            this.disClusData = disClusData;
            this.tcdis = tcdis; this.nxcdis = nxcdis; this.ncdis = ncdis;
            this.mhdis = mhdis; this.nijTable = nijTable; this.kbdis = kbdis;
            this.phaseClusterData = phaseClusterData;
            this.ordClusData = ordClusData;
            this.lc = lc; this.tc = tc; this.nxc = nxc; this.nc = nc;
            this.mh = mh; this.rc = rc;
            this.disCFData = disCFData; this.tcfdis = tcfdis;
            this.phaseCFDataRaw = phaseCFDataRaw;
            this.ordCFData = ordCFData;
            this.cfData = cfData;
            this.lcf = lcf; this.tcf = tcf; this.nxcf = nxcf; this.ncf = ncf;
            this.cfBasisIndices = cfBasisIndices;
            this.numComponents = numComponents;
            this.matrixData = matrixData;
            this.uList = uList;
            this.vList = vList;
            this.eOList = eOList;
            this.eList = eList;
            this.equiatomicCVCF = equiatomicCVCF;
        }

        // Stage 1a getters
        public ClusCoordListData getDisClusData() { return disClusData; }
        public int getTcdis() { return tcdis; }
        public int getNxcdis() { return nxcdis; }
        public int getNcdis() { return ncdis; }
        public double[] getMhdis() { return mhdis; }
        public int[][] getNijTable() { return nijTable; }
        public double[] getKbdis() { return kbdis; }

        // Stage 1b getters
        public ClusCoordListData getPhaseClusterData() { return phaseClusterData; }
        public ClassifiedData getOrdClusData() { return ordClusData; }
        public int[] getLc() { return lc; }
        public int getTc() { return tc; }
        public int getNxc() { return nxc; }
        public int getNc() { return nc; }
        public double[][] getMh() { return mh; }
        public List<List<List<Integer>>> getRc() { return rc; }

        // Stage 2a getters
        public ClusCoordListData getDisCFData() { return disCFData; }
        public int getTcfdis() { return tcfdis; }

        // Stage 2b getters
        public ClusCoordListData getPhaseCFDataRaw() { return phaseCFDataRaw; }
        public ClassifiedData getOrdCFData() { return ordCFData; }
        public GroupedCFData getCfData() { return cfData; }
        public int[][] getLcf() { return lcf; }
        public int getTcf() { return tcf; }
        public int getNxcf() { return nxcf; }
        public int getNcf() { return ncf; }
        public int[][] getCfBasisIndices() { return cfBasisIndices; }
        public int getNumComponents() { return numComponents; }

        public CMatrixPipeline.CMatrixData getMatrixData() { return matrixData; }
        public List<String> getUList() { return uList; }
        public List<String> getVList() { return vList; }
        public List<String> getEOList() { return eOList; }
        public List<String> getEList() { return eList; }
        public double[] getEquiatomicCVCF() { return equiatomicCVCF; }

        // =================================================================
        //  Random CF computation — CFs as functions of composition
        //
        //  At the disordered (random) state, each multi-site CF factors as
        //  a product of point CFs:
        //    CF_random[col] = Π_{b ∈ basisIndices[col]} pointCF[b - 1]
        //  where pointCF[k] = ⟨σ^(k+1)⟩ = Σ_i x_i · basis_i^(k+1)
        //  and basis_i comes from the Inden (1992) symmetric integer basis.
        // =================================================================

        /**
         * Computes the full CF vector at the disordered (random) state.
         *
         * <p>The returned vector has length {@code ncf + K} (= {@code tcf + 1}),
         * matching the C-matrix column count:</p>
         * <pre>
         *   [0 .. ncf-1]        non-point CFs (multi-site, optimization variables)
         *   [ncf .. ncf+K-2]    orthogonal point CFs: ⟨σ^k⟩ for k = 1..K-1
         *   [ncf + K - 1]       empty cluster constant = 1.0
         * </pre>
         *
         * @param moleFractions  mole fractions (length K, Σ = 1)
         * @return full random CF vector (length ncf + K)
         */
        public double[] computeRandomCFs(double[] moleFractions) {
            int K = moleFractions.length;
            int pointCfCount = K - 1; // = nxcf for orthogonal basis
            int fullLength = ncf + K; // = tcf + 1

            // Step 1: Compute K-1 point CFs from composition
            //   pointCF[k] = ⟨σ^(k+1)⟩ = Σ_i x_i · basis_i^(k+1)
            double[] basis = ClusterMath.buildBasis(K);
            double[] pointCF = new double[pointCfCount];
            for (int k = 0; k < pointCfCount; k++) {
                for (int i = 0; i < K; i++) {
                    pointCF[k] += moleFractions[i] * Math.pow(basis[i], k + 1);
                }
            }

            // Step 2: Build the full vector
            double[] uFull = new double[fullLength];

            // [0..ncf-1]: non-point CFs = product of point CFs per basis indices
            for (int col = 0; col < ncf; col++) {
                int[] indices = cfBasisIndices[col];
                double val = 1.0;
                for (int b : indices) {
                    val *= pointCF[b - 1]; // basisIndex is 1-based
                }
                uFull[col] = val;
            }

            // [ncf..ncf+K-2]: point CFs, placed using cfBasisIndices
            for (int k = 0; k < pointCfCount; k++) {
                int col = ncf + k;
                int power = cfBasisIndices[col][0]; // single decoration for point CF
                uFull[col] = pointCF[power - 1];
            }

            // [ncf+K-1]: empty cluster constant = 1.0
            uFull[fullLength - 1] = 1.0;

            return uFull;
        }

        /**
         * Computes only the non-point (optimization variable) CF values
         * at the disordered state.
         *
         * @param moleFractions  mole fractions (length K, Σ = 1)
         * @return random non-point CF values (length ncf)
         */
        public double[] computeRandomNonPointCFs(double[] moleFractions) {
            double[] all = computeRandomCFs(moleFractions);
            return java.util.Arrays.copyOf(all, ncf);
        }

        // =====================================================================
        //  ADAPTERS TO LEGACY TYPES
        // =====================================================================

        public ClusterIdentificationResult toClusterIdentificationResult() {
            return new ClusterIdentificationResult(
                    disClusData,
                    nijTable,
                    kbdis,
                    phaseClusterData,
                    ordClusData,
                    lc,
                    mh,
                    tcdis,
                    nxcdis,
                    tc,
                    nxc
            );
        }

        public CFIdentificationResult toCFIdentificationResult() {
            // Generate symbolic names (matches old CFIdentifier logic)
            List<String> uNames = new ArrayList<>();
            List<String> eoNames = new ArrayList<>();
            for (int t = 0; t < tcdis; t++) {
                for (int j = 0; j < lcf[t].length; j++) {
                    int numCfs = lcf[t][j];
                    for (int p = 0; p < numCfs; p++) {
                        uNames.add(String.format("u[%d][%d][%d]", t + 1, j + 1, p + 1));
                        if (t < tcdis - 1) {
                            eoNames.add(String.format("e[%d][%d][%d]", t + 1, j + 1, p + 1));
                        }
                    }
                }
            }
            uNames.add(String.format("u[%d][1][1]", tcdis + 1));

            return new CFIdentificationResult(
                    disCFData,
                    tcfdis,
                    phaseCFDataRaw,
                    ordCFData,
                    cfData,
                    lcf,
                    tcf,
                    nxcf,
                    ncf,
                    uNames,
                    eoNames
            );
        }

        // Compatibility aliases matching AllClusterData
        public ClusterIdentificationResult getDisorderedClusterResult() { return toClusterIdentificationResult(); }
        public ClusterIdentificationResult getOrderedClusterResult() { return toClusterIdentificationResult(); }
        public CFIdentificationResult getDisorderedCFResult() { return toCFIdentificationResult(); }
        public CFIdentificationResult getOrderedCFResult() { return toCFIdentificationResult(); }

        // =========================================================================
        // Print Methods (Moved from AllClusterData)
        // =========================================================================

        /**
         * Prints a detailed summary of all identification results.
         */
        public void printSummary(Consumer<String> sink) {
            if (sink == null) return;
            sink.accept("================================================================================");
            sink.accept("                       CLUSTER IDENTIFICATION RESULT");
            sink.accept("================================================================================");

            sink.accept(String.format("\nIDENTIFICATION PIPELINE: ncf=%d, total-cfs=%d", getNcf(), getTcf()));

            if (matrixData != null) {
                sink.accept(String.format("C-MATRIX PIPELINE: %d types, %d sites",
                        matrixData.getCmat().size(), matrixData.getSiteList().size()));
            }

            printStage3Diagnostics(sink);
            printStage4Diagnostics(sink);

            sink.accept("================================================================================");
        }

        private void printStage3Diagnostics(Consumer<String> sink) {
            if (sink == null || matrixData == null) return;

            sink.accept("\n  [STAGE 3 DIAGNOSTICS: Random state at equiatomic composition]");
            int K = getNumComponents();
            double[] x = new double[K];
            java.util.Arrays.fill(x, 1.0 / K);
            CMatrixPipeline.verifyRandomCVs(x, this, matrixData, sink);
        }

        private void printStage4Diagnostics(Consumer<String> sink) {
            if (sink == null || matrixData == null || equiatomicCVCF == null) return;

            sink.accept("\n  [STAGE 4 DIAGNOSTICS: Random state equiatomic (CVCF Basis)]");
            double[][][] cv = CMatrixPipeline.evaluateCVs(
                    equiatomicCVCF,
                    matrixData.getCmat(),
                    matrixData.getLcv(),
                    getTcdis(),
                    getLc());

            sink.accept("    Resulting CVs (rho_v):");
            for (int t = 0; t < cv.length; t++) {
                for (int j = 0; j < cv[t].length; j++) {
                    sink.accept(String.format("      Type %d, Group %d: %s", t, j, java.util.Arrays.toString(cv[t][j])));
                }
            }
        }
    }



    // =====================================================================
    //  PIPELINE ENTRY POINT
    // =====================================================================

    /**
     * Runs the complete Stages 1a → 1b → 2a → 2b pipeline.
     *
     * @param disMaxClusCoord HSP maximal clusters
     * @param disSymOpList    HSP space-group operations
     * @param maxClusCoord    ordered-phase maximal clusters
     * @param symOpList       ordered-phase space-group operations
     * @param rotateMat       3×3 rotation (ordered → HSP frame)
     * @param translateMat    translation vector (ordered → HSP frame)
     * @param numComp         number of chemical components (≥ 2)
     * @param sink            optional progress sink (may be null)
     * @return fully populated {@link PipelineResult}
     */
    public static PipelineResult run(
            List<Cluster> disMaxClusCoord,
            List<SymmetryOperation> disSymOpList,
            List<Cluster> maxClusCoord,
            List<SymmetryOperation> symOpList,
            double[][] rotateMat,
            double[] translateMat,
            int numComp,
            Consumer<String> sink) {

        return run(disMaxClusCoord, disSymOpList, maxClusCoord, symOpList,
                   rotateMat, translateMat, numComp,
                   null, null, null, null, null, sink);
    }

    /** Inner core run with expansion fields for bundling. */
    private static PipelineResult run(
            List<Cluster> disMaxClusCoord,
            List<SymmetryOperation> disSymOpList,
            List<Cluster> maxClusCoord,
            List<SymmetryOperation> symOpList,
            double[][] rotateMat,
            double[] translateMat,
            int numComp,
            CMatrixPipeline.CMatrixData matrixData,
            List<String> uList, List<String> vList,
            List<String> eOList, List<String> eList,
            Consumer<String> sink) {

        emit(sink, "=== ClusterCFIdentificationPipeline: START ===");

        // ---- Stage 1a: HSP clusters (binary basis) ----
        int numCompBin = 2;
        List<String> basisSymbolListBin = genBasisSymbolList(numCompBin);
        emit(sink, "[1a] basisSymbolListBin: " + basisSymbolListBin);

        // Mathematica: disClusData = genClusCoordList[disMaxClusCoord, disSymOpList, basisBin]
        // genClusCoordList now excludes empty clusters → tc is always Tier (ii)
        ClusCoordListData disClusData =
                genClusCoordList(disMaxClusCoord, disSymOpList, basisSymbolListBin, sink);

        int tcdis = disClusData.getTc();
        int nxcdis = 1;
        int ncdis = tcdis - nxcdis;

        List<Cluster> disClusList = disClusData.getClusCoordList();
        List<Double> mhdisList = disClusData.getMultiplicities();
        List<List<Cluster>> disOrbitList = disClusData.getOrbitList();

        double[] mhdis = new double[tcdis];
        for (int i = 0; i < tcdis; i++) mhdis[i] = mhdisList.get(i);

        // Mathematica: nijTable = getNijTable[disClusCoordList, mhdis, disClusOrbitList]
        int[][] nijTable = getNijTable(disClusList, mhdisList, disOrbitList);

        // Mathematica: kbdis = generateKikuchiBakerCoefficients[mhdis, nijTable]
        double[] kbdis = generateKikuchiBakerCoefficients(mhdis, nijTable);

        emit(sink, "[1a] tcdis=" + tcdis + ", nxcdis=" + nxcdis + ", ncdis=" + ncdis);
        emit(sink, "[1a] mhdis=" + Arrays.toString(mhdis));
        emit(sink, "[1a] kbdis=" + Arrays.toString(kbdis));

        // ---- Stage 1b: Phase clusters ----
        // Mathematica: clusData = genClusCoordList[maxClusCoord, symOpList, basisBin]
        ClusCoordListData phaseClusterData =
                genClusCoordList(maxClusCoord, symOpList, basisSymbolListBin, sink);
        int tc = phaseClusterData.getTc();

        // Mathematica: clusCoordList = ordToDisordCoord[rotateMat, translateMat, clusData[[1]]]
        List<Cluster> transformedClusList =
                ordToDisordCoord(rotateMat, translateMat, phaseClusterData.getClusCoordList());

        // Mathematica: ordClusData = transClusCoordList[disClusData, clusData, clusCoordList]
        ClassifiedData ordClusData =
                transClusCoordList(disClusData, phaseClusterData, transformedClusList);

        // lc = Map[Length, ordClusCoordList]
        int[] lc = new int[tcdis];
        for (int t = 0; t < tcdis; t++) {
            lc[t] = ordClusData.getCoordList().get(t).size();
        }
        int nxc = lc[tcdis - 1];
        int nc = tc - nxc;

        // rc = ordClusData[[4]]
        List<List<List<Integer>>> rc = ordClusData.getRcList();

        // mh = ordClusMList / mhdis
        double[][] mh = new double[tcdis][];
        for (int t = 0; t < tcdis; t++) {
            mh[t] = new double[lc[t]];
            for (int j = 0; j < lc[t]; j++) {
                mh[t][j] = ordClusData.getMultiplicityList().get(t).get(j) / mhdis[t];
            }
        }

        emit(sink, "[1b] tc=" + tc + ", nxc=" + nxc + ", nc=" + nc);
        emit(sink, "[1b] lc=" + Arrays.toString(lc));

        // ---- Stage 2a: HSP CFs (n-component basis) ----
        List<String> basisSymbolList = genBasisSymbolList(numComp);
        emit(sink, "[2a] basisSymbolList: " + basisSymbolList);

        // Mathematica: disCFData = genClusCoordList[disMaxClusCoord, disSymOpList, basisN]
        ClusCoordListData disCFData =
                genClusCoordList(disMaxClusCoord, disSymOpList, basisSymbolList, sink);
        int tcfdis = disCFData.getTc();
        emit(sink, "[2a] tcfdis=" + tcfdis);

        // ---- Stage 2b: Phase CFs ----
        // Mathematica: CFData = genClusCoordList[maxClusCoord, symOpList, basisN]
        ClusCoordListData phaseCFDataRaw =
                genClusCoordList(maxClusCoord, symOpList, basisSymbolList, sink);
        int tcfRaw = phaseCFDataRaw.getTc();

        // Mathematica: CFCoordList = ordToDisordCoord[rotateMat, translateMat, CFData[[1]]]
        List<Cluster> transformedCFList =
                ordToDisordCoord(rotateMat, translateMat, phaseCFDataRaw.getClusCoordList());

        // Mathematica: ordCFData = transClusCoordList[disCFData, CFData, CFCoordList]
        ClassifiedData ordCFData =
                transClusCoordList(disCFData, phaseCFDataRaw, transformedCFList);

        // Mathematica: cfData = groupCFData[disClusData, disCFData, ordCFData, ..., basisBin]
        GroupedCFData cfData =
                groupCFData(disClusData, disCFData, ordCFData, basisSymbolListBin);

        // lcf = readLength[cfData[[1]]]
        int[][] lcf = readLength(cfData.getCoordData());

        // tcf = Sum[Sum[lcf[i][j]]]
        int tcf = 0;
        for (int i = 0; i < lcf.length; i++)
            for (int j = 0; j < lcf[i].length; j++)
                tcf += lcf[i][j];

        // nxcf = Sum[lcf[tcdis-1][j]]  (last HSP type = point cluster type)
        int nxcf = 0;
        if (tcdis - 1 < lcf.length)
            for (int j = 0; j < lcf[tcdis - 1].length; j++)
                nxcf += lcf[tcdis - 1][j];

        int ncf = tcf - nxcf;

        // ---- Derive CF basis indices from grouped CF data ----
        // For each CF column, extract the 1-based σ-power decorations
        // from the representative cluster's site symbols ("s1"→1, "s2"→2, ...).
        int[][] cfBasisIndices = deriveCfBasisIndices(cfData, lcf);

        emit(sink, "[2b] tcf=" + tcf + ", nxcf=" + nxcf + ", ncf=" + ncf);
        emit(sink, "=== ClusterCFIdentificationPipeline: COMPLETE ===");

        return new PipelineResult(
                disClusData, tcdis, nxcdis, ncdis, mhdis, nijTable, kbdis,
                phaseClusterData, ordClusData, lc, tc, nxc, nc, mh, rc,
                disCFData, tcfdis,
                phaseCFDataRaw, ordCFData, cfData, lcf, tcf, nxcf, ncf,
                cfBasisIndices, numComp,
                matrixData,
                uList, vList,
                eOList, eList,
                null // equiatomicCVCF filled by runFullWorkflow
        );
    }

    // =====================================================================
    //  genBasisSymbolList[numComp, siteOpSymbol]
    //
    //  Return[Table[siteOpSymbol[i], {i, 1, numComp-1}]]
    // =====================================================================

    public static List<String> genBasisSymbolList(int numComp) {
        List<String> result = new ArrayList<>();
        for (int i = 1; i <= numComp - 1; i++) {
            result.add("s" + i);
        }
        return result;
    }

    // =====================================================================
    //  sortClusCoord[clusCoord]
    //
    //  Insertion sort on sites by coordinate (x, then y, then z).
    //  Mathematica accesses site[[1]][[1..3]] for the 3D coordinate.
    // =====================================================================

    public static List<Site> sortClusCoord(List<Site> sites) {
        List<Site> sorted = new ArrayList<>(sites);
        for (int i = 1; i < sorted.size(); i++) {
            Site x = sorted.get(i);
            int j = i - 1;
            while (j >= 0 && compareSites(sorted.get(j), x) > 0) {
                sorted.set(j + 1, sorted.get(j));
                j--;
            }
            sorted.set(j + 1, x);
        }
        return sorted;
    }

    /**
     * Site comparator matching Mathematica sortClusCoord ordering:
     * primary: x-coordinate, secondary: y-coordinate, tertiary: z-coordinate.
     */
    static int compareSites(Site a, Site b) {
        Position pa = a.getPosition();
        Position pb = b.getPosition();
        int cx = Double.compare(pa.getX(), pb.getX());
        if (cx != 0) return cx;
        int cy = Double.compare(pa.getY(), pb.getY());
        if (cy != 0) return cy;
        return Double.compare(pa.getZ(), pb.getZ());
    }

    // =====================================================================
    //  applySymOpPoint[symOp, coord]
    //
    //  matrixPart = Drop[symOP, None, {4}]
    //  transPart  = symOP[[All, 4]]
    //  transCoord = matrixPart . coord[[1]] + transPart
    //  Return[{transCoord, coord[[2]]}]
    // =====================================================================

    static Site applySymOpPoint(SymmetryOperation symOp, Site site) {
        Position r = site.getPosition();
        double[][] rot = symOp.getRotation();
        double[] trans = symOp.getTranslation();
        double x = rot[0][0] * r.getX() + rot[0][1] * r.getY() + rot[0][2] * r.getZ() + trans[0];
        double y = rot[1][0] * r.getX() + rot[1][1] * r.getY() + rot[1][2] * r.getZ() + trans[1];
        double z = rot[2][0] * r.getX() + rot[2][1] * r.getY() + rot[2][2] * r.getZ() + trans[2];
        return new Site(new Position(x, y, z), site.getSymbol());
    }

    // =====================================================================
    //  applySymOpClus[symOp, clusCoord]
    //
    //  For each sublattice: transform each site, then sortClusCoord.
    // =====================================================================

    static Cluster applySymOpClus(SymmetryOperation symOp, Cluster cluster) {
        List<Sublattice> newSublattices = new ArrayList<>();
        for (Sublattice sub : cluster.getSublattices()) {
            List<Site> newSites = new ArrayList<>();
            for (Site site : sub.getSites()) {
                newSites.add(applySymOpPoint(symOp, site));
            }
            newSites = sortClusCoord(newSites);
            newSublattices.add(new Sublattice(newSites));
        }
        return new Cluster(newSublattices);
    }

    // =====================================================================
    //  genOrbit[clusCoord, spaceGroup]
    //
    //  For each symmetry op, transform cluster.
    //  If not already in orbit, append.
    // =====================================================================

    public static List<Cluster> genOrbit(Cluster cluster, List<SymmetryOperation> spaceGroup) {
        List<Cluster> orbit = new ArrayList<>();
        for (SymmetryOperation op : spaceGroup) {
            Cluster transformed = applySymOpClus(op, cluster);
            if (!isContained(orbit, transformed)) {
                orbit.add(transformed);
            }
        }
        return orbit;
    }

    // =====================================================================
    //  isTranslated[clusCoord, clusCoord2]
    //
    //  Computes site-wise differences {coordDiff, symbolDiff}.
    //  DeleteDuplicates. Checks length, symbol match, integer shift.
    // =====================================================================

    public static boolean isTranslated(Cluster c1, Cluster c2) {
        if (c1.getSublattices().size() != c2.getSublattices().size()) return false;

        // Collect all site-pair differences as {dx, dy, dz, symbolMatch}
        List<double[]> diffs = new ArrayList<>();

        for (int i = 0; i < c1.getSublattices().size(); i++) {
            List<Site> s1 = c1.getSublattices().get(i).getSites();
            List<Site> s2 = c2.getSublattices().get(i).getSites();
            if (s1.size() != s2.size()) return false;

            for (int j = 0; j < s1.size(); j++) {
                Position p1 = s1.get(j).getPosition();
                Position p2 = s2.get(j).getPosition();
                double dx = p2.getX() - p1.getX();
                double dy = p2.getY() - p1.getY();
                double dz = p2.getZ() - p1.getZ();
                // Mathematica: symbolDiff = sym2 - sym1; equals 0 if symbols match
                double symMatch = s1.get(j).getSymbol().equals(s2.get(j).getSymbol()) ? 1.0 : 0.0;
                diffs.add(new double[]{dx, dy, dz, symMatch});
            }
        }

        // DeleteDuplicates
        List<double[]> unique = new ArrayList<>();
        for (double[] d : diffs) {
            boolean isDup = false;
            for (double[] u : unique) {
                if (Math.abs(d[0] - u[0]) < DELTA && Math.abs(d[1] - u[1]) < DELTA
                        && Math.abs(d[2] - u[2]) < DELTA && d[3] == u[3]) {
                    isDup = true;
                    break;
                }
            }
            if (!isDup) unique.add(d);
        }

        // Mathematica: If[Length[clusCoordDiff] > 1, False, ...]
        if (unique.size() > 1) return false;

        // Mathematica: empty cluster → True
        if (unique.isEmpty()) return true;

        double[] diff = unique.get(0);

        // Mathematica: If[clusCoordDiff[[1]][[2]] != 0, False, ...]
        // symbolDiff != 0 means symbols don't match → symMatch == 0.0
        if (diff[3] == 0.0) return false;

        // Mathematica: FractionalPart[Abs[clusCoordDiff[[1]][[1]][[j]]]] >= delta
        for (int j = 0; j < 3; j++) {
            double absVal = Math.abs(diff[j]);
            double fracPart = absVal - Math.floor(absVal);
            if (fracPart >= DELTA && fracPart <= (1.0 - DELTA)) {
                return false;
            }
        }
        return true;
    }

    // =====================================================================
    //  isContained[orbit, clusCoord]
    //
    //  Checks if cluster is in orbit (size check + isTranslated).
    // =====================================================================

    public static boolean isContained(List<Cluster> orbit, Cluster cluster) {
        List<Integer> candidateSizes = sublatticeSizes(cluster);
        for (Cluster existing : orbit) {
            if (sublatticeSizes(existing).equals(candidateSizes)) {
                if (isTranslated(existing, cluster)) {
                    return true;
                }
            }
        }
        return false;
    }

    // =====================================================================
    //  genSubClusCoord[clusCoord]  (undecorated — Subsets version)
    //
    //  Flatten sites, generate all subsets, redistribute to sublattices.
    // =====================================================================

    public static List<Cluster> genSubClusCoord(Cluster cluster) {
        List<Sublattice> originalSubs = cluster.getSublattices();
        int numSubLattice = originalSubs.size();

        // Flatten and sort: disClusCoord = sortClusCoord[Flatten[clusCoord, 1]]
        List<Site> allSites = new ArrayList<>();
        for (Sublattice sub : originalSubs) {
            allSites.addAll(sub.getSites());
        }
        List<Site> sorted = sortClusCoord(allSites);

        // Subsets[disClusCoord]
        List<List<Site>> subsets = generateAllSubsets(sorted);

        // Redistribute each subset back to original sublattices
        List<Cluster> result = new ArrayList<>();
        for (List<Site> subset : subsets) {
            List<Sublattice> subClus = new ArrayList<>();
            for (int k = 0; k < numSubLattice; k++) {
                subClus.add(new Sublattice(new ArrayList<>()));
            }
            for (Site site : subset) {
                for (int k = 0; k < numSubLattice; k++) {
                    if (containsSite(originalSubs.get(k).getSites(), site)) {
                        subClus.get(k).getSites().add(site);
                        break;
                    }
                }
            }
            result.add(new Cluster(subClus));
        }
        return result;
    }

    // =====================================================================
    //  genSubClusCoord[clusCoord, basisSymbolList]  (decorated — Tuples)
    //
    //  For each flattened site position, create options:
    //    empty (absent) + one option per basis symbol.
    //  Cartesian product (Tuples), then redistribute to sublattices.
    // =====================================================================

    public static List<Cluster> genSubClusCoord(Cluster cluster, List<String> basisSymbolList) {
        List<Sublattice> originalSubs = cluster.getSublattices();
        int numSubLattice = originalSubs.size();

        // Flatten and sort positions
        List<Site> allSites = new ArrayList<>();
        for (Sublattice sub : originalSubs) {
            allSites.addAll(sub.getSites());
        }
        allSites = sortClusCoord(allSites);

        // Build options per site: {null(absent), (pos,s1), (pos,s2), ...}
        List<List<Site>> disClus = new ArrayList<>();
        for (Site site : allSites) {
            List<Site> options = new ArrayList<>();
            options.add(null); // empty option: {{}}
            for (String symbol : basisSymbolList) {
                options.add(new Site(site.getPosition(), symbol));
            }
            disClus.add(options);
        }

        // Tuples[disClus]
        List<List<Site>> tuples = cartesianProduct(disClus);

        // Redistribute to sublattices
        List<Cluster> result = new ArrayList<>();
        for (List<Site> tuple : tuples) {
            List<Sublattice> subClus = new ArrayList<>();
            for (int k = 0; k < numSubLattice; k++) {
                subClus.add(new Sublattice(new ArrayList<>()));
            }
            for (Site decoratedSite : tuple) {
                if (decoratedSite == null) continue; // absent option
                for (int k = 0; k < numSubLattice; k++) {
                    if (containsPosition(originalSubs.get(k).getSites(),
                            decoratedSite.getPosition())) {
                        subClus.get(k).getSites().add(decoratedSite);
                        break;
                    }
                }
            }
            result.add(new Cluster(subClus));
        }
        return result;
    }

    // =====================================================================
    //  genClusCoordList[maxClusCoord, spaceGroup, basisSymbolList]
    //
    //  Main cluster enumeration.
    //  1. Generate decorated subclusters for each maximal cluster
    //  2. Sort by descending size, iterate smallest→largest
    //  3. For each new subcluster, record orbit, multiplicity, rc
    //  4. Normalize multiplicities by point-cluster orbit count
    //  5. Final sort by descending cluster size
    // =====================================================================

    public static ClusCoordListData genClusCoordList(
            List<Cluster> maxClusCoord,
            List<SymmetryOperation> spaceGroup,
            List<String> basisSymbolList,
            Consumer<String> sink) {

        emit(sink, "  genClusCoordList called: maxClus=" + maxClusCoord.size()
                + ", symOps=" + spaceGroup.size() + ", basis=" + basisSymbolList);

        List<Cluster> clusCoordList = new ArrayList<>();
        List<List<Cluster>> subClusOrbitList = new ArrayList<>();
        List<Integer> subClusMList = new ArrayList<>();
        List<List<Integer>> rc = new ArrayList<>();

        int numMaxClus = maxClusCoord.size();

        // STEP 1: Iterate maximal clusters
        for (int k = 0; k < numMaxClus; k++) {
            // genSubClusCoord[maxClusCoord[[k]], basisSymbolList]
            List<Cluster> subClusCoord = genSubClusCoord(maxClusCoord.get(k), basisSymbolList);

            // Sort descending by total site count
            subClusCoord.sort((a, b) ->
                    Integer.compare(b.getAllSites().size(), a.getAllSites().size()));

            int numSubClus = subClusCoord.size();

            // Iterate from smallest to largest (Mathematica: For[i=numSubClus, i>=1, i--])
            for (int i = numSubClus - 1; i >= 0; i--) {
                Cluster candidate = subClusCoord.get(i);

                // Skip empty clusters — ensures tc is always Tier (ii)
                if (candidate.getAllSites().isEmpty()) continue;

                boolean foundNewCluster = true;

                for (List<Cluster> existingOrbit : subClusOrbitList) {
                    if (isContained(existingOrbit, candidate)) {
                        foundNewCluster = false;
                        break;
                    }
                }

                if (foundNewCluster) {
                    clusCoordList.add(candidate);
                    List<Cluster> orbit = genOrbit(candidate, spaceGroup);
                    subClusOrbitList.add(orbit);
                    subClusMList.add(orbit.size());

                    // rc = Map[Length, subClusCoord[[i]]]
                    List<Integer> rcEntry = new ArrayList<>();
                    for (Sublattice sub : candidate.getSublattices()) {
                        rcEntry.add(sub.getSites().size());
                    }
                    rc.add(rcEntry);
                }
            }
        }

        int tcCount = clusCoordList.size();

        // STEP 2: Count point clusters for normalization
        int numPointSubClusFound = 0;
        double pointM = 0;
        List<Position> pointPositions = new ArrayList<>();

        for (int i = 0; i < tcCount; i++) {
            List<Site> flatList = clusCoordList.get(i).getAllSites();
            if (flatList.size() == 1) {
                // Mathematica: flatList[[1]][[1]] = coordinate (ignoring decoration)
                Position pos = flatList.get(0).getPosition();
                boolean alreadyCounted = false;
                for (Position existing : pointPositions) {
                    if (positionsMatch(existing, pos)) {
                        alreadyCounted = true;
                        break;
                    }
                }
                if (!alreadyCounted) {
                    pointPositions.add(pos);
                    pointM += subClusOrbitList.get(i).size();
                }
                numPointSubClusFound++;
            }
        }

        // STEP 3: Normalize multiplicities
        List<Double> normalizedM = new ArrayList<>();
        for (int m : subClusMList) {
            normalizedM.add(m / pointM);
        }

        // STEP 4: Sort by descending total site count
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < tcCount; i++) indices.add(i);
        indices.sort((i1, i2) -> Integer.compare(
                clusCoordList.get(i2).getAllSites().size(),
                clusCoordList.get(i1).getAllSites().size()));

        List<Cluster> finalClusList = new ArrayList<>();
        List<Double> finalM = new ArrayList<>();
        List<List<Cluster>> finalOrbit = new ArrayList<>();
        List<List<Integer>> finalRc = new ArrayList<>();

        for (int idx : indices) {
            finalClusList.add(clusCoordList.get(idx));
            finalM.add(normalizedM.get(idx));
            finalOrbit.add(subClusOrbitList.get(idx));
            finalRc.add(rc.get(idx));
        }

        emit(sink, "  genClusCoordList ended: tc=" + tcCount);

        return new ClusCoordListData(finalClusList, finalM, finalOrbit, finalRc,
                tcCount, numPointSubClusFound);
    }

    // =====================================================================
    //  getNijTable[clusCoordList, clusMList, clusOrbitList]
    //
    //  For each cluster i, generate its subclusters.
    //  For each subcluster k, check containment in orbit j (j >= i).
    // =====================================================================

    public static int[][] getNijTable(
            List<Cluster> clusCoordList,
            List<Double> clusMList,
            List<List<Cluster>> clusOrbitList) {

        int numClus = clusCoordList.size();
        int[][] nijTable = new int[numClus][numClus];

        for (int i = 0; i < numClus; i++) {
            // genSubClusCoord[clusCoordList[[i]]] (undecorated)
            List<Cluster> subClusCoord = genSubClusCoord(clusCoordList.get(i));
            int numSubClus = subClusCoord.size();

            // Mathematica: For[k=numSubClus, k>=1, k--]
            for (int k = numSubClus - 1; k >= 0; k--) {
                Cluster subCluster = subClusCoord.get(k);
                List<Integer> subSizes = sublatticeSizes(subCluster);

                for (int j = 0; j < numClus; j++) {
                    // Mathematica: If[j >= i, ...]
                    if (j >= i) {
                        // Mathematica: Map[Length, subClusCoord[[k]]] == Map[Length, clusCoordList[[j]]]
                        if (subSizes.equals(sublatticeSizes(clusCoordList.get(j)))) {
                            if (isContained(clusOrbitList.get(j), subCluster)) {
                                nijTable[i][j]++;
                            }
                        }
                    }
                }
            }
        }

        return nijTable;
    }

    // =====================================================================
    //  generateKikuchiBakerCoefficients[mList, nijTable]
    //
    //  kb[j] = (m[j] - sum_{i<j} m[i]*nij[i][j]*kb[i]) / m[j]
    // =====================================================================

    public static double[] generateKikuchiBakerCoefficients(double[] mList, int[][] nijTable) {
        int n = mList.length;
        double[] kb = new double[n];
        for (int j = 0; j < n; j++) {
            double tempSum = 0.0;
            for (int i = 0; i < j; i++) {
                tempSum += mList[i] * nijTable[i][j] * kb[i];
            }
            kb[j] = (mList[j] - tempSum) / mList[j];
        }
        return kb;
    }

    // =====================================================================
    //  ordToDisordCoord[rotateMat, translateMat, clusCoordList]
    //
    //  Transforms each site: newCoord = R · coord + t, preserves symbol.
    // =====================================================================

    public static List<Cluster> ordToDisordCoord(
            double[][] rotateMat, double[] translateMat, List<Cluster> clusters) {

        List<Cluster> result = new ArrayList<>();
        for (Cluster cluster : clusters) {
            List<Sublattice> newSublattices = new ArrayList<>();
            for (Sublattice sub : cluster.getSublattices()) {
                List<Site> newSites = new ArrayList<>();
                for (Site site : sub.getSites()) {
                    Position r = site.getPosition();
                    double x = rotateMat[0][0] * r.getX() + rotateMat[0][1] * r.getY()
                             + rotateMat[0][2] * r.getZ() + translateMat[0];
                    double y = rotateMat[1][0] * r.getX() + rotateMat[1][1] * r.getY()
                             + rotateMat[1][2] * r.getZ() + translateMat[1];
                    double z = rotateMat[2][0] * r.getX() + rotateMat[2][1] * r.getY()
                             + rotateMat[2][2] * r.getZ() + translateMat[2];
                    newSites.add(new Site(new Position(x, y, z), site.getSymbol()));
                }
                newSublattices.add(new Sublattice(newSites));
            }
            result.add(new Cluster(newSublattices));
        }
        return result;
    }

    // =====================================================================
    //  transClusCoordList[disClusData, clusData, transformedCoords]
    //
    //  Classifies ordered-phase clusters into HSP types by flattening
    //  transformed coordinates and checking orbit containment.
    // =====================================================================

    public static ClassifiedData transClusCoordList(
            ClusCoordListData disClusData,
            ClusCoordListData clusData,
            List<Cluster> transformedCoords) {

        List<Cluster> disClusCoordList = disClusData.getClusCoordList();
        List<List<Cluster>> disClusOrbitList = disClusData.getOrbitList();
        int tcdis = disClusCoordList.size(); // includes empty cluster

        List<Cluster> clusCoordList1 = clusData.getClusCoordList();
        List<Double> clusMList = clusData.getMultiplicities();
        List<List<Cluster>> clusOrbitList = clusData.getOrbitList();
        List<List<Integer>> clusRcList = clusData.getRcList();

        int tc = transformedCoords.size();

        // Flatten each transformed cluster into single-sublattice form
        // Mathematica: flattenClusCoord = sortClusCoord[Flatten[clusCoordList[[i]], 1]]
        //              flattenClusCoordList = Append[..., {flattenClusCoord}]
        List<Cluster> flattenClusCoordList = new ArrayList<>();
        for (int i = 0; i < tc; i++) {
            List<Site> flat = new ArrayList<>(transformedCoords.get(i).getAllSites());
            flat = sortClusCoord(flat);
            // Wrap in single sublattice: {flattenClusCoord}
            List<Sublattice> singleSub = new ArrayList<>();
            singleSub.add(new Sublattice(flat));
            flattenClusCoordList.add(new Cluster(singleSub));
        }

        // Classification
        List<List<Cluster>> classifiedCoord = new ArrayList<>();
        List<List<Double>> classifiedM = new ArrayList<>();
        List<List<List<Cluster>>> classifiedOrbit = new ArrayList<>();
        List<List<List<Integer>>> classifiedRc = new ArrayList<>();

        for (int j = 0; j < tcdis; j++) {
            classifiedCoord.add(new ArrayList<>());
            classifiedM.add(new ArrayList<>());
            classifiedOrbit.add(new ArrayList<>());
            classifiedRc.add(new ArrayList<>());

            for (int i = 0; i < tc; i++) {
                if (isContained(disClusOrbitList.get(j), flattenClusCoordList.get(i))) {
                    classifiedCoord.get(j).add(clusCoordList1.get(i));
                    classifiedM.get(j).add(clusMList.get(i));
                    classifiedOrbit.get(j).add(clusOrbitList.get(i));
                    classifiedRc.get(j).add(clusRcList.get(i));
                }
            }
        }

        return new ClassifiedData(classifiedCoord, classifiedM, classifiedOrbit, classifiedRc);
    }

    // =====================================================================
    //  groupCFData[disClusData, disCFData, ordCFData, basisSymbolListBin]
    //
    //  Groups CFs into HSP clusters by stripping decorations from
    //  disordered CFs and checking orbit containment.
    //  Note: Mathematica signature includes rotateMat/translateMat but
    //  they are NOT used in the function body (vestigial parameters).
    // =====================================================================

    public static GroupedCFData groupCFData(
            ClusCoordListData disClusData,
            ClusCoordListData disCFData,
            ClassifiedData ordCFData,
            List<String> basisSymbolListBin) {

        List<Cluster> disClusCoordData = disClusData.getClusCoordList();
        List<List<Cluster>> disClusOrbitList = disClusData.getOrbitList();
        List<Cluster> disCFCoordData = disCFData.getClusCoordList();

        // ----- Remove decorations from disordered CFs -----
        // Replace each site's symbol with basisSymbolListBin[[1]] (= "s1")
        List<Cluster> transDisCFCoordData = new ArrayList<>();
        for (Cluster cf : disCFCoordData) {
            List<Sublattice> newSubs = new ArrayList<>();
            for (Sublattice sub : cf.getSublattices()) {
                List<Site> newSites = new ArrayList<>();
                for (Site s : sub.getSites()) {
                    newSites.add(new Site(s.getPosition(), basisSymbolListBin.get(0)));
                }
                newSubs.add(new Sublattice(newSites));
            }
            transDisCFCoordData.add(new Cluster(newSubs));
        }

        // ----- Group by orbit containment -----
        List<List<List<Cluster>>> groupedCoord = new ArrayList<>();
        List<List<List<Double>>> groupedM = new ArrayList<>();
        List<List<List<List<Cluster>>>> groupedOrbit = new ArrayList<>();
        List<List<List<List<Integer>>>> groupedRc = new ArrayList<>();

        for (int i = 0; i < disClusCoordData.size(); i++) {
            groupedCoord.add(new ArrayList<>());
            groupedM.add(new ArrayList<>());
            groupedOrbit.add(new ArrayList<>());
            groupedRc.add(new ArrayList<>());

            for (int j = 0; j < transDisCFCoordData.size(); j++) {
                if (isContained(disClusOrbitList.get(i), transDisCFCoordData.get(j))) {
                    groupedCoord.get(i).add(ordCFData.getCoordList().get(j));
                    groupedM.get(i).add(ordCFData.getMultiplicityList().get(j));
                    groupedOrbit.get(i).add(ordCFData.getOrbitList().get(j));
                    groupedRc.get(i).add(ordCFData.getRcList().get(j));
                }
            }
        }

        return new GroupedCFData(groupedCoord, groupedM, groupedOrbit, groupedRc);
    }

    // =====================================================================
    //  readLength[array]
    //
    //  Returns nested length: lcf[i][j] = array[i][j].size()
    // =====================================================================

    public static int[][] readLength(List<List<List<Cluster>>> array) {
        int[][] lenArray = new int[array.size()][];
        for (int i = 0; i < array.size(); i++) {
            lenArray[i] = new int[array.get(i).size()];
            for (int j = 0; j < array.get(i).size(); j++) {
                lenArray[i][j] = array.get(i).get(j).size();
            }
        }
        return lenArray;
    }

    // =====================================================================
    //  deriveCfBasisIndices
    //
    //  Extracts per-CF-column basis indices from the grouped CF data.
    //  Each CF cluster's sites carry symbols like "s1", "s2", etc.
    //  The integer suffix is the 1-based σ-power for that site.
    //
    //  The iteration order over (t, j, p) must match the CF column
    //  ordering used throughout the pipeline (same as lcf enumeration).
    //
    //  Result: cfBasisIndices[col] = int[] of 1-based σ-powers
    //  Example: a pair CF decorated (s1, s2) → {1, 2}
    //           a point CF decorated (s1)    → {1}
    // =====================================================================

    static int[][] deriveCfBasisIndices(GroupedCFData cfData, int[][] lcf) {
        // Total CFs = sum of lcf
        int totalCfs = 0;
        for (int[] row : lcf)
            for (int val : row)
                totalCfs += val;

        int[][] indices = new int[totalCfs][];
        int col = 0;

        List<List<List<Cluster>>> coordData = cfData.getCoordData();

        for (int t = 0; t < coordData.size(); t++) {
            List<List<Cluster>> groupsT = coordData.get(t);
            for (int j = 0; j < groupsT.size(); j++) {
                List<Cluster> cfList = groupsT.get(j);
                for (int p = 0; p < cfList.size(); p++) {
                    Cluster cfCluster = cfList.get(p);
                    List<Site> sites = cfCluster.getAllSites();
                    int[] basisIdx = new int[sites.size()];
                    for (int s = 0; s < sites.size(); s++) {
                        String symbol = sites.get(s).getSymbol();
                        // Extract integer from symbol "s1" → 1, "s2" → 2, etc.
                        basisIdx[s] = Integer.parseInt(symbol.substring(1));
                    }
                    indices[col++] = basisIdx;
                }
            }
        }
        return indices;
    }

    // =====================================================================
    //  genConfig[clusCoord, numElements]
    //
    //  Returns all configurations (Tuples of element indices 1..numElements)
    //  for the given cluster's sites.
    // =====================================================================

    public static List<int[]> genConfig(Cluster cluster, int numElements) {
        int numClusSite = cluster.getAllSites().size();
        // transBasis = Table[i, {i, 1, numElements}]
        int[] transBasis = new int[numElements];
        for (int i = 0; i < numElements; i++) transBasis[i] = i + 1;

        // Tuples[transBasis, numClusSite]
        List<int[]> configList = new ArrayList<>();
        generateTuples(transBasis, numClusSite, new int[numClusSite], 0, configList);
        return configList;
    }

    // =====================================================================
    //  calMultiplicity[clusterList, spaceGroup]
    //
    //  Computes multiplicity of each cluster, normalized by last entry.
    // =====================================================================

    public static List<Double> calMultiplicity(
            List<Cluster> clusterList, List<SymmetryOperation> spaceGroup) {

        List<Integer> mList = new ArrayList<>();
        for (Cluster cluster : clusterList) {
            List<Cluster> orbit = genOrbit(cluster, spaceGroup);
            mList.add(orbit.size());
        }
        int lastM = mList.get(mList.size() - 1);
        List<Double> result = new ArrayList<>();
        for (int m : mList) {
            result.add((double) m / lastM);
        }
        return result;
    }

    // =====================================================================
    //  groupSubClus[maxClusCoord, cfData, basisSymbolList]
    //
    //  Generates all decorated subclusters of maximal clusters and
    //  classifies them into CF orbits.
    //  Returns classifiedSubClusList[i][j][k] = list of matching subclusters
    // =====================================================================

    public static List<List<List<List<Cluster>>>> groupSubClus(
            List<Cluster> maxClusCoord,
            GroupedCFData cfData,
            List<String> basisSymbolList) {

        List<List<List<List<Cluster>>>> cfOrbitList = cfData.getOrbitData();
        List<List<List<List<Cluster>>>> classifiedSubClusList = new ArrayList<>();

        for (int i = 0; i < cfOrbitList.size(); i++) {
            classifiedSubClusList.add(new ArrayList<>());
            for (int j = 0; j < cfOrbitList.get(i).size(); j++) {
                classifiedSubClusList.get(i).add(new ArrayList<>());
                for (int k = 0; k < cfOrbitList.get(i).get(j).size(); k++) {
                    classifiedSubClusList.get(i).get(j).add(new ArrayList<>());

                    List<Cluster> cfOrbit = cfOrbitList.get(i).get(j).get(k);
                    for (int m = 0; m < maxClusCoord.size(); m++) {
                        List<Cluster> subClusCoordList =
                                genSubClusCoord(maxClusCoord.get(m), basisSymbolList);
                        for (Cluster subClus : subClusCoordList) {
                            if (isContained(cfOrbit, subClus)) {
                                classifiedSubClusList.get(i).get(j).get(k).add(subClus);
                            }
                        }
                    }
                }
            }
        }
        return classifiedSubClusList;
    }

    // =====================================================================
    //  groupCMatData[cMatData, lc]
    //
    //  Groups flat C-matrix data into HSP cluster-type groups.
    //  cMatData = {cMatList, lcv, wcv}
    //  Returns grouped {classifiedCMatList, classifiedLcv, classifiedWcv}
    // =====================================================================

    /**
     * Groups flat lists into sub-lists according to {@code lc} sizes.
     *
     * @param flatList the flat list to group
     * @param lc       per-HSP-type counts (excluding empty cluster)
     * @param <T>      element type
     * @return grouped list
     */
    public static <T> List<List<T>> groupByLc(List<T> flatList, int[] lc) {
        List<List<T>> result = new ArrayList<>();
        int index = 0;
        for (int count : lc) {
            List<T> group = new ArrayList<>();
            for (int j = 0; j < count; j++) {
                if (index < flatList.size()) {
                    group.add(flatList.get(index++));
                }
            }
            result.add(group);
        }
        return result;
    }

    // =====================================================================
    //  Helper: generateAllSubsets (Subsets equivalent)
    // =====================================================================

    private static List<List<Site>> generateAllSubsets(List<Site> sites) {
        List<List<Site>> subsets = new ArrayList<>();
        int n = sites.size();
        int total = 1 << n;
        for (int mask = 0; mask < total; mask++) {
            List<Site> subset = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    subset.add(sites.get(i));
                }
            }
            subsets.add(subset);
        }
        return subsets;
    }

    // =====================================================================
    //  Helper: cartesianProduct (Tuples equivalent)
    // =====================================================================

    private static List<List<Site>> cartesianProduct(List<List<Site>> lists) {
        List<List<Site>> result = new ArrayList<>();
        cartesianRecursive(lists, 0, new ArrayList<>(), result);
        return result;
    }

    private static void cartesianRecursive(
            List<List<Site>> lists, int depth,
            List<Site> current, List<List<Site>> result) {
        if (depth == lists.size()) {
            result.add(new ArrayList<>(current));
            return;
        }
        for (Site s : lists.get(depth)) {
            current.add(s);
            cartesianRecursive(lists, depth + 1, current, result);
            current.remove(current.size() - 1);
        }
    }

    // =====================================================================
    //  Helper: integer-tuple generation for genConfig
    // =====================================================================

    private static void generateTuples(
            int[] basis, int depth, int[] current, int pos, List<int[]> result) {
        if (pos == depth) {
            result.add(Arrays.copyOf(current, depth));
            return;
        }
        for (int b : basis) {
            current[pos] = b;
            generateTuples(basis, depth, current, pos + 1, result);
        }
    }

    // =====================================================================
    //  Helper: position / site matching
    // =====================================================================

    private static boolean positionsMatch(Position p1, Position p2) {
        return Math.abs(p1.getX() - p2.getX()) < DELTA
            && Math.abs(p1.getY() - p2.getY()) < DELTA
            && Math.abs(p1.getZ() - p2.getZ()) < DELTA;
    }

    private static boolean containsPosition(List<Site> sites, Position pos) {
        for (Site s : sites) {
            if (positionsMatch(s.getPosition(), pos)) return true;
        }
        return false;
    }

    private static boolean containsSite(List<Site> sites, Site target) {
        for (Site s : sites) {
            if (positionsMatch(s.getPosition(), target.getPosition())
                    && s.getSymbol().equals(target.getSymbol())) {
                return true;
            }
        }
        return false;
    }

    private static List<Integer> sublatticeSizes(Cluster c) {
        List<Integer> sizes = new ArrayList<>();
        for (Sublattice sub : c.getSublattices()) {
            sizes.add(sub.getSites().size());
        }
        return sizes;
    }


    private static void emit(Consumer<String> sink, String message) {
        if (sink != null) {
            sink.accept(message);
        }
    }

    /**
     * Orchestrates the full Stages 0-4 workflow, returning a consolidated PipelineResult.
     * (Logic moved from AllClusterData.identify)
     */
    public static PipelineResult runFullWorkflow(ClusterIdentificationRequest config, Consumer<String> progressSink) {
        java.util.logging.Logger LOG = java.util.logging.Logger.getLogger(ClusterCFIdentificationPipeline.class.getName());
        LOG.info("ClusterCFIdentificationPipeline.runFullWorkflow — START");

        int numComponents = config.getNumComponents();
        String structurePhase = config.getStructurePhase();
        String model = config.getModel();

        // 1. Load resources
        emit(progressSink, "\n[STAGE 0]: Loading Inputs...");
        List<Cluster> disorderedClusters = org.ce.model.storage.InputLoader.parseClusterFile(config.getDisorderedClusterFile());
        disorderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup disorderedSpaceGroup = org.ce.model.storage.InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());

        List<Cluster> orderedClusters = org.ce.model.storage.InputLoader.parseClusterFile(config.getOrderedClusterFile());
        orderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup orderedSpaceGroup = org.ce.model.storage.InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());

        // 2. Stage 1 & 2: Cluster + CF Identification
        emit(progressSink, "\n[STAGE 1/2]: Running Identification Pipeline...");
        PipelineResult partialResult = run(
                disorderedClusters,
                disorderedSpaceGroup.getOperations(),
                orderedClusters,
                orderedSpaceGroup.getOperations(),
                config.getTransformationMatrix(),
                new double[] { config.getTranslationVector().getX(),
                        config.getTranslationVector().getY(),
                        config.getTranslationVector().getZ() },
                numComponents,
                progressSink);

        // 3. Stage 3: C-Matrix foundation
        emit(progressSink, "\n[STAGE 3]: Running C-Matrix Pipeline...");
        CMatrixPipeline.CMatrixData matrixData = CMatrixPipeline.run(
                partialResult.toClusterIdentificationResult(),
                partialResult.toCFIdentificationResult(),
                disorderedClusters,
                numComponents,
                progressSink);

        // 4. Stage 4: CVCF Transformation
        emit(progressSink, "\n[STAGE 4]: Basis Transformation...");
        String parentStructure = resolveParentStructure(structurePhase);
        org.ce.model.cvm.CvCfBasis cvcfBasis = org.ce.model.cvm.CvCfBasis.generate(
                parentStructure,
                partialResult,
                matrixData,
                model, 
                progressSink);

        // Final bundle
        double[] equiX = new double[numComponents];
        java.util.Arrays.fill(equiX, 1.0 / numComponents);
        double[] vRandEqui = cvcfBasis.computeRandomCvcfCFs(equiX, partialResult);

        PipelineResult finalResult = new PipelineResult(
                partialResult.getDisClusData(),
                partialResult.getTcdis(), partialResult.getNxcdis(), partialResult.getNcdis(),
                partialResult.getMhdis(), partialResult.getNijTable(), partialResult.getKbdis(),
                partialResult.getPhaseClusterData(), partialResult.getOrdClusData(),
                partialResult.getLc(), partialResult.getTc(), partialResult.getNxc(), partialResult.getNc(),
                partialResult.getMh(), partialResult.getRc(),
                partialResult.getDisCFData(), partialResult.getTcfdis(),
                partialResult.getPhaseCFDataRaw(), partialResult.getOrdCFData(),
                partialResult.getCfData(), partialResult.getLcf(), partialResult.getTcf(), partialResult.getNxcf(), partialResult.getNcf(),
                partialResult.getCfBasisIndices(), numComponents,
                matrixData,
                partialResult.toCFIdentificationResult().getUNames(),
                cvcfBasis.cfNames,
                partialResult.toCFIdentificationResult().getEONames(),
                cvcfBasis.eciNames,
                vRandEqui
        );

        LOG.info("ClusterCFIdentificationPipeline.runFullWorkflow — EXIT");
        return finalResult;
    }
    private static String resolveParentStructure(String structure) {
        if (structure == null) return null;
        String base = structure.replace("_CVCF", "");
        if (base.equals("BCC_B2")) return "BCC_A2";
        if (base.equals("FCC_L12")) return "FCC_A1";
        return base;
    }
}
