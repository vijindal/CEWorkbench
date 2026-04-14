package org.ce.model.cluster;

import java.util.*;
import java.util.function.Consumer;

import static org.ce.model.cluster.ClusterPrimitives.*;
import static org.ce.model.cluster.ClusterResults.*;
import static org.ce.model.cluster.ClusterKeys.*;

/**
 * CMatrixPipeline — Self-contained Java translation of the Mathematica
 * C-matrix generation pipeline.
 *
 * <p>Each Mathematica function is translated 1-to-1 as a Java method.
 * Symbolic algebra ({@code Expand}, {@code /. rules}, {@code CoefficientArrays})
 * is replaced by numerical polynomial maps ({@code Map<SiteOpProductKey, Double>}).
 * No logic has been added or removed relative to the Mathematica source.</p>
 *
 * <h2>Mathematica functions → Java methods</h2>
 * <pre>
 *   genSiteList          → {@link #buildSiteList(List)}
 *   genRMat              → {@link #buildRMatrix(int)}
 *   genPRules            → R-matrix coefficients, applied in {@link #expandConfiguration}
 *   groupSubClus         → {@link #groupSubClusters}
 *   genCfSiteOpList      → {@link #buildCfSiteOpList}
 *   genSubstituteRules   → {@link #buildSubstituteRules}
 *   tranClus             → {@link #translateCluster(Cluster, List)}
 *   genConfig            → {@link #generateConfigurations(int, int)}
 *   genCV                → {@link #generateCMatrix}
 * </pre>
 */
public final class CMatrixPipeline {

    private CMatrixPipeline() {}

    // =====================================================================
    //  Output container — equivalent to Return[{cmat, lcv, wcv}] in genCV
    // =====================================================================

    public static final class CMatrixData {
        /** cmat[t][j] = coefficient matrix; rows = distinct CVs, cols = totalCFs + 1 */
        public final List<List<double[][]>> cmat;
        /** lcv[t][j] = number of distinct CV rows for cluster type t, group j */
        public final int[][] lcv;
        /** wcv[t][j] = weight array; wcv[t][j][v] = count of configs producing CV row v */
        public final List<List<int[]>> wcv;

        // --- Structural metadata for Stage 4 reuse ---
        /** siteList = unique coordinates of all maximal clusters */
        public final List<Position> siteList;
        /** pRules = expansion coefficients for site operators */
        public final ClusterMath.PRules pRules;
        /** substituteRules = mapping from site-op products to CF indices */
        public final SubstituteRules substituteRules;
        /** cfBasisIndices[col] = basis indices for each CF column */
        public final int[][] cfBasisIndices;

        public List<List<double[][]>> getCmat() { return cmat; }
        public int[][] getLcv() { return lcv; }
        public List<List<int[]>> getWcv() { return wcv; }
        public List<Position> getSiteList() { return siteList; }
        public ClusterMath.PRules getPRules() { return pRules; }
        public SubstituteRules getSubstituteRules() { return substituteRules; }
        public int[][] getCfBasisIndices() { return cfBasisIndices; }

        /**
         * Transforms all C-matrix blocks into a new basis using transformation matrix T.
         * C_new = C_orth Â· T
         */
        public CMatrixData transform(double[][] T) {
            List<List<double[][]>> newCmat = new ArrayList<>();
            for (List<double[][]> typeBlock : cmat) {
                List<double[][]> newTypeBlock = new ArrayList<>();
                for (double[][] groupBlock : typeBlock) {
                    newTypeBlock.add(LinearAlgebra.multiply(groupBlock, T));
                }
                newCmat.add(newTypeBlock);
            }
            return new CMatrixData(newCmat, lcv, wcv, siteList, pRules, substituteRules, cfBasisIndices);
        }

        /**
         * Converts a generic product of physical site-atom probabilities into an
         * equivalent row vector expressed in the Orthogonal CF basis.
         *
         * @param siteIndices List of physical site indices.
         * @param config Sequence of elements/atoms corresponding strictly to siteIndices.
         * @param cfColumnMap Pre-compiled column mapping from CFIndex to orthogonal column.
         * @param totalCfs Total number of correlation functions (width of returned array is totalCfs + 1).
         * @return Orthogonal row vector describing the coefficients.
         */
        public double[] expandProbabilityExpression(
                List<Integer> siteIndices,
                int[] config,
                Map<org.ce.model.cluster.ClusterKeys.CFIndex, Integer> cfColumnMap,
                int totalCfs) {
            
            return CMatrixPipeline.expandConfiguration(
                    siteIndices, config,
                    this.pRules.getRMatrix(),
                    this.substituteRules.getRules(),
                    cfColumnMap, totalCfs
            );
        }

        CMatrixData(
                List<List<double[][]>> cmat,
                int[][] lcv,
                List<List<int[]>> wcv,
                List<Position> siteList,
                ClusterMath.PRules pRules,
                SubstituteRules substituteRules,
                int[][] cfBasisIndices) {
            this.cmat = cmat;
            this.lcv  = lcv;
            this.wcv  = wcv;
            this.siteList = siteList;
            this.pRules = pRules;
            this.substituteRules = substituteRules;
            this.cfBasisIndices = cfBasisIndices;
        }
    }

    // =====================================================================
    //  MAIN ENTRY POINT — matches Mathematica script block line-for-line
    //
    //  maxClusSiteList    = genSiteList[maxClusCoord];
    //  groupCfCoordList   = groupSubClus[maxClusCoord, cfData, basisSymbolList];
    //  cfSiteOpList       = genCfSiteOpList[groupCfCoordList, maxClusSiteList];
    //  substituteRules    = genSubstituteRules[cfSiteOpList, cfSymbol];
    //  numSites           = Length[maxClusSiteList];
    //  pRules             = genPRules[numSites, numComp, siteOcSymbol, siteOpSymbol];
    //  cMatData           = genCV[maxClusSiteList, ordClusData, ordCFData,
    //                             substituteRules, pRules, numComp, ...];
    // =====================================================================

    public static CMatrixData run(
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            List<Cluster> maxClusters,
            int numElements,
            Consumer<String> sink) {

        emit(sink, "========================================================");
        emit(sink, "  C-MATRIX PIPELINE (Mathematica translation)");
        emit(sink, "========================================================");

        // maxClusSiteList = genSiteList[maxClusCoord]
        List<Position> siteList = buildSiteList(maxClusters);
        emit(sink, "\n[1] genSiteList -> " + siteList.size() + " unique sites");

        // genRMat[numComp]
        double[][] rMat = buildRMatrix(numElements);
        emit(sink, "[2] genRMat -> " + numElements + "x" + numElements + " R-matrix");

        // basisSymbolList = genBasisSymbolList[numComp, siteOpSymbol]
        List<String> basisSymbolList = new ArrayList<>();
        for (int s = 1; s < numElements; s++) {
            basisSymbolList.add("s" + s);
        }

        // groupCfCoordList = groupSubClus[maxClusCoord, cfData, basisSymbolList]
        GroupedCFResult groupedCF = cfResult.getGroupedCFData();
        List<List<List<List<Cluster>>>> cfOrbitList = groupedCF.getOrbitData();
        List<List<List<List<Cluster>>>> groupCfCoordList =
                groupSubClusters(maxClusters, cfOrbitList, basisSymbolList, sink);

        // cfSiteOpList = genCfSiteOpList[groupCfCoordList, maxClusSiteList]
        List<List<List<List<List<SiteOp>>>>> cfSiteOpList =
                buildCfSiteOpList(groupCfCoordList, siteList, sink);

        // substituteRules = genSubstituteRules[cfSiteOpList, cfSymbol]
        SubstituteRules substituteRules =
                new SubstituteRules(buildSubstituteRules(cfSiteOpList, sink));
        emit(sink, "[3] substituteRules -> " + substituteRules.getRules().size() + " rules");

        // uList / CF column map
        int[][] lcf = cfResult.getLcf();
        int totalCfs = cfResult.getTcf();
        Map<CFIndex, Integer> cfColumnMap = buildCfColumnMap(lcf);
        int[][] cfBasisIndices = deriveCfBasisIndices(cfSiteOpList, cfColumnMap, totalCfs);
        emit(sink, "[4] CF column map -> " + totalCfs + " total CFs");

        // pRules equivalent: ClusterMath.PRules (Stage 4 reuse)
        ClusterMath.PRules pRules = ClusterMath.PRules.build(siteList.size(), numElements);

        // cMatData = genCV[...]
        List<List<Cluster>> ordClusCoordList =
                clusterResult.getOrdClusterData().getCoordList();
        CMatrixData result = generateCMatrix(
                ordClusCoordList, siteList, rMat, substituteRules.getRules(),
                cfColumnMap, totalCfs, numElements, sink);

        CMatrixData finalResult = new CMatrixData(
                result.cmat, result.lcv, result.wcv,
                siteList, pRules, substituteRules, cfBasisIndices);

        printResults(finalResult, sink);
        emit(sink, "\n========================================================");
        emit(sink, "  PIPELINE COMPLETE");
        emit(sink, "========================================================");
        return finalResult;
    }

    // =====================================================================
    //  genSiteList[maxClusCoord]
    //
    //  For[i, maxClusters]
    //    For[j, sublattices]
    //      For[k, sites]
    //        If[!MemberQ[siteList, pos], siteList = Append[siteList, pos]]
    //  Return[siteList]
    // =====================================================================

    static List<Position> buildSiteList(List<Cluster> maxClusters) {
        List<Position> siteList = new ArrayList<>();
        for (Cluster cluster : maxClusters) {
            for (Sublattice sub : cluster.getSublattices()) {
                for (Site site : sub.getSites()) {
                    Position pos = site.getPosition();
                    if (!containsPosition(siteList, pos)) {
                        siteList.add(pos);
                    }
                }
            }
        }
        return siteList;
    }

    // =====================================================================
    //  genRMat[numElements]
    //
    //  Even K: basis = {-K/2, ..., -1, 1, ..., K/2}
    //  Odd  K: basis = {-(K-1)/2, ..., 0, ..., (K-1)/2}
    //  matM[i][j] = If[i==0, 1, basis[[j]]^i]
    //  RMat = Inverse[matM]
    // =====================================================================

    static double[][] buildRMatrix(int numElements) {
        return ClusterMath.buildRMatrix(numElements);
    }

    // =====================================================================
    //  groupSubClus[maxClusCoord, cfData, basisSymbolList]
    //
    //  cfOrbitList = cfData[[3]];
    //  For[i, cfOrbitList]                  (* cluster types *)
    //    For[j, cfOrbitList[[i]]]            (* CF groups  *)
    //      For[k, cfOrbitList[[i]][[j]]]      (* CFs in group *)
    //        For[m, maxClusCoord]              (* maximal clusters *)
    //          subClusCoordList = genSubClusCoord[maxClusCoord[[m]], basisSymbolList]
    //          For[l, subClusCoordList]
    //            If[isContained[cfOrbitList[[i]][[j]][[k]], subClusCoordList[[l]]],
    //              Append[..., subClus]]
    //  Return[classifiedSubClusList]
    // =====================================================================

    static List<List<List<List<Cluster>>>> groupSubClusters(
            List<Cluster> maxClusters,
            List<List<List<List<Cluster>>>> cfOrbitList,
            List<String> basisSymbolList,
            Consumer<String> sink) {

        emit(sink, "\n  [3a] groupSubClus");
        List<List<List<List<Cluster>>>> classifiedSubClusList = new ArrayList<>();

        for (int i = 0; i < cfOrbitList.size(); i++) {
            List<List<List<Cluster>>> typeLevel = new ArrayList<>();
            for (int j = 0; j < cfOrbitList.get(i).size(); j++) {
                List<List<Cluster>> groupLevel = new ArrayList<>();
                for (int k = 0; k < cfOrbitList.get(i).get(j).size(); k++) {
                    List<Cluster> cfOrbit = cfOrbitList.get(i).get(j).get(k);
                    List<Cluster> matched = new ArrayList<>();

                    for (Cluster maxClus : maxClusters) {
                        List<Cluster> subClusCoordList =
                                ClusterCFIdentificationPipeline.genSubClusCoord(maxClus, basisSymbolList);
                        for (Cluster subClus : subClusCoordList) {
                            if (ClusterUtils.isContained(cfOrbit, subClus)) {
                                matched.add(subClus);
                            }
                        }
                    }

                    groupLevel.add(matched);
                    emit(sink, "    classified[" + i + "][" + j + "][" + k + "] = "
                            + matched.size() + " subclusters");
                }
                typeLevel.add(groupLevel);
            }
            classifiedSubClusList.add(typeLevel);
        }

        return classifiedSubClusList;
    }

    // =====================================================================
    //  genCfSiteOpList[groupClusCoordList, siteList]
    //
    //  For[i=1, i<Length[groupClusCoordList], i++,   (* NOTE: strict < *)
    //    For[j, groups]
    //      For[k, cluster groups]
    //        For[l, individual clusters in group]
    //          siteOp = {};
    //          For[n, sublattices]
    //            For[m, sites]
    //              siteOp = Append[siteOp,
    //                tempClusCoord[[n]][[m]][[2]][
    //                  Position[siteList, tempClusCoord[[n]][[m]][[1]]][[1]][[1]]
    //                ]]
    //          rules[[i]][[j]][[k]] = Append[..., siteOp]
    //  Return[rules]
    // =====================================================================

    static List<List<List<List<List<SiteOp>>>>> buildCfSiteOpList(
            List<List<List<List<Cluster>>>> groupClusCoordList,
            List<Position> siteList,
            Consumer<String> sink) {

        emit(sink, "  [3b] genCfSiteOpList");
        List<List<List<List<List<SiteOp>>>>> rules = new ArrayList<>();

        for (int i = 0; i < groupClusCoordList.size(); i++) {
            List<List<List<List<SiteOp>>>> typeLevel = new ArrayList<>();
            for (int j = 0; j < groupClusCoordList.get(i).size(); j++) {
                List<List<List<SiteOp>>> groupLevel = new ArrayList<>();
                for (int k = 0; k < groupClusCoordList.get(i).get(j).size(); k++) {
                    List<List<SiteOp>> cfLevel = new ArrayList<>();

                    for (Cluster subCluster : groupClusCoordList.get(i).get(j).get(k)) {
                        List<SiteOp> siteOp = new ArrayList<>();
                        for (Sublattice sub : subCluster.getSublattices()) {
                            for (Site site : sub.getSites()) {
                                int siteIdx = positionIndexOf(siteList, site.getPosition());
                                if (siteIdx < 0) {
                                    throw new IllegalStateException(
                                            "Site not found in siteList: " + site.getPosition());
                                }
                                int basisIdx = parseBasisIndex(site.getSymbol());
                                siteOp.add(new SiteOp(siteIdx, basisIdx));
                            }
                        }
                        cfLevel.add(siteOp);
                    }

                    groupLevel.add(cfLevel);
                }
                typeLevel.add(groupLevel);
            }
            rules.add(typeLevel);
        }

        emit(sink, "    cfSiteOpList built: " + rules.size() + " types");
        return rules;
    }

    // =====================================================================
    //  genSubstituteRules[cfSiteOpList, cfSymbol]
    // =====================================================================

    /**
     * Mapping from site-operator products to CF indices.
     * Groups structural mapping logic directly with the pipeline.
     */
    public static final class SubstituteRules {
        private final Map<SiteOpProductKey, CFIndex> rules;

        public SubstituteRules(Map<SiteOpProductKey, CFIndex> rules) {
            this.rules = Collections.unmodifiableMap(new LinkedHashMap<>(rules));
        }

        public CFIndex lookup(List<SiteOp> ops) {
            return rules.get(new SiteOpProductKey(ops));
        }

        public CFIndex lookup(SiteOpProductKey key) {
            return rules.get(key);
        }

        public Map<SiteOpProductKey, CFIndex> getRules() {
            return rules;
        }
    }

    static Map<SiteOpProductKey, CFIndex> buildSubstituteRules(
            List<List<List<List<List<SiteOp>>>>> cfSiteOpList,
            Consumer<String> sink) {

        emit(sink, "  [3c] genSubstituteRules");
        Map<SiteOpProductKey, CFIndex> rules = new LinkedHashMap<>();

        for (int i = 0; i < cfSiteOpList.size(); i++) {
            for (int j = 0; j < cfSiteOpList.get(i).size(); j++) {
                for (int k = 0; k < cfSiteOpList.get(i).get(j).size(); k++) {
                    for (List<SiteOp> tempClusCoord : cfSiteOpList.get(i).get(j).get(k)) {
                        // siteOp = 1; For[n, siteOp = siteOp * tempClusCoord[[n]]]
                        // → canonical SiteOpProductKey
                        SiteOpProductKey productKey = new SiteOpProductKey(tempClusCoord);
                        CFIndex cfIndex = new CFIndex(i, j, k);
                        rules.putIfAbsent(productKey, cfIndex);
                    }
                }
            }
        }

        emit(sink, "    substitute rules: " + rules.size() + " unique entries");
        return rules;
    }

    // =====================================================================
    //  uList = Flatten[Table[..cfSymbol[i][j][k]..], 2]
    //  Maps each CF triplet (t,j,k) to a flat column index.
    // =====================================================================

    private static Map<CFIndex, Integer> buildCfColumnMap(int[][] lcf) {
        Map<CFIndex, Integer> map = new LinkedHashMap<>();
        int col = 0;
        for (int t = 0; t < lcf.length; t++) {
            for (int j = 0; j < lcf[t].length; j++) {
                for (int k = 0; k < lcf[t][j]; k++) {
                    map.put(new CFIndex(t, j, k), col++);
                }
            }
        }
        return map;
    }

    // =====================================================================
    //  tranClus[clusCoord, maxClusSiteList]
    //
    //  For[i, sublattices]
    //    For[j, sites]
    //      Position[maxClusSiteList, clusCoord[[i]][[j]][[1]]][[1]][[1]]
    //  (Returns flat list after Flatten[..., 1] in genCV)
    // =====================================================================

    static List<Integer> translateCluster(Cluster cluster, List<Position> siteList) {
        List<Integer> indices = new ArrayList<>();
        for (Sublattice sub : cluster.getSublattices()) {
            for (Site site : sub.getSites()) {
                int idx = positionIndexOf(siteList, site.getPosition());
                if (idx < 0) {
                    throw new IllegalStateException(
                            "Cluster site not found in siteList: " + site.getPosition());
                }
                indices.add(idx);
            }
        }
        return indices;
    }

    // =====================================================================
    //  genConfig[clusCoord, numElements]
    //  All K^N element assignments for an N-site cluster.
    // =====================================================================

    static List<int[]> generateConfigurations(int numSites, int K) {
        int total = 1;
        for (int i = 0; i < numSites; i++) total *= K;
        List<int[]> configs = new ArrayList<>(total);
        for (int i = 0; i < total; i++) {
            int[] cfg = new int[numSites];
            int x = i;
            for (int s = numSites - 1; s >= 0; s--) {
                cfg[s] = x % K;
                x /= K;
            }
            configs.add(cfg);
        }
        return configs;
    }

    // =====================================================================
    //  Numerical equivalent of:
    //    tempCV = 1
    //    For[l, tempCV = tempCV * siteOcSymbol[clusSiteList[[l]]][config[[l]]]]
    //    tempCV = Expand[tempCV /. pRules] /. substituteRules
    //
    //  Polynomial is Map<SiteOpProductKey, Double>, grown site-by-site using
    //  R-matrix coefficients (= genPRules inline). After all sites, each
    //  monomial key is looked up in substituteRules for its CF column.
    // =====================================================================

    static double[] expandConfiguration(
            List<Integer> siteIndices,
            int[] config,
            double[][] rMat,
            Map<SiteOpProductKey, CFIndex> substituteRules,
            Map<CFIndex, Integer> cfColumnMap,
            int totalCfs) {

        int K = rMat.length;
        Map<SiteOpProductKey, Double> poly = new LinkedHashMap<>();
        poly.put(new SiteOpProductKey(List.of()), 1.0);

        for (int s = 0; s < siteIndices.size(); s++) {
            int globalSite = siteIndices.get(s);
            int element    = config[s];
            double[] R_row = rMat[element];
            Map<SiteOpProductKey, Double> nextPoly = new LinkedHashMap<>();

            for (Map.Entry<SiteOpProductKey, Double> term : poly.entrySet()) {
                List<SiteOp> existingOps = term.getKey().getOps();
                double existingCoeff     = term.getValue();

                for (int a = 0; a < K; a++) {
                    double c = R_row[a];
                    if (Math.abs(c) < 1e-14) continue;
                    List<SiteOp> newOps = new ArrayList<>(existingOps);
                    if (a > 0) {
                        newOps.add(new SiteOp(globalSite, a));
                    }
                    SiteOpProductKey newKey = new SiteOpProductKey(newOps);
                    nextPoly.merge(newKey, existingCoeff * c, Double::sum);
                }
            }
            poly = nextPoly;
        }

        // /. substituteRules + CoefficientArrays[cvList, uList]
        double[] row = new double[totalCfs + 1];   // +1 for constant column
        for (Map.Entry<SiteOpProductKey, Double> term : poly.entrySet()) {
            SiteOpProductKey key = term.getKey();
            double coeff         = term.getValue();
            if (Math.abs(coeff) < 1e-14) continue;

            if (key.getOps().isEmpty()) {
                row[totalCfs] += coeff;          // constant term -> last column
            } else {
                CFIndex cfIdx = substituteRules.get(key);
                if (cfIdx == null) {
                    throw new IllegalStateException(
                            "No substitute rule for site-op product: " + key);
                }
                Integer col = cfColumnMap.get(cfIdx);
                if (col == null) {
                    throw new IllegalStateException(
                            "No column mapping for CF index: " + cfIdx);
                }
                row[col] += coeff;
            }
        }
        return row;
    }

    // =====================================================================
    //  genCV[maxClusSiteList, ordClusData, ordCFData, substituteRules,
    //        pRules, numElements, ...]
    //
    //  For[i=1, i<=Length[ordClusCoordList]-1, i++,   (* cluster types *)
    //    For[j, ordClusCoordList[[i]]]                 (* clusters *)
    //      clusSiteList = Flatten[tranClus[...], 1]
    //      configList   = genConfig[..., numElements]
    //      For[k, configList]
    //        tempCV = product of occupation ops
    //        tempCV = Expand[tempCV /. pRules] /. substituteRules
    //      tallyCVList = Tally[tempCVList]
    //      tempArray   = CoefficientArrays[cvList, uList]
    //      crow        = Join[tempArray[[2]], Partition[tempArray[[1]],{1}], 2]
    //  Return[{cmat, lcv, wcv}]
    // =====================================================================

    private static CMatrixData generateCMatrix(
            List<List<Cluster>> ordClusCoordList,
            List<Position> siteList,
            double[][] rMat,
            Map<SiteOpProductKey, CFIndex> substituteRules,
            Map<CFIndex, Integer> cfColumnMap,
            int totalCfs,
            int numElements,
            Consumer<String> sink) {

        emit(sink, "\n[5] genCV -- " + ordClusCoordList.size() + " cluster types, "
                + totalCfs + " CFs, K=" + numElements);

        List<List<double[][]>> cmat = new ArrayList<>();
        List<List<int[]>>      wcv  = new ArrayList<>();
        int[][] lcv = new int[ordClusCoordList.size()][];

        for (int t = 0; t < ordClusCoordList.size(); t++) {
            List<Cluster> groups = ordClusCoordList.get(t);
            List<double[][]> cmatType = new ArrayList<>();
            List<int[]>      wcvType  = new ArrayList<>();
            lcv[t] = new int[groups.size()];

            for (int j = 0; j < groups.size(); j++) {
                Cluster cluster = groups.get(j);
                List<Integer> siteIndices = translateCluster(cluster, siteList);
                List<int[]> configs = generateConfigurations(siteIndices.size(), numElements);

                // Expand all configs and tally
                Map<RowKey, double[]> seenRows   = new LinkedHashMap<>();
                Map<RowKey, Integer>  seenCounts = new LinkedHashMap<>();
                for (int[] config : configs) {
                    double[] row = expandConfiguration(
                            siteIndices, config, rMat,
                            substituteRules, cfColumnMap, totalCfs);
                    RowKey key = new RowKey(row);
                    seenRows.putIfAbsent(key, row);
                    seenCounts.merge(key, 1, Integer::sum);
                }

                // Build cmat block and wcv
                lcv[t][j] = seenRows.size();
                double[][] cmatBlock = new double[seenRows.size()][];
                int[]      wcvBlock  = new int[seenRows.size()];
                int idx = 0;
                for (Map.Entry<RowKey, double[]> entry : seenRows.entrySet()) {
                    cmatBlock[idx] = entry.getValue();
                    wcvBlock[idx]  = seenCounts.get(entry.getKey());
                    idx++;
                }

                cmatType.add(cmatBlock);
                wcvType.add(wcvBlock);

                emit(sink, "  t=" + t + " j=" + j + ": "
                        + siteIndices.size() + " sites, ncv=" + seenRows.size()
                        + ", wcv=" + Arrays.toString(wcvBlock));
            }

            cmat.add(cmatType);
            wcv.add(wcvType);
        }

        return new CMatrixData(cmat, lcv, wcv, null, null, null, null);
    }

    // =====================================================================
    //  Tolerance-aware row key for Tally[tempCVList]
    // =====================================================================

    private static final class RowKey {
        private final long[] rounded;

        RowKey(double[] row) {
            rounded = new long[row.length];
            for (int i = 0; i < row.length; i++) {
                rounded[i] = Math.round(row[i] * 1e10);
            }
        }

        @Override
        public boolean equals(Object o) {
            return o instanceof RowKey && Arrays.equals(rounded, ((RowKey) o).rounded);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(rounded);
        }
    }

    // =====================================================================
    //  Inverse[matM] — Gaussian elimination with partial pivoting
    // =====================================================================

    private static double[][] invertMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] aug = new double[n][2 * n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, aug[i], 0, n);
            aug[i][n + i] = 1.0;
        }

        for (int col = 0; col < n; col++) {
            int maxRow = col;
            for (int row = col + 1; row < n; row++) {
                if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) {
                    maxRow = row;
                }
            }
            double[] tmp = aug[col]; aug[col] = aug[maxRow]; aug[maxRow] = tmp;

            double pivot = aug[col][col];
            if (Math.abs(pivot) < 1e-15) {
                throw new ArithmeticException("Singular matrix in R-matrix computation");
            }
            for (int j = 0; j < 2 * n; j++) aug[col][j] /= pivot;

            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double factor = aug[row][col];
                for (int j = 0; j < 2 * n; j++) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        double[][] inv = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(aug[i], n, inv[i], 0, n);
        }
        return inv;
    }

    // =====================================================================
    //  Position and symbol helpers
    // =====================================================================

    private static int positionIndexOf(List<Position> siteList, Position pos) {
        for (int i = 0; i < siteList.size(); i++) {
            if (siteList.get(i).equals(pos)) return i;
        }
        return -1;
    }

    private static boolean containsPosition(List<Position> list, Position pos) {
        for (Position p : list) {
            if (p.equals(pos)) return true;
        }
        return false;
    }

    private static int parseBasisIndex(String symbol) {
        if (symbol == null || !symbol.startsWith("s")) {
            throw new IllegalArgumentException("Invalid basis symbol: " + symbol);
        }
        return Integer.parseInt(symbol.substring(1));
    }

    // =====================================================================
    //  Result printing
    // =====================================================================

    private static void printResults(CMatrixData result, Consumer<String> sink) {
        if (sink == null) return;

        emit(sink, "\n[RESULT] lcv:");
        for (int t = 0; t < result.lcv.length; t++) {
            emit(sink, "  t=" + t + ": " + Arrays.toString(result.lcv[t]));
        }

        emit(sink, "\n[RESULT] wcv:");
        for (int t = 0; t < result.wcv.size(); t++) {
            for (int j = 0; j < result.wcv.get(t).size(); j++) {
                emit(sink, "  t=" + t + " j=" + j + ": "
                        + Arrays.toString(result.wcv.get(t).get(j)));
            }
        }

        emit(sink, "\n[RESULT] cmat:");
        for (int t = 0; t < result.cmat.size(); t++) {
            for (int j = 0; j < result.cmat.get(t).size(); j++) {
                double[][] block = result.cmat.get(t).get(j);
                emit(sink, "  cmat[" + t + "][" + j + "] = "
                        + block.length + " rows x "
                        + (block.length > 0 ? block[0].length : 0) + " cols");
                for (int r = 0; r < block.length; r++) {
                    StringBuilder sb = new StringBuilder("    [");
                    for (int c = 0; c < block[r].length; c++) {
                        sb.append(String.format(" %8.5f", block[r][c]));
                    }
                    emit(sink, sb.append(" ]").toString());
                }
            }
        }
    }

    static int[][] deriveCfBasisIndices(
            List<List<List<List<List<SiteOp>>>>> cfSiteOpList,
            Map<CFIndex, Integer> cfColumnMap,
            int totalCfs) {

        int[][] result = new int[totalCfs][];
        for (int t = 0; t < cfSiteOpList.size(); t++) {
            for (int j = 0; j < cfSiteOpList.get(t).size(); j++) {
                for (int k = 0; k < cfSiteOpList.get(t).get(j).size(); k++) {
                    Integer col = cfColumnMap.get(new CFIndex(t, j, k));
                    if (col != null) {
                        // Use the first variant (index 0) of the subcluster
                        // membership as representative for basis indices.
                        List<SiteOp> ops = cfSiteOpList.get(t).get(j).get(k).get(0);
                        result[col] = new int[ops.size()];
                        for (int m = 0; m < ops.size(); m++) {
                            result[col][m] = ops.get(m).getBasisIndex();
                        }
                    }
                }
            }
        }
        return result;
    }

    // =====================================================================
    //  CV Verification Logic (Self-Sufficient)
    // =====================================================================

    /**
     * Verifies the generated C-matrix by evaluating cluster variables at a
     * given composition under the disordered (random) state assumption.
     *
     * <p>At equiatomic composition (x_i = 1/K), for a cluster of size N,
     * every configuration is equally likely with probability 1/(K^N).
     * This method verifies that the C-matrix correctly reproduces these
     * statistical weights.</p>
     *
     * @param moleFractions  mole fractions (length K, Σ = 1)
     * @param pipelineResult results from Stage 2 (provides random CFs)
     * @param cmatData       results from Stage 3 (provides C-matrix)
     * @param sink           diagnostic output
     */
    public static void verifyRandomCVs(
            double[] moleFractions,
            ClusterCFIdentificationPipeline.PipelineResult pipelineResult,
            CMatrixData cmatData,
            Consumer<String> sink) {

        emit(sink, "\n========================================================");
        emit(sink, "  CV VERIFICATION: Random State at Composition");
        emit(sink, "========================================================");
        emit(sink, "  Composition: " + Arrays.toString(moleFractions));

        // 1. Get random CF vector (full length ncf + K)
        double[] uFull = pipelineResult.computeRandomCFs(moleFractions);

        // 2. Evaluate CVs using local logic (self-sufficient)
        double[][][] cv = evaluateCVs(uFull, cmatData.cmat, cmatData.lcv,
                pipelineResult.getTcdis(), pipelineResult.getLc());

        int K = moleFractions.length;
        boolean isEquiatomic = true;
        for (double x : moleFractions) {
            if (Math.abs(x - 1.0 / K) > 1e-6) {
                isEquiatomic = false;
                break;
            }
        }

        for (int t = 0; t < cv.length; t++) {
            // Representative cluster for this type provides the site count
            Cluster representative = pipelineResult.getDisClusData().getClusCoordList().get(t);
            int n = representative.getAllSites().size();
            double expected = Math.pow(1.0 / K, n);

            for (int j = 0; j < cv[t].length; j++) {
                emit(sink, String.format("\n  Cluster Type t=%d, Group j=%d (n=%d sites):", t, j, n));
                for (int v = 0; v < cv[t][j].length; v++) {
                    double val = cv[t][j][v];
                    String msg = String.format("    CV[%2d] = %12.8f", v, val);
                    if (isEquiatomic) {
                        double error = Math.abs(val - expected);
                        msg += String.format(" (Expected: %12.8f, Diff: %.2e)", expected, error);
                        if (error > 1e-9) msg += " [!] DISCREPANCY";
                    }
                    emit(sink, msg);
                }
            }
        }
        emit(sink, "\n=== CV Verification: COMPLETE ===");
    }

    /**
     * Builds the full CF vector (Tier iii) for C-matrix multiplication in orthogonal basis.
     * Length = ncf + K. 
     * [0..ncf-1] = multi-site CFs, [ncf..ncf+K-2] = point CFs, [ncf+K-1] = constant 1.0.
     */
    public static double[] buildFullCFVector(double[] u, double[] moleFractions, 
                                            int numComponents, int[][] cfBasisIndices, int ncf) {
        int K = numComponents;
        double[] uFull = new double[ncf + K];
        System.arraycopy(u, 0, uFull, 0, ncf);

        double[] basis = ClusterMath.buildBasis(K);
        for (int k = 0; k < K - 1; k++) {
            double pcf = 0;
            for (int i = 0; i < K; i++) pcf += moleFractions[i] * Math.pow(basis[i], k + 1);
            int col = ncf + k;
            uFull[col] = pcf;
        }
        uFull[ncf + K - 1] = 1.0;
        return uFull;
    }

    /**
     * Builds the full CVCF vector (Tier iii) for C-matrix multiplication.
     * Length = ncf + K.
     * [0..ncf-1] = non-point variables, [ncf..ncf+K-1] = mole fractions.
     */
    public static double[] buildFullCVCFVector(double[] v, double[] x, int ncf) {
        int K = x.length;
        double[] vFull = new double[ncf + K];
        System.arraycopy(v, 0, vFull, 0, ncf);
        System.arraycopy(x, 0, vFull, ncf, K);
        return vFull;
    }

    /**
     * Internal CV evaluator to ensure self-sufficiency.
     * CV[t][j][v] = Σ_k cmat[t][j][v][k] * uFull[k]
     */
    public static double[][][] evaluateCVs(
            double[] uFull,
            List<List<double[][]>> cmat,
            int[][] lcv,
            int tcdis,
            int[] lc) {

        int width = uFull.length;
        double[][][] cv = new double[tcdis][][];

        for (int t = 0; t < tcdis; t++) {
            cv[t] = new double[lc[t]][];
            for (int j = 0; j < lc[t]; j++) {
                double[][] block = cmat.get(t).get(j);
                int nv = lcv[t][j];
                cv[t][j] = new double[nv];
                for (int v = 0; v < nv; v++) {
                    double val = 0.0;
                    for (int k = 0; k < width; k++) {
                        val += block[v][k] * uFull[k];
                    }
                    cv[t][j][v] = val;
                }
            }
        }
        return cv;
    }

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }
}
