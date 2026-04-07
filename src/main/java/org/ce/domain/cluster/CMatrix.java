package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterKeys.*;
import static org.ce.domain.cluster.ClusterPrimitives.*;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.cluster.cvcf.CvCfBasisTransformer;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * High-level API for C-matrix operations: building and data container.
 */
public final class CMatrix {

    private static final Logger LOG = Logger.getLogger(CMatrix.class.getName());

    private CMatrix() {
        // Utility container
    }

    /**
     * Holds C-matrix data: coefficients, counts, weights, and random CF generation.
     */
    public static final class Result {

        private final List<List<double[][]>> cmat;
        private final int[][] lcv;
        private final List<List<int[]>> wcv;
        private final int[][] cfBasisIndices;
        private final List<String> cmatCfNames;

        @JsonCreator
        public Result(
                @JsonProperty("cmat")            List<List<double[][]>> cmat,
                @JsonProperty("lcv")             int[][]                lcv,
                @JsonProperty("wcv")             List<List<int[]>>      wcv,
                @JsonProperty("cfBasisIndices")  int[][]                cfBasisIndices,
                @JsonProperty("cmatCfNames")     List<String>           cmatCfNames) {
            this.cmat = cmat;
            this.lcv = lcv;
            this.wcv = wcv;
            this.cfBasisIndices = cfBasisIndices;
            this.cmatCfNames = cmatCfNames;
        }

        public Result() {
            this.cmat = null;
            this.lcv = null;
            this.wcv = null;
            this.cfBasisIndices = null;
            this.cmatCfNames = null;
        }

        public List<List<double[][]>> getCmat() { return cmat; }
        public int[][] getLcv() { return lcv; }
        public List<List<int[]>> getWcv() { return wcv; }
        public int[][] getCfBasisIndices() { return cfBasisIndices; }
        public List<String> getCmatCfNames() { return cmatCfNames; }

        public void validateCols(int expectedCols, String contextMsg) {
            if (cmat == null || cmat.isEmpty()) return;
            var firstType = cmat.get(0);
            if (firstType == null || firstType.isEmpty()) return;
            double[][] firstBlock = firstType.get(0);
            if (firstBlock == null || firstBlock.length == 0) return;
            int actual = firstBlock[0].length;
            if (actual != expectedCols) {
                throw new IllegalStateException(contextMsg
                        + ": expected " + expectedCols + " columns, got " + actual);
            }
        }

        public double[] evaluateRandomCFs(double[] moleFractions, int numElements) {
            if (moleFractions == null || moleFractions.length != numElements) {
                throw new IllegalArgumentException("Invalid moleFractions");
            }
            if (cfBasisIndices == null) {
                throw new IllegalStateException("cfBasisIndices is null");
            }

            double[] basis = ClusterMath.buildBasis(numElements);
            int nxcf = numElements - 1;
            double[] pointCFs = new double[nxcf];
            for (int k = 0; k < nxcf; k++) {
                int power = k + 1;
                for (int i = 0; i < numElements; i++) {
                    pointCFs[k] += moleFractions[i] * Math.pow(basis[i], power);
                }
            }

            int ncf = cfBasisIndices.length - nxcf;
            double[] uRandom = new double[ncf];
            for (int icf = 0; icf < ncf; icf++) {
                int[] indices = cfBasisIndices[icf];
                double val = 1.0;
                for (int b : indices) val *= pointCFs[b - 1];
                uRandom[icf] = val;
            }
            return uRandom;
        }

        public void printSummary(String title, java.util.function.Consumer<String> sink) {
            sink.accept("--------------------------------------------------------------------------------");
            sink.accept("  C-MATRIX: " + title);
            sink.accept("--------------------------------------------------------------------------------");

            if (cmat == null) {
                sink.accept("  (No C-Matrix data present)");
                return;
            }

            for (int t = 0; t < cmat.size(); t++) {
                for (int j = 0; j < cmat.get(t).size(); j++) {
                    double[][] block = cmat.get(t).get(j);
                    sink.accept(String.format("  - Block (t=%-2d, j=%-2d): %d rows x %d cols",
                        t, j, block.length, block[0].length));
                    for (int r = 0; r < block.length; r++) {
                        StringBuilder row = new StringBuilder("    [ ");
                        for (int c = 0; c < block[r].length; c++) {
                            row.append(String.format("%7.3f ", block[r][c]));
                        }
                        row.append("]");
                        sink.accept(row.toString());
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Builder logic (formerly CMatrixBuilder)
    // -------------------------------------------------------------------------

    public static Result build(
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            List<Cluster> maxClusters,
            int numElements,
            CvCfBasis basis) {
        return build(clusterResult, cfResult, maxClusters, numElements, basis, null);
    }

    public static Result build(
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            List<Cluster> maxClusters,
            int numElements,
            CvCfBasis basis,
            Consumer<String> progressSink) {

        Result orthResult = buildOrthogonal(clusterResult, cfResult, maxClusters, numElements);

        if (progressSink != null) {
            emit(progressSink, "  [Orthogonal Basis C-Matrix Data]");
            dumpCmat("ORTHO", orthResult, progressSink);
        }

        Result cvcfResult = CvCfBasisTransformer.transform(orthResult, basis);

        if (progressSink != null) {
            emit(progressSink, "  [CVCF Basis C-Matrix Data]");
            dumpCmat("CVCF", cvcfResult, progressSink);
        }

        return cvcfResult;
    }

    /**
     * Internal logic for building the raw orthogonal C-matrix.
     */
    public static Result buildOrthogonal(
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            List<Cluster> maxClusters,
            int numElements) {

        if (clusterResult == null || cfResult == null || maxClusters == null) {
            throw new IllegalArgumentException("Inputs must not be null");
        }

        LOG.fine("CMatrix.buildOrthogonal — ENTER");

        List<Position> siteList = ClusterBuilders.buildSiteList(maxClusters);
        ClusterMath.PRules pRules = ClusterMath.PRules.build(siteList.size(), numElements);

        List<List<List<List<SiteOp>>>> cfSiteOpList =
                ClusterBuilders.buildCFSiteOpList(cfResult.getGroupedCFData(), siteList);
        SubstituteRules substituteRules = SubstituteRules.build(cfSiteOpList, siteList);

        int totalCfs = cfResult.getTcf();
        int[][] lcf = cfResult.getLcf();
        Map<CFIndex, Integer> cfColumn = buildCfColumnMap(lcf);

        int[][] cfBasisIndices = extractCfBasisIndices(cfSiteOpList, cfColumn, totalCfs);
        List<List<Cluster>> ordClustersByType = clusterResult.getOrdClusterData().getCoordList();

        List<List<double[][]>> cmat = new ArrayList<>();
        List<List<int[]>> wcv = new ArrayList<>();
        int[][] lcv = new int[ordClustersByType.size()][];

        for (int t = 0; t < ordClustersByType.size(); t++) {
            List<Cluster> groups = ordClustersByType.get(t);
            List<double[][]> cmatType = new ArrayList<>();
            List<int[]> wcvType = new ArrayList<>();
            lcv[t] = new int[groups.size()];

            for (int j = 0; j < groups.size(); j++) {
                Cluster cluster = groups.get(j);
                List<Integer> siteIndices = flattenSiteIndices(cluster, siteList);

                List<int[]> configs = generateConfigs(siteIndices.size(), numElements);
                Map<PolynomialKey, Integer> countMap = new LinkedHashMap<>();
                Map<PolynomialKey, double[]> rowMap = new LinkedHashMap<>();

                for (int[] config : configs) {
                    double[] row = computeCvRow(siteIndices, config, pRules, substituteRules,
                            cfColumn, totalCfs);
                    PolynomialKey key = new PolynomialKey(row);
                    countMap.put(key, countMap.getOrDefault(key, 0) + 1);
                    rowMap.putIfAbsent(key, row);
                }

                int ncv = countMap.size();
                lcv[t][j] = ncv;

                double[][] cmatGroup = new double[ncv][totalCfs + 1];
                int[] wcvGroup = new int[ncv];

                int idx = 0;
                for (Map.Entry<PolynomialKey, Integer> entry : countMap.entrySet()) {
                    cmatGroup[idx] = rowMap.get(entry.getKey());
                    wcvGroup[idx] = entry.getValue();
                    idx++;
                }
                cmatType.add(cmatGroup);
                wcvType.add(wcvGroup);
            }
            cmat.add(cmatType);
            wcv.add(wcvType);
        }

        return new Result(cmat, lcv, wcv, cfBasisIndices, null);
    }

    private static void dumpCmat(String label, Result result, Consumer<String> sink) {
        List<List<double[][]>> cmat = result.getCmat();
        for (int t = 0; t < cmat.size(); t++) {
            List<double[][]> groups = cmat.get(t);
            for (int j = 0; j < groups.size(); j++) {
                double[][] block = groups.get(j);
                int rows = block.length;
                int cols = rows > 0 ? block[0].length : 0;
                emit(sink, String.format("    %s cmat[t=%d][j=%d] dims=%dx%d", label, t, j, rows, cols));
            }
        }
    }

    private static void emit(Consumer<String> sink, String line) {
        sink.accept(line);
    }

    private static int[][] extractCfBasisIndices(
            List<List<List<List<SiteOp>>>> cfSiteOpList,
            Map<CFIndex, Integer> cfColumn,
            int totalCfs) {
        int[][] result = new int[totalCfs][];
        for (int t = 0; t < cfSiteOpList.size(); t++) {
            List<List<List<SiteOp>>> typeGroups = cfSiteOpList.get(t);
            for (int j = 0; j < typeGroups.size(); j++) {
                List<List<SiteOp>> group = typeGroups.get(j);
                for (int k = 0; k < group.size(); k++) {
                    int col = cfColumn.get(new CFIndex(t, j, k));
                    List<SiteOp> ops = group.get(k);
                    result[col] = new int[ops.size()];
                    for (int s = 0; s < ops.size(); s++) {
                        result[col][s] = ops.get(s).getBasisIndex();
                    }
                }
            }
        }
        return result;
    }

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

    private static List<Integer> flattenSiteIndices(Cluster cluster, List<Position> siteList) {
        List<Integer> indices = new ArrayList<>();
        for (Sublattice sub : cluster.getSublattices()) {
            for (Site site : sub.getSites()) {
                int idx = indexOf(siteList, site.getPosition());
                if (idx < 0) throw new IllegalStateException("Site position not found");
                indices.add(idx);
            }
        }
        return indices;
    }

    private static int indexOf(List<Position> siteList, Position pos) {
        for (int i = 0; i < siteList.size(); i++) {
            if (siteList.get(i).equals(pos)) return i;
        }
        return -1;
    }

    private static List<int[]> generateConfigs(int numSites, int numElements) {
        List<int[]> configs = new ArrayList<>();
        int total = (int) Math.pow(numElements, numSites);
        for (int i = 0; i < total; i++) {
            int[] cfg = new int[numSites];
            int x = i;
            for (int s = numSites - 1; s >= 0; s--) {
                cfg[s] = x % numElements;
                x /= numElements;
            }
            configs.add(cfg);
        }
        return configs;
    }

    private static double[] computeCvRow(
            List<Integer> siteIndices,
            int[] config,
            ClusterMath.PRules pRules,
            SubstituteRules substituteRules,
            Map<CFIndex, Integer> cfColumn,
            int totalCfs) {

        Map<SiteOpProductKey, Double> poly = new LinkedHashMap<>();
        poly.put(new SiteOpProductKey(List.of()), 1.0);

        for (int s = 0; s < siteIndices.size(); s++) {
            int siteIndex = siteIndices.get(s);
            int elementIndex = config[s];
            double[] coeffs = pRules.coefficientsFor(siteIndex, elementIndex);

            Map<SiteOpProductKey, Double> next = new LinkedHashMap<>();
            for (Map.Entry<SiteOpProductKey, Double> entry : poly.entrySet()) {
                List<SiteOp> baseOps = entry.getKey().getOps();
                double baseCoeff = entry.getValue();
                for (int a = 0; a < coeffs.length; a++) {
                    double c = coeffs[a];
                    if (Math.abs(c) < 1e-12) continue;
                    List<SiteOp> newOps = new ArrayList<>(baseOps);
                    if (a > 0) newOps.add(new SiteOp(siteIndex, a));
                    SiteOpProductKey key = new SiteOpProductKey(newOps);
                    next.put(key, next.getOrDefault(key, 0.0) + baseCoeff * c);
                }
            }
            poly = next;
        }

        double[] row = new double[totalCfs + 1];
        for (Map.Entry<SiteOpProductKey, Double> entry : poly.entrySet()) {
            SiteOpProductKey key = entry.getKey();
            double coeff = entry.getValue();
            List<SiteOp> ops = key.getOps();
            if (ops.isEmpty()) {
                row[totalCfs] += coeff;
                continue;
            }
            CFIndex cfIndex = substituteRules.lookup(ops);
            if (cfIndex == null) throw new IllegalStateException("No CF mapping for " + key);
            Integer col = cfColumn.get(cfIndex);
            if (col == null) throw new IllegalStateException("No CF column for " + cfIndex);
            row[col] += coeff;
        }
        return row;
    }
}
