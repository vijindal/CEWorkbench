package org.ce.model.cluster.cvcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.logging.Logger;

import org.ce.model.cluster.CFIdentificationResult;
import org.ce.model.cluster.ClusterKeys.CFIndex;
import org.ce.model.cluster.CMatrix;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterIdentificationResult;
import org.ce.model.cluster.ClusterMath;
import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cluster.ClusterPrimitives.Position;
import org.ce.model.cluster.ClusterPrimitives.Sublattice;
import org.ce.model.cluster.ClusterPrimitives.Site;
import org.ce.model.cluster.ClusterKeys.SiteOp;
import org.ce.model.cluster.ClusterKeys.SiteOpProductKey;
import org.ce.model.cluster.SubstituteRules;

/**
 * Generates CVCF transformation matrices dynamically from physical definitions.
 *
 * Implements the Mathematica pipeline:
 *   CVS  = Expand[cvs /. pRules] /. substituteRules
 *   M    = CoefficientArrays[CVS, uList]   (v = M * u)
 *   T    = M^-1                             (u = T * v)
 *
 * The physical v-function definitions (which sites and atoms define each CV)
 * are looked up from CvCfDefinition -- nothing structure-specific is hardcoded
 * here. Add new structures/models/components by registering them in CvCfDefinition.
 *
 * pRules and substituteRules are reused from the CMatrix.Result produced by
 * Stage 3 (CMatrix.buildOrthogonal) rather than being rebuilt from scratch.
 */
public final class CvCfMatrixGenerator {

    private static final Logger LOG = Logger.getLogger(CvCfMatrixGenerator.class.getName());

    private CvCfMatrixGenerator() {}

    // =========================================================================
    // Public entry point
    // =========================================================================

    /**
     * Generates a CVCF basis for the given (structurePhase, model, numComponents)
     * combination.
     *
     * @param structurePhase  e.g. "BCC_A2"
     * @param clusterResult   Stage 2a output
     * @param cfResult        Stage 2b output
     * @param orthMatrix      Stage 3 output -- provides pRules, substituteRules, siteList
     * @param model           e.g. "T"
     * @param sink            optional progress sink (may be null)
     * @return a fully constructed CvCfBasis with T, Tinv, cfNames, eciNames
     */
    public static CvCfBasis generate(
            String structurePhase,
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            CMatrix.Result orthMatrix,
            String model,
            Consumer<String> sink) {

        if (!"T".equalsIgnoreCase(model)) {
            throw new UnsupportedOperationException(
                "Dynamic generation only supported for T-model; got: " + model);
        }

        int numComponents = cfResult.getNxcf() + 1;
        // 1. Look up the physical v-function definition for this combination
        CvCfDefinition def = CvCfDefinition.get(structurePhase, model, numComponents);

        // 2. Resolve logical site indices -> physical site indices.
        // Use the maximal ordered-phase cluster (type 0, representative 0) to obtain
        // site positions and their physical indices. This guarantees that the
        // site-operator products we generate are exactly the same ones registered in
        // cfSiteOpList by Stage 2b, so SubstituteRules.lookup() never returns null.
        List<Position> siteList = orthMatrix.getSiteList();
        Cluster maxCluster = clusterResult.getOrdClusterData().getCoordList().get(0).get(0);
        Map<Integer, Integer> siteMap = resolveSiteMap(def.logicalSiteCoords, maxCluster, siteList, sink);

        // 3. Reuse pRules and substituteRules from Stage 3
        ClusterMath.PRules pRules   = orthMatrix.getPRules();
        SubstituteRules    subRules = orthMatrix.getSubstituteRules();

        // [STEP P] pRules
        emit(sink, "\n[STEP P] pRules (R-Matrix coefficients):");
        double[][] rMat = pRules.getRMatrix();
        for (int i = 0; i < rMat.length; i++) {
            emit(sink, String.format("  atom %d -> %s", i, Arrays.toString(rMat[i])));
        }

        // [STEP S] substitutionRules
        emit(sink, "\n[STEP S] substitutionRules:");
        for (SiteOpProductKey key : subRules.getRules().keySet()) {
            CFIndex cfIdx = subRules.lookup(key.getOps());
            emit(sink, String.format(
                "  ops=%s  -> CFIndex=%s",
                key.getOps(), cfIdx
            ));
        }

        // [STEP S2] substitutionRules consistency check
        emit(sink, "\n[STEP S2] substitutionRules consistency check:");
        for (SiteOpProductKey key : subRules.getRules().keySet()) {
            CFIndex cfIdx = subRules.lookup(key.getOps());
            if (cfIdx == null) {
                emit(sink, "  [ERROR] Missing mapping for ops=" + key.getOps());
            }
        }

        // 4. Expand each v-function into the orthogonal basis
        List<LinearForm> vFunctions = buildVFunctions(def, siteMap, pRules);
        /*
        emit(sink, "\n[DIAG] V-FUNCTION EXPANSION:");
        for (int i = 0; i < vFunctions.size(); i++) {
            LinearForm f = vFunctions.get(i);
            emit(sink, "\n  v[" + i + "]");
            for (Map.Entry<SiteOpProductKey, Double> e : f.terms.entrySet()) {
                emit(sink, String.format(
                    "    coeff=% .6f  ops=%s",
                    e.getValue(),
                    e.getKey().getOps()
                ));
            }
        }
        */

        // [STEP 1] RAW CV DEFINITIONS (PROBABILITY FORM)
        emit(sink, "\n[STEP 1] RAW CV DEFINITIONS (PROBABILITY FORM):");
        for (int i = 0; i < def.vSpecs.size(); i++) {
            CvCfDefinition.VSpec spec = def.vSpecs.get(i);
            String name = def.cfNames.get(i);
            emit(sink, String.format("  v[%d] (%s) = %s", i, name, formatVSpec(spec)));
        }

        // [STEP 2] AFTER pRules (SITE OPERATOR EXPANSION)
        emit(sink, "\n[STEP 2] AFTER pRules (SITE OPERATOR EXPANSION):");
        for (int i = 0; i < vFunctions.size(); i++) {
            LinearForm f = vFunctions.get(i);
            emit(sink, "\n  v[" + i + "] expanded:");
            for (Map.Entry<SiteOpProductKey, Double> e : f.terms.entrySet()) {
                emit(sink, String.format(
                    "    coeff=% .6f  ops=%s",
                    e.getValue(),
                    e.getKey().getOps()
                ));
            }
        }

        // 6. Build M matrix: v[i] = sum_j M[i][j] * u[j]
        int totalCfs  = cfResult.getTcf();
        int basisSize = vFunctions.size();
        if (basisSize != totalCfs + 1) {
            LOG.warning(String.format(
                "Basis size mismatch: CVCF has %d, orthogonal has %d (tcf+1). Inversion may fail.",
                basisSize, totalCfs + 1));
            emit(sink, String.format("  [WARN] Basis size mismatch: CVCF=%d, orthogonal tcf+1=%d",
                    basisSize, totalCfs + 1));
        }
        Map<CFIndex, Integer> cfColMap = buildCfColumnMap(cfResult.getLcf());

        emit(sink, "\n[DEBUG] CFIndex â†’ column mapping:");
        for (Map.Entry<CFIndex, Integer> e : cfColMap.entrySet()) {
            emit(sink, String.format("  %s -> col %d", e.getKey(), e.getValue()));
        }

        emit(sink, "\n[DEBUG] Column meanings:");
        List<String> uNames = cfResult.getUNames();
        for (int i = 0; i < uNames.size(); i++) {
            emit(sink, String.format("  col %d -> %s", i, uNames.get(i)));
        }
        emit(sink, String.format("  col %d -> (Empty Cluster / 1.0)", basisSize - 1));

        emit(sink, "\n[DEBUG] Point CFIndex mapping check:");
        for (Map.Entry<CFIndex, Integer> e : cfColMap.entrySet()) {
            CFIndex idx = e.getKey();
            // In BCC_A2 ternary, Point CFs are cluster type 1
            if (idx.getTypeIndex() == 1) {
                emit(sink, String.format("  POINT CF: %s -> col %d", idx, e.getValue()));
            }
        }

        double[][] M = buildMMatrix(vFunctions, subRules, cfColMap, basisSize, numComponents, sink);

        /*
        emit(sink, "\n[DIAG] M MATRIX:");
        for (int i = 0; i < M.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("  row %d: ", i));
            for (int j = 0; j < M[i].length; j++) {
                sb.append(String.format("% .6f ", M[i][j]));
            }
            emit(sink, sb.toString());
        }
        */

        // 7. T = M^-1  (u = T * v)
        double[][] T    = LinearAlgebra.invert(M);

        /*
        emit(sink, "\n[DIAG] Checking T * M:");
        for (int i = 0; i < T.length; i++) {
            for (int j = 0; j < T.length; j++) {
                double sum = 0;
                for (int k = 0; k < T.length; k++) {
                    sum += T[i][k] * M[k][j];
                }
                emit(sink, String.format("I[%d][%d] = %.6f", i, j, sum));
            }
        }
        */

        double[][] Tinv = M;

        // 8. Derive numNonPointCfs and eciNames
        int numNonPointCfs = def.cfNames.size() - numComponents;
        List<String> eciNames = cfResult.getEONames().subList(0, numNonPointCfs);

        return new CvCfBasis(structurePhase, model, numComponents,
                def.cfNames, eciNames, numNonPointCfs, T, Tinv);
    }

    private static String formatVSpec(CvCfDefinition.VSpec spec) {
        String plus = formatTerm(spec.plusTerm);
        if (spec.isDiff()) {
            return plus + " - " + formatTerm(spec.minusTerm);
        }
        return plus;
    }

    private static String formatTerm(int[] pairs) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < pairs.length; i += 2) {
            if (i > 0) sb.append("*");
            sb.append("p[").append(pairs[i]).append("][").append(pairs[i+1]).append("]");
        }
        return sb.toString();
    }

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }

    // =========================================================================
    // Site coordinate matching
    // =========================================================================

    private static Map<Integer, Integer> resolveSiteMap(
            double[][] logicalSiteCoords,
            Cluster maxCluster,
            List<Position> siteList,
            Consumer<String> sink) {

        // Collect the positions of sites in the maximal cluster (in cluster order).
        // These are exactly the sites whose products appear in cfSiteOpList.
        List<Position> clusterPositions = new ArrayList<>();
        for (Sublattice sub : maxCluster.getSublattices()) {
            for (Site site : sub.getSites()) {
                clusterPositions.add(site.getPosition());
            }
        }

        Map<Integer, Integer> siteMap = new LinkedHashMap<>();
        emit(sink, "  [4.2] Site coordinate matching (against maximal cluster sites):");
        for (int logIdx = 0; logIdx < logicalSiteCoords.length; logIdx++) {
            double[] coord = logicalSiteCoords[logIdx];
            // Match against cluster positions, then look up the global siteList index.
            Position matched = findMatchingPosition(coord, clusterPositions);
            int physIdx = indexOf(matched, siteList);
            if (physIdx < 0) {
                throw new IllegalStateException(String.format(
                    "Matched cluster position {%.4f,%.4f,%.4f} not found in siteList.",
                    matched.getX(), matched.getY(), matched.getZ()));
            }
            siteMap.put(logIdx + 1, physIdx);
        }
        return siteMap;
    }

    private static Position findMatchingPosition(double[] coord, List<Position> positions) {
        for (Position p : positions) {
            double dx = p.getX() - coord[0];
            double dy = p.getY() - coord[1];
            double dz = p.getZ() - coord[2];
            if (Math.sqrt(dx*dx + dy*dy + dz*dz) < 1e-4) {
                return p;
            }
        }
        throw new IllegalStateException(String.format(
            "No site in maximal cluster matches logical site coordinate {%.4f, %.4f, %.4f}. "
            + "Check that the definition coordinates match the cluster input file.",
            coord[0], coord[1], coord[2]));
    }

    private static int indexOf(Position pos, List<Position> siteList) {
        for (int i = 0; i < siteList.size(); i++) {
            if (siteList.get(i).equals(pos)) return i;
        }
        return -1;
    }

    // =========================================================================
    // V-function expansion
    // =========================================================================

    private static List<LinearForm> buildVFunctions(
            CvCfDefinition def,
            Map<Integer, Integer> siteMap,
            ClusterMath.PRules pRules) {

        List<LinearForm> result = new ArrayList<>();
        for (CvCfDefinition.VSpec spec : def.vSpecs) {
            LinearForm plus = expandProduct(spec.plusTerm, siteMap, pRules);
            if (spec.isDiff()) {
                LinearForm minus = expandProduct(spec.minusTerm, siteMap, pRules);
                result.add(plus.subtract(minus));
            } else {
                result.add(plus);
            }
        }
        return result;
    }

    private static LinearForm expandProduct(
            int[] siteAtomPairs,
            Map<Integer, Integer> siteMap,
            ClusterMath.PRules pRules) {

        LinearForm res = new LinearForm();
        res.add(new SiteOpProductKey(List.of()), 1.0);
        for (int i = 0; i < siteAtomPairs.length; i += 2) {
            LinearForm p = probability(siteMap, pRules, siteAtomPairs[i], siteAtomPairs[i + 1]);
            res = res.multiply(p);
        }
        return res;
    }

    private static LinearForm probability(
            Map<Integer, Integer> siteMap,
            ClusterMath.PRules pRules,
            int logicalSite,
            int atom) {

        int physicalSite = siteMap.get(logicalSite);
        double[] coeffs = pRules.coefficientsFor(physicalSite, atom);
        LinearForm f = new LinearForm();
        f.add(new SiteOpProductKey(List.of()), coeffs[0]);
        for (int i = 1; i < coeffs.length; i++) {
            f.add(new SiteOpProductKey(List.of(new SiteOp(physicalSite, i))), coeffs[i]);
        }
        return f;
    }

    // =========================================================================
    // M matrix assembly
    // =========================================================================

    private static double[][] buildMMatrix(
            List<LinearForm> vFunctions,
            SubstituteRules subRules,
            Map<CFIndex, Integer> cfColMap,
            int basisSize,
            int numComponents,
            Consumer<String> sink) {

        double[][] M = new double[basisSize][basisSize];
        emit(sink, "\n[STEP 3] SUBSTITUTE RULES (OPS → CF INDEX):");

        for (int i = 0; i < basisSize; i++) {
            for (Map.Entry<SiteOpProductKey, Double> entry : vFunctions.get(i).terms.entrySet()) {
                double coeff = entry.getValue();
                List<SiteOp> ops = entry.getKey().getOps();
                if (ops.isEmpty()) {
                    M[i][basisSize - 1] += coeff;
                } else {
                    CFIndex cfIdx = subRules.lookup(ops);

                    /*
                    if (cfIdx == null) {
                        emit(sink, "  [ERROR] Missing CF mapping for ops: " + ops);
                    } else {
                        Integer col = cfColMap.get(cfIdx);
                        emit(sink, String.format(
                            "  [MAP] ops=%s -> cfIdx=%s -> col=%d",
                            ops, cfIdx, col
                        ));
                        if (col != null) {
                            M[i][col] += coeff;
                        }
                    }
                    */

                    emit(sink, String.format("  ops=%s  -> cfIdx=%s", ops, cfIdx));

                    if (cfIdx != null) {
                        Integer col = cfColMap.get(cfIdx);
                        emit(sink, String.format("[DEBUG] Using CFIndex %s -> col %d for ops %s", cfIdx, col, ops));
                        if (col != null) {
                            M[i][col] += coeff;
                        }
                    }
                }
            }

            // [STEP 4] FINAL EQUATION ROW
            emit(sink, "\n[STEP 4] FINAL EQUATION ROW " + i);
            for (int j = 0; j < M[i].length; j++) {
                if (Math.abs(M[i][j]) > 1e-10) {
                    emit(sink, String.format(
                        "    col %d  coeff=% .6f",
                        j, M[i][j]
                    ));
                }
            }
        }
        return M;
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

    // =========================================================================
    // LinearForm - symbolic polynomial in site-operator basis
    // =========================================================================

    static final class LinearForm {
        final Map<SiteOpProductKey, Double> terms = new LinkedHashMap<>();

        void add(SiteOpProductKey key, double val) {
            if (Math.abs(val) < 1e-15) return;
            terms.put(key, terms.getOrDefault(key, 0.0) + val);
        }

        LinearForm subtract(LinearForm other) {
            LinearForm res = new LinearForm();
            res.terms.putAll(this.terms);
            for (Map.Entry<SiteOpProductKey, Double> e : other.terms.entrySet()) {
                res.add(e.getKey(), -e.getValue());
            }
            return res;
        }

        LinearForm multiply(LinearForm other) {
            LinearForm res = new LinearForm();
            for (Map.Entry<SiteOpProductKey, Double> e1 : this.terms.entrySet()) {
                for (Map.Entry<SiteOpProductKey, Double> e2 : other.terms.entrySet()) {
                    List<SiteOp> ops = new ArrayList<>(e1.getKey().getOps());
                    ops.addAll(e2.getKey().getOps());
                    res.add(new SiteOpProductKey(ops), e1.getValue() * e2.getValue());
                }
            }
            return res;
        }
    }
}
