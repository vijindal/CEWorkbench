package org.ce.model.cluster.cvcf;

import java.util.ArrayList;
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
        emit(sink, String.format("  [4.0] Generating CVCF basis: structure=%s, model=%s, K=%d",
                structurePhase, model, numComponents));

        // 1. Look up the physical v-function definition for this combination
        CvCfDefinition def = CvCfDefinition.get(structurePhase, model, numComponents);
        emit(sink, String.format("  [4.1] Definition loaded: %d CFs (%d non-point + %d point)",
                def.cfNames.size(), def.cfNames.size() - numComponents, numComponents));

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
        emit(sink, String.format(
                "  [4.3] PRules reused from Stage 3 (K=%d, %d sites); SubstituteRules reused (%d rules)",
                pRules.getNumElements(), pRules.getNumSites(), subRules.getRules().size()));
        // Debug: print PRules
        // emit(sink, "  [4.3] PRules: p[site][atom] coefficients (c0, c1, c2, ...):");
        // int K = pRules.getNumElements();
        // int S = pRules.getNumSites();
        // for (int site = 0; site < S; site++) {
        //     for (int atom = 0; atom < K; atom++) {
        //         double[] coeffs = pRules.coefficientsFor(site, atom);
        //         StringBuilder sb = new StringBuilder(String.format("    p[%d][%d] =", site, atom));
        //         for (int i = 0; i < coeffs.length; i++) {
        //             if (i == 0) sb.append(String.format(" %.4f", coeffs[i]));
        //             else sb.append(String.format(" + %.4f*s%d", coeffs[i], i));
        //         }
        //         emit(sink, sb.toString());
        //     }
        // }

        // 4. Expand each v-function into the orthogonal basis
        List<LinearForm> vFunctions = buildVFunctions(def, siteMap, pRules);
        emit(sink, String.format("  [4.4] V-functions expanded: %d LinearForms built",
                vFunctions.size()));
        // Debug: print expanded LinearForms (site-operator products after pRules substitution)
        // emit(sink, "  [4.4] Expanded CVs (site-operator products after pRules):");
        // for (int i = 0; i < vFunctions.size(); i++) {
        //     StringBuilder sb = new StringBuilder(String.format("    %-10s:", def.cfNames.get(i)));
        //     for (Map.Entry<SiteOpProductKey, Double> e : vFunctions.get(i).terms.entrySet()) {
        //         double coeff = e.getValue();
        //         List<SiteOp> ops = e.getKey().getOps();
        //         sb.append(String.format("  %.4f*%s", coeff, ops.isEmpty() ? "1" : ops.toString()));
        //     }
        //     emit(sink, sb.toString());
        // }

        // 5. Print symbolic expressions for verification
        printSymbolicExpressions(vFunctions, def.cfNames, cfResult, subRules, sink);

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
        double[][] M = buildMMatrix(vFunctions, subRules, cfColMap, basisSize);
        emit(sink, String.format("  [4.6] M matrix built: %dx%d", M.length, M[0].length));

        // 7. T = M^-1  (u = T * v)
        double[][] T    = LinearAlgebra.invert(M);
        double[][] Tinv = M;
        emit(sink, String.format("  [4.7] T matrix computed (M^-1): %dx%d", T.length, T[0].length));
        printMatrix("T (orthogonal->CVCF)", T, def.cfNames, cfResult.getUNames(), sink);
        printMatrix("Tinv (CVCF->orthogonal)", Tinv, cfResult.getUNames(), def.cfNames, sink);

        // 8. Derive numNonPointCfs and eciNames
        int numNonPointCfs = def.cfNames.size() - numComponents;
        List<String> eciNames = cfResult.getEONames().subList(0, numNonPointCfs);
        emit(sink, String.format("  [4.8] Basis complete: cfNames=%s, eciNames=%s",
                def.cfNames, eciNames));

        return new CvCfBasis(structurePhase, model, numComponents,
                def.cfNames, eciNames, numNonPointCfs, T, Tinv);
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
            emit(sink, String.format(
                    "        p%d {%.4f,%.4f,%.4f} -> physicalSite[%d] {%.4f,%.4f,%.4f} OK",
                    logIdx + 1, coord[0], coord[1], coord[2],
                    physIdx, matched.getX(), matched.getY(), matched.getZ()));
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
            int basisSize) {

        double[][] M = new double[basisSize][basisSize];
        for (int i = 0; i < basisSize; i++) {
            for (Map.Entry<SiteOpProductKey, Double> entry : vFunctions.get(i).terms.entrySet()) {
                double coeff = entry.getValue();
                List<SiteOp> ops = entry.getKey().getOps();
                if (ops.isEmpty()) {
                    M[i][basisSize - 1] += coeff;
                } else {
                    CFIndex cfIdx = subRules.lookup(ops);
                    if (cfIdx != null) {
                        Integer col = cfColMap.get(cfIdx);
                        if (col != null && col < basisSize - 1) {
                            M[i][col] += coeff;
                        }
                    }
                }
            }
        }
        return M;
    }

    // =========================================================================
    // Diagnostic: symbolic expressions
    // =========================================================================

    private static void printSymbolicExpressions(
            List<LinearForm> vFunctions,
            List<String> vNames,
            CFIdentificationResult cfResult,
            SubstituteRules subRules,
            Consumer<String> sink) {

        if (sink == null) return;
        sink.accept("\n  [DIAGNOSTIC] CVCF Basis Symbolic Expressions:");
        List<String> uNames = cfResult.getUNames();
        Map<CFIndex, Integer> cfColMap = buildCfColumnMap(cfResult.getLcf());

        for (int i = 0; i < vFunctions.size(); i++) {
            LinearForm f = vFunctions.get(i);
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("    %-8s = ", vNames.get(i)));

            Map<Integer, Double> groupedTerms = new java.util.TreeMap<>();
            double constantPart = 0.0;

            for (Map.Entry<SiteOpProductKey, Double> entry : f.terms.entrySet()) {
                List<SiteOp> ops = entry.getKey().getOps();
                if (ops.isEmpty()) {
                    constantPart += entry.getValue();
                } else {
                    CFIndex cfIdx = subRules.lookup(ops);
                    if (cfIdx != null) {
                        Integer col = cfColMap.get(cfIdx);
                        if (col != null) groupedTerms.merge(col, entry.getValue(), Double::sum);
                    }
                }
            }

            boolean first = true;
            for (Map.Entry<Integer, Double> entry : groupedTerms.entrySet()) {
                double val = entry.getValue();
                if (Math.abs(val) < 1e-10) continue;
                if (!first && val > 0) sb.append(" + ");
                else if (val < 0) sb.append(" - ");
                double absVal = Math.abs(val);
                if (absVal != 1.0) sb.append(String.format("%.4f", absVal)).append("*");
                sb.append(uNames.get(entry.getKey()));
                first = false;
            }

            if (Math.abs(constantPart) > 1e-10) {
                if (!first && constantPart > 0) { sb.append(" + "); }
                else if (constantPart < 0) { sb.append(" - "); }
                sb.append(String.format("%.4f", Math.abs(constantPart)));
            } else if (first) {
                sb.append("0.0");
            }

            sink.accept(sb.toString());
        }
    }

    // =========================================================================
    // Utilities
    // =========================================================================

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }

    private static void printMatrix(String title, double[][] mat,
                                    List<String> rowLabels, List<String> colLabels,
                                    Consumer<String> sink) {
        if (sink == null) return;
        sink.accept("  ---- " + title + " ----");
        int cols = mat[0].length;
        StringBuilder hdr = new StringBuilder("         ");
        for (int j = 0; j < cols; j++) {
            hdr.append(String.format(" %9s", j < colLabels.size() ? colLabels.get(j) : ("c" + j)));
        }
        sink.accept(hdr.toString());
        for (int i = 0; i < mat.length; i++) {
            String label = i < rowLabels.size() ? rowLabels.get(i) : ("r" + i);
            StringBuilder row = new StringBuilder(String.format("  %-8s", label));
            for (int j = 0; j < cols; j++) {
                row.append(String.format(" %9.4f", mat[i][j]));
            }
            sink.accept(row.toString());
        }
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
