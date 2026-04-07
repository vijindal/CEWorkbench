package org.ce.domain.cluster.cvcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.logging.Logger;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.ClusterKeys.CFIndex;
import org.ce.domain.cluster.CMatrix;
import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.ClusterBuilders;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.ClusterMath;
import org.ce.domain.cluster.LinearAlgebra;
import org.ce.domain.cluster.ClusterPrimitives.Position;
import org.ce.domain.cluster.ClusterPrimitives.Site;
import org.ce.domain.cluster.ClusterKeys.SiteOp;
import org.ce.domain.cluster.ClusterKeys.SiteOpProductKey;
import org.ce.domain.cluster.SubstituteRules;

/**
 * Generates CVCF transformation matrices dynamically from physical definitions.
 * 
 * <p>Standardizes the basis transition by expressing CVCF functions (v) as 
 * products of site probabilities (p), and expanding them into the orthogonal 
 * basis (u) using identification rules.</p>
 */
public final class CvCfMatrixGenerator {

    private static final Logger LOG = Logger.getLogger(CvCfMatrixGenerator.class.getName());

    private CvCfMatrixGenerator() {}

    /**
     * Generates a CVCF basis for the T-model.
     */
    public static CvCfBasis generate(
            String structurePhase,
            ClusterIdentificationResult clusterResult,
            CFIdentificationResult cfResult,
            CMatrix.Result orthMatrix,
            String model) {
        if (!"T".equalsIgnoreCase(model)) {
            throw new UnsupportedOperationException("Dynamic generation only supported for T-model");
        }

        int numComponents = cfResult.getNxcf() + 1;
        
        // 1. Identify representative tetrahedron sites
        Cluster tetrahedron = findTetrahedron(clusterResult);
        List<Position> siteList = ClusterBuilders.buildSiteList(List.of(tetrahedron));
        Map<Integer, Integer> siteMap = resolveSiteMapping(tetrahedron, siteList);

        // 2. Build SubstituteRules and PRules for expansion
        ClusterMath.PRules pRules = ClusterMath.PRules.build(siteList.size(), numComponents);
        SubstituteRules substituteRules = rebuildSubstituteRules(cfResult, siteList);

        // 3. Define v-functions (Binary BCC_A2)
        List<String> cfNames = new ArrayList<>();
        List<LinearForm> vFunctions = new ArrayList<>();

        if (numComponents == 2) {
            cfNames.addAll(List.of("v4AB", "v3AB", "v22AB", "v21AB", "xA", "xB"));
            
            // v4AB = p[1][A]*p[2][B]*p[3][B]*p[4][A]
            vFunctions.add(product(siteMap, pRules, List.of(1, 0, 2, 1, 3, 1, 4, 0)));
            
            // v3AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
            LinearForm v3Term1 = product(siteMap, pRules, List.of(1, 0, 2, 1, 3, 1));
            LinearForm v3Term2 = product(siteMap, pRules, List.of(1, 1, 2, 0, 3, 0));
            vFunctions.add(v3Term1.subtract(v3Term2));

            vFunctions.add(product(siteMap, pRules, List.of(1, 0, 4, 1))); // v22AB
            vFunctions.add(product(siteMap, pRules, List.of(1, 0, 2, 1))); // v21AB
            vFunctions.add(probability(siteMap, pRules, 1, 0)); // xA
            vFunctions.add(probability(siteMap, pRules, 1, 1)); // xB
        } else {
            throw new UnsupportedOperationException("Dynamic generation for component count " + numComponents + " not yet implemented.");
        }

        // Print symbolic expressions for user verification
        printSymbolicExpressions(vFunctions, cfNames, cfResult, substituteRules);

        // 4. Extract M matrix (v = M * u)
        int totalCfs = cfResult.getTcf();
        int expectedSize = vFunctions.size();
        if (expectedSize != (totalCfs + 1)) {
             LOG.warning(String.format("Basis size mismatch: CVCF has %d, Orthogonal has %d. Inversion may fail.", 
                  expectedSize, totalCfs + 1));
        }

        double[][] M = new double[expectedSize][expectedSize];
        Map<CFIndex, Integer> cfColMap = buildCfColumnMap(cfResult.getLcf());

        for (int i = 0; i < expectedSize; i++) {
            LinearForm f = vFunctions.get(i);
            for (Map.Entry<SiteOpProductKey, Double> entry : f.terms.entrySet()) {
                double coeff = entry.getValue();
                List<SiteOp> ops = entry.getKey().getOps();
                if (ops.isEmpty()) {
                    M[i][expectedSize - 1] += coeff; // Dummy is last column
                } else {
                    CFIndex cfIdx = substituteRules.lookup(ops);
                    if (cfIdx != null) {
                        Integer col = cfColMap.get(cfIdx);
                        if (col != null && col < expectedSize - 1) {
                            M[i][col] += coeff;
                        }
                    }
                }
            }
        }

        // T satisfies u = T * v. My M is rows=v, cols=u.
        // So v = M * u  =>  u = M^-1 * v.
        // Thus T = M^-1.
        double[][] T = LinearAlgebra.invert(M);
        double[][] Tinv = M;

        return new CvCfBasis(structurePhase, model, numComponents,
                cfNames, List.of(), cfResult.getNcf(), T, Tinv);
    }

    private static void printSymbolicExpressions(
            List<LinearForm> vFunctions, 
            List<String> vNames, 
            CFIdentificationResult cfResult,
            SubstituteRules substituteRules) {
        
        System.out.println("\n  [DIAGNOSTIC] CVCF Basis Symbolic Expressions:");
        List<String> uNames = cfResult.getUNames();
        Map<CFIndex, Integer> cfColMap = buildCfColumnMap(cfResult.getLcf());
        int dummyIdx = uNames.size();

        for (int i = 0; i < vFunctions.size(); i++) {
            LinearForm f = vFunctions.get(i);
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("    %-8s = ", vNames.get(i)));
            
            boolean first = true;
            // Collect terms by u-index to group them
            Map<Integer, Double> groupedTerms = new java.util.TreeMap<>();
            double constantPart = 0.0;

            for (Map.Entry<SiteOpProductKey, Double> entry : f.terms.entrySet()) {
                List<SiteOp> ops = entry.getKey().getOps();
                if (ops.isEmpty()) {
                    constantPart += entry.getValue();
                } else {
                    CFIndex cfIdx = substituteRules.lookup(ops);
                    if (cfIdx != null) {
                        Integer col = cfColMap.get(cfIdx);
                        if (col != null) {
                            groupedTerms.merge(col, entry.getValue(), Double::sum);
                        }
                    }
                }
            }

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
                if (!first && constantPart > 0) sb.append(" + ");
                else if (constantPart < 0) sb.append(" - ");
                sb.append(String.format("%.4f", Math.abs(constantPart)));
            } else if (first) {
                sb.append("0.0");
            }
            
            System.out.println(sb.toString());
        }
    }

    private static Cluster findTetrahedron(ClusterIdentificationResult clusterResult) {
        for (List<Cluster> groups : clusterResult.getOrdClusterData().getCoordList()) {
            for (Cluster c : groups) {
                if (c.getAllSites().size() == 4) return c;
            }
        }
        throw new IllegalStateException("No tetrahedron cluster found in Stage 3 results.");
    }

    private static Map<Integer, Integer> resolveSiteMapping(Cluster cluster, List<Position> siteList) {
        List<Position> sites = new ArrayList<>();
        for (int i = 0; i < siteList.size(); i++) sites.add(siteList.get(i));

        int n = sites.size();
        double minPair = Double.MAX_VALUE;
        double maxPair = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double d = sites.get(i).distance(sites.get(j));
                minPair = Math.min(minPair, d);
                maxPair = Math.max(maxPair, d);
            }
        }

        int p1 = 0;
        int p2 = -1;
        for (int i = 0; i < n; i++) {
            if (i == p1) continue;
            if (Math.abs(sites.get(p1).distance(sites.get(i)) - minPair) < 1e-4) {
                p2 = i; break;
            }
        }
        int p4 = -1;
        for (int i = 0; i < n; i++) {
            if (i == p1) continue;
            if (Math.abs(sites.get(p1).distance(sites.get(i)) - maxPair) < 1e-4) {
                p4 = i; break;
            }
        }
        int p3 = -1;
        for (int i = 0; i < n; i++) {
            if (i == p1 || i == p2 || i == p4) continue;
            p3 = i; break;
        }

        Map<Integer, Integer> map = new LinkedHashMap<>();
        map.put(1, p1); map.put(2, p2); map.put(3, p3); map.put(4, p4);
        return map;
    }

    private static LinearForm probability(Map<Integer, Integer> siteMap, ClusterMath.PRules pRules, int logicalSite, int atom) {
        int physicalSite = siteMap.get(logicalSite);
        double[] coeffs = pRules.coefficientsFor(physicalSite, atom);
        LinearForm f = new LinearForm();
        f.add(new SiteOpProductKey(List.of()), coeffs[0]);
        for (int i = 1; i < coeffs.length; i++) {
            f.add(new SiteOpProductKey(List.of(new SiteOp(physicalSite, i))), coeffs[i]);
        }
        return f;
    }

    private static LinearForm product(Map<Integer, Integer> siteMap, ClusterMath.PRules pRules, List<Integer> siteAtomPairs) {
        LinearForm res = new LinearForm();
        res.add(new SiteOpProductKey(List.of()), 1.0);
        for (int i = 0; i < siteAtomPairs.size(); i += 2) {
            LinearForm p = probability(siteMap, pRules, siteAtomPairs.get(i), siteAtomPairs.get(i + 1));
            res = res.multiply(p);
        }
        return res;
    }

    private static SubstituteRules rebuildSubstituteRules(CFIdentificationResult cfResult, List<Position> siteList) {
        return SubstituteRules.build(ClusterBuilders.buildCFSiteOpList(
                cfResult.getGroupedCFData(), siteList), siteList);
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

    private static class LinearForm {
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
