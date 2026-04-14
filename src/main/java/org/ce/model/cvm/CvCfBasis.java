package org.ce.model.cvm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.logging.Logger;

import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterMath;
import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cluster.CFIdentificationResult;
import org.ce.model.cluster.ClusterIdentificationResult;
import org.ce.model.cluster.ClusterKeys.CFIndex;
import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterPrimitives.Position;
import org.ce.model.cluster.ClusterPrimitives.Sublattice;
import org.ce.model.cluster.ClusterPrimitives.Site;
import org.ce.model.cluster.ClusterKeys.SiteOp;
import org.ce.model.cluster.ClusterKeys.SiteOpProductKey;
import org.ce.model.cluster.CMatrixPipeline.SubstituteRules;
import java.io.Serializable;

public final class CvCfBasis {
    private static final Logger LOG = Logger.getLogger(CvCfBasis.class.getName());

    public final String structurePhase;
    public final String model;
    public final int numComponents;
    public final List<String> cfNames;
    public final int numNonPointCfs;
    public final List<String> eciNames;
    public final double[][] T;
    public final double[][] Tinv;
    public final CMatrixPipeline.CMatrixData cvcfCMatrixData;
    private final Map<String, Integer> cfNameIndex;

    private CvCfBasis(String structurePhase, String model, int numComponents,
                      List<String> cfNames, List<String> eciNames, int numNonPointCfs,
                      double[][] T, double[][] Tinv, CMatrixPipeline.CMatrixData cvcfCMatrixData) {
        this.structurePhase  = structurePhase;
        this.model           = (model == null ? "" : model.toUpperCase());
        this.numComponents   = numComponents;
        this.cfNames         = cfNames;
        this.eciNames        = eciNames;
        this.numNonPointCfs  = numNonPointCfs;
        this.T               = T;
        this.Tinv            = Tinv;
        this.cvcfCMatrixData = cvcfCMatrixData;
        Map<String, Integer> idx = new LinkedHashMap<>(cfNames.size() * 2);
        for (int i = 0; i < cfNames.size(); i++) idx.put(cfNames.get(i), i);
        this.cfNameIndex = idx;
    }

    public static boolean isSupported(String structurePhase, String model, int numComponents) {
        return REGISTRY.containsKey(structurePhase + "_" + model.toUpperCase() + "_" + numComponents);
    }

    public static int getNumNonPointCfs(String structurePhase, String model, int numComponents) {
        Definition def = REGISTRY.get(structurePhase + "_" + model.toUpperCase() + "_" + numComponents);
        if (def == null) throw new IllegalArgumentException("Unregistered CVCF combination.");
        return def.cfNames.size() - numComponents;
    }

    public static String supportedSummary() {
        return "Supported combinations: " + String.join(", ", REGISTRY.keySet());
    }

    public int indexOfCf(String name) { return cfNameIndex.getOrDefault(name, -1); }
    public int totalCfs() { return cfNames.size(); }

    public double[] computeRandomState(double[] moleFractions, int[][] orthCfBasisIndices) {
        if (orthCfBasisIndices == null || orthCfBasisIndices.length == 0) {
            throw new IllegalArgumentException("orthCfBasisIndices must not be null or empty.");
        }
        int K    = moleFractions.length;
        int nxcf = K - 1;
        int orthTcf = orthCfBasisIndices.length;
        int orthNcf = orthTcf - nxcf;

        double[] basisVec = ClusterMath.buildBasis(K);
        double[] pointCF  = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            for (int i = 0; i < K; i++) {
                pointCF[k] += moleFractions[i] * Math.pow(basisVec[i], k + 1);
            }
        }

        double[] uNonPoint = new double[orthNcf];
        for (int col = 0; col < orthNcf; col++) {
            double val = 1.0;
            for (int b : orthCfBasisIndices[col]) val *= pointCF[b - 1];
            uNonPoint[col] = val;
        }

        double[] uPoint = new double[nxcf];
        for (int k = 0; k < nxcf; k++) {
            uPoint[k] = pointCF[k];
        }

        double[][] tInv = resolvedTinv();
        int tRows = tInv[0].length;

        double[] uOrthFull = new double[tRows];
        System.arraycopy(uNonPoint, 0, uOrthFull, 0, orthNcf);
        System.arraycopy(uPoint,    0, uOrthFull, orthNcf, nxcf);
        if (tRows == orthTcf + 1) {
            uOrthFull[tRows - 1] = 1.0;
        }

        int ncf  = numNonPointCfs;
        int tcf  = totalCfs();
        double[] vFull = new double[tcf];
        for (int j = 0; j < tcf; j++) {
            double sum = 0.0;
            double[] row = tInv[j];
            for (int i = 0; i < row.length; i++) sum += row[i] * uOrthFull[i];
            vFull[j] = sum;
        }

        for (int i = 0; i < K; i++) vFull[ncf + i] = moleFractions[i];

        return vFull;
    }

    private double[][] resolvedTinv() {
        if (Tinv != null) return Tinv;
        return LinearAlgebra.invert(T);
    }

    /**
     * Computes the full random CF vector in the CVCF basis.
     * v_random = T^-1 * u_random
     */
    public double[] computeRandomCvcfCFs(double[] moleFractions, 
                                        ClusterCFIdentificationPipeline.PipelineResult pr) {
        return computeRandomCvcfCFs(moleFractions, pr, resolvedTinv());
    }

    private static double[] computeRandomCvcfCFs(double[] moleFractions, 
                                               ClusterCFIdentificationPipeline.PipelineResult pr,
                                               double[][] mInv) {
        double[] uOrth = pr.computeRandomCFs(moleFractions);
        int dim = mInv.length;
        double[] vFull = new double[dim];
        for (int i = 0; i < dim; i++) {
            double sum = 0;
            for (int k = 0; k < dim; k++) {
                sum += mInv[i][k] * uOrth[k];
            }
            vFull[i] = sum;
        }
        return vFull;
    }

    /**
     * Diagnostic verification of the C-matrix in this basis at the disordered state.
     */
    public void verifyCvcfCMatrix(double[] moleFractions, 
                                 ClusterCFIdentificationPipeline.PipelineResult pr, 
                                 Consumer<String> sink) {
        
        emit(sink, "\n========================================================");
        emit(sink, "  CVCF VERIFICATION: Random State at Composition");
        emit(sink, "========================================================");
        emit(sink, "  Composition: " + Arrays.toString(moleFractions));

        // 1. Get random CFs in CVCF basis
        double[] vFull = computeRandomCvcfCFs(moleFractions, pr);

        emit(sink, "\n  RANDOM CVCF CFs:");
        for (int i = 0; i < vFull.length; i++) {
            emit(sink, String.format("    %-10s = %12.8f", cfNames.get(i), vFull[i]));
        }

        // 2. Evaluate CVs using CMatrixPipeline helper
        double[][][] cv = CMatrixPipeline.evaluateCVs(
                vFull, 
                this.cvcfCMatrixData.cmat, 
                this.cvcfCMatrixData.lcv,
                pr.getTcdis(), 
                pr.getLc()
        );

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
            Cluster representative = pr.getDisClusData().getClusCoordList().get(t);
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
        emit(sink, "\n=== CVCF Verification: COMPLETE ===");
    }

    public static final class VSpec {

        /** First (or only) product term: flat [site1, atom1, site2, atom2, ...]. */
        public final int[] plusTerm;

        /**
         * Subtracted product term, or {@code null} for a pure product.
         * If non-null, the v-function = product(plusTerm) - product(minusTerm).
         */
        public final int[] minusTerm;

        private VSpec(int[] plusTerm, int[] minusTerm) {
            this.plusTerm  = plusTerm;
            this.minusTerm = minusTerm;
        }

        /** v = p[s1][a1] * p[s2][a2] * ... */
        public static VSpec product(int... siteAtomPairs) {
            if (siteAtomPairs.length % 2 != 0)
                throw new IllegalArgumentException("siteAtomPairs must be even-length");
            return new VSpec(siteAtomPairs, null);
        }

        /** v = product(plusPairs) - product(minusPairs). */
        public static VSpec diff(int[] plusPairs, int[] minusPairs) {
            if (plusPairs.length % 2 != 0 || minusPairs.length % 2 != 0)
                throw new IllegalArgumentException("siteAtomPairs must be even-length");
            return new VSpec(plusPairs, minusPairs);
        }

        /** Convenience: single-site probability, i.e. point CF. v = p[site][atom]. */
        public static VSpec point(int logicalSite, int atom) {
            return product(logicalSite, atom);
        }

        public boolean isDiff() { return minusTerm != null; }
    }

    private static final class Definition {
        final double[][] logicalSiteCoords;
        final List<String> cfNames;
        final List<VSpec> vSpecs;

        Definition(double[][] logicalSiteCoords, List<String> cfNames, List<VSpec> vSpecs) {
            this.logicalSiteCoords = logicalSiteCoords;
            this.cfNames = Collections.unmodifiableList(new ArrayList<>(cfNames));
            this.vSpecs = Collections.unmodifiableList(new ArrayList<>(vSpecs));
        }
    }


    private static final Map<String, Definition> REGISTRY = new LinkedHashMap<>();

    private static void register(String structurePhase, String model, int numComponents,
                                 double[][] coords, List<String> cfNames, List<VSpec> vSpecs) {
        REGISTRY.put(structurePhase + "_" + model.toUpperCase() + "_" + numComponents,
                     new Definition(coords, cfNames, vSpecs));
    }

    static {
        // -----------------------------------------------------------------
        // BCC_A2  |  T-model  |  binary (K=2)
        //
        // Logical site coordinates (fractional):
        //   p1 = {0.0,  0.0,  0.0}
        //   p2 = {0.5, -0.5,  0.5}
        //   p3 = {0.5,  0.5,  0.5}
        //   p4 = {1.0,  0.0,  0.0}
        //
        // Pair types:
        //   p1-p4 and p2-p3 are II-n pairs  (used in v22AB)
        //   p1-p2 (and others) are I-n pairs (used in v21AB)
        //
        // Atom indices: A=0, B=1
        // -----------------------------------------------------------------
        register("BCC_A2", "T", 2,
            new double[][] {
                { 0.0,  0.0,  0.0 },   // p1
                { 0.5, -0.5,  0.5 },   // p2
                { 0.5,  0.5,  0.5 },   // p3
                { 1.0,  0.0,  0.0 }    // p4
            },
            List.of("v4AB", "v3AB", "v22AB", "v21AB", "xA", "xB"),
            List.of(
                // v4AB = p[1][A]*p[2][B]*p[3][B]*p[4][A]
                VSpec.product(1,0, 2,1, 3,1, 4,0),

                // v3AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                VSpec.diff(
                    new int[]{1,0, 2,1, 3,1},
                    new int[]{1,1, 2,0, 3,0}
                ),

                // v22AB = p[1][A]*p[4][B]  (II-n pair)
                VSpec.product(1,0, 4,1),

                // v21AB = p[1][A]*p[2][B]  (I-n pair)
                VSpec.product(1,0, 2,1),

                // xA = p[1][A]
                VSpec.point(1, 0),

                // xB = p[1][B]
                VSpec.point(1, 1)
            )
        );

        // -----------------------------------------------------------------
        // BCC_A2  |  T-model  |  ternary (K=3)
        //
        // Same logical site coordinates as binary:
        //   p1 = {0.0,  0.0,  0.0}
        //   p2 = {0.5, -0.5,  0.5}
        //   p3 = {0.5,  0.5,  0.5}
        //   p4 = {1.0,  0.0,  0.0}
        //
        // Atom indices: A=0, B=1, C=2
        //
        // 21 CVs total: 6 tetr + 6 tri + 6 pair + 3 point
        //   Tetrahedron (6):
        //     3 binary:  v4AB, v4AC, v4BC
        //     3 ternary: v4ABC1, v4ABC2, v4ABC3
        //   Triangle (6):
        //     3 binary:  v3AB, v3AC, v3BC
        //     3 ternary: v3ABC1, v3ABC2, v3ABC3
        //   Pair (6):
        //     3 II-n:    v22AB, v22AC, v22BC
        //     3 I-n:     v21AB, v21AC, v21BC
        //   Point (3):  xA, xB, xC
        // -----------------------------------------------------------------
        register("BCC_A2", "T", 3,
            new double[][] {
                { 0.0,  0.0,  0.0 },   // p1
                { 0.5, -0.5,  0.5 },   // p2
                { 0.5,  0.5,  0.5 },   // p3
                { 1.0,  0.0,  0.0 }    // p4
            },
            List.of(
                "v4AB",   "v4AC",   "v4BC",
                "v4ABC1", "v4ABC2", "v4ABC3",
                "v3AB",   "v3AC",   "v3BC",
                "v3ABC1", "v3ABC2", "v3ABC3",
                "v22AB",  "v22AC",  "v22BC",
                "v21AB",  "v21AC",  "v21BC",
                "xA", "xB", "xC"
            ),
            List.of(
                // ---- tetrahedra ----

                // v4AB  = p[1][A]*p[2][B]*p[3][B]*p[4][A]
                VSpec.product(1,0, 2,1, 3,1, 4,0),

                // v4AC  = p[1][A]*p[2][C]*p[3][C]*p[4][A]
                VSpec.product(1,0, 2,2, 3,2, 4,0),

                // v4BC  = p[1][B]*p[2][C]*p[3][C]*p[4][B]
                VSpec.product(1,1, 2,2, 3,2, 4,1),

                // v4ABC1 = p[1][B]*p[2][A]*p[3][A]*p[4][C]
                VSpec.product(1,1, 2,0, 3,0, 4,2),

                // v4ABC2 = p[1][A]*p[2][B]*p[3][B]*p[4][C]
                VSpec.product(1,0, 2,1, 3,1, 4,2),

                // v4ABC3 = p[1][A]*p[2][C]*p[3][C]*p[4][B]
                VSpec.product(1,0, 2,2, 3,2, 4,1),

                // ---- triangles ----

                // v3AB   = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                VSpec.diff(
                    new int[]{1,0, 2,1, 3,1},
                    new int[]{1,1, 2,0, 3,0}
                ),

                // v3AC   = p[1][A]*p[2][C]*p[3][C] - p[1][C]*p[2][A]*p[3][A]
                VSpec.diff(
                    new int[]{1,0, 2,2, 3,2},
                    new int[]{1,2, 2,0, 3,0}
                ),

                // v3BC   = p[1][B]*p[2][C]*p[3][C] - p[1][C]*p[2][B]*p[3][B]
                VSpec.diff(
                    new int[]{1,1, 2,2, 3,2},
                    new int[]{1,2, 2,1, 3,1}
                ),

                // v3ABC1 = p[1][C]*p[2][A]*p[3][B]
                VSpec.product(1,2, 2,0, 3,1),

                // v3ABC2 = p[1][B]*p[2][A]*p[3][C]
                VSpec.product(1,1, 2,0, 3,2),

                // v3ABC3 = p[1][A]*p[2][B]*p[3][C]
                VSpec.product(1,0, 2,1, 3,2),

                // ---- pairs ----

                // v22AB = p[1][A]*p[4][B]  (II-n pair)
                VSpec.product(1,0, 4,1),

                // v22AC = p[1][A]*p[4][C]  (II-n pair)
                VSpec.product(1,0, 4,2),

                // v22BC = p[1][B]*p[4][C]  (II-n pair)
                VSpec.product(1,1, 4,2),

                // v21AB = p[1][A]*p[2][B]  (I-n pair)
                VSpec.product(1,0, 2,1),

                // v21AC = p[1][A]*p[2][C]  (I-n pair)
                VSpec.product(1,0, 2,2),

                // v21BC = p[1][B]*p[2][C]  (I-n pair)
                VSpec.product(1,1, 2,2),

                // ---- points ----

                // xA = p[1][A]
                VSpec.point(1, 0),

                // xB = p[1][B]
                VSpec.point(1, 1),

                // xC = p[1][C]
                VSpec.point(1, 2)
            )
        );

        // -----------------------------------------------------------------
        // BCC_A2  |  T-model  |  quaternary (K=4)
        //
        // Same logical site coordinates as binary/ternary:
        //   p1 = {0.0,  0.0,  0.0}
        //   p2 = {0.5, -0.5,  0.5}
        //   p3 = {0.5,  0.5,  0.5}
        //   p4 = {1.0,  0.0,  0.0}
        //
        // Atom indices: A=0, B=1, C=2, D=3
        //
        // 55 CVs total: 21 tetr + 18 tri + 12 pair + 4 point
        //   Tetrahedron (21):
        //     6 binary:      v4AB, v4AC, v4AD, v4BC, v4BD, v4CD
        //     12 ternary:    v4ABC1/2/3, v4ABD1/2/3, v4ACD1/2/3, v4BCD1/2/3
        //     3 quaternary:  v4ABCD1/2/3
        //   Triangle (18):
        //     NOTE: quaternary uses sites p1,p2,p4 (not p1,p2,p3)
        //     6 binary:      v3AB, v3AC, v3AD, v3BC, v3BD, v3CD
        //     12 ternary:    v3ABC1/2/3, v3ABD1/2/3, v3ACD1/2/3, v3BCD1/2/3
        //   Pair (12):
        //     6 II-n (p1,p4): v22AB, v22AC, v22AD, v22BC, v22BD, v22CD
        //     6 I-n  (p1,p2): v21AB, v21AC, v21AD, v21BC, v21BD, v21CD
        //   Point (4): xA, xB, xC, xD
        // -----------------------------------------------------------------
        register("BCC_A2", "T", 4,
            new double[][] {
                { 0.0,  0.0,  0.0 },   // p1
                { 0.5, -0.5,  0.5 },   // p2
                { 0.5,  0.5,  0.5 },   // p3
                { 1.0,  0.0,  0.0 }    // p4
            },
            List.of(
                "v4AB",    "v4AC",    "v4AD",    "v4BC",    "v4BD",    "v4CD",
                "v4ABC1",  "v4ABC2",  "v4ABC3",
                "v4ABD1",  "v4ABD2",  "v4ABD3",
                "v4ACD1",  "v4ACD2",  "v4ACD3",
                "v4BCD1",  "v4BCD2",  "v4BCD3",
                "v4ABCD1", "v4ABCD2", "v4ABCD3",
                "v3AB",    "v3AC",    "v3AD",    "v3BC",    "v3BD",    "v3CD",
                "v3ABC1",  "v3ABC2",  "v3ABC3",
                "v3ABD1",  "v3ABD2",  "v3ABD3",
                "v3ACD1",  "v3ACD2",  "v3ACD3",
                "v3BCD1",  "v3BCD2",  "v3BCD3",
                "v22AB",   "v22AC",   "v22AD",   "v22BC",   "v22BD",   "v22CD",
                "v21AB",   "v21AC",   "v21AD",   "v21BC",   "v21BD",   "v21CD",
                "xA", "xB", "xC", "xD"
            ),
            List.of(
                // ---- tetrahedra: 6 binary ----

                // v4AB  = p[1][A]*p[2][B]*p[3][B]*p[4][A]
                VSpec.product(1,0, 2,1, 3,1, 4,0),
                // v4AC  = p[1][A]*p[2][C]*p[3][C]*p[4][A]
                VSpec.product(1,0, 2,2, 3,2, 4,0),
                // v4AD  = p[1][A]*p[2][D]*p[3][D]*p[4][A]
                VSpec.product(1,0, 2,3, 3,3, 4,0),
                // v4BC  = p[1][B]*p[2][C]*p[3][C]*p[4][B]
                VSpec.product(1,1, 2,2, 3,2, 4,1),
                // v4BD  = p[1][B]*p[2][D]*p[3][D]*p[4][B]
                VSpec.product(1,1, 2,3, 3,3, 4,1),
                // v4CD  = p[1][C]*p[2][D]*p[3][D]*p[4][C]
                VSpec.product(1,2, 2,3, 3,3, 4,2),

                // ---- tetrahedra: 12 ternary ----

                // v4ABC1 = p[1][B]*p[2][A]*p[3][A]*p[4][C]
                VSpec.product(1,1, 2,0, 3,0, 4,2),
                // v4ABC2 = p[1][A]*p[2][B]*p[3][B]*p[4][C]
                VSpec.product(1,0, 2,1, 3,1, 4,2),
                // v4ABC3 = p[1][A]*p[2][C]*p[3][C]*p[4][B]
                VSpec.product(1,0, 2,2, 3,2, 4,1),

                // v4ABD1 = p[1][B]*p[2][A]*p[3][A]*p[4][D]
                VSpec.product(1,1, 2,0, 3,0, 4,3),
                // v4ABD2 = p[1][A]*p[2][B]*p[3][B]*p[4][D]
                VSpec.product(1,0, 2,1, 3,1, 4,3),
                // v4ABD3 = p[1][A]*p[2][D]*p[3][D]*p[4][B]
                VSpec.product(1,0, 2,3, 3,3, 4,1),

                // v4ACD1 = p[1][C]*p[2][A]*p[3][A]*p[4][D]
                VSpec.product(1,2, 2,0, 3,0, 4,3),
                // v4ACD2 = p[1][A]*p[2][C]*p[3][C]*p[4][D]
                VSpec.product(1,0, 2,2, 3,2, 4,3),
                // v4ACD3 = p[1][A]*p[2][D]*p[3][D]*p[4][C]
                VSpec.product(1,0, 2,3, 3,3, 4,2),

                // v4BCD1 = p[1][C]*p[2][B]*p[3][B]*p[4][D]
                VSpec.product(1,2, 2,1, 3,1, 4,3),
                // v4BCD2 = p[1][B]*p[2][C]*p[3][C]*p[4][D]
                VSpec.product(1,1, 2,2, 3,2, 4,3),
                // v4BCD3 = p[1][B]*p[2][D]*p[3][D]*p[4][C]
                VSpec.product(1,1, 2,3, 3,3, 4,2),

                // ---- tetrahedra: 3 quaternary ----

                // v4ABCD1 = p[1][A]*p[2][C]*p[3][D]*p[4][B]
                VSpec.product(1,0, 2,2, 3,3, 4,1),
                // v4ABCD2 = p[1][A]*p[2][B]*p[3][D]*p[4][C]
                VSpec.product(1,0, 2,1, 3,3, 4,2),
                // v4ABCD3 = p[1][A]*p[2][B]*p[3][C]*p[4][D]
                VSpec.product(1,0, 2,1, 3,2, 4,3),

                // ---- triangles: 6 binary (sites p1,p2,p4) ----

                // v3AB  = p[1][B]*p[2][A]*p[4][B] - p[1][A]*p[2][B]*p[4][A]
                VSpec.diff(new int[]{1,1, 2,0, 4,1}, new int[]{1,0, 2,1, 4,0}),
                // v3AC  = p[1][C]*p[2][A]*p[4][C] - p[1][A]*p[2][C]*p[4][A]
                VSpec.diff(new int[]{1,2, 2,0, 4,2}, new int[]{1,0, 2,2, 4,0}),
                // v3AD  = p[1][D]*p[2][A]*p[4][D] - p[1][A]*p[2][D]*p[4][A]
                VSpec.diff(new int[]{1,3, 2,0, 4,3}, new int[]{1,0, 2,3, 4,0}),
                // v3BC  = p[1][C]*p[2][B]*p[4][C] - p[1][B]*p[2][C]*p[4][B]
                VSpec.diff(new int[]{1,2, 2,1, 4,2}, new int[]{1,1, 2,2, 4,1}),
                // v3BD  = p[1][D]*p[2][B]*p[4][D] - p[1][B]*p[2][D]*p[4][B]
                VSpec.diff(new int[]{1,3, 2,1, 4,3}, new int[]{1,1, 2,3, 4,1}),
                // v3CD  = p[1][D]*p[2][C]*p[4][D] - p[1][C]*p[2][D]*p[4][C]
                VSpec.diff(new int[]{1,3, 2,2, 4,3}, new int[]{1,2, 2,3, 4,2}),

                // ---- triangles: 12 ternary (sites p1,p2,p4) ----

                // v3ABC1 = p[1][B]*p[2][A]*p[4][C]
                VSpec.product(1,1, 2,0, 4,2),
                // v3ABC2 = p[1][A]*p[2][B]*p[4][C]
                VSpec.product(1,0, 2,1, 4,2),
                // v3ABC3 = p[1][A]*p[2][C]*p[4][B]
                VSpec.product(1,0, 2,2, 4,1),

                // v3ABD1 = p[1][B]*p[2][A]*p[4][D]
                VSpec.product(1,1, 2,0, 4,3),
                // v3ABD2 = p[1][A]*p[2][B]*p[4][D]
                VSpec.product(1,0, 2,1, 4,3),
                // v3ABD3 = p[1][A]*p[2][D]*p[4][B]
                VSpec.product(1,0, 2,3, 4,1),

                // v3ACD1 = p[1][C]*p[2][A]*p[4][D]
                VSpec.product(1,2, 2,0, 4,3),
                // v3ACD2 = p[1][A]*p[2][C]*p[4][D]
                VSpec.product(1,0, 2,2, 4,3),
                // v3ACD3 = p[1][A]*p[2][D]*p[4][C]
                VSpec.product(1,0, 2,3, 4,2),

                // v3BCD1 = p[1][C]*p[2][B]*p[4][D]
                VSpec.product(1,2, 2,1, 4,3),
                // v3BCD2 = p[1][B]*p[2][C]*p[4][D]
                VSpec.product(1,1, 2,2, 4,3),
                // v3BCD3 = p[1][B]*p[2][D]*p[4][C]
                VSpec.product(1,1, 2,3, 4,2),

                // ---- pairs: 6 II-n (p1,p4) ----

                // v22AB = p[1][A]*p[4][B]
                VSpec.product(1,0, 4,1),
                // v22AC = p[1][A]*p[4][C]
                VSpec.product(1,0, 4,2),
                // v22AD = p[1][A]*p[4][D]
                VSpec.product(1,0, 4,3),
                // v22BC = p[1][B]*p[4][C]
                VSpec.product(1,1, 4,2),
                // v22BD = p[1][B]*p[4][D]
                VSpec.product(1,1, 4,3),
                // v22CD = p[1][C]*p[4][D]
                VSpec.product(1,2, 4,3),

                // ---- pairs: 6 I-n (p1,p2) ----

                // v21AB = p[1][A]*p[2][B]
                VSpec.product(1,0, 2,1),
                // v21AC = p[1][A]*p[2][C]
                VSpec.product(1,0, 2,2),
                // v21AD = p[1][A]*p[2][D]
                VSpec.product(1,0, 2,3),
                // v21BC = p[1][B]*p[2][C]
                VSpec.product(1,1, 2,2),
                // v21BD = p[1][B]*p[2][D]
                VSpec.product(1,1, 2,3),
                // v21CD = p[1][C]*p[2][D]
                VSpec.product(1,2, 2,3),

                // ---- points ----

                // xA = p[1][A]
                VSpec.point(1, 0),
                // xB = p[1][B]
                VSpec.point(1, 1),
                // xC = p[1][C]
                VSpec.point(1, 2),
                // xD = p[1][D]
                VSpec.point(1, 3)
            )
        );

        // -----------------------------------------------------------------
        // Add FCC_A1, HCP_A3 etc. here when ready.
        // -----------------------------------------------------------------
    }
    public static CvCfBasis generate(
            String structurePhase,
            ClusterCFIdentificationPipeline.PipelineResult pr,
            CMatrixPipeline.CMatrixData matrixData,
            String model,
            Consumer<String> sink) {

        if (!"T".equalsIgnoreCase(model)) {
            throw new UnsupportedOperationException(
                "Dynamic generation only supported for T-model; got: " + model);
        }

        emit(sink, "========================================================");
        emit(sink, "  CVCF BASIS GENERATION");
        emit(sink, "========================================================");

        int numComponents = pr.getNumComponents();
        ClusterIdentificationResult clusterResult = pr.toClusterIdentificationResult();
        CFIdentificationResult cfResult = pr.toCFIdentificationResult();
        
        // ---------------------------------------------------------
        // Stage 1: Load Specification
        // ---------------------------------------------------------
        emit(sink, "\n[STAGE 1] Loading CVCF Specification");
        Definition def = REGISTRY.get(structurePhase + "_" + model.toUpperCase() + "_" + numComponents);
        if (def == null) throw new IllegalArgumentException("Unregistered CVCF combination.");
        emit(sink, "  Loaded " + def.cfNames.size() + " CVs for " + structurePhase + " K=" + numComponents);

        for (int i = 0; i < def.vSpecs.size(); i++) {
            emit(sink, String.format("  v[%d] (%s) = %s", i, def.cfNames.get(i), formatVSpec(def.vSpecs.get(i))));
        }

        // ---------------------------------------------------------
        // Stage 2: Resolve Coordinates
        // ---------------------------------------------------------
        emit(sink, "\n[STAGE 2] Resolving Fractional Coordinates to Physical Indices");
        List<Position> siteList = matrixData.getSiteList();
        Cluster maxCluster = clusterResult.getOrdClusterData().getCoordList().get(0).get(0);
        Map<Integer, Integer> siteMap = resolveSiteMap(def.logicalSiteCoords, maxCluster, siteList, sink);

        for (Map.Entry<Integer, Integer> entry : siteMap.entrySet()) {
            emit(sink, "  Logical p[" + entry.getKey() + "] mapping to physical internal site #" + entry.getValue());
        }

        // ---------------------------------------------------------
        // Stage 3: Orthogonal CF Data prep
        // ---------------------------------------------------------
        emit(sink, "\n[STAGE 3] Extracting Orthogonal Basis Metadata");
        int totalCfs  = cfResult.getTcf();
        int basisSize = def.vSpecs.size();
        
        if (basisSize != totalCfs + 1) {
            LOG.warning(String.format("Basis size mismatch: CVCF=%d, Orthogonal=%d (tcf+1).", basisSize, totalCfs + 1));
            emit(sink, String.format("  [WARN] Size mismatch! CVCF=%d, Orth=(%d+1)", basisSize, totalCfs));
        }

        Map<CFIndex, Integer> cfColMap = buildCfColumnMap(cfResult.getLcf());
        emit(sink, "  Built CF column map for " + cfColMap.size() + " unique orthogonal physical correlation functions.");

        // ---------------------------------------------------------
        // Stage 4: Matrix M construction via Pipeline 
        // ---------------------------------------------------------
        emit(sink, "\n[STAGE 4] Assembling Transformation M-Matrix by querying Pipeline");
        double[][] M = buildMMatrix(def.vSpecs, siteMap, matrixData, cfColMap, totalCfs, basisSize, sink);

        // ---------------------------------------------------------
        // Stage 5: Inverting to T-matrix
        // ---------------------------------------------------------
        emit(sink, "\n[STAGE 5] Inverting M-Matrix to T-Matrix Basis");
        double[][] T    = LinearAlgebra.invert(M);
        double[][] Tinv = M;

        emit(sink, "  Geometric Inversion complete, T-matrix bound (" + T.length + "x" + T[0].length + ")");
        
        emit(sink, "\n  T-MATRIX (Orthogonal â†’ CVCF MAPPING):");
        for (int i = 0; i < T.length; i++) {
             StringBuilder line = new StringBuilder();
            line.append(String.format("  Row %2d -> ", i));
            boolean hasTerms = false;
            for (int j = 0; j < T[i].length; j++) {
                 if (Math.abs(T[i][j]) > 1e-10) {
                     line.append(String.format("[V%d=%.3f] ", j, T[i][j]));
                     hasTerms = true;
                 }
             }
             if (!hasTerms) line.append("[0.000]");
             emit(sink, line.toString());
        }
        emit(sink, "  Geometric Inversion complete, T-matrix bound (" + T.length + "x" + T[0].length + ")");
        
        emit(sink, "\n[STAGE 6] Transforming Orthogonal C-Matrix to CVCF Basis");
        CMatrixPipeline.CMatrixData cvcfData = matrixData.transform(T);
        emit(sink, "  C-Matrix natively bound inside Basis mapping.");

        int numNonPointCfs = def.cfNames.size() - numComponents;
        List<String> eciNames = cfResult.getEONames().subList(0, numNonPointCfs);
        emit(sink, "  Binding " + eciNames.size() + " CVCF non-point components to Hamiltonian map.");

        // --- Integrated Verification ---
        double[] xEqui = new double[numComponents];
        java.util.Arrays.fill(xEqui, 1.0 / numComponents);
        double[] vFull = computeRandomCvcfCFs(xEqui, pr, Tinv);

        emit(sink, "\n  [SELF-TEST] RANDOM CVCF CFs (Equiatomic):");
        for (int i = 0; i < vFull.length; i++) {
            emit(sink, String.format("    %-10s = %12.8f", def.cfNames.get(i), vFull[i]));
        }

        emit(sink, "\n  [SELF-TEST] CV VERIFICATION (Equiatomic Disorder):");
        double[][][] cv = CMatrixPipeline.evaluateCVs(vFull, cvcfData.cmat, cvcfData.lcv, pr.getTcdis(), pr.getLc());
        for (int t = 0; t < cv.length; t++) {
            Cluster representative = pr.getDisClusData().getClusCoordList().get(t);
            int n = representative.getAllSites().size();
            double expected = Math.pow(1.0 / numComponents, n);
            for (int j = 0; j < cv[t].length; j++) {
                emit(sink, String.format("    Type t=%d, Group j=%d (n=%d): %12.8f (Diff: %.2e)", 
                    t, j, n, cv[t][j][0], Math.abs(cv[t][j][0] - expected)));
            }
        }

        emit(sink, "\n========================================================");
        emit(sink, "  CVCF PIPELINE COMPLETE");
        emit(sink, "========================================================");

        return new CvCfBasis(structurePhase, model, numComponents,
                def.cfNames, eciNames, numNonPointCfs, T, Tinv, cvcfData);
    }

    private static String formatVSpec(VSpec spec) {
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
        for (int logIdx = 0; logIdx < logicalSiteCoords.length; logIdx++) {
            double[] coord = logicalSiteCoords[logIdx];
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
    // M matrix assembly natively derived via CMatrixPipeline
    // =========================================================================

    private static double[][] buildMMatrix(
            List<VSpec> vSpecs,
            Map<Integer, Integer> siteMap,
            CMatrixPipeline.CMatrixData matrixData,
            Map<CFIndex, Integer> cfColMap,
            int totalCfs,
            int basisSize,
            Consumer<String> sink) {

        double[][] M = new double[basisSize][basisSize];

        for (int i = 0; i < basisSize; i++) {
            VSpec spec = vSpecs.get(i);
            double[] plusRow = evaluateSpecTerm(spec.plusTerm, siteMap, matrixData, cfColMap, totalCfs);
            
            if (spec.isDiff()) {
                double[] minusRow = evaluateSpecTerm(spec.minusTerm, siteMap, matrixData, cfColMap, totalCfs);
                for (int j = 0; j < M[i].length; j++) {
                    M[i][j] = plusRow[j] - minusRow[j];
                }
            } else {
                M[i] = plusRow;
            }

            // Optional trace of the array coefficients
            if (i < 3 || i > basisSize - 4) { // keep logs somewhat clipped
                 StringBuilder line = new StringBuilder();
                 line.append(String.format("  Row %2d -> ", i));
                 for (int j = 0; j < M[i].length; j++) {
                     if (Math.abs(M[i][j]) > 1e-10) {
                         line.append(String.format("[C%d=%.3f] ", j, M[i][j]));
                     }
                 }
                 emit(sink, line.toString());
            }
            if (i == 3 && basisSize > 6) {
                 emit(sink, "  ... (rows clipped for brevity)");
            }
        }
        return M;
    }

    private static double[] evaluateSpecTerm(
            int[] termPairs, 
            Map<Integer, Integer> siteMap, 
            CMatrixPipeline.CMatrixData matrixData, 
            Map<CFIndex, Integer> cfColMap, 
            int totalCfs) {

        List<Integer> siteIndices = new ArrayList<>();
        int[] config = new int[termPairs.length / 2];
        for (int i = 0; i < termPairs.length; i += 2) {
            int logicalSite = termPairs[i];
            int atom = termPairs[i + 1];
            siteIndices.add(siteMap.get(logicalSite));
            config[i / 2] = atom;
        }
        return matrixData.expandProbabilityExpression(siteIndices, config, cfColMap, totalCfs);
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
}
