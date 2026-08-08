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

import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterMath;
import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.cluster.CFIdentificationResult;
import org.ce.model.cluster.ClusterKeys.CFIndex;
import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterPrimitives.Position;
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
    public final List<String> cfDefinitions;
    private final Map<String, Integer> cfNameIndex;

    private CvCfBasis(String structurePhase, String model, int numComponents,
            List<String> cfNames, List<String> eciNames, int numNonPointCfs,
            double[][] T, double[][] Tinv, CMatrixPipeline.CMatrixData cvcfCMatrixData,
            List<String> cfDefinitions) {
        this.structurePhase = structurePhase;
        this.model = (model == null ? "" : model.toUpperCase());
        this.numComponents = numComponents;
        this.cfNames = cfNames;
        this.eciNames = eciNames;
        this.numNonPointCfs = numNonPointCfs;
        this.T = T;
        this.Tinv = Tinv;
        this.cvcfCMatrixData = cvcfCMatrixData;
        this.cfDefinitions = cfDefinitions;
        Map<String, Integer> idx = new LinkedHashMap<>(cfNames.size() * 2);
        for (int i = 0; i < cfNames.size(); i++)
            idx.put(cfNames.get(i), i);
        this.cfNameIndex = idx;
    }

    public static boolean isSupported(String structurePhase, String model, int numComponents) {
        return REGISTRY.containsKey(structurePhase + "_" + model.toUpperCase() + "_" + numComponents);
    }

    public static int getNumNonPointCfs(String structurePhase, String model, int numComponents) {
        Definition def = REGISTRY.get(structurePhase + "_" + model.toUpperCase() + "_" + numComponents);
        if (def == null)
            throw new IllegalArgumentException("Unregistered CVCF combination.");
        return def.cfNames.size() - def.numPointCfs;
    }

    public static String supportedSummary() {
        return "Supported combinations: " + String.join(", ", REGISTRY.keySet());
    }

    /**
     * Canonical names of the non-point CFs for a registered combination — i.e. exactly
     * the ECI terms a Hamiltonian must supply, in basis order.
     *
     * <p>Point CFs (mole fractions, and any order parameter such as B2's {@code eta})
     * trail the list and carry no ECI, so they are excluded.</p>
     *
     * <p>Exposed for API capability discovery, so external callers can ask which ECI
     * names a system expects instead of hard-coding them.</p>
     *
     * @throws IllegalArgumentException if the combination is not registered
     */
    public static List<String> getNonPointCfNames(String structurePhase, String model, int numComponents) {
        Definition def = REGISTRY.get(structurePhase + "_" + model.toUpperCase() + "_" + numComponents);
        if (def == null)
            throw new IllegalArgumentException("Unregistered CVCF combination. " + supportedSummary());
        return List.copyOf(def.cfNames.subList(0, def.cfNames.size() - def.numPointCfs));
    }

    /** Registry keys of all supported {structure}_{MODEL}_{ncomp} combinations. */
    public static List<String> supportedKeys() {
        return List.copyOf(REGISTRY.keySet());
    }

    public int indexOfCf(String name) {
        return cfNameIndex.getOrDefault(name, -1);
    }

    public int totalCfs() {
        return cfNames.size();
    }

    /**
     * All intermediate vectors produced when computing the random disordered state.
     */
    public static final class RandomStateVectors {
        /** Point basis functions φₖ(x) evaluated at the composition (length = K-1). */
        public final double[] pointCF;
        /**
         * Orthogonal non-point CFs at random state — product of pointCFs (length =
         * orthNcf).
         */
        public final double[] uOrthNonPoint;
        /** Full orthogonal vector passed to T⁻¹: [uOrthNonPoint | pointCF | 1?] */
        public final double[] uOrthFull;
        /** CVCF-transformed non-point CFs (length = ncf) — what N-R minimizes over. */
        public final double[] uCvcfNonPoint;
        /**
         * Full CVCF vector including point vars appended: [uCvcfNonPoint |
         * moleFractions]
         */
        public final double[] uCvcfFull;

        RandomStateVectors(double[] pointCF, double[] uOrthNonPoint, double[] uOrthFull,
                double[] uCvcfNonPoint, double[] uCvcfFull) {
            this.pointCF = pointCF;
            this.uOrthNonPoint = uOrthNonPoint;
            this.uOrthFull = uOrthFull;
            this.uCvcfNonPoint = uCvcfNonPoint;
            this.uCvcfFull = uCvcfFull;
        }
    }

    /**
     * Computes the random disordered state and returns all intermediate vectors.
     * Useful for debugging the CVCF transformation.
     */
    public RandomStateVectors computeRandomStateVectors(double[] moleFractions, int[][] orthCfBasisIndices) {
        if (orthCfBasisIndices == null || orthCfBasisIndices.length == 0) {
            throw new IllegalArgumentException("orthCfBasisIndices must not be null or empty.");
        }
        int K = moleFractions.length;
        int nxcf = K - 1;
        int orthTcf = orthCfBasisIndices.length;
        int orthNcf = orthTcf - nxcf;

        double[] basisVec = ClusterMath.buildBasis(K);
        double[] pointCF = new double[nxcf];
        for (int k = 0; k < nxcf; k++)
            for (int i = 0; i < K; i++)
                pointCF[k] += moleFractions[i] * Math.pow(basisVec[i], k + 1);

        double[] uOrthNonPoint = new double[orthNcf];
        for (int col = 0; col < orthNcf; col++) {
            double val = 1.0;
            for (int b : orthCfBasisIndices[col]) {
                if (b < 1 || b > nxcf) {
                    throw new IllegalArgumentException(
                            "orthCfBasisIndices[" + col + "] contains out-of-range basis index " + b
                                    + " (valid range 1.." + nxcf + " for " + K + " components).");
                }
                val *= pointCF[b - 1];
            }
            uOrthNonPoint[col] = val;
        }

        double[][] tInv = resolvedTinv();
        if (tInv.length == 0 || tInv[0].length == 0) {
            throw new IllegalStateException("resolvedTinv() returned an empty T-inverse matrix.");
        }
        int tRows = tInv[0].length;
        if (tRows != orthTcf && tRows != orthTcf + 1) {
            throw new IllegalStateException(
                    "T-matrix dimension mismatch: T-inverse has " + tRows + " columns but orthCfBasisIndices "
                            + "implies orthTcf=" + orthTcf + " (expected " + orthTcf + " or " + (orthTcf + 1) + ").");
        }
        double[] uOrthFull = new double[tRows];
        System.arraycopy(uOrthNonPoint, 0, uOrthFull, 0, orthNcf);
        System.arraycopy(pointCF, 0, uOrthFull, orthNcf, nxcf);
        if (tRows == orthTcf + 1)
            uOrthFull[tRows - 1] = 1.0;

        int ncf = numNonPointCfs;
        int tcf = totalCfs();
        double[] vFull = new double[tcf];
        for (int j = 0; j < tcf; j++) {
            double sum = 0.0;
            for (int i = 0; i < tInv[j].length; i++)
                sum += tInv[j][i] * uOrthFull[i];
            vFull[j] = sum;
        }
        for (int i = 0; i < K; i++)
            vFull[ncf + i] = moleFractions[i];

        double[] uCvcfNonPoint = Arrays.copyOf(vFull, ncf);
        return new RandomStateVectors(pointCF, uOrthNonPoint, uOrthFull, uCvcfNonPoint, vFull);
    }

    public double[] computeRandomState(double[] moleFractions, int[][] orthCfBasisIndices) {
        return computeRandomStateVectors(moleFractions, orthCfBasisIndices).uCvcfFull;
    }

    private double[][] resolvedTinv() {
        if (Tinv != null)
            return Tinv;
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
     * Diagnostic verification of the C-matrix in this basis at the disordered
     * state.
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
                pr.getLc());

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
                        if (error > 1e-9)
                            msg += " [!] DISCREPANCY";
                    }
                    emit(sink, msg);
                }
            }
        }
        emit(sink, "\n=== CVCF Verification: COMPLETE ===");
    }

    public static final class VSpec {

        /** One weighted product term in the linear combination. */
        public static final class Term {
            public final double coefficient;
            /** Flat [site1, atom1, site2, atom2, ...]. */
            public final int[] siteAtomPairs;

            Term(double coefficient, int[] siteAtomPairs) {
                this.coefficient = coefficient;
                this.siteAtomPairs = siteAtomPairs;
            }
        }

        /**
         * v = Σ coefficient_i * p[site1][atom1] * p[site2][atom2] * ...
         * General weighted linear combination of product terms. product(),
         * diff(), and point() are convenience constructors for the common
         * single- and two-term cases.
         */
        public final List<Term> terms;

        private VSpec(List<Term> terms) {
            this.terms = Collections.unmodifiableList(terms);
        }

        /** v = p[s1][a1] * p[s2][a2] * ... */
        public static VSpec product(int... siteAtomPairs) {
            return combo(1.0, siteAtomPairs);
        }

        /** v = product(plusPairs) - product(minusPairs). */
        public static VSpec diff(int[] plusPairs, int[] minusPairs) {
            List<Term> terms = new ArrayList<>();
            terms.add(makeTerm(1.0, plusPairs));
            terms.add(makeTerm(-1.0, minusPairs));
            return new VSpec(terms);
        }

        /** Convenience: single-site probability, i.e. point CF. v = p[site][atom]. */
        public static VSpec point(int logicalSite, int atom) {
            return product(logicalSite, atom);
        }

        /** v = coefficient * product(siteAtomPairs). Single-term convenience for {@link #combo(Object...)}. */
        public static VSpec combo(double coefficient, int[] siteAtomPairs) {
            List<Term> terms = new ArrayList<>();
            terms.add(makeTerm(coefficient, siteAtomPairs));
            return new VSpec(terms);
        }

        /**
         * v = Σ coefficient_i * product(siteAtomPairs_i). General weighted
         * linear combination, for CVs that are neither a pure product nor a
         * simple plus/minus difference — e.g. the FCC ternary/quaternary
         * tetrahedron CVs, which are {@code (t1 + 2*t2 - t3 - t4)/3}.
         *
         * @param coeffAndPairs alternating {@code double coefficient, int[] siteAtomPairs}
         *                      arguments, one pair per term.
         */
        public static VSpec combo(Object... coeffAndPairs) {
            if (coeffAndPairs.length % 2 != 0)
                throw new IllegalArgumentException("combo() takes (coefficient, pairs) pairs");
            List<Term> terms = new ArrayList<>();
            for (int i = 0; i < coeffAndPairs.length; i += 2) {
                double coeff = ((Number) coeffAndPairs[i]).doubleValue();
                int[] pairs = (int[]) coeffAndPairs[i + 1];
                terms.add(makeTerm(coeff, pairs));
            }
            return new VSpec(terms);
        }

        private static Term makeTerm(double coefficient, int[] siteAtomPairs) {
            if (siteAtomPairs.length % 2 != 0)
                throw new IllegalArgumentException("siteAtomPairs must be even-length");
            return new Term(coefficient, siteAtomPairs);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < terms.size(); i++) {
                Term t = terms.get(i);
                if (i > 0) {
                    sb.append(t.coefficient < 0 ? " - " : " + ");
                } else if (t.coefficient < 0) {
                    sb.append("-");
                }
                double absCoeff = Math.abs(t.coefficient);
                if (absCoeff != 1.0)
                    sb.append(absCoeff).append("*");
                sb.append(formatTerm(t.siteAtomPairs));
            }
            return sb.toString();
        }

        private String formatTerm(int[] pairs) {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < pairs.length; i += 2) {
                if (i > 0)
                    sb.append("*");
                sb.append("p[").append(pairs[i]).append("][").append(atomSymbol(pairs[i + 1])).append("]");
            }
            return sb.toString();
        }

        private String atomSymbol(int atom) {
            switch (atom) {
                case 0: return "A";
                case 1: return "B";
                case 2: return "C";
                case 3: return "D";
                default: return String.valueOf(atom);
            }
        }
    }

    private static final class Definition {
        final double[][] logicalSiteCoords;
        final List<String> cfNames;
        final List<VSpec> vSpecs;
        /**
         * Number of trailing entries in cfNames/vSpecs that are point-like
         * (composition/order-parameter) CFs carrying no ECI — the rest
         * (leading entries) are non-point CFs bound to Hamiltonian ECI
         * terms. Defaults to numComponents (one point CF per species), but
         * ordered structures may need more — e.g. BCC_B2 has xA, xB, plus
         * a long-range-order parameter eta, so numPointCfs=3 with
         * numComponents=2.
         */
        final int numPointCfs;

        Definition(double[][] logicalSiteCoords, List<String> cfNames, List<VSpec> vSpecs, int numPointCfs) {
            this.logicalSiteCoords = logicalSiteCoords;
            this.cfNames = Collections.unmodifiableList(new ArrayList<>(cfNames));
            this.vSpecs = Collections.unmodifiableList(new ArrayList<>(vSpecs));
            this.numPointCfs = numPointCfs;
        }
    }

    private static final Map<String, Definition> REGISTRY = new LinkedHashMap<>();

    private static void register(String structurePhase, String model, int numComponents,
            double[][] coords, List<String> cfNames, List<VSpec> vSpecs) {
        register(structurePhase, model, numComponents, coords, cfNames, vSpecs, numComponents);
    }

    private static void register(String structurePhase, String model, int numComponents,
            double[][] coords, List<String> cfNames, List<VSpec> vSpecs, int numPointCfs) {
        REGISTRY.put(structurePhase + "_" + model.toUpperCase() + "_" + numComponents,
                new Definition(coords, cfNames, vSpecs, numPointCfs));
    }

    static {
        // -----------------------------------------------------------------
        // BCC_A2 | T-model | binary (K=2)
        //
        // Logical site coordinates (fractional):
        // p1 = {0.0, 0.0, 0.0}
        // p2 = {0.5, -0.5, 0.5}
        // p3 = {0.5, 0.5, 0.5}
        // p4 = {1.0, 0.0, 0.0}
        //
        // Pair types:
        // p1-p4 and p2-p3 are II-n pairs (used in v22AB)
        // p1-p2 (and others) are I-n pairs (used in v21AB)
        //
        // Atom indices: A=0, B=1
        // -----------------------------------------------------------------
        register("BCC_A2", "T", 2,
                new double[][] {
                        { 0.0, 0.0, 0.0 }, // p1
                        { 0.5, -0.5, 0.5 }, // p2
                        { 0.5, 0.5, 0.5 }, // p3
                        { 1.0, 0.0, 0.0 } // p4
                },
                List.of("v4AB", "v3AB", "v22AB", "v21AB", "xA", "xB"),
                List.of(
                        // v4AB = p[1][A]*p[2][B]*p[3][B]*p[4][A]
                        VSpec.product(1, 0, 2, 1, 3, 1, 4, 0),

                        // v3AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 1, 2, 0, 3, 0 }),

                        // v22AB = p[1][A]*p[4][B] (II-n pair)
                        VSpec.product(1, 0, 4, 1),

                        // v21AB = p[1][A]*p[2][B] (I-n pair)
                        VSpec.product(1, 0, 2, 1),

                        // xA = p[1][A]
                        VSpec.point(1, 0),

                        // xB = p[1][B]
                        VSpec.point(1, 1)));

        // -----------------------------------------------------------------
        // BCC_A2 | T-model | ternary (K=3)
        //
        // Same logical site coordinates as binary:
        // p1 = {0.0, 0.0, 0.0}
        // p2 = {0.5, -0.5, 0.5}
        // p3 = {0.5, 0.5, 0.5}
        // p4 = {1.0, 0.0, 0.0}
        //
        // Atom indices: A=0, B=1, C=2
        //
        // 21 CVs total: 6 tetr + 6 tri + 6 pair + 3 point
        // Tetrahedron (6):
        // 3 binary: v4AB, v4AC, v4BC
        // 3 ternary: v4ABC1, v4ABC2, v4ABC3
        // Triangle (6):
        // 3 binary: v3AB, v3AC, v3BC
        // 3 ternary: v3ABC1, v3ABC2, v3ABC3
        // Pair (6):
        // 3 II-n: v22AB, v22AC, v22BC
        // 3 I-n: v21AB, v21AC, v21BC
        // Point (3): xA, xB, xC
        // -----------------------------------------------------------------
        register("BCC_A2", "T", 3,
                new double[][] {
                        { 0.0, 0.0, 0.0 }, // p1
                        { 0.5, -0.5, 0.5 }, // p2
                        { 0.5, 0.5, 0.5 }, // p3
                        { 1.0, 0.0, 0.0 } // p4
                },
                List.of(
                        "v4AB", "v4AC", "v4BC",
                        "v4ABC1", "v4ABC2", "v4ABC3",
                        "v3AB", "v3AC", "v3BC",
                        "v3ABC1", "v3ABC2", "v3ABC3",
                        "v22AB", "v22AC", "v22BC",
                        "v21AB", "v21AC", "v21BC",
                        "xA", "xB", "xC"),
                List.of(
                        // ---- tetrahedra ----

                        // v4AB = p[1][A]*p[2][B]*p[3][B]*p[4][A]
                        VSpec.product(1, 0, 2, 1, 3, 1, 4, 0),

                        // v4AC = p[1][A]*p[2][C]*p[3][C]*p[4][A]
                        VSpec.product(1, 0, 2, 2, 3, 2, 4, 0),

                        // v4BC = p[1][B]*p[2][C]*p[3][C]*p[4][B]
                        VSpec.product(1, 1, 2, 2, 3, 2, 4, 1),

                        // v4ABC1 = p[1][B]*p[2][A]*p[3][A]*p[4][C]
                        VSpec.product(1, 1, 2, 0, 3, 0, 4, 2),

                        // v4ABC2 = p[1][A]*p[2][B]*p[3][B]*p[4][C]
                        VSpec.product(1, 0, 2, 1, 3, 1, 4, 2),

                        // v4ABC3 = p[1][A]*p[2][C]*p[3][C]*p[4][B]
                        VSpec.product(1, 0, 2, 2, 3, 2, 4, 1),

                        // ---- triangles ----

                        // v3AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 1, 2, 0, 3, 0 }),

                        // v3AC = p[1][A]*p[2][C]*p[3][C] - p[1][C]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 0, 3, 0 }),

                        // v3BC = p[1][B]*p[2][C]*p[3][C] - p[1][C]*p[2][B]*p[3][B]
                        VSpec.diff(
                                new int[] { 1, 1, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 1, 3, 1 }),

                        // v3ABC1 = p[1][C]*p[2][A]*p[3][B]
                        VSpec.product(1, 2, 2, 0, 3, 1),

                        // v3ABC2 = p[1][B]*p[2][A]*p[3][C]
                        VSpec.product(1, 1, 2, 0, 3, 2),

                        // v3ABC3 = p[1][A]*p[2][B]*p[3][C]
                        VSpec.product(1, 0, 2, 1, 3, 2),

                        // ---- pairs ----

                        // v22AB = p[1][A]*p[4][B] (II-n pair)
                        VSpec.product(1, 0, 4, 1),

                        // v22AC = p[1][A]*p[4][C] (II-n pair)
                        VSpec.product(1, 0, 4, 2),

                        // v22BC = p[1][B]*p[4][C] (II-n pair)
                        VSpec.product(1, 1, 4, 2),

                        // v21AB = p[1][A]*p[2][B] (I-n pair)
                        VSpec.product(1, 0, 2, 1),

                        // v21AC = p[1][A]*p[2][C] (I-n pair)
                        VSpec.product(1, 0, 2, 2),

                        // v21BC = p[1][B]*p[2][C] (I-n pair)
                        VSpec.product(1, 1, 2, 2),

                        // ---- points ----

                        // xA = p[1][A]
                        VSpec.point(1, 0),

                        // xB = p[1][B]
                        VSpec.point(1, 1),

                        // xC = p[1][C]
                        VSpec.point(1, 2)));

        // -----------------------------------------------------------------
        // BCC_A2 | T-model | quaternary (K=4)
        //
        // Same logical site coordinates as binary/ternary:
        // p1 = {0.0, 0.0, 0.0}
        // p2 = {0.5, -0.5, 0.5}
        // p3 = {0.5, 0.5, 0.5}
        // p4 = {1.0, 0.0, 0.0}
        //
        // Atom indices: A=0, B=1, C=2, D=3
        //
        // 55 CVs total: 21 tetr + 18 tri + 12 pair + 4 point
        // Tetrahedron (21):
        // 6 binary: v4AB, v4AC, v4AD, v4BC, v4BD, v4CD
        // 12 ternary: v4ABC1/2/3, v4ABD1/2/3, v4ACD1/2/3, v4BCD1/2/3
        // 3 quaternary: v4ABCD1/2/3
        // Triangle (18):
        // NOTE: quaternary uses sites p1,p2,p4 (not p1,p2,p3)
        // 6 binary: v3AB, v3AC, v3AD, v3BC, v3BD, v3CD
        // 12 ternary: v3ABC1/2/3, v3ABD1/2/3, v3ACD1/2/3, v3BCD1/2/3
        // Pair (12):
        // 6 II-n (p1,p4): v22AB, v22AC, v22AD, v22BC, v22BD, v22CD
        // 6 I-n (p1,p2): v21AB, v21AC, v21AD, v21BC, v21BD, v21CD
        // Point (4): xA, xB, xC, xD
        // -----------------------------------------------------------------
        register("BCC_A2", "T", 4,
                new double[][] {
                        { 0.0, 0.0, 0.0 }, // p1
                        { 0.5, -0.5, 0.5 }, // p2
                        { 0.5, 0.5, 0.5 }, // p3
                        { 1.0, 0.0, 0.0 } // p4
                },
                List.of(
                        "v4AB", "v4AC", "v4AD", "v4BC", "v4BD", "v4CD",
                        "v4ABC1", "v4ABC2", "v4ABC3",
                        "v4ABD1", "v4ABD2", "v4ABD3",
                        "v4ACD1", "v4ACD2", "v4ACD3",
                        "v4BCD1", "v4BCD2", "v4BCD3",
                        "v4ABCD1", "v4ABCD2", "v4ABCD3",
                        "v3AB", "v3AC", "v3AD", "v3BC", "v3BD", "v3CD",
                        "v3ABC1", "v3ABC2", "v3ABC3",
                        "v3ABD1", "v3ABD2", "v3ABD3",
                        "v3ACD1", "v3ACD2", "v3ACD3",
                        "v3BCD1", "v3BCD2", "v3BCD3",
                        "v22AB", "v22AC", "v22AD", "v22BC", "v22BD", "v22CD",
                        "v21AB", "v21AC", "v21AD", "v21BC", "v21BD", "v21CD",
                        "xA", "xB", "xC", "xD"),
                List.of(
                        // ---- tetrahedra: 6 binary ----

                        // v4AB = p[1][A]*p[2][B]*p[3][B]*p[4][A]
                        VSpec.product(1, 0, 2, 1, 3, 1, 4, 0),
                        // v4AC = p[1][A]*p[2][C]*p[3][C]*p[4][A]
                        VSpec.product(1, 0, 2, 2, 3, 2, 4, 0),
                        // v4AD = p[1][A]*p[2][D]*p[3][D]*p[4][A]
                        VSpec.product(1, 0, 2, 3, 3, 3, 4, 0),
                        // v4BC = p[1][B]*p[2][C]*p[3][C]*p[4][B]
                        VSpec.product(1, 1, 2, 2, 3, 2, 4, 1),
                        // v4BD = p[1][B]*p[2][D]*p[3][D]*p[4][B]
                        VSpec.product(1, 1, 2, 3, 3, 3, 4, 1),
                        // v4CD = p[1][C]*p[2][D]*p[3][D]*p[4][C]
                        VSpec.product(1, 2, 2, 3, 3, 3, 4, 2),

                        // ---- tetrahedra: 12 ternary ----

                        // v4ABC1 = p[1][B]*p[2][A]*p[3][A]*p[4][C]
                        VSpec.product(1, 1, 2, 0, 3, 0, 4, 2),
                        // v4ABC2 = p[1][A]*p[2][B]*p[3][B]*p[4][C]
                        VSpec.product(1, 0, 2, 1, 3, 1, 4, 2),
                        // v4ABC3 = p[1][A]*p[2][C]*p[3][C]*p[4][B]
                        VSpec.product(1, 0, 2, 2, 3, 2, 4, 1),

                        // v4ABD1 = p[1][B]*p[2][A]*p[3][A]*p[4][D]
                        VSpec.product(1, 1, 2, 0, 3, 0, 4, 3),
                        // v4ABD2 = p[1][A]*p[2][B]*p[3][B]*p[4][D]
                        VSpec.product(1, 0, 2, 1, 3, 1, 4, 3),
                        // v4ABD3 = p[1][A]*p[2][D]*p[3][D]*p[4][B]
                        VSpec.product(1, 0, 2, 3, 3, 3, 4, 1),

                        // v4ACD1 = p[1][C]*p[2][A]*p[3][A]*p[4][D]
                        VSpec.product(1, 2, 2, 0, 3, 0, 4, 3),
                        // v4ACD2 = p[1][A]*p[2][C]*p[3][C]*p[4][D]
                        VSpec.product(1, 0, 2, 2, 3, 2, 4, 3),
                        // v4ACD3 = p[1][A]*p[2][D]*p[3][D]*p[4][C]
                        VSpec.product(1, 0, 2, 3, 3, 3, 4, 2),

                        // v4BCD1 = p[1][C]*p[2][B]*p[3][B]*p[4][D]
                        VSpec.product(1, 2, 2, 1, 3, 1, 4, 3),
                        // v4BCD2 = p[1][B]*p[2][C]*p[3][C]*p[4][D]
                        VSpec.product(1, 1, 2, 2, 3, 2, 4, 3),
                        // v4BCD3 = p[1][B]*p[2][D]*p[3][D]*p[4][C]
                        VSpec.product(1, 1, 2, 3, 3, 3, 4, 2),

                        // ---- tetrahedra: 3 quaternary ----

                        // v4ABCD1 = p[1][A]*p[2][C]*p[3][D]*p[4][B]
                        VSpec.product(1, 0, 2, 2, 3, 3, 4, 1),
                        // v4ABCD2 = p[1][A]*p[2][B]*p[3][D]*p[4][C]
                        VSpec.product(1, 0, 2, 1, 3, 3, 4, 2),
                        // v4ABCD3 = p[1][A]*p[2][B]*p[3][C]*p[4][D]
                        VSpec.product(1, 0, 2, 1, 3, 2, 4, 3),

                        // ---- triangles: 6 binary (sites p1,p2,p4) ----

                        // v3AB = p[1][B]*p[2][A]*p[4][B] - p[1][A]*p[2][B]*p[4][A]
                        VSpec.diff(new int[] { 1, 1, 2, 0, 4, 1 }, new int[] { 1, 0, 2, 1, 4, 0 }),
                        // v3AC = p[1][C]*p[2][A]*p[4][C] - p[1][A]*p[2][C]*p[4][A]
                        VSpec.diff(new int[] { 1, 2, 2, 0, 4, 2 }, new int[] { 1, 0, 2, 2, 4, 0 }),
                        // v3AD = p[1][D]*p[2][A]*p[4][D] - p[1][A]*p[2][D]*p[4][A]
                        VSpec.diff(new int[] { 1, 3, 2, 0, 4, 3 }, new int[] { 1, 0, 2, 3, 4, 0 }),
                        // v3BC = p[1][C]*p[2][B]*p[4][C] - p[1][B]*p[2][C]*p[4][B]
                        VSpec.diff(new int[] { 1, 2, 2, 1, 4, 2 }, new int[] { 1, 1, 2, 2, 4, 1 }),
                        // v3BD = p[1][D]*p[2][B]*p[4][D] - p[1][B]*p[2][D]*p[4][B]
                        VSpec.diff(new int[] { 1, 3, 2, 1, 4, 3 }, new int[] { 1, 1, 2, 3, 4, 1 }),
                        // v3CD = p[1][D]*p[2][C]*p[4][D] - p[1][C]*p[2][D]*p[4][C]
                        VSpec.diff(new int[] { 1, 3, 2, 2, 4, 3 }, new int[] { 1, 2, 2, 3, 4, 2 }),

                        // ---- triangles: 12 ternary (sites p1,p2,p4) ----

                        // v3ABC1 = p[1][B]*p[2][A]*p[4][C]
                        VSpec.product(1, 1, 2, 0, 4, 2),
                        // v3ABC2 = p[1][A]*p[2][B]*p[4][C]
                        VSpec.product(1, 0, 2, 1, 4, 2),
                        // v3ABC3 = p[1][A]*p[2][C]*p[4][B]
                        VSpec.product(1, 0, 2, 2, 4, 1),

                        // v3ABD1 = p[1][B]*p[2][A]*p[4][D]
                        VSpec.product(1, 1, 2, 0, 4, 3),
                        // v3ABD2 = p[1][A]*p[2][B]*p[4][D]
                        VSpec.product(1, 0, 2, 1, 4, 3),
                        // v3ABD3 = p[1][A]*p[2][D]*p[4][B]
                        VSpec.product(1, 0, 2, 3, 4, 1),

                        // v3ACD1 = p[1][C]*p[2][A]*p[4][D]
                        VSpec.product(1, 2, 2, 0, 4, 3),
                        // v3ACD2 = p[1][A]*p[2][C]*p[4][D]
                        VSpec.product(1, 0, 2, 2, 4, 3),
                        // v3ACD3 = p[1][A]*p[2][D]*p[4][C]
                        VSpec.product(1, 0, 2, 3, 4, 2),

                        // v3BCD1 = p[1][C]*p[2][B]*p[4][D]
                        VSpec.product(1, 2, 2, 1, 4, 3),
                        // v3BCD2 = p[1][B]*p[2][C]*p[4][D]
                        VSpec.product(1, 1, 2, 2, 4, 3),
                        // v3BCD3 = p[1][B]*p[2][D]*p[4][C]
                        VSpec.product(1, 1, 2, 3, 4, 2),

                        // ---- pairs: 6 II-n (p1,p4) ----

                        // v22AB = p[1][A]*p[4][B]
                        VSpec.product(1, 0, 4, 1),
                        // v22AC = p[1][A]*p[4][C]
                        VSpec.product(1, 0, 4, 2),
                        // v22AD = p[1][A]*p[4][D]
                        VSpec.product(1, 0, 4, 3),
                        // v22BC = p[1][B]*p[4][C]
                        VSpec.product(1, 1, 4, 2),
                        // v22BD = p[1][B]*p[4][D]
                        VSpec.product(1, 1, 4, 3),
                        // v22CD = p[1][C]*p[4][D]
                        VSpec.product(1, 2, 4, 3),

                        // ---- pairs: 6 I-n (p1,p2) ----

                        // v21AB = p[1][A]*p[2][B]
                        VSpec.product(1, 0, 2, 1),
                        // v21AC = p[1][A]*p[2][C]
                        VSpec.product(1, 0, 2, 2),
                        // v21AD = p[1][A]*p[2][D]
                        VSpec.product(1, 0, 2, 3),
                        // v21BC = p[1][B]*p[2][C]
                        VSpec.product(1, 1, 2, 2),
                        // v21BD = p[1][B]*p[2][D]
                        VSpec.product(1, 1, 2, 3),
                        // v21CD = p[1][C]*p[2][D]
                        VSpec.product(1, 2, 2, 3),

                        // ---- points ----

                        // xA = p[1][A]
                        VSpec.point(1, 0),
                        // xB = p[1][B]
                        VSpec.point(1, 1),
                        // xC = p[1][C]
                        VSpec.point(1, 2),
                        // xD = p[1][D]
                        VSpec.point(1, 3)));

        // -----------------------------------------------------------------
        // BCC_B2 | T-model | binary (K=2)
        //
        // Same 4-site tetrahedron as BCC_A2, but the two sublattices are
        // now physically distinct (B2 ordering): alpha = {p1, p4},
        // beta = {p2, p3} (per BCC_B2-T.txt's maxClusCoord grouping
        // {{p1,p4}},{{p2,p3}}).
        //
        // Logical site coordinates (fractional), same tetrahedron as A2:
        // p1 = {0.0, 0.0, 0.0}   (alpha)
        // p2 = {0.5, -0.5, 0.5}  (beta)
        // p3 = {0.5, 0.5, 0.5}   (beta)
        // p4 = {1.0, 0.0, 0.0}   (alpha)
        //
        // Orthogonal-CF orbit structure (tcf=8, from Stage 2b):
        //   t=0 tetrahedron          -- unsplit (1 orbit, ncv=9)
        //   t=1 triangle             -- splits: {1,2,3} and {1,2,4}
        //   t=2 II-n pair            -- unsplit as ONE orbit by the
        //                               pipeline's space-group symmetry,
        //                               but {1,4} (alpha-alpha) and {2,3}
        //                               (beta-beta) are physically distinct
        //                               pairs once sublattice occupation is
        //                               tracked, so we pick two independent
        //                               CVCF variables from it.
        //   t=3 I-n pair             -- two orbit entries ({1,2},{1,3}) but
        //                               same physical pair type; ONE CVCF
        //                               variable (unchanged from A2).
        //   t=4 point                -- splits: alpha (site 1), beta (site 2)
        //
        // 9 CVCF variables total (matches tcf+1 = 9):
        //   V4AB    tetrahedron, antisymmetric (same formula as A2's v4AB)
        //   V31AB   triangle {1,2,3} (alpha vertex = site 1)
        //   V32AB   triangle {1,2,4} (beta vertex = site 2)
        //   V221AB  II-n pair, alpha-alpha {1,4}, antisymmetric
        //   V222AB  II-n pair, beta-beta {2,3}, antisymmetric
        //   V21AB   I-n pair {1,2} (unchanged from A2)
        //   xA, xB  bulk composition (alpha-sublattice convention, site 1)
        //   eta     LRO parameter: xA(alpha) - xA(beta) = p[1][A] - p[2][A]
        //
        // Atom indices: A=0, B=1
        // -----------------------------------------------------------------
        register("BCC_B2", "T", 2,
                new double[][] {
                        { 0.0, 0.0, 0.0 }, // p1 (alpha)
                        { 0.5, -0.5, 0.5 }, // p2 (beta)
                        { 0.5, 0.5, 0.5 }, // p3 (beta)
                        { 1.0, 0.0, 0.0 } // p4 (alpha)
                },
                List.of("V4AB", "V31AB", "V32AB", "V221AB", "V222AB", "V21AB", "xA", "xB", "eta"),
                List.of(
                        // V4AB = p[1][A]*p[2][A]*p[3][B]*p[4][B]  (self-paired orbit,
                        // AABB == ABAB == BABA == BBAA under B2's (1<->4)(2<->3)
                        // symmetry; no antisymmetric partner exists for this orbit)
                        VSpec.product(1, 0, 2, 0, 3, 1, 4, 1),

                        // V31AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]  (triangle {1,2,3})
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 1, 2, 0, 3, 0 }),

                        // V32AB = p[2][A]*p[1][B]*p[4][B] - p[2][B]*p[1][A]*p[4][A]  (triangle {1,2,4})
                        VSpec.diff(
                                new int[] { 2, 0, 1, 1, 4, 1 },
                                new int[] { 2, 1, 1, 0, 4, 0 }),

                        // V221AB = p[1][A]*p[4][B]  (II-n pair, alpha-alpha)
                        VSpec.product(1, 0, 4, 1),

                        // V222AB = p[2][A]*p[3][B]  (II-n pair, beta-beta)
                        VSpec.product(2, 0, 3, 1),

                        // V21AB = p[1][A]*p[2][B]  (I-n pair, unchanged from A2)
                        VSpec.product(1, 0, 2, 1),

                        // xA = p[1][A]
                        VSpec.point(1, 0),

                        // xB = p[1][B]
                        VSpec.point(1, 1),

                        // eta = p[1][A] - p[2][A]  (LRO parameter)
                        VSpec.diff(
                                new int[] { 1, 0 },
                                new int[] { 2, 0 })),
                // numPointCfs=3: xA, xB, eta are all point-like (no ECI);
                // eta is a long-range-order parameter, not a species
                // composition, so it doesn't fit the numComponents default.
                3);

        // -----------------------------------------------------------------
        // FCC_A1 | T-model | binary (K=2)
        //
        // Logical site coordinates (fractional), matching FCC_A1-T.txt's
        // tetrahedron maximal cluster:
        // p1 = {0.0, 0.0, 0.0}
        // p2 = {0.0, 0.5, 0.5}
        // p3 = {0.5, 0.0, 0.5}
        // p4 = {0.5, 0.5, 0.0}
        //
        // Unlike BCC_A2's tetrahedron, all 6 pairs of the FCC tetrahedron
        // are symmetry-equivalent (single pair type, no I-n/II-n split).
        //
        // Atom indices: A=0, B=1
        // -----------------------------------------------------------------
        register("FCC_A1", "T", 2,
                new double[][] {
                        { 0.0, 0.0, 0.0 }, // p1
                        { 0.0, 0.5, 0.5 }, // p2
                        { 0.5, 0.0, 0.5 }, // p3
                        { 0.5, 0.5, 0.0 } // p4
                },
                List.of("v4AB", "v3AB", "v2AB", "xA", "xB"),
                List.of(
                        // v4AB = p[1][A]*p[2][A]*p[3][B]*p[4][B]
                        VSpec.product(1, 0, 2, 0, 3, 1, 4, 1),

                        // v3AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 1, 2, 0, 3, 0 }),

                        // v2AB = p[1][A]*p[2][B]
                        VSpec.product(1, 0, 2, 1),

                        // xA = p[1][A]
                        VSpec.point(1, 0),

                        // xB = p[1][B]
                        VSpec.point(1, 1)));

        // -----------------------------------------------------------------
        // FCC_A1 | T-model | ternary (K=3)
        //
        // Same logical site coordinates as binary:
        // p1 = {0.0, 0.0, 0.0}
        // p2 = {0.0, 0.5, 0.5}
        // p3 = {0.5, 0.0, 0.5}
        // p4 = {0.5, 0.5, 0.0}
        //
        // Atom indices: A=0, B=1, C=2
        //
        // 15 CVs total: 6 tetr + 3 tri + 3 pair + 3 point
        // Tetrahedron (6): 3 binary (v4AB,v4AC,v4BC) + 3 ternary
        // (v4ABC1/2/3, each a weighted combination of 4 product terms —
        // NOT a pure product/diff; see VSpec.combo)
        // Triangle (3): v3AB, v3AC, v3BC (single pair type; FCC's
        // tetrahedron has no I-n/II-n split, unlike BCC_A2)
        // Pair (3): v2AB, v2AC, v2BC
        // Point (3): xA, xB, xC
        //
        // v4ABC1/2/3 transcribed verbatim from source, including the first
        // term of each (p[1]*p[2]*p[3], 3 sites only — no p[4] factor,
        // unlike the other 3 terms in the same sum which all include p[4]).
        // -----------------------------------------------------------------
        register("FCC_A1", "T", 3,
                new double[][] {
                        { 0.0, 0.0, 0.0 }, // p1
                        { 0.0, 0.5, 0.5 }, // p2
                        { 0.5, 0.0, 0.5 }, // p3
                        { 0.5, 0.5, 0.0 } // p4
                },
                List.of(
                        "v4AB", "v4AC", "v4BC",
                        "v4ABC1", "v4ABC2", "v4ABC3",
                        "v3AB", "v3AC", "v3BC",
                        "v2AB", "v2AC", "v2BC",
                        "xA", "xB", "xC"),
                List.of(
                        // ---- tetrahedra (binary) ----

                        // v4AB = p[1][A]*p[2][A]*p[3][B]*p[4][B]
                        VSpec.product(1, 0, 2, 0, 3, 1, 4, 1),

                        // v4AC = p[1][A]*p[2][A]*p[3][C]*p[4][C]
                        VSpec.product(1, 0, 2, 0, 3, 2, 4, 2),

                        // v4BC = p[1][B]*p[2][B]*p[3][C]*p[4][C]
                        VSpec.product(1, 1, 2, 1, 3, 2, 4, 2),

                        // ---- tetrahedra (ternary) ----

                        // v4ABC1 = (p1[A]p2[B]p3[C]
                        //           + 2*p1[A]p2[A]p3[B]p4[C]
                        //           - p1[A]p2[B]p3[B]p4[C]
                        //           - p1[A]p2[B]p3[C]p4[C]) / 3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2 },
                                2.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2, 4, 2 }),

                        // v4ABC2 = (p1[A]p2[B]p3[C]
                        //           + 2*p1[A]p2[B]p3[B]p4[C]
                        //           - p1[A]p2[A]p3[B]p4[C]
                        //           - p1[A]p2[B]p3[C]p4[C]) / 3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2 },
                                2.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2, 4, 2 }),

                        // v4ABC3 = (p1[A]p2[B]p3[C]
                        //           + 2*p1[A]p2[B]p3[C]p4[C]
                        //           - p1[A]p2[A]p3[B]p4[C]
                        //           - p1[A]p2[B]p3[B]p4[C]) / 3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2 },
                                2.0 / 3, new int[] { 1, 0, 2, 1, 3, 2, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 2 }),

                        // ---- triangles ----

                        // v3AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 1, 2, 0, 3, 0 }),

                        // v3AC = p[1][A]*p[2][C]*p[3][C] - p[1][C]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 0, 3, 0 }),

                        // v3BC = p[1][B]*p[2][C]*p[3][C] - p[1][C]*p[2][B]*p[3][B]
                        VSpec.diff(
                                new int[] { 1, 1, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 1, 3, 1 }),

                        // ---- pairs ----

                        // v2AB = p[1][A]*p[2][B]
                        VSpec.product(1, 0, 2, 1),

                        // v2AC = p[1][A]*p[2][C]
                        VSpec.product(1, 0, 2, 2),

                        // v2BC = p[1][B]*p[2][C]
                        VSpec.product(1, 1, 2, 2),

                        // ---- points ----

                        // xA = p[1][A]
                        VSpec.point(1, 0),

                        // xB = p[1][B]
                        VSpec.point(1, 1),

                        // xC = p[1][C]
                        VSpec.point(1, 2)));

        // -----------------------------------------------------------------
        // -----------------------------------------------------------------
        // FCC_A1 | T-model | quaternary (K=4)
        //
        // Same logical site coordinates as binary/ternary:
        // p1 = {0.0, 0.0, 0.0}
        // p2 = {0.0, 0.5, 0.5}
        // p3 = {0.5, 0.0, 0.5}
        // p4 = {0.5, 0.5, 0.0}
        //
        // Atom indices: A=0, B=1, C=2, D=3
        //
        // 35 CVs total: 19 tetr (6 bin + 12 ter + 1 quat) + 6 tri + 6 pair + 4 point
        // -----------------------------------------------------------------
        register("FCC_A1", "T", 4,
                new double[][] {
                        { 0.0, 0.0, 0.0 }, // p1
                        { 0.0, 0.5, 0.5 }, // p2
                        { 0.5, 0.0, 0.5 }, // p3
                        { 0.5, 0.5, 0.0 } // p4
                },
                List.of(
                        "v4AB", "v4AC", "v4AD", "v4BC", "v4BD", "v4CD",
                        "v4ABC1", "v4ABC2", "v4ABC3",
                        "v4ABD1", "v4ABD2", "v4ABD3",
                        "v4ACD1", "v4ACD2", "v4ACD3",
                        "v4BCD1", "v4BCD2", "v4BCD3",
                        "v4ABCD",
                        "v3AB", "v3AC", "v3AD", "v3BC", "v3BD", "v3CD",
                        "v2AB", "v2AC", "v2AD", "v2BC", "v2BD", "v2CD",
                        "xA", "xB", "xC", "xD"),
                List.of(
                        // ---- tetrahedra (binary) ----

                        // v4AB = p[1][A]*p[2][A]*p[3][B]*p[4][B]
                        VSpec.product(1, 0, 2, 0, 3, 1, 4, 1),
                        // v4AC = p[1][A]*p[2][A]*p[3][C]*p[4][C]
                        VSpec.product(1, 0, 2, 0, 3, 2, 4, 2),
                        // v4AD = p[1][A]*p[2][A]*p[3][D]*p[4][D]
                        VSpec.product(1, 0, 2, 0, 3, 3, 4, 3),
                        // v4BC = p[1][B]*p[2][B]*p[3][C]*p[4][C]
                        VSpec.product(1, 1, 2, 1, 3, 2, 4, 2),
                        // v4BD = p[1][B]*p[2][B]*p[3][D]*p[4][D]
                        VSpec.product(1, 1, 2, 1, 3, 3, 4, 3),
                        // v4CD = p[1][C]*p[2][C]*p[3][D]*p[4][D]
                        VSpec.product(1, 2, 2, 2, 3, 3, 4, 3),

                        // ---- tetrahedra (ternary: ABC) ----

                        // v4ABC1 = (p1[A]p2[B]p3[C] + 2*p1[A]p2[A]p3[B]p4[C]
                        //           - p1[A]p2[B]p3[B]p4[C] - p1[A]p2[B]p3[C]p4[C]) / 3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2 },
                                2.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2, 4, 2 }),
                        // v4ABC2
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2 },
                                2.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2, 4, 2 }),
                        // v4ABC3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 2 },
                                2.0 / 3, new int[] { 1, 0, 2, 1, 3, 2, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 2 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 2 }),

                        // ---- tetrahedra (ternary: ABD) ----

                        // v4ABD1 = (p1[A]p2[B]p3[D] + 2*p1[A]p2[A]p3[B]p4[D]
                        //           - p1[A]p2[B]p3[B]p4[D] - p1[A]p2[B]p3[D]p4[D]) / 3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 3 },
                                2.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 3, 4, 3 }),
                        // v4ABD2
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 3 },
                                2.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 3, 4, 3 }),
                        // v4ABD3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 1, 3, 3 },
                                2.0 / 3, new int[] { 1, 0, 2, 1, 3, 3, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 1, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 1, 3, 1, 4, 3 }),

                        // ---- tetrahedra (ternary: ACD) ----

                        // v4ACD1 = (p1[A]p2[C]p3[D] + 2*p1[A]p2[A]p3[C]p4[D]
                        //           - p1[A]p2[C]p3[C]p4[D] - p1[A]p2[C]p3[D]p4[D]) / 3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 2, 3, 3 },
                                2.0 / 3, new int[] { 1, 0, 2, 0, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 2, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 2, 3, 3, 4, 3 }),
                        // v4ACD2
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 2, 3, 3 },
                                2.0 / 3, new int[] { 1, 0, 2, 2, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 2, 3, 3, 4, 3 }),
                        // v4ACD3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 0, 2, 2, 3, 3 },
                                2.0 / 3, new int[] { 1, 0, 2, 2, 3, 3, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 0, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 0, 2, 2, 3, 2, 4, 3 }),

                        // ---- tetrahedra (ternary: BCD) ----

                        // v4BCD1 = (p1[B]p2[C]p3[D] + 2*p1[B]p2[B]p3[C]p4[D]
                        //           - p1[B]p2[C]p3[C]p4[D] - p1[B]p2[C]p3[D]p4[D]) / 3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 1, 2, 2, 3, 3 },
                                2.0 / 3, new int[] { 1, 1, 2, 1, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 1, 2, 2, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 1, 2, 2, 3, 3, 4, 3 }),
                        // v4BCD2
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 1, 2, 2, 3, 3 },
                                2.0 / 3, new int[] { 1, 1, 2, 2, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 1, 2, 1, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 1, 2, 2, 3, 3, 4, 3 }),
                        // v4BCD3
                        VSpec.combo(
                                1.0 / 3, new int[] { 1, 1, 2, 2, 3, 3 },
                                2.0 / 3, new int[] { 1, 1, 2, 2, 3, 3, 4, 3 },
                                -1.0 / 3, new int[] { 1, 1, 2, 1, 3, 2, 4, 3 },
                                -1.0 / 3, new int[] { 1, 1, 2, 2, 3, 2, 4, 3 }),

                        // ---- tetrahedron (quaternary) ----

                        // v4ABCD = p[1][A]*p[2][B]*p[3][C]*p[4][D]
                        VSpec.product(1, 0, 2, 1, 3, 2, 4, 3),

                        // ---- triangles ----

                        // v3AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 1, 2, 0, 3, 0 }),
                        // v3AC = p[1][A]*p[2][C]*p[3][C] - p[1][C]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 0, 3, 0 }),
                        // v3AD = p[1][A]*p[2][D]*p[3][D] - p[1][D]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 3, 3, 3 },
                                new int[] { 1, 3, 2, 0, 3, 0 }),
                        // v3BC = p[1][B]*p[2][C]*p[3][C] - p[1][C]*p[2][B]*p[3][B]
                        VSpec.diff(
                                new int[] { 1, 1, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 1, 3, 1 }),
                        // v3BD = p[1][B]*p[2][D]*p[3][D] - p[1][D]*p[2][B]*p[3][B]
                        VSpec.diff(
                                new int[] { 1, 1, 2, 3, 3, 3 },
                                new int[] { 1, 3, 2, 1, 3, 1 }),
                        // v3CD = p[1][C]*p[2][D]*p[3][D] - p[1][D]*p[2][C]*p[3][C]
                        VSpec.diff(
                                new int[] { 1, 2, 2, 3, 3, 3 },
                                new int[] { 1, 3, 2, 2, 3, 2 }),

                        // ---- pairs ----

                        // v2AB = p[1][A]*p[2][B]
                        VSpec.product(1, 0, 2, 1),
                        // v2AC = p[1][A]*p[2][C]
                        VSpec.product(1, 0, 2, 2),
                        // v2AD = p[1][A]*p[2][D]
                        VSpec.product(1, 0, 2, 3),
                        // v2BC = p[1][B]*p[2][C]
                        VSpec.product(1, 1, 2, 2),
                        // v2BD = p[1][B]*p[2][D]
                        VSpec.product(1, 1, 2, 3),
                        // v2CD = p[1][C]*p[2][D]
                        VSpec.product(1, 2, 2, 3),

                        // ---- points ----

                        // xA = p[1][A]
                        VSpec.point(1, 0),
                        // xB = p[1][B]
                        VSpec.point(1, 1),
                        // xC = p[1][C]
                        VSpec.point(1, 2),
                        // xD = p[1][D]
                        VSpec.point(1, 3)));

        // -----------------------------------------------------------------
        // HCP_A3 | T-model | binary (K=2)
        //
        // Site numbering follows Jindal & Lele, Calphad 89 (2025) 102825,
        // Appendix 2 / Fig. 27 exactly:
        //   1,2,3 = tetrahedral basal triangle (z=3/4 plane)
        //   4     = tetrahedron apex (z=1/4 plane) -- SHARED with the
        //           octahedral basal triangle (this is the site common to
        //           both of HCP_A3-T.txt's maximal clusters)
        //   5,6   = remaining octahedral basal triangle sites (z=1/4 plane)
        // Clusters: point(1), in-plane I-n pair(1,2), out-of-plane I-n
        // pair(1,4), tetrahedral basal triangle(1,2,3), octahedral basal
        // triangle(4,5,6), out-of-plane pyramidal triangle(1,2,4),
        // tetrahedron(1,2,3,4).
        //
        // p1 = { 2/3, 4/3, 3/4}
        // p2 = { 2/3, 1/3, 3/4}
        // p3 = {-1/3, 1/3, 3/4}
        // p4 = { 1/3, 2/3, 1/4}
        // p5 = {-2/3, 2/3, 1/4}
        // p6 = { 1/3, 5/3, 1/4}
        //
        // Atom indices: A=0, B=1
        //
        // 7 CVs total (binary, per paper Table 19): 1 tetr + 3 tri + 2 pair
        // + 1 point set (xA, xB, both retained per point-CF convention).
        // v2PR1 = in-plane pair (1,2); v2PR2 = out-of-plane pair (1,4);
        // v3PR1 = tetrahedral basal triangle (1,2,3);
        // v3PR2 = octahedral basal triangle (4,5,6);
        // v3PR3 = out-of-plane pyramidal triangle (1,2,4);
        // v4PR  = tetrahedron (1,2,3,4).
        // -----------------------------------------------------------------
        register("HCP_A3", "T", 2,
                new double[][] {
                        { 2.0 / 3, 4.0 / 3, 3.0 / 4 }, // p1
                        { 2.0 / 3, 1.0 / 3, 3.0 / 4 }, // p2
                        { -1.0 / 3, 1.0 / 3, 3.0 / 4 }, // p3
                        { 1.0 / 3, 2.0 / 3, 1.0 / 4 }, // p4
                        { -2.0 / 3, 2.0 / 3, 1.0 / 4 }, // p5
                        { 1.0 / 3, 5.0 / 3, 1.0 / 4 } // p6
                },
                List.of("v4AB", "v3AB1", "v3AB2", "v3AB3", "v2AB1", "v2AB2", "xA", "xB"),
                List.of(
                        // v4AB = (p[1][A]p[2][B]p[3][B]p[4][A] + p[1][B]p[2][A]p[3][A]p[4][B]) / 2
                        //   (eq. 64: v4PR = (yPRRP + yPPRR)/2, sites 1,2,3,4)
                        VSpec.combo(
                                0.5, new int[] { 1, 0, 2, 1, 3, 1, 4, 0 },
                                0.5, new int[] { 1, 1, 2, 1, 3, 0, 4, 0 }),

                        // v3AB1 (tetrahedral basal, 1-2-3) = yPRR1 - yPPR1
                        //   = p[1][A]p[2][B]p[3][B] - p[1][A]p[2][A]p[3][B]  (eq. 61)
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 0, 2, 0, 3, 1 }),

                        // v3AB2 (octahedral basal, 4-5-6) = yPRR2 - yPPR2
                        //   = p[4][A]p[5][B]p[6][B] - p[4][A]p[5][A]p[6][B]  (eq. 62)
                        VSpec.diff(
                                new int[] { 4, 0, 5, 1, 6, 1 },
                                new int[] { 4, 0, 5, 0, 6, 1 }),

                        // v3AB3 (out-of-plane pyramidal, 1-2-4) = yRRP3 - yPPR3
                        //   = p[1][B]p[2][B]p[4][A] - p[1][A]p[2][A]p[4][B]  (eq. 63)
                        VSpec.diff(
                                new int[] { 1, 1, 2, 1, 4, 0 },
                                new int[] { 1, 0, 2, 0, 4, 1 }),

                        // v2AB1 (in-plane pair, 1-2) = p[1][A]*p[2][B]
                        VSpec.product(1, 0, 2, 1),

                        // v2AB2 (out-of-plane pair, 1-4) = p[1][A]*p[4][B]
                        VSpec.product(1, 0, 4, 1),

                        // xA = p[1][A]
                        VSpec.point(1, 0),

                        // xB = p[1][B]
                        VSpec.point(1, 1)));

        // -----------------------------------------------------------------
        // HCP_A3 | T-model | ternary (K=3)
        //
        // Same 6 logical sites/coordinates as binary (see the K=2 block
        // above for the maximal-cluster geometry and site numbering).
        //
        // Atom indices: A=0, B=1, C=2
        //
        // 31 CVs total: 8 tetr + 14 tri + 6 pair + 3 point
        // -----------------------------------------------------------------
        register("HCP_A3", "T", 3,
                new double[][] {
                        { 2.0 / 3, 4.0 / 3, 3.0 / 4 }, // p1
                        { 2.0 / 3, 1.0 / 3, 3.0 / 4 }, // p2
                        { -1.0 / 3, 1.0 / 3, 3.0 / 4 }, // p3
                        { 1.0 / 3, 2.0 / 3, 1.0 / 4 }, // p4
                        { -2.0 / 3, 2.0 / 3, 1.0 / 4 }, // p5
                        { 1.0 / 3, 5.0 / 3, 1.0 / 4 } // p6
                },
                List.of(
                        "v4AB", "v4AC", "v4BC",
                        "v4ABC1", "v4ABC2", "v4ABC3", "v4ABC4", "v4ABC5",
                        "v31AB", "v31AC", "v31BC", "v31ABC",
                        "v32AB", "v32AC", "v32BC", "v32ABC",
                        "v33AB", "v33AC", "v33BC", "v33ABC1", "v33ABC2", "v33ABC3",
                        "v21AB", "v21AC", "v21BC",
                        "v22AB", "v22AC", "v22BC",
                        "xA", "xB", "xC"),
                List.of(
                        // ---- tetrahedra ----
                        // Literal translation of Mathematica tetrcvs (1=A,2=B,3=C
                        // there -> 0=A,1=B,2=C here). Plain products, not combos.
                        // Mathematica site labels remapped to our binary-verified
                        // convention: math 1->1, 2->4, 3->2, 4->3, 5->5, 6->6.

                        // v4AB = p[1][B]*p[4][A]*p[2][B]*p[3][B]
                        VSpec.product(1, 1, 4, 0, 2, 1, 3, 1),
                        // v4AC = p[1][C]*p[4][A]*p[2][C]*p[3][C]
                        VSpec.product(1, 2, 4, 0, 2, 2, 3, 2),
                        // v4BC = p[1][C]*p[4][B]*p[2][C]*p[3][C]
                        VSpec.product(1, 2, 4, 1, 2, 2, 3, 2),
                        // v4ABC1 = p[1][B]*p[4][A]*p[2][B]*p[3][C]
                        VSpec.product(1, 1, 4, 0, 2, 1, 3, 2),
                        // v4ABC2 = p[1][A]*p[4][B]*p[2][A]*p[3][C]
                        VSpec.product(1, 0, 4, 1, 2, 0, 3, 2),
                        // v4ABC3 = p[1][C]*p[4][B]*p[2][A]*p[3][C]
                        VSpec.product(1, 2, 4, 1, 2, 0, 3, 2),
                        // v4ABC4 = p[1][A]*p[4][C]*p[2][A]*p[3][B]
                        VSpec.product(1, 0, 4, 2, 2, 0, 3, 1),
                        // v4ABC5 = p[1][A]*p[4][C]*p[2][B]*p[3][B]
                        VSpec.product(1, 0, 4, 2, 2, 1, 3, 1),

                        // ---- triangles: basal-octahedral (sites 4,5,6) ----

                        // v31AB = p[4][A]*p[5][B]*p[6][B] - p[4][B]*p[5][A]*p[6][A]
                        VSpec.diff(
                                new int[] { 4, 0, 5, 1, 6, 1 },
                                new int[] { 4, 1, 5, 0, 6, 0 }),
                        // v31AC = p[4][A]*p[5][C]*p[6][C] - p[4][C]*p[5][A]*p[6][A]
                        VSpec.diff(
                                new int[] { 4, 0, 5, 2, 6, 2 },
                                new int[] { 4, 2, 5, 0, 6, 0 }),
                        // v31BC = p[4][B]*p[5][C]*p[6][C] - p[4][C]*p[5][B]*p[6][B]
                        VSpec.diff(
                                new int[] { 4, 1, 5, 2, 6, 2 },
                                new int[] { 4, 2, 5, 1, 6, 1 }),
                        // v31ABC = p[4][A]*p[5][B]*p[6][C]
                        VSpec.product(4, 0, 5, 1, 6, 2),

                        // ---- triangles: basal-tetrahedral (sites 1,2,3) ----

                        // v32AB = p[1][A]*p[2][B]*p[3][B] - p[1][B]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 1, 3, 1 },
                                new int[] { 1, 1, 2, 0, 3, 0 }),
                        // v32AC = p[1][A]*p[2][C]*p[3][C] - p[1][C]*p[2][A]*p[3][A]
                        VSpec.diff(
                                new int[] { 1, 0, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 0, 3, 0 }),
                        // v32BC = p[1][B]*p[2][C]*p[3][C] - p[1][C]*p[2][B]*p[3][B]
                        VSpec.diff(
                                new int[] { 1, 1, 2, 2, 3, 2 },
                                new int[] { 1, 2, 2, 1, 3, 1 }),
                        // v32ABC = p[1][A]*p[2][B]*p[3][C]
                        VSpec.product(1, 0, 2, 1, 3, 2),

                        // ---- triangles: lateral (sites 4,1,2) ----

                        // v33AB = p[4][A]*p[1][B]*p[2][B] - p[4][B]*p[1][A]*p[2][A]
                        VSpec.diff(
                                new int[] { 4, 0, 1, 1, 2, 1 },
                                new int[] { 4, 1, 1, 0, 2, 0 }),
                        // v33AC = p[4][A]*p[1][C]*p[2][C] - p[4][C]*p[1][A]*p[2][A]
                        VSpec.diff(
                                new int[] { 4, 0, 1, 2, 2, 2 },
                                new int[] { 4, 2, 1, 0, 2, 0 }),
                        // v33BC = p[4][B]*p[1][C]*p[2][C] - p[4][C]*p[1][B]*p[2][B]
                        VSpec.diff(
                                new int[] { 4, 1, 1, 2, 2, 2 },
                                new int[] { 4, 2, 1, 1, 2, 1 }),
                        // v33ABC1 = p[4][A]*p[1][B]*p[2][C]
                        VSpec.product(4, 0, 1, 1, 2, 2),
                        // v33ABC2 = p[4][B]*p[1][A]*p[2][C]
                        VSpec.product(4, 1, 1, 0, 2, 2),
                        // v33ABC3 = p[4][C]*p[1][A]*p[2][B]
                        VSpec.product(4, 2, 1, 0, 2, 1),

                        // ---- pairs: in-plane (sites 1,2) ----

                        // v21AB = p[1][A]*p[2][B]
                        VSpec.product(1, 0, 2, 1),
                        // v21AC = p[1][A]*p[2][C]
                        VSpec.product(1, 0, 2, 2),
                        // v21BC = p[1][B]*p[2][C]
                        VSpec.product(1, 1, 2, 2),

                        // ---- pairs: out-of-plane (sites 1,4) ----

                        // v22AB = p[1][A]*p[4][B]
                        VSpec.product(1, 0, 4, 1),
                        // v22AC = p[1][A]*p[4][C]
                        VSpec.product(1, 0, 4, 2),
                        // v22BC = p[1][B]*p[4][C]
                        VSpec.product(1, 1, 4, 2),

                        // ---- points ----

                        // xA = p[1][A]
                        VSpec.point(1, 0),
                        // xB = p[1][B]
                        VSpec.point(1, 1),
                        // xC = p[1][C]
                        VSpec.point(1, 2)));

        // -----------------------------------------------------------------
        // HCP_A3 | T-model | quaternary (K=4)
        //
        // Same 6 logical sites/coordinates as binary/ternary (see the K=2
        // block above for the maximal-cluster geometry and site numbering).
        //
        // Atom indices: A=0, B=1, C=2, D=3
        //
        // 84 CVs total: 30 tetr (6 bin + 20 tern + 4 quat) + 38 tri
        // (10 basal-octa + 10 basal-tetr + 18 lateral) + 12 pair + 4 point
        // -----------------------------------------------------------------
        register("HCP_A3", "T", 4,
                new double[][] {
                        { 2.0 / 3, 4.0 / 3, 3.0 / 4 }, // p1
                        { 2.0 / 3, 1.0 / 3, 3.0 / 4 }, // p2
                        { -1.0 / 3, 1.0 / 3, 3.0 / 4 }, // p3
                        { 1.0 / 3, 2.0 / 3, 1.0 / 4 }, // p4
                        { -2.0 / 3, 2.0 / 3, 1.0 / 4 }, // p5
                        { 1.0 / 3, 5.0 / 3, 1.0 / 4 } // p6
                },
                List.of(
                        "v4AB", "v4AC", "v4AD", "v4BC", "v4BD", "v4CD",
                        "v4ABC1", "v4ABC2", "v4ABC3", "v4ABC4", "v4ABC5",
                        "v4ABD1", "v4ABD2", "v4ABD3", "v4ABD4", "v4ABD5",
                        "v4ACD1", "v4ACD2", "v4ACD3", "v4ACD4", "v4ACD5",
                        "v4BCD1", "v4BCD2", "v4BCD3", "v4BCD4", "v4BCD5",
                        "v4ABCD1", "v4ABCD2", "v4ABCD3", "v4ABCD4",
                        "v31AB", "v31AC", "v31AD", "v31BC", "v31BD", "v31CD",
                        "v31ABC", "v31ABD", "v31ACD", "v31BCD",
                        "v32AB", "v32AC", "v32AD", "v32BC", "v32BD", "v32CD",
                        "v32ABC", "v32ABD", "v32ACD", "v32BCD",
                        "v33AB", "v33AC", "v33AD", "v33BC", "v33BD", "v33CD",
                        "v33ABC1", "v33ABC2", "v33ABC3",
                        "v33ABD1", "v33ABD2", "v33ABD3",
                        "v33ACD1", "v33ACD2", "v33ACD3",
                        "v33BCD1", "v33BCD2", "v33BCD3",
                        "v21AB", "v21AC", "v21AD", "v21BC", "v21BD", "v21CD",
                        "v22AB", "v22AC", "v22AD", "v22BC", "v22BD", "v22CD",
                        "xA", "xB", "xC", "xD"),
                List.of(
                        // ---- tetrahedra (binary) ----

                        // v4AB = p[1][B]*p[2][A]*p[3][B]*p[4][B]
                        VSpec.product(1, 1, 4, 0, 2, 1, 3, 1),
                        // v4AC = p[1][C]*p[2][A]*p[3][C]*p[4][C]
                        VSpec.product(1, 2, 4, 0, 2, 2, 3, 2),
                        // v4AD = p[1][D]*p[2][A]*p[3][D]*p[4][D]
                        VSpec.product(1, 3, 4, 0, 2, 3, 3, 3),
                        // v4BC = p[1][C]*p[2][B]*p[3][C]*p[4][C]
                        VSpec.product(1, 2, 4, 1, 2, 2, 3, 2),
                        // v4BD = p[1][D]*p[2][B]*p[3][D]*p[4][D]
                        VSpec.product(1, 3, 4, 1, 2, 3, 3, 3),
                        // v4CD = p[1][D]*p[2][C]*p[3][D]*p[4][D]
                        VSpec.product(1, 3, 4, 2, 2, 3, 3, 3),

                        // ---- tetrahedra (ternary: ABC) ----

                        // v4ABC1 = p[1][B]*p[2][A]*p[3][B]*p[4][C]
                        VSpec.product(1, 1, 4, 0, 2, 1, 3, 2),
                        // v4ABC2 = p[1][A]*p[2][B]*p[3][A]*p[4][C]
                        VSpec.product(1, 0, 4, 1, 2, 0, 3, 2),
                        // v4ABC3 = p[1][C]*p[2][B]*p[3][A]*p[4][C]
                        VSpec.product(1, 2, 4, 1, 2, 0, 3, 2),
                        // v4ABC4 = p[1][A]*p[2][C]*p[3][A]*p[4][B]
                        VSpec.product(1, 0, 4, 2, 2, 0, 3, 1),
                        // v4ABC5 = p[1][A]*p[2][C]*p[3][B]*p[4][B]
                        VSpec.product(1, 0, 4, 2, 2, 1, 3, 1),

                        // ---- tetrahedra (ternary: ABD) ----

                        // v4ABD1 = p[1][B]*p[2][A]*p[3][B]*p[4][D]
                        VSpec.product(1, 1, 4, 0, 2, 1, 3, 3),
                        // v4ABD2 = p[1][A]*p[2][B]*p[3][A]*p[4][D]
                        VSpec.product(1, 0, 4, 1, 2, 0, 3, 3),
                        // v4ABD3 = p[1][D]*p[2][B]*p[3][A]*p[4][D]
                        VSpec.product(1, 3, 4, 1, 2, 0, 3, 3),
                        // v4ABD4 = p[1][A]*p[2][D]*p[3][A]*p[4][B]
                        VSpec.product(1, 0, 4, 3, 2, 0, 3, 1),
                        // v4ABD5 = p[1][A]*p[2][D]*p[3][B]*p[4][B]
                        VSpec.product(1, 0, 4, 3, 2, 1, 3, 1),

                        // ---- tetrahedra (ternary: ACD) ----

                        // v4ACD1 = p[1][C]*p[2][A]*p[3][C]*p[4][D]
                        VSpec.product(1, 2, 4, 0, 2, 2, 3, 3),
                        // v4ACD2 = p[1][A]*p[2][C]*p[3][A]*p[4][D]
                        VSpec.product(1, 0, 4, 2, 2, 0, 3, 3),
                        // v4ACD3 = p[1][D]*p[2][C]*p[3][A]*p[4][D]
                        VSpec.product(1, 3, 4, 2, 2, 0, 3, 3),
                        // v4ACD4 = p[1][A]*p[2][D]*p[3][A]*p[4][C]
                        VSpec.product(1, 0, 4, 3, 2, 0, 3, 2),
                        // v4ACD5 = p[1][A]*p[2][D]*p[3][C]*p[4][C]
                        VSpec.product(1, 0, 4, 3, 2, 2, 3, 2),

                        // ---- tetrahedra (ternary: BCD) ----

                        // v4BCD1 = p[1][C]*p[2][B]*p[3][C]*p[4][D]
                        VSpec.product(1, 2, 4, 1, 2, 2, 3, 3),
                        // v4BCD2 = p[1][B]*p[2][C]*p[3][B]*p[4][D]
                        VSpec.product(1, 1, 4, 2, 2, 1, 3, 3),
                        // v4BCD3 = p[1][D]*p[2][C]*p[3][B]*p[4][D]
                        VSpec.product(1, 3, 4, 2, 2, 1, 3, 3),
                        // v4BCD4 = p[1][B]*p[2][D]*p[3][B]*p[4][C]
                        VSpec.product(1, 1, 4, 3, 2, 1, 3, 2),
                        // v4BCD5 = p[1][B]*p[2][D]*p[3][C]*p[4][C]
                        VSpec.product(1, 1, 4, 3, 2, 2, 3, 2),

                        // ---- tetrahedra (quaternary) ----

                        // v4ABCD1 = p[1][B]*p[2][A]*p[3][C]*p[4][D]
                        VSpec.product(1, 1, 4, 0, 2, 2, 3, 3),
                        // v4ABCD2 = p[1][A]*p[2][B]*p[3][C]*p[4][D]
                        VSpec.product(1, 0, 4, 1, 2, 2, 3, 3),
                        // v4ABCD3 = p[1][A]*p[2][C]*p[3][B]*p[4][D]
                        VSpec.product(1, 0, 4, 2, 2, 1, 3, 3),
                        // v4ABCD4 = p[1][A]*p[2][D]*p[3][B]*p[4][C]
                        VSpec.product(1, 0, 4, 3, 2, 1, 3, 2),

                        // ---- triangles: basal-octahedral (sites 2,5,6) — binary ----

                        // v31AB = p[2][A]*p[5][B]*p[6][B] - p[2][B]*p[5][A]*p[6][A]
                        VSpec.diff(new int[] { 4, 0, 5, 1, 6, 1 }, new int[] { 4, 1, 5, 0, 6, 0 }),
                        // v31AC = p[2][A]*p[5][C]*p[6][C] - p[2][C]*p[5][A]*p[6][A]
                        VSpec.diff(new int[] { 4, 0, 5, 2, 6, 2 }, new int[] { 4, 2, 5, 0, 6, 0 }),
                        // v31AD = p[2][A]*p[5][D]*p[6][D] - p[2][D]*p[5][A]*p[6][A]
                        VSpec.diff(new int[] { 4, 0, 5, 3, 6, 3 }, new int[] { 4, 3, 5, 0, 6, 0 }),
                        // v31BC = p[2][B]*p[5][C]*p[6][C] - p[2][C]*p[5][B]*p[6][B]
                        VSpec.diff(new int[] { 4, 1, 5, 2, 6, 2 }, new int[] { 4, 2, 5, 1, 6, 1 }),
                        // v31BD = p[2][B]*p[5][D]*p[6][D] - p[2][D]*p[5][B]*p[6][B]
                        VSpec.diff(new int[] { 4, 1, 5, 3, 6, 3 }, new int[] { 4, 3, 5, 1, 6, 1 }),
                        // v31CD = p[2][C]*p[5][D]*p[6][D] - p[2][D]*p[5][C]*p[6][C]
                        VSpec.diff(new int[] { 4, 2, 5, 3, 6, 3 }, new int[] { 4, 3, 5, 2, 6, 2 }),

                        // ---- triangles: basal-octahedral — ternary ----

                        // v31ABC = p[2][A]*p[5][B]*p[6][C]
                        VSpec.product(4, 0, 5, 1, 6, 2),
                        // v31ABD = p[2][A]*p[5][B]*p[6][D]
                        VSpec.product(4, 0, 5, 1, 6, 3),
                        // v31ACD = p[2][A]*p[5][C]*p[6][D]
                        VSpec.product(4, 0, 5, 2, 6, 3),
                        // v31BCD = p[2][B]*p[5][C]*p[6][D]
                        VSpec.product(4, 1, 5, 2, 6, 3),

                        // ---- triangles: basal-tetrahedral (sites 1,3,4) — binary ----

                        // v32AB = p[1][A]*p[3][B]*p[4][B] - p[1][B]*p[3][A]*p[4][A]
                        VSpec.diff(new int[] { 1, 0, 2, 1, 3, 1 }, new int[] { 1, 1, 2, 0, 3, 0 }),
                        // v32AC = p[1][A]*p[3][C]*p[4][C] - p[1][C]*p[3][A]*p[4][A]
                        VSpec.diff(new int[] { 1, 0, 2, 2, 3, 2 }, new int[] { 1, 2, 2, 0, 3, 0 }),
                        // v32AD = p[1][A]*p[3][D]*p[4][D] - p[1][D]*p[3][A]*p[4][A]
                        VSpec.diff(new int[] { 1, 0, 2, 3, 3, 3 }, new int[] { 1, 3, 2, 0, 3, 0 }),
                        // v32BC = p[1][B]*p[3][C]*p[4][C] - p[1][C]*p[3][B]*p[4][B]
                        VSpec.diff(new int[] { 1, 1, 2, 2, 3, 2 }, new int[] { 1, 2, 2, 1, 3, 1 }),
                        // v32BD = p[1][B]*p[3][D]*p[4][D] - p[1][D]*p[3][B]*p[4][B]
                        VSpec.diff(new int[] { 1, 1, 2, 3, 3, 3 }, new int[] { 1, 3, 2, 1, 3, 1 }),
                        // v32CD = p[1][C]*p[3][D]*p[4][D] - p[1][D]*p[3][C]*p[4][C]
                        VSpec.diff(new int[] { 1, 2, 2, 3, 3, 3 }, new int[] { 1, 3, 2, 2, 3, 2 }),

                        // ---- triangles: basal-tetrahedral — ternary ----

                        // v32ABC = p[1][A]*p[3][B]*p[4][C]
                        VSpec.product(1, 0, 2, 1, 3, 2),
                        // v32ABD = p[1][A]*p[3][B]*p[4][D]
                        VSpec.product(1, 0, 2, 1, 3, 3),
                        // v32ACD = p[1][A]*p[3][C]*p[4][D]
                        VSpec.product(1, 0, 2, 2, 3, 3),
                        // v32BCD = p[1][B]*p[3][C]*p[4][D]
                        VSpec.product(1, 1, 2, 2, 3, 3),

                        // ---- triangles: lateral (sites 1,2,3) — binary ----

                        // v33AB = p[2][A]*p[1][B]*p[3][B] - p[2][B]*p[1][A]*p[3][A]
                        VSpec.diff(new int[] { 4, 0, 1, 1, 2, 1 }, new int[] { 4, 1, 1, 0, 2, 0 }),
                        // v33AC = p[2][A]*p[1][C]*p[3][C] - p[2][C]*p[1][A]*p[3][A]
                        VSpec.diff(new int[] { 4, 0, 1, 2, 2, 2 }, new int[] { 4, 2, 1, 0, 2, 0 }),
                        // v33AD = p[2][A]*p[1][D]*p[3][D] - p[2][D]*p[1][A]*p[3][A]
                        VSpec.diff(new int[] { 4, 0, 1, 3, 2, 3 }, new int[] { 4, 3, 1, 0, 2, 0 }),
                        // v33BC = p[2][B]*p[1][C]*p[3][C] - p[2][C]*p[1][B]*p[3][B]
                        VSpec.diff(new int[] { 4, 1, 1, 2, 2, 2 }, new int[] { 4, 2, 1, 1, 2, 1 }),
                        // v33BD = p[2][B]*p[1][D]*p[3][D] - p[2][D]*p[1][B]*p[3][B]
                        VSpec.diff(new int[] { 4, 1, 1, 3, 2, 3 }, new int[] { 4, 3, 1, 1, 2, 1 }),
                        // v33CD = p[2][C]*p[1][D]*p[3][D] - p[2][D]*p[1][C]*p[3][C]
                        VSpec.diff(new int[] { 4, 2, 1, 3, 2, 3 }, new int[] { 4, 3, 1, 2, 2, 2 }),

                        // ---- triangles: lateral — ternary ----

                        // v33ABC1 = p[2][A]*p[1][B]*p[3][C]
                        VSpec.product(4, 0, 1, 1, 2, 2),
                        // v33ABC2 = p[2][B]*p[1][A]*p[3][C]
                        VSpec.product(4, 1, 1, 0, 2, 2),
                        // v33ABC3 = p[2][C]*p[1][A]*p[3][B]
                        VSpec.product(4, 2, 1, 0, 2, 1),
                        // v33ABD1 = p[2][A]*p[1][B]*p[3][D]
                        VSpec.product(4, 0, 1, 1, 2, 3),
                        // v33ABD2 = p[2][B]*p[1][A]*p[3][D]
                        VSpec.product(4, 1, 1, 0, 2, 3),
                        // v33ABD3 = p[2][D]*p[1][A]*p[3][B]
                        VSpec.product(4, 3, 1, 0, 2, 1),
                        // v33ACD1 = p[2][A]*p[1][C]*p[3][D]
                        VSpec.product(4, 0, 1, 2, 2, 3),
                        // v33ACD2 = p[2][C]*p[1][A]*p[3][D]
                        VSpec.product(4, 2, 1, 0, 2, 3),
                        // v33ACD3 = p[2][D]*p[1][A]*p[3][C]
                        VSpec.product(4, 3, 1, 0, 2, 2),
                        // v33BCD1 = p[2][B]*p[1][C]*p[3][D]
                        VSpec.product(4, 1, 1, 2, 2, 3),
                        // v33BCD2 = p[2][C]*p[1][B]*p[3][D]
                        VSpec.product(4, 2, 1, 1, 2, 3),
                        // v33BCD3 = p[2][D]*p[1][B]*p[3][C]
                        VSpec.product(4, 3, 1, 1, 2, 2),

                        // ---- pairs: in-plane (sites 1,3) ----

                        // v21AB = p[1][A]*p[3][B]
                        VSpec.product(1, 0, 2, 1),
                        // v21AC = p[1][A]*p[3][C]
                        VSpec.product(1, 0, 2, 2),
                        // v21AD = p[1][A]*p[3][D]
                        VSpec.product(1, 0, 2, 3),
                        // v21BC = p[1][B]*p[3][C]
                        VSpec.product(1, 1, 2, 2),
                        // v21BD = p[1][B]*p[3][D]
                        VSpec.product(1, 1, 2, 3),
                        // v21CD = p[1][C]*p[3][D]
                        VSpec.product(1, 2, 2, 3),

                        // ---- pairs: out-of-plane (sites 1,2) ----

                        // v22AB = p[1][A]*p[2][B]
                        VSpec.product(1, 0, 4, 1),
                        // v22AC = p[1][A]*p[2][C]
                        VSpec.product(1, 0, 4, 2),
                        // v22AD = p[1][A]*p[2][D]
                        VSpec.product(1, 0, 4, 3),
                        // v22BC = p[1][B]*p[2][C]
                        VSpec.product(1, 1, 4, 2),
                        // v22BD = p[1][B]*p[2][D]
                        VSpec.product(1, 1, 4, 3),
                        // v22CD = p[1][C]*p[2][D]
                        VSpec.product(1, 2, 4, 3),

                        // ---- points ----

                        // xA = p[1][A]
                        VSpec.point(1, 0),
                        // xB = p[1][B]
                        VSpec.point(1, 1),
                        // xC = p[1][C]
                        VSpec.point(1, 2),
                        // xD = p[1][D]
                        VSpec.point(1, 3)));
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

        emit(sink, "\n  [Basis-Gen] Building CVCF transformation...");

        int numComponents = pr.getNumComponents();
        CFIdentificationResult cfResult = pr.toCFIdentificationResult();

        // ---------------------------------------------------------
        // Phase 1: Load Specification
        // ---------------------------------------------------------
        emit(sink, "    - Phase 1: Load Specification");
        Definition def = REGISTRY.get(structurePhase + "_" + model.toUpperCase() + "_" + numComponents);
        if (def == null)
            throw new IllegalArgumentException("Unregistered CVCF combination.");
        emit(sink, "  Loaded " + def.cfNames.size() + " CVs for " + structurePhase + " K=" + numComponents);

        for (int i = 0; i < def.vSpecs.size(); i++) {
            emit(sink, String.format("  v[%d] (%s) = %s", i, def.cfNames.get(i), def.vSpecs.get(i).toString()));
        }

        // ---------------------------------------------------------
        // Phase 2: Resolve Coordinates
        // ---------------------------------------------------------
        emit(sink, "    - Phase 2: Resolve Coordinates");
        List<Position> siteList = matrixData.getSiteList();
        Map<Integer, Integer> siteMap = resolveSiteMap(def.logicalSiteCoords, siteList);

        for (Map.Entry<Integer, Integer> entry : siteMap.entrySet()) {
            emit(sink, "  Logical p[" + entry.getKey() + "] mapping to physical internal site #" + entry.getValue());
        }

        // ---------------------------------------------------------
        // Phase 3: Orthogonal CF Data prep
        // ---------------------------------------------------------
        emit(sink, "    - Phase 3: Orthogonal CF Metadata");
        int totalCfs = cfResult.getTcf();
        int basisSize = def.vSpecs.size();

        if (basisSize != totalCfs + 1) {
            LOG.warning(String.format("Basis size mismatch: CVCF=%d, Orthogonal=%d (tcf+1).", basisSize, totalCfs + 1));
            emit(sink, String.format("  [WARN] Size mismatch! CVCF=%d, Orth=(%d+1)", basisSize, totalCfs));
        }

        Map<CFIndex, Integer> cfColMap = buildCfColumnMap(cfResult.getLcf());
        emit(sink,
                "  Built CF column map for " + cfColMap.size() + " unique orthogonal physical correlation functions.");

        // ---------------------------------------------------------
        // Phase 4: Matrix M construction via Pipeline
        // ---------------------------------------------------------
        emit(sink, "    - Phase 4: Build M-Matrix");
        double[][] M = buildMMatrix(def.vSpecs, siteMap, matrixData, cfColMap, totalCfs, basisSize, sink);

        // ---------------------------------------------------------
        // Phase 5: Inverting to T-matrix
        // ---------------------------------------------------------
        emit(sink, "    - Phase 5: Inverting to T-Matrix");
        double[][] T = LinearAlgebra.invert(M);
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
            if (!hasTerms)
                line.append("[0.000]");
            emit(sink, line.toString());
        }
        emit(sink, "  Geometric Inversion complete, T-matrix bound (" + T.length + "x" + T[0].length + ")");

        // ---------------------------------------------------------
        // Phase 6: Transforming Orthogonal C-Matrix to CVCF Basis
        // ---------------------------------------------------------
        if (matrixData.getSiteList().isEmpty()) {
            throw new IllegalStateException("Cannot transform C-Matrix: site list is empty.");
        }
        emit(sink, "    - Phase 6: Coordinate Transformation");
        CMatrixPipeline.CMatrixData cvcfData = matrixData.transform(T);
        emit(sink, "  C-Matrix natively bound inside Basis mapping.");

        int numNonPointCfs = def.cfNames.size() - def.numPointCfs;
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

        List<String> cfDefinitions = new ArrayList<>();
        for (VSpec spec : def.vSpecs) {
            cfDefinitions.add(spec.toString());
        }

        emit(sink, "  [Basis-Gen] ✓ Basis generation complete.");

        return new CvCfBasis(structurePhase, model, numComponents,
                def.cfNames, eciNames, numNonPointCfs, T, Tinv, cvcfData, cfDefinitions);
    }

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null)
            sink.accept(msg);
    }

    // =========================================================================
    // Site coordinate matching
    // =========================================================================

    /**
     * Maps each logical site index (1-based, as used in {@code vSpecs}) to
     * its physical index in {@code siteList} — the global set of unique
     * site positions across ALL of the structure's maximal clusters (see
     * {@link CMatrixPipeline#buildSiteList}). Matching against the full
     * global list (rather than a single representative ordered cluster)
     * is required for structures like HCP_A3 whose logical sites span more
     * than one maximal cluster (e.g. a 4-site tetrahedron plus a disjoint
     * 3-site triangle) — no single cluster contains every logical site.
     */
    private static Map<Integer, Integer> resolveSiteMap(
            double[][] logicalSiteCoords,
            List<Position> siteList) {

        Map<Integer, Integer> siteMap = new LinkedHashMap<>();
        for (int logIdx = 0; logIdx < logicalSiteCoords.length; logIdx++) {
            double[] coord = logicalSiteCoords[logIdx];
            Position matched = findMatchingPosition(coord, siteList);
            int physIdx = indexOf(matched, siteList);
            siteMap.put(logIdx + 1, physIdx);
        }
        return siteMap;
    }

    private static Position findMatchingPosition(double[] coord, List<Position> positions) {
        for (Position p : positions) {
            double dx = p.getX() - coord[0];
            double dy = p.getY() - coord[1];
            double dz = p.getZ() - coord[2];
            if (Math.sqrt(dx * dx + dy * dy + dz * dz) < 1e-4) {
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
            if (siteList.get(i).equals(pos))
                return i;
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
            double[] row = new double[basisSize];
            for (VSpec.Term term : spec.terms) {
                double[] termRow = evaluateSpecTerm(term.siteAtomPairs, siteMap, matrixData, cfColMap, totalCfs);
                for (int j = 0; j < basisSize; j++) {
                    row[j] += term.coefficient * termRow[j];
                }
            }
            M[i] = row;
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
            Integer physicalSite = siteMap.get(logicalSite);
            if (physicalSite == null) {
                throw new IllegalStateException(
                        "No physical site mapped for logical site " + logicalSite
                                + " — resolveSiteMap() did not cover a site referenced by this CvCfBasis "
                                + "Definition's vSpecs; the registered logicalSiteCoords likely doesn't include "
                                + "every site index used in vSpecs.");
            }
            siteIndices.add(physicalSite);
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
