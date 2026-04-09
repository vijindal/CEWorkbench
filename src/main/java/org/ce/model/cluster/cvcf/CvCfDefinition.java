package org.ce.domain.cluster.cvcf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Registry of physical CVCF basis definitions keyed by (structurePhase, model, numComponents).
 *
 * <p>Each definition specifies:
 * <ul>
 *   <li>The fractional-coordinate positions of the logical sites (p1, p2, ...) used in the
 *       v-function expressions. These are matched at runtime against the physical site list
 *       from Stage 3 ({@code CMatrix.Result.getSiteList()}) to resolve logical→physical
 *       site indices without relying on fragile geometry heuristics.</li>
 *   <li>The ordered CF names (non-point CFs first, then one per component).</li>
 *   <li>The v-function specifications ({@link VSpec}): products or differences of site
 *       occupation probabilities p[logicalSite][atomIndex].</li>
 * </ul>
 *
 * <p>To add a new combination, add a {@code register(...)} call in the static block.
 * Nothing else needs to change in {@link CvCfMatrixGenerator}.</p>
 *
 * <h3>Coordinate convention</h3>
 * <p>Site coordinates are in <em>fractional lattice coordinates</em>, matching the
 * {@code Position} class used throughout the cluster pipeline.</p>
 *
 * <h3>Index conventions</h3>
 * <ul>
 *   <li>Logical site indices are 1-based (p1, p2, ...).</li>
 *   <li>Atom indices are 0-based (A=0, B=1, C=2, D=3).</li>
 * </ul>
 */
public final class CvCfDefinition {

    // =========================================================================
    // Public data
    // =========================================================================

    /**
     * Fractional coordinates of each logical site.
     * {@code logicalSiteCoords[i]} = {x, y, z} of logical site (i+1).
     */
    public final double[][] logicalSiteCoords;

    /** Ordered CF names: non-point CFs first, then one per component (xA, xB, …). */
    public final List<String> cfNames;

    /** One {@link VSpec} per entry in {@code cfNames}. */
    public final List<VSpec> vSpecs;

    // =========================================================================
    // Constructor (package-private — use registry)
    // =========================================================================

    private CvCfDefinition(double[][] logicalSiteCoords,
                            List<String> cfNames,
                            List<VSpec> vSpecs) {
        this.logicalSiteCoords = logicalSiteCoords;
        this.cfNames           = Collections.unmodifiableList(new ArrayList<>(cfNames));
        this.vSpecs            = Collections.unmodifiableList(new ArrayList<>(vSpecs));
    }

    // =========================================================================
    // VSpec — symbolic v-function specification
    // =========================================================================

    /**
     * Symbolic specification of one v-function as a product of site-occupation
     * probabilities, or a difference of two such products.
     *
     * <p>A term is encoded as a flat int array of (logicalSite, atomIndex) pairs:
     * {@code [s1, a1, s2, a2, ...]} where each pair represents p[s][a].</p>
     */
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

    // =========================================================================
    // Registry
    // =========================================================================

    private static final Map<Key, CvCfDefinition> REGISTRY = new LinkedHashMap<>();

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

    // =========================================================================
    // Public API
    // =========================================================================

    /**
     * Returns the definition for the given combination, or throws
     * {@link UnsupportedOperationException} with a helpful message listing
     * all registered keys.
     */
    public static CvCfDefinition get(String structurePhase, String model, int numComponents) {
        CvCfDefinition def = REGISTRY.get(new Key(structurePhase, model, numComponents));
        if (def == null) {
            throw new UnsupportedOperationException(
                "No CvCfDefinition registered for ("
                + structurePhase + ", " + model + ", K=" + numComponents + ").\n"
                + "Registered combinations: " + registeredSummary());
        }
        return def;
    }

    /** Returns {@code true} if a definition exists for the given combination. */
    public static boolean isSupported(String structurePhase, String model, int numComponents) {
        return REGISTRY.containsKey(new Key(structurePhase, model, numComponents));
    }

    // =========================================================================
    // Internal helpers
    // =========================================================================

    private static void register(String structurePhase,
                                  String model,
                                  int numComponents,
                                  double[][] logicalSiteCoords,
                                  List<String> cfNames,
                                  List<VSpec> vSpecs) {
        if (cfNames.size() != vSpecs.size()) {
            throw new IllegalArgumentException(
                "cfNames.size() != vSpecs.size() for ("
                + structurePhase + ", " + model + ", K=" + numComponents + ")");
        }
        REGISTRY.put(new Key(structurePhase, model, numComponents),
                     new CvCfDefinition(logicalSiteCoords, cfNames, vSpecs));
    }

    private static String registeredSummary() {
        StringBuilder sb = new StringBuilder();
        for (Key k : REGISTRY.keySet()) {
            sb.append("\n  (").append(k.structurePhase)
              .append(", ").append(k.model)
              .append(", K=").append(k.numComponents).append(")");
        }
        return sb.toString().isEmpty() ? " (none)" : sb.toString();
    }

    // =========================================================================
    // Key
    // =========================================================================

    private static final class Key {
        final String structurePhase;
        final String model;
        final int    numComponents;

        Key(String structurePhase, String model, int numComponents) {
            this.structurePhase = structurePhase;
            this.model          = model == null ? "" : model.toUpperCase();
            this.numComponents  = numComponents;
        }

        @Override public boolean equals(Object o) {
            if (!(o instanceof Key)) return false;
            Key k = (Key) o;
            return numComponents == k.numComponents
                && structurePhase.equals(k.structurePhase)
                && model.equals(k.model);
        }

        @Override public int hashCode() {
            return 31 * (31 * structurePhase.hashCode() + model.hashCode()) + numComponents;
        }
    }
}
