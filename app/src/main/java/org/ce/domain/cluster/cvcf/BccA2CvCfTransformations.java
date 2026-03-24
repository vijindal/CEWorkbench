package org.ce.domain.cluster.cvcf;

import java.util.List;

/**
 * Hardcoded CVCF basis transformation matrices for BCC_A2.
 *
 * <p>Each matrix T satisfies: {@code u_old[i] = Σ_j T[i][j] · v_new[j]}</p>
 *
 * <p>Row order = old orthogonal CF order (u[type][group][1] in Mathematica notation).
 * Column order = new CVCF CF order as listed in {@code *_CF_NAMES}.</p>
 *
 * <p>Point variables are listed explicitly (xA, xB for binary; xA, xB, xC for ternary)
 * so that expressions like {@code xA + xC} and {@code -xA + xC} are captured
 * by T with no need for a separate constant vector.</p>
 *
 * <p>Transcribed verbatim from user-provided Mathematica transformation rules.</p>
 */
public final class BccA2CvCfTransformations {

    private BccA2CvCfTransformations() {}

    // =========================================================================
    // Binary BCC_A2  (2 components: A, B)
    //
    // Old CFs (rows, 6 total):
    //   u[1][1][1], u[2][1][1], u[3][1][1], u[4][1][1], u[5][1][1], u[6][1][1]
    //
    // New CFs (cols, 6 total):
    //   v21AB  v22AB  v3AB  v4AB  xA  xB
    //   col 0    1     2     3    4   5
    //
    // Transformation rules:
    //   u[1][1][1] = -16·v21AB + 8·v22AB + 16·v4AB + xA + xB
    //   u[2][1][1] =            - 4·v3AB            - xA + xB
    //   u[3][1][1] =  -4·v22AB                      + xA + xB
    //   u[4][1][1] =  -4·v21AB                      + xA + xB
    //   u[5][1][1] =                                - xA + xB
    //   u[6][1][1] =                                + xA + xB
    // =========================================================================

    /** CF names for binary BCC_A2 (non-point first, then point variables). */
    public static final List<String> BINARY_CF_NAMES = List.of(
            "v21AB", "v22AB", "v3AB", "v4AB", "xA", "xB"
    );

    /** Number of non-point CFs (optimization variables) for binary BCC_A2. */
    public static final int BINARY_NUM_NON_POINT_CFS = 4;

    /**
     * Transformation matrix for binary BCC_A2.
     * Dimensions: 6 rows (old CFs) × 6 columns (new CVCF CFs).
     */
    public static final double[][] BINARY_T = {
        // v21AB   v22AB   v3AB    v4AB    xA      xB
        { -16.0,    8.0,   0.0,   16.0,   1.0,    1.0 },  // u[1][1][1]
        {   0.0,    0.0,  -4.0,    0.0,  -1.0,    1.0 },  // u[2][1][1]
        {   0.0,   -4.0,   0.0,    0.0,   1.0,    1.0 },  // u[3][1][1]
        {  -4.0,    0.0,   0.0,    0.0,   1.0,    1.0 },  // u[4][1][1]
        {   0.0,    0.0,   0.0,    0.0,  -1.0,    1.0 },  // u[5][1][1]
        {   0.0,    0.0,   0.0,    0.0,   1.0,    1.0 },  // u[6][1][1]
    };

    // =========================================================================
    // Ternary BCC_A2  (3 components: A, B, C)
    //
    // Old CFs (rows, 21 total):
    //   u[1][1..6][1]  (6 rows)
    //   u[2][1..6][1]  (6 rows)
    //   u[3][1..3][1]  (3 rows)
    //   u[4][1..3][1]  (3 rows)
    //   u[5][1..2][1]  (2 rows)
    //   u[6][1][1]     (1 row)
    //
    // New CFs (cols, 21 total):
    //   v21AB  v21BC  v21AC  v22AB  v22BC  v22AC  v3AB  v3BC  v3AC
    //   col 0    1      2      3      4      5      6     7     8
    //
    //   v3ABC1  v3ABC2  v3ABC3  v4AB  v4BC  v4AC  v4ABC1  v4ABC2  v4ABC3  xA  xB  xC
    //   col 9    10      11     12    13    14     15      16      17     18  19  20
    //
    // Transformation rules (from user-provided Mathematica expressions):
    //   u[1][1][1] -> -2 v21AB - 2 v21BC + v3AB - 2 v3ABC2 - v3BC
    //                  + v4AB + 2 v4ABC2 + v4BC + xA + xC
    //   u[1][2][1] -> 2 v21AB - 2 v21BC - v3AB + 2 v3ABC1 - 2 v3ABC3 - v3BC
    //                  - v4AB + v4BC - xA + xC
    //   u[1][3][1] -> -2 v21AB - 2 v21BC - 4 v22AC + v3AB + 6 v3ABC2 - v3BC
    //                  + v4AB - 2 v4ABC2 + v4BC + xA + xC
    //   u[1][4][1] -> -2 v21AB - 4 v21AC - 2 v21BC + v3AB + 2 v3ABC1 + 2 v3ABC3 - v3BC
    //                  + v4AB + 2 v4ABC1 + 2 v4ABC3 + v4BC + xA + xC
    //   u[1][5][1] -> 2 v21AB - 2 v21BC - v3AB + 2 v3ABC1 - 2 v3ABC3 - 4 v3AC - v3BC
    //                  - v4AB - 4 v4ABC1 + 4 v4ABC3 + v4BC - xA + xC
    //   u[1][6][1] -> -2 v21AB - 16 v21AC - 2 v21BC + 8 v22AC + v3AB
    //                  + 8 v3ABC1 - 10 v3ABC2 + 8 v3ABC3 - v3BC + v4AB
    //                  + 8 v4ABC1 + 2 v4ABC2 + 8 v4ABC3 + 16 v4AC + v4BC + xA + xC
    //   u[2][1][1] -> 1/2 (-2 v21AB - 2 v21BC - v22AB - v22BC + v3AB - 2 v3ABC2 - v3BC
    //                  + 2 xA + 2 xC)
    //   u[2][2][1] -> 1/2 (2 v21AB - 2 v21BC + v22AB - v22BC - v3AB - 2 v3ABC1
    //                  + 2 v3ABC3 - v3BC - 2 xA + 2 xC)
    //   u[2][3][1] -> 1/2 (2 v21AB - 2 v21BC + v22AB - v22BC - v3AB + 2 v3ABC1
    //                  - 2 v3ABC3 - v3BC - 2 xA + 2 xC)
    //   u[2][4][1] -> 1/2 (-2 v21AB - 8 v21AC - 2 v21BC - v22AB - v22BC + v3AB
    //                  + 4 v3ABC1 - 2 v3ABC2 + 4 v3ABC3 - v3BC + 2 xA + 2 xC)
    //   u[2][5][1] -> 1/2 (-2 v21AB - 2 v21BC - v22AB - 8 v22AC - v22BC + v3AB
    //                  + 6 v3ABC2 - v3BC + 2 xA + 2 xC)
    //   u[2][6][1] -> 1/2 (2 v21AB - 2 v21BC + v22AB - v22BC - v3AB + 6 v3ABC1
    //                  - 6 v3ABC3 - 8 v3AC - v3BC - 2 xA + 2 xC)
    //   u[3][1][1] -> -v22AB - v22BC + xA + xC
    //   u[3][2][1] -> v22AB - v22BC - xA + xC
    //   u[3][3][1] -> -v22AB - 4 v22AC - v22BC + xA + xC
    //   u[4][1][1] -> -v21AB - v21BC + xA + xC
    //   u[4][2][1] -> v21AB - v21BC - xA + xC
    //   u[4][3][1] -> -v21AB - 4 v21AC - v21BC + xA + xC
    //   u[5][1][1] -> xA + xC
    //   u[5][2][1] -> -xA + xC
    //   u[6][1][1] -> xA + xB + xC
    // =========================================================================

    /**
     * CF names for ternary BCC_A2 (non-point first, then point variables).
     *
     * <p>Non-point CFs (18): binary-subset pairs, II-neighbor pairs, triangles,
     * ternary triangles (3 types), tetrahedra, ternary tetrahedra (3 types).</p>
     * <p>Point variables (3): xA, xB, xC.</p>
     */
    public static final List<String> TERNARY_CF_NAMES = List.of(
            // I-neighbor pairs (binary subsets + ternary cross-pair)
            "v21AB", "v21BC", "v21AC",
            // II-neighbor pairs
            "v22AB", "v22BC", "v22AC",
            // Triangles (binary subsets)
            "v3AB", "v3BC", "v3AC",
            // Triangles (ternary types 1, 2, 3)
            "v3ABC1", "v3ABC2", "v3ABC3",
            // Tetrahedra (binary subsets)
            "v4AB", "v4BC", "v4AC",
            // Tetrahedra (ternary types 1, 2, 3)
            "v4ABC1", "v4ABC2", "v4ABC3",
            // Point variables
            "xA", "xB", "xC"
    );

    /** Number of non-point CFs (optimization variables) for ternary BCC_A2. */
    public static final int TERNARY_NUM_NON_POINT_CFS = 18;

    /**
     * Transformation matrix for ternary BCC_A2.
     * Dimensions: 21 rows (old CFs) × 21 columns (new CVCF CFs).
     */
    public static final double[][] TERNARY_T = {
        // col:  v21AB  v21BC  v21AC  v22AB  v22BC  v22AC  v3AB   v3BC   v3AC   v3ABC1 v3ABC2 v3ABC3 v4AB   v4BC   v4AC   v4ABC1 v4ABC2 v4ABC3  xA     xB     xC
        {        -2.0,  -2.0,   0.0,   0.0,   0.0,   0.0,  1.0,  -1.0,   0.0,   0.0,  -2.0,   0.0,  1.0,   1.0,   0.0,   0.0,   2.0,   0.0,  1.0,   0.0,   1.0 }, // u[1][1][1]
        {         2.0,  -2.0,   0.0,   0.0,   0.0,   0.0, -1.0,  -1.0,   0.0,   2.0,   0.0,  -2.0, -1.0,   1.0,   0.0,   0.0,   0.0,   0.0, -1.0,   0.0,   1.0 }, // u[1][2][1]
        {        -2.0,  -2.0,   0.0,   0.0,   0.0,  -4.0,  1.0,  -1.0,   0.0,   0.0,   6.0,   0.0,  1.0,   1.0,   0.0,   0.0,  -2.0,   0.0,  1.0,   0.0,   1.0 }, // u[1][3][1]
        {        -2.0,  -2.0,  -4.0,   0.0,   0.0,   0.0,  1.0,  -1.0,   0.0,   2.0,   0.0,   2.0,  1.0,   1.0,   0.0,   2.0,   0.0,   2.0,  1.0,   0.0,   1.0 }, // u[1][4][1]
        {         2.0,  -2.0,   0.0,   0.0,   0.0,   0.0, -1.0,  -1.0,  -4.0,   2.0,   0.0,  -2.0, -1.0,   1.0,   0.0,  -4.0,   0.0,   4.0, -1.0,   0.0,   1.0 }, // u[1][5][1]
        {        -2.0,  -2.0, -16.0,   0.0,   0.0,   8.0,  1.0,  -1.0,   0.0,   8.0, -10.0,   8.0,  1.0,   1.0,  16.0,   8.0,   2.0,   8.0,  1.0,   0.0,   1.0 }, // u[1][6][1]
        {        -1.0,  -1.0,   0.0,  -0.5,  -0.5,   0.0,  0.5,  -0.5,   0.0,   0.0,  -1.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[2][1][1]
        {         1.0,  -1.0,   0.0,   0.5,  -0.5,   0.0, -0.5,  -0.5,   0.0,  -1.0,   0.0,   1.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.0,   0.0,   1.0 }, // u[2][2][1]
        {         1.0,  -1.0,   0.0,   0.5,  -0.5,   0.0, -0.5,  -0.5,   0.0,   1.0,   0.0,  -1.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.0,   0.0,   1.0 }, // u[2][3][1]
        {        -1.0,  -1.0,  -4.0,  -0.5,  -0.5,   0.0,  0.5,  -0.5,   0.0,   2.0,  -1.0,   2.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[2][4][1]
        {        -1.0,  -1.0,   0.0,  -0.5,  -0.5,  -4.0,  0.5,  -0.5,   0.0,   0.0,   3.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[2][5][1]
        {         1.0,  -1.0,   0.0,   0.5,  -0.5,   0.0, -0.5,  -0.5,  -4.0,   3.0,   0.0,  -3.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.0,   0.0,   1.0 }, // u[2][6][1]
        {         0.0,   0.0,   0.0,  -1.0,  -1.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[3][1][1]
        {         0.0,   0.0,   0.0,   1.0,  -1.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.0,   0.0,   1.0 }, // u[3][2][1]
        {         0.0,   0.0,   0.0,  -1.0,  -1.0,  -4.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[3][3][1]
        {        -1.0,  -1.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[4][1][1]
        {         1.0,  -1.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.0,   0.0,   1.0 }, // u[4][2][1]
        {        -1.0,  -1.0,  -4.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[4][3][1]
        {         0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   0.0,   1.0 }, // u[5][1][1]
        {         0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.0,   0.0,   1.0 }, // u[5][2][1]
        {         0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,   1.0,   1.0 }, // u[6][1][1]
    };

    // =========================================================================
    // Quaternary T matrix will be added once the user provides the rules.
    // =========================================================================

    // =========================================================================
    // Factory methods
    // =========================================================================

    /** Returns the {@link CvCfBasis} for binary BCC_A2. */
    public static CvCfBasis binaryBasis() {
        return new CvCfBasis("BCC_A2", 2, BINARY_CF_NAMES, BINARY_NUM_NON_POINT_CFS, BINARY_T);
    }

    /** Returns the {@link CvCfBasis} for ternary BCC_A2. */
    public static CvCfBasis ternaryBasis() {
        return new CvCfBasis("BCC_A2", 3, TERNARY_CF_NAMES, TERNARY_NUM_NON_POINT_CFS, TERNARY_T);
    }
}
