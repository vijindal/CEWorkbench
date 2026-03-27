package org.ce.domain.cluster.cvcf;

import org.ce.domain.cluster.CMatrixResult;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Transforms an old orthogonal-basis {@link CMatrixResult} to the CVCF basis.
 *
 * <p>Given the relationship {@code u_old[i] = Σ_j T[i][j] · v_new[j]},
 * the cluster variable evaluation:
 * <pre>
 *   CV[v] = Σ_{k_old} cmat_old[v][k_old] · u_old[k_old] + cmat_old[v][tcf_old]
 * </pre>
 * becomes:
 * <pre>
 *   CV[v] = Σ_{k_new} cmat_new[v][k_new] · v_new[k_new] + cmat_new[v][tcf_new]
 * </pre>
 * where:
 * <pre>
 *   cmat_new[v][k_new] = Σ_{k_old} cmat_old[v][k_old] · T[k_old][k_new]
 *   cmat_new[v][tcf_new] = cmat_old[v][tcf_old]   (constant term unchanged)
 * </pre>
 * i.e., {@code cmat_new[:, :tcf_new] = cmat_old[:, :tcf_old] · T}.
 */
public final class CvCfBasisTransformer {

    private static final Logger LOG = Logger.getLogger(CvCfBasisTransformer.class.getName());
    private static final double VERIFY_TOL = 1e-10;

    // Mathematica reference for BCC_A2 binary (orthogonal basis cmat blocks).
    private static final double[][][] ORTHO_MATH_BCC_A2_BIN = {
            {
                    {1.0 / 16.0, -1.0 / 4.0, 1.0 / 8.0, 1.0 / 4.0, -1.0 / 4.0, 1.0 / 16.0},
                    {-1.0 / 16.0, 1.0 / 8.0, 0.0, 0.0, -1.0 / 8.0, 1.0 / 16.0},
                    {1.0 / 16.0, 0.0, -1.0 / 8.0, 0.0, 0.0, 1.0 / 16.0},
                    {1.0 / 16.0, 0.0, 1.0 / 8.0, -1.0 / 4.0, 0.0, 1.0 / 16.0},
                    {-1.0 / 16.0, -1.0 / 8.0, 0.0, 0.0, 1.0 / 8.0, 1.0 / 16.0},
                    {1.0 / 16.0, 1.0 / 4.0, 1.0 / 8.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 16.0}
            },
            {
                    {0.0, -1.0 / 8.0, 1.0 / 8.0, 1.0 / 4.0, -3.0 / 8.0, 1.0 / 8.0},
                    {0.0, 1.0 / 8.0, 1.0 / 8.0, -1.0 / 4.0, -1.0 / 8.0, 1.0 / 8.0},
                    {0.0, 1.0 / 8.0, -1.0 / 8.0, 0.0, -1.0 / 8.0, 1.0 / 8.0},
                    {0.0, -1.0 / 8.0, -1.0 / 8.0, 0.0, 1.0 / 8.0, 1.0 / 8.0},
                    {0.0, -1.0 / 8.0, 1.0 / 8.0, -1.0 / 4.0, 1.0 / 8.0, 1.0 / 8.0},
                    {0.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 8.0}
            },
            {
                    {0.0, 0.0, 1.0 / 4.0, 0.0, -1.0 / 2.0, 1.0 / 4.0},
                    {0.0, 0.0, -1.0 / 4.0, 0.0, 0.0, 1.0 / 4.0},
                    {0.0, 0.0, 1.0 / 4.0, 0.0, 1.0 / 2.0, 1.0 / 4.0}
            },
            {
                    {0.0, 0.0, 0.0, 1.0 / 4.0, -1.0 / 2.0, 1.0 / 4.0},
                    {0.0, 0.0, 0.0, -1.0 / 4.0, 0.0, 1.0 / 4.0},
                    {0.0, 0.0, 0.0, 1.0 / 4.0, 1.0 / 2.0, 1.0 / 4.0}
            },
            {
                    {0.0, 0.0, 0.0, 0.0, -1.0 / 2.0, 1.0 / 2.0},
                    {0.0, 0.0, 0.0, 0.0, 1.0 / 2.0, 1.0 / 2.0}
            }
    };

    // Mathematica reference for BCC_A2 binary (CVCF basis cmat blocks).
    private static final double[][][] CVCF_MATH_BCC_A2_BIN = {
            {
                    {1.0, 1.0, 0.0, -2.0, 1.0, 0.0},
                    {-1.0, -1.0 / 2.0, -1.0 / 2.0, 1.0, 0.0, 0.0},
                    {1.0, 0.0, 1.0, -1.0, 0.0, 0.0},
                    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {-1.0, 1.0 / 2.0, -1.0 / 2.0, 1.0, 0.0, 0.0},
                    {1.0, -1.0, 0.0, -2.0, 0.0, 1.0}
            },
            {
                    {0.0, 1.0 / 2.0, -1.0 / 2.0, -1.0, 1.0, 0.0},
                    {0.0, -1.0 / 2.0, -1.0 / 2.0, 1.0, 0.0, 0.0},
                    {0.0, -1.0 / 2.0, 1.0 / 2.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0 / 2.0, 1.0 / 2.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0 / 2.0, -1.0 / 2.0, 1.0, 0.0, 0.0},
                    {0.0, -1.0 / 2.0, -1.0 / 2.0, -1.0, 0.0, 1.0}
            },
            {
                    {0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, -1.0, 0.0, 0.0, 1.0}
            },
            {
                    {0.0, 0.0, 0.0, -1.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, -1.0, 0.0, 1.0}
            },
            {
                    {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
            }
    };

    private CvCfBasisTransformer() {}

    /**
     * Transforms a CMatrixResult from the old orthogonal basis to the CVCF basis.
     *
     * @param oldResult  the CMatrixResult in the old orthogonal basis
     * @param basis      the CVCF basis description, including T matrix and CF names
     * @return a new CMatrixResult in the CVCF basis
     */
    public static CMatrixResult transform(CMatrixResult oldResult, CvCfBasis basis) {
        if (oldResult == null) {
            throw new IllegalArgumentException("oldResult must not be null");
        }
        if (basis == null) {
            throw new IllegalArgumentException("basis must not be null");
        }

        List<List<double[][]>> oldCmat = oldResult.getCmat();
        int[][] oldLcv = oldResult.getLcv();
        List<List<int[]>> wcv = oldResult.getWcv();

        if (oldCmat == null || oldLcv == null || wcv == null) {
            throw new IllegalArgumentException("oldResult is missing required cmat/lcv/wcv fields");
        }
        if (oldCmat.isEmpty()) {
            throw new IllegalArgumentException("oldResult.cmat is empty");
        }
        if (basis.T == null || basis.T.length == 0 || basis.T[0].length == 0) {
            throw new IllegalArgumentException("basis.T must be a non-empty matrix");
        }

        int tcdis = oldCmat.size();
        double[][] T = basis.T;
        int tcfOld = T.length;          // number of old CF columns (excluding constant column)
        int tcfNew = basis.totalCfs();  // number of new CF columns
        int tCols = T[0].length;
        if (basis.numNonPointCfs + basis.numComponents != tcfNew) {
            throw new IllegalArgumentException("Invalid CVCF basis cardinality: numNonPointCfs("
                    + basis.numNonPointCfs + ") + numComponents(" + basis.numComponents
                    + ") must equal totalCfs(" + tcfNew + ")");
        }

        if (tCols != tcfNew) {
            throw new IllegalArgumentException("Basis shape mismatch: T has " + tCols
                    + " columns but basis.totalCfs() = " + tcfNew);
        }
        for (int r = 1; r < T.length; r++) {
            if (T[r] == null || T[r].length != tCols) {
                throw new IllegalArgumentException("Basis T must be rectangular. Row " + r
                        + " has length " + (T[r] == null ? 0 : T[r].length) + ", expected " + tCols);
            }
        }

        List<List<double[][]>> newCmat = new ArrayList<>(tcdis);

        for (int t = 0; t < tcdis; t++) {
            List<double[][]> oldGroup = oldCmat.get(t);
            int lc = oldGroup.size();
            List<double[][]> newGroup = new ArrayList<>(lc);

            for (int j = 0; j < lc; j++) {
                double[][] cmOld = oldGroup.get(j);  // [nv][numOldCols]
                int nv = cmOld.length;
                int numOldCols = (nv > 0) ? cmOld[0].length : 0;
                if (numOldCols < 1) {
                    throw new IllegalArgumentException("Invalid cmat row width at (t=" + t + ", j=" + j
                            + "): expected at least 1 column for constant term");
                }
                int oldTcfFromCmat = numOldCols - 1;
                // Supported shape modes:
                // 0) numOldCols == T.rows:
                //    C-matrix already stores exactly the columns represented by T rows
                //    (e.g. Mathematica-style {u1..u6} including the old "empty"/sum term).
                //    In this mode we transform all columns by T and DO NOT append a
                //    separate constant column.
                // 1) oldTcfFromCmat == T.rows:
                //    T rows map exactly to old CF columns (no empty-cluster row in T).
                // 2) oldTcfFromCmat + 1 == T.rows:
                //    T includes an extra empty-cluster row. Since old C-matrix CF columns
                //    exclude the empty cluster, we use only the first oldTcfFromCmat rows.
                int transformRows;
                boolean transformIncludesLastOldColumn = false;
                if (numOldCols == tcfOld) {
                    transformRows = tcfOld;
                    transformIncludesLastOldColumn = true;
                } else if (oldTcfFromCmat == tcfOld) {
                    transformRows = tcfOld;
                } else if (oldTcfFromCmat + 1 == tcfOld) {
                    transformRows = oldTcfFromCmat;
                } else {
                    throw new IllegalArgumentException("Dimension mismatch at (t=" + t + ", j=" + j + "): "
                            + "cmat old CF columns = " + oldTcfFromCmat + ", but basis.T rows = " + tcfOld
                            + ". Supported modes are rows==oldCFs or rows==oldCFs+1 (extra empty-cluster row).");
                }
                for (int v = 1; v < nv; v++) {
                    if (cmOld[v] == null || cmOld[v].length != numOldCols) {
                        throw new IllegalArgumentException("Non-rectangular cmat block at (t=" + t + ", j=" + j
                                + "): row " + v + " has length "
                                + (cmOld[v] == null ? 0 : cmOld[v].length) + ", expected " + numOldCols);
                    }
                }

                // The last column of old C-matrix is the constant term
                // Structure: [CF0, CF1, ..., CFn, constant]
                // where constant column contains the fixed term (usually same value across rows)

                // New cmat shape:
                // - full-column transform mode: exactly tcfNew columns (no separate constant column)
                // - legacy mode: tcfNew CF columns + 1 constant column
                int newCols = transformIncludesLastOldColumn ? tcfNew : (tcfNew + 1);
                double[][] cmNew = new double[nv][newCols];

                for (int v = 0; v < nv; v++) {
                    // Multiply: cmat_new[v][k_new] = Σ_{k_old} cmat_old[v][k_old] · T[k_old][k_new]
                    for (int kNew = 0; kNew < tcfNew; kNew++) {
                        double sum = 0.0;
                        for (int kOld = 0; kOld < transformRows; kOld++) {
                            sum += cmOld[v][kOld] * T[kOld][kNew];
                        }
                        cmNew[v][kNew] = sum;
                    }
                    if (!transformIncludesLastOldColumn) {
                        // Legacy mode: preserve constant term from last column of old matrix.
                        cmNew[v][tcfNew] = cmOld[v][numOldCols - 1];
                    }
                }
                newGroup.add(cmNew);
            }
            newCmat.add(newGroup);
        }

        if (isBinaryBccA2(basis) && hasOneGroupPerType(oldCmat) && hasOneGroupPerType(newCmat)) {
            logFullCmat("ORTHOGONAL CMAT (pre-transform)", oldCmat);
            logFullCmat("CVCF CMAT (post-transform)", newCmat);
            verifyAgainstMathematica("ORTHOGONAL", oldCmat, ORTHO_MATH_BCC_A2_BIN);
            verifyAgainstMathematica("CVCF", newCmat, CVCF_MATH_BCC_A2_BIN);
        }

        // Preserve orthogonal-basis decoration metadata for exact random-state
        // initialization (orthogonal random CFs -> CVCF via T^{-1}).
        return new CMatrixResult(newCmat, oldLcv, wcv, oldResult.getCfBasisIndices());
    }

    private static boolean isBinaryBccA2(CvCfBasis basis) {
        return basis != null
                && "BCC_A2".equalsIgnoreCase(basis.structurePhase)
                && basis.numComponents == 2;
    }

    private static boolean hasOneGroupPerType(List<List<double[][]>> cmat) {
        if (cmat == null || cmat.isEmpty()) {
            return false;
        }
        for (List<double[][]> groups : cmat) {
            if (groups == null || groups.size() != 1 || groups.get(0) == null) {
                return false;
            }
        }
        return true;
    }

    private static void logFullCmat(String title, List<List<double[][]>> cmat) {
        LOG.info("========== " + title + " ==========");
        for (int t = 0; t < cmat.size(); t++) {
            List<double[][]> groups = cmat.get(t);
            for (int j = 0; j < groups.size(); j++) {
                double[][] block = groups.get(j);
                int rows = block.length;
                int cols = rows > 0 ? block[0].length : 0;
                LOG.info("cmat[t=" + t + "][j=" + j + "] dims=" + rows + "x" + cols);
                for (int r = 0; r < rows; r++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append("  [");
                    for (int c = 0; c < cols; c++) {
                        if (c > 0) {
                            sb.append(", ");
                        }
                        sb.append(String.format("% .12f", block[r][c]));
                    }
                    sb.append("]");
                    LOG.info(sb.toString());
                }
            }
        }
    }

    private static void verifyAgainstMathematica(
            String label,
            List<List<double[][]>> actualCmat,
            double[][][] expectedBlocks) {

        if (actualCmat.size() != expectedBlocks.length) {
            LOG.warning(label + " verification: type-count mismatch. actual=" + actualCmat.size()
                    + ", expected=" + expectedBlocks.length);
            return;
        }

        double maxAbsDiff = 0.0;
        String maxWhere = "(none)";

        for (int t = 0; t < expectedBlocks.length; t++) {
            List<double[][]> groups = actualCmat.get(t);
            if (groups.size() != 1) {
                LOG.warning(label + " verification: t=" + t + " has " + groups.size()
                        + " groups; expected exactly 1");
                return;
            }
            double[][] actual = groups.get(0);
            double[][] expected = expectedBlocks[t];

            if (actual.length != expected.length) {
                LOG.warning(label + " verification: row mismatch at t=" + t
                        + ". actual=" + actual.length + ", expected=" + expected.length);
                return;
            }
            for (int r = 0; r < expected.length; r++) {
                if (actual[r].length != expected[r].length) {
                    LOG.warning(label + " verification: col mismatch at t=" + t + ", r=" + r
                            + ". actual=" + actual[r].length + ", expected=" + expected[r].length);
                    return;
                }
                for (int c = 0; c < expected[r].length; c++) {
                    double d = Math.abs(actual[r][c] - expected[r][c]);
                    if (d > maxAbsDiff) {
                        maxAbsDiff = d;
                        maxWhere = "t=" + t + ", j=0, r=" + r + ", c=" + c
                                + " (actual=" + actual[r][c] + ", expected=" + expected[r][c] + ")";
                    }
                }
            }
        }

        if (maxAbsDiff <= VERIFY_TOL) {
            LOG.info(label + " verification against Mathematica PASSED. max |Δ| = "
                    + String.format("%.3e", maxAbsDiff));
        } else {
            LOG.warning(label + " verification against Mathematica FAILED. max |Δ| = "
                    + String.format("%.3e", maxAbsDiff) + " at " + maxWhere);
        }
    }
}
