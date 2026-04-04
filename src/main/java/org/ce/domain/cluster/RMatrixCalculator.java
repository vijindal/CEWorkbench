package org.ce.domain.cluster;

/**
 * Builds the Inden (1992) R-matrix used to convert site-occupation operators
 * to site-operator basis functions.
 */
public final class RMatrixCalculator {

    private static final double EPS = 1e-12;

    private RMatrixCalculator() {}

    /**
     * Builds the R-matrix for the given number of elements.
     *
     * <p>The matrix is the inverse of a Vandermonde-like matrix constructed
     * from a symmetric integer basis sequence (Mathematica genRMat).</p>
     *
     * @param numElements number of components (>= 2)
     * @return R-matrix of size numElements x numElements
     */
    public static double[][] buildRMatrix(int numElements) {
        if (numElements < 2) {
            throw new IllegalArgumentException("numElements must be >= 2");
        }

        double[] basis = buildBasis(numElements);
        double[][] matM = new double[numElements][numElements];

        for (int i = 0; i < numElements; i++) {
            for (int j = 0; j < numElements; j++) {
                if (i == 0) {
                    matM[i][j] = 1.0;
                } else {
                    matM[i][j] = Math.pow(basis[j], i);
                }
            }
        }

        return invert(matM);
    }

    /**
     * Returns the symmetric integer basis sequence for K components.
     *
     * <p>For even K: {-K/2, ..., -1, 1, ..., K/2}.<br>
     * For odd K: {-(K-1)/2, ..., 0, ..., (K-1)/2}.</p>
     *
     * <p>Binary (K=2): {-1, 1}. Ternary (K=3): {-1, 0, 1}.</p>
     *
     * @param numElements number of components K (â‰¥ 2)
     * @return basis array of length K
     */
    public static double[] buildBasis(int numElements) {
        double[] basis = new double[numElements];
        if (numElements % 2 == 0) {
            int half = numElements / 2;
            for (int i = 0; i < half; i++) {
                basis[i] = -half + i;
            }
            for (int i = 0; i < half; i++) {
                basis[half + i] = 1 + i;
            }
        } else {
            int start = -((numElements - 1) / 2);
            for (int i = 0; i < numElements; i++) {
                basis[i] = start + i;
            }
        }
        return basis;
    }

    private static double[][] invert(double[][] matrix) {
        int n = matrix.length;
        double[][] a = new double[n][n];
        double[][] inv = new double[n][n];

        for (int i = 0; i < n; i++) {
            if (matrix[i].length != n) {
                throw new IllegalArgumentException("Matrix must be square");
            }
            for (int j = 0; j < n; j++) {
                a[i][j] = matrix[i][j];
            }
            inv[i][i] = 1.0;
        }

        for (int i = 0; i < n; i++) {
            int pivot = i;
            double max = Math.abs(a[i][i]);
            for (int r = i + 1; r < n; r++) {
                double v = Math.abs(a[r][i]);
                if (v > max) {
                    max = v;
                    pivot = r;
                }
            }
            if (max < EPS) {
                throw new IllegalStateException("Matrix is singular or ill-conditioned");
            }
            if (pivot != i) {
                double[] tmp = a[i];
                a[i] = a[pivot];
                a[pivot] = tmp;
                double[] tmpInv = inv[i];
                inv[i] = inv[pivot];
                inv[pivot] = tmpInv;
            }

            double diag = a[i][i];
            for (int j = 0; j < n; j++) {
                a[i][j] /= diag;
                inv[i][j] /= diag;
            }

            for (int r = 0; r < n; r++) {
                if (r == i) {
                    continue;
                }
                double factor = a[r][i];
                if (Math.abs(factor) < EPS) {
                    continue;
                }
                for (int c = 0; c < n; c++) {
                    a[r][c] -= factor * a[i][c];
                    inv[r][c] -= factor * inv[i][c];
                }
            }
        }

        return inv;
    }
}

