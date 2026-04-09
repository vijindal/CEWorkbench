package org.ce.domain.cluster;

/**
 * Minimal linear algebra utilities for thermodynamic solvers.
 *
 * <p>Provides Gaussian elimination with partial pivoting for solving
 * {@code A Â· x = b}. This avoids a dependency on external linear algebra
 * libraries for the small (ncf Ã— ncf, typically 4Ã—4) systems encountered
 * in CVM calculations.</p>
 *
 * <p>Thread Safety: All methods are stateless and thread-safe.</p>
 *
 * @since 2.0
 */
public final class LinearAlgebra {

    private LinearAlgebra() { /* utility class */ }

    /**
     * Solves the linear system {@code A Â· x = b} using Gaussian elimination
     * with partial pivoting.
     *
     * <p>The input arrays are <em>not</em> modified (copies are made internally).</p>
     *
     * @param A  coefficient matrix (n Ã— n)
     * @param b  right-hand side vector (length n)
     * @return   solution vector x (length n)
     * @throws IllegalArgumentException if A is singular or dimensions mismatch
     */
    public static double[] solve(double[][] A, double[] b) {
        int n = A.length;
        if (n == 0) throw new IllegalArgumentException("Empty matrix");
        if (A[0].length != n) throw new IllegalArgumentException("Matrix must be square");
        if (b.length != n) throw new IllegalArgumentException("RHS length must match matrix size");

        // Work on copies
        double[][] M = new double[n][n];
        double[] rhs = new double[n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            rhs[i] = b[i];
        }

        // Forward elimination with partial pivoting
        for (int col = 0; col < n; col++) {
            // Find pivot
            int maxRow = col;
            double maxVal = Math.abs(M[col][col]);
            for (int row = col + 1; row < n; row++) {
                double v = Math.abs(M[row][col]);
                if (v > maxVal) {
                    maxVal = v;
                    maxRow = row;
                }
            }

            if (maxVal < 1e-30) {
                throw new IllegalArgumentException(
                        "Singular or near-singular matrix (pivot = " + maxVal + " at column " + col + ")");
            }

            // Swap rows
            if (maxRow != col) {
                double[] tmp = M[col]; M[col] = M[maxRow]; M[maxRow] = tmp;
                double t = rhs[col]; rhs[col] = rhs[maxRow]; rhs[maxRow] = t;
            }

            // Eliminate below
            double pivot = M[col][col];
            for (int row = col + 1; row < n; row++) {
                double factor = M[row][col] / pivot;
                for (int j = col; j < n; j++) {
                    M[row][j] -= factor * M[col][j];
                }
                rhs[row] -= factor * rhs[col];
            }
        }

        // Back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = rhs[i];
            for (int j = i + 1; j < n; j++) {
                sum -= M[i][j] * x[j];
            }
            x[i] = sum / M[i][i];
        }

        return x;
    }

    /**
     * Computes the L2 norm (Euclidean length) of a vector.
     *
     * @param v the vector
     * @return â€–vâ€–â‚‚ = âˆš(Î£ v[i]Â²)
     */
    public static double norm(double[] v) {
        double sum = 0.0;
        for (double x : v) {
            sum += x * x;
        }
        return Math.sqrt(sum);
    }

    /**
     * Computes the maximum absolute element of a vector (infinity norm).
     *
     * @param v the vector
     * @return â€–vâ€–âˆž = max|v[i]|
     */
    public static double normInf(double[] v) {
        double max = 0.0;
        for (double x : v) {
            double abs = Math.abs(x);
            if (abs > max) max = abs;
        }
        return max;
    }

    /**
     * Computes the inverse of a square matrix using Gaussian elimination with
     * partial pivoting.
     *
     * @param A the n Ã— n matrix to invert
     * @return the n Ã— n inverse matrix
     * @throws IllegalArgumentException if the matrix is singular or not square
     */
    public static double[][] invert(double[][] A) {
        int n = A.length;
        if (n == 0) throw new IllegalArgumentException("Empty matrix");
        if (A[0].length != n) throw new IllegalArgumentException("Matrix must be square");

        // Augmented matrix [A | I]
        double[][] aug = new double[n][2 * n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, aug[i], 0, n);
            aug[i][n + i] = 1.0;
        }

        // Forward elimination
        for (int col = 0; col < n; col++) {
            // Find pivot
            int pivot = col;
            for (int row = col + 1; row < n; row++) {
                if (Math.abs(aug[row][col]) > Math.abs(aug[pivot][col])) pivot = row;
            }

            // Check for singularity
            if (Math.abs(aug[pivot][col]) < 1e-18) {
                throw new IllegalArgumentException("Matrix is singular or near-singular at column " + col);
            }

            // Swap rows
            double[] tmp = aug[col];
            aug[col] = aug[pivot];
            aug[pivot] = tmp;

            // Normalize pivot row
            double scale = aug[col][col];
            for (int j = 0; j < 2 * n; j++) {
                aug[col][j] /= scale;
            }

            // Eliminate other rows
            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double f = aug[row][col];
                for (int j = 0; j < 2 * n; j++) {
                    aug[row][j] -= f * aug[col][j];
                }
            }
        }

        // Extract inverse [I | A^-1]
        double[][] inv = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(aug[i], n, inv[i], 0, n);
        }
        return inv;
    }

    /**
     * Computes the dot product of two vectors.
     *
     * @param a first vector
     * @param b second vector
     * @return a Â· b = Î£ a[i]Â·b[i]
     * @throws IllegalArgumentException if vectors have different lengths
     */
    public static double dot(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("Vector lengths must match");
        }
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }
}

