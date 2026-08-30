package org.ce.model.cluster;

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
     * with partial pivoting, after symmetric (Jacobi) diagonal scaling.
     *
     * <p>The CVM Hessian this is normally called on ({@code Guu} in
     * {@code CVMGibbsModel.minimize}) can be extremely badly scaled near a
     * dilute composition -- diagonal entries spanning ~10 orders of
     * magnitude between a majority-element cluster variable and a
     * near-zero rare-pair one. Plain partial-pivoting elimination on such a
     * matrix lets large-magnitude rows dominate the elimination of small
     * ones, degrading the solved direction for the small-magnitude
     * variables even though the matrix is not actually singular. Scaling
     * each row/column by {@code 1/sqrt(|A[i][i]|)} before elimination (and
     * undoing it on the solution) brings every diagonal entry to
     * O(1) magnitude, which is the standard remedy for this failure mode
     * and does not change the solution of a consistent system -- only its
     * numerical conditioning during elimination.</p>
     *
     * <p>The input arrays are <em>not</em> modified (copies are made internally).</p>
     *
     * @param A  coefficient matrix (n Ã— n)
     * @param b  right-hand side vector (length n)
     * @return   solution vector x (length n)
     * @throws IllegalArgumentException if A is singular or dimensions mismatch
     */
    public static double[] solve(double[][] A, double[] b) {
        return solveChecked(A, b).x();
    }

    /**
     * One solved system together with a cheap ill-conditioning signal.
     *
     * @param x                 the solution vector
     * @param relativeResidual  {@code ||A*x - b|| / ||b||} (0 if {@code b} is
     *                          the zero vector), computed in the original
     *                          (unscaled) system so it reflects the answer's
     *                          actual accuracy, not the scaled intermediate
     *                          one. A well-conditioned solve gives a value at
     *                          or near machine epsilon; a value orders of
     *                          magnitude larger means the diagonal rescaling
     *                          in {@link #solve} was not enough to keep
     *                          round-off from corrupting the result -- the
     *                          matrix itself is fundamentally hard for
     *                          floating-point elimination on this right-hand
     *                          side, not merely differently scaled.
     */
    public record Solution(double[] x, double relativeResidual) {}

    /**
     * Same elimination as {@link #solve}, but also reports the relative
     * residual of the returned solution -- see {@link Solution}.
     *
     * <p>Exists because {@link #solve} alone gives no signal when a matrix is
     * merely <em>badly</em> conditioned rather than exactly singular: the
     * {@code 1e-30} pivot guard only catches the latter. A CVM widened
     * Hessian near a dilute/near-boundary composition can have diagonal
     * entries spanning 10+ orders of magnitude (see this class's main doc)
     * and still clear that pivot threshold at every step while still losing
     * several digits of accuracy to round-off -- exactly the kind of failure
     * this residual check is meant to surface to a caller that cares (e.g.
     * {@code HillertSolver}, which can then warn via its progress sink rather
     * than silently accepting a degraded Newton step).</p>
     *
     * <p><b>Also applies iterative refinement</b> when the initial solve's
     * relative residual exceeds {@link #REFINEMENT_THRESHOLD}: it solves the
     * same matrix again against the residual {@code b - A*x} and corrects
     * {@code x} with that answer, up to {@link #MAX_REFINEMENT_STEPS} times
     * or until the residual stops improving. This is the standard fix for
     * exactly the failure mode {@link Solution#relativeResidual()} detects --
     * a matrix that is not singular but is badly scaled loses accuracy during
     * elimination in a way a second pass against the leftover residual
     * recovers, without changing what system is being solved or introducing
     * any new physics/algorithm. The returned residual reflects the
     * <em>refined</em> solution, so a caller reading it after this method
     * returns sees whether refinement actually helped.</p>
     */
    public static Solution solveChecked(double[][] A, double[] b) {
        double[] x = eliminate(A, b);
        double relative = relativeResidual(A, b, x);

        for (int step = 0; step < MAX_REFINEMENT_STEPS && relative > REFINEMENT_THRESHOLD; step++) {
            int n = b.length;
            double[] residual = new double[n];
            for (int i = 0; i < n; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < n; j++) {
                    rowSum += A[i][j] * x[j];
                }
                residual[i] = b[i] - rowSum;
            }
            double[] correction = eliminate(A, residual);
            double[] refined = new double[n];
            for (int i = 0; i < n; i++) {
                refined[i] = x[i] + correction[i];
            }
            double refinedRelative = relativeResidual(A, b, refined);
            if (refinedRelative >= relative) {
                // No further improvement (or it got worse, e.g. the
                // correction solve is itself hitting the same ill-conditioning) --
                // stop rather than churn or drift away from the best answer found.
                break;
            }
            x = refined;
            relative = refinedRelative;
        }

        return new Solution(x, relative);
    }

    /** A relative residual at or below this is treated as already accurate enough. */
    private static final double REFINEMENT_THRESHOLD = 1e-10;

    /** Cap on refinement passes -- each is one more full elimination on the same-size matrix. */
    private static final int MAX_REFINEMENT_STEPS = 3;

    private static double relativeResidual(double[][] A, double[] b, double[] x) {
        double residualNorm = 0.0;
        double bNorm = 0.0;
        int n = b.length;
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += A[i][j] * x[j];
            }
            double r = rowSum - b[i];
            residualNorm += r * r;
            bNorm += b[i] * b[i];
        }
        residualNorm = Math.sqrt(residualNorm);
        bNorm = Math.sqrt(bNorm);
        return (bNorm > 0.0) ? residualNorm / bNorm : residualNorm;
    }

    private static double[] eliminate(double[][] A, double[] b) {
        int n = A.length;
        if (n == 0) throw new IllegalArgumentException("Empty matrix");
        if (A[0].length != n) throw new IllegalArgumentException("Matrix must be square");
        if (b.length != n) throw new IllegalArgumentException("RHS length must match matrix size");

        // Diagonal scaling factors: s[i] = 1/sqrt(|A[i][i]|), falling back to 1
        // for a zero/negligible diagonal entry (elimination's own pivoting and
        // singularity check still apply to the scaled matrix).
        double[] s = new double[n];
        for (int i = 0; i < n; i++) {
            double d = Math.abs(A[i][i]);
            s[i] = (d > 1e-300) ? 1.0 / Math.sqrt(d) : 1.0;
        }

        // Work on scaled copies: M[i][j] = s[i]*A[i][j]*s[j], rhs[i] = s[i]*b[i].
        // Solving M*y = rhs then gives the true solution as x[i] = s[i]*y[i]
        // (substitute x = S*y into A*x=b to get S*A*S*y = S*b, i.e. M*y=rhs).
        double[][] M = new double[n][n];
        double[] rhs = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] = s[i] * A[i][j] * s[j];
            }
            rhs[i] = s[i] * b[i];
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

        // Back substitution (in the scaled variable y)
        double[] y = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = rhs[i];
            for (int j = i + 1; j < n; j++) {
                sum -= M[i][j] * y[j];
            }
            y[i] = sum / M[i][i];
        }

        // Undo the scaling: x[i] = s[i]*y[i]
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = s[i] * y[i];
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
    /**
     * Computes the matrix product C = A Â· B.
     *
     * @param A  left matrix (m Ã— n)
     * @param B  right matrix (n Ã— p)
     * @return   product matrix C (m Ã— p)
     * @throws IllegalArgumentException if dimensions mismatch
     */
    public static double[][] multiply(double[][] A, double[][] B) {
        int m = A.length;
        int n = A[0].length;
        int p = B[0].length;
        if (B.length != n) {
            throw new IllegalArgumentException("Matrix dimension mismatch: " + n + " != " + B.length);
        }

        double[][] C = new double[m][p];
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < n; k++) {
                double aik = A[i][k];
                if (aik == 0.0) continue; // Optimization for sparse C-matrix
                for (int j = 0; j < p; j++) {
                    C[i][j] += aik * B[k][j];
                }
            }
        }
        return C;
    }

    /**
     * Computes the matrix-vector product y = A Â· x.
     */
    public static double[] multiply(double[][] A, double[] x) {
        int m = A.length;
        int n = A[0].length;
        if (x.length != n) {
            throw new IllegalArgumentException("Vector length mismatch: " + n + " != " + x.length);
        }
        double[] y = new double[m];
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            y[i] = sum;
        }
        return y;
    }

    /**
     * Returns the transpose of a matrix.
     */
    public static double[][] transpose(double[][] A) {
        int m = A.length;
        int n = A[0].length;
        double[][] AT = new double[n][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                AT[j][i] = A[i][j];
            }
        }
        return AT;
    }
}

