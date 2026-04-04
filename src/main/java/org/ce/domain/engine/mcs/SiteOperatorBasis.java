package org.ce.domain.engine.mcs;

/**
 * Orthogonal point-function basis for an n-component system.
 * Evaluates basis functions phi_alpha at a given site occupation.
 *
 * Binary: phi_1 = [-1, +1]  (A=-1, B=+1)
 */
public class SiteOperatorBasis {

    private final int       numComp;
    private final double[][] basisMatrix; // basisMatrix[alpha-1][sigma]

    public SiteOperatorBasis(int numComp) {
        if (numComp < 2) throw new IllegalArgumentException("numComp must be >= 2");
        this.numComp     = numComp;
        this.basisMatrix = buildBasis(numComp);
    }

    public int getNumComp() { return numComp; }
    public int getNumBasisFunctions() { return numComp - 1; }

    public double evaluate(int alpha, int sigma) {
        if (alpha < 1 || alpha > numComp - 1)
            throw new IllegalArgumentException("alpha must be in [1," + (numComp-1) + "], got " + alpha);
        if (sigma < 0 || sigma >= numComp)
            throw new IllegalArgumentException("sigma must be in [0," + (numComp-1) + "], got " + sigma);
        return basisMatrix[alpha - 1][sigma];
    }

    public static int alphaFromSymbol(String symbol) {
        if (symbol == null || !symbol.startsWith("s"))
            throw new IllegalArgumentException("Site symbol must start with 's', got: " + symbol);
        try {
            return Integer.parseInt(symbol.substring(1));
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Cannot parse alpha from symbol: " + symbol, e);
        }
    }

    private static double[][] buildBasis(int n) {
        double[][] basis = new double[n - 1][n];
        double[][] ortho = new double[n][n];

        for (int s = 0; s < n; s++) ortho[0][s] = 1.0 / Math.sqrt(n);

        for (int k = 1; k < n; k++) {
            double[] vk = new double[n];
            for (int s = 0; s < n; s++) vk[s] = Math.pow(s, k);

            double[][] prev = new double[k][];
            for (int j = 0; j < k; j++) prev[j] = ortho[j];

            for (int j = 0; j < k; j++) {
                double ip    = 0.0;
                for (int s = 0; s < n; s++) ip += vk[s] * prev[j][s];
                ip /= n;
                double norm2 = 0.0;
                for (int s = 0; s < n; s++) norm2 += prev[j][s] * prev[j][s];
                norm2 /= n;
                for (int s = 0; s < n; s++) vk[s] -= (ip / norm2) * prev[j][s];
            }

            double norm = 0.0;
            for (int s = 0; s < n; s++) norm += vk[s] * vk[s];
            norm = Math.sqrt(norm / n);
            if (norm < 1e-14)
                throw new IllegalStateException("Near-zero norm at k=" + k);
            for (int s = 0; s < n; s++) vk[s] /= norm;

            ortho[k]     = vk;
            basis[k - 1] = vk;
        }
        return basis;
    }
}
