package org.ce.model.mcs;

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
        double[] sequence = org.ce.model.cluster.ClusterMath.buildBasis(n);
        double[][] basis = new double[n - 1][n];
        for (int alpha = 1; alpha <= n - 1; alpha++) {
            for (int s = 0; s < n; s++) {
                basis[alpha - 1][s] = Math.pow(sequence[s], alpha);
            }
        }
        return basis;
    }
}
