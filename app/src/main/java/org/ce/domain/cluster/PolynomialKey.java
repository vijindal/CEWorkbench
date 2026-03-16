package org.ce.domain.cluster;

import java.util.Arrays;

/**
 * Tolerance-aware key for polynomial coefficient vectors.
 */
public final class PolynomialKey {

    private static final double ROUND = 1e10;

    private final double[] values;

    public PolynomialKey(double[] values) {
        this.values = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            this.values[i] = round(values[i]);
        }
    }

    public double[] getValues() {
        return Arrays.copyOf(values, values.length);
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof PolynomialKey)) {
            return false;
        }
        PolynomialKey other = (PolynomialKey) obj;
        return Arrays.equals(values, other.values);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(values);
    }

    private static double round(double v) {
        return Math.round(v * ROUND) / ROUND;
    }
}

