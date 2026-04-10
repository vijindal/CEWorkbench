package org.ce.model.cluster;

import static org.ce.model.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

/**
 * Unified container for cluster-related keys and value types:
 * SiteOp, SiteOpProductKey, CFIndex, and PolynomialKey.
 */
public final class ClusterKeys {

    private ClusterKeys() {
        // Utility container
    }

    /**
     * Identifies a site-operator basis function at a specific site.
     */
    public static final class SiteOp {
        private final int siteIndex;
        private final int basisIndex;

        public SiteOp(int siteIndex, int basisIndex) {
            if (siteIndex < 0) {
                throw new IllegalArgumentException("siteIndex must be >= 0");
            }
            if (basisIndex < 1) {
                throw new IllegalArgumentException("basisIndex must be >= 1");
            }
            this.siteIndex = siteIndex;
            this.basisIndex = basisIndex;
        }

        public int getSiteIndex() { return siteIndex; }
        public int getBasisIndex() { return basisIndex; }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof SiteOp)) return false;
            SiteOp other = (SiteOp) obj;
            return siteIndex == other.siteIndex && basisIndex == other.basisIndex;
        }

        @Override
        public int hashCode() {
            return Objects.hash(siteIndex, basisIndex);
        }

        @Override
        public String toString() {
            return "SiteOp{site=" + siteIndex + ", basis=" + basisIndex + "}";
        }
    }

    /**
     * Order-independent key for a product of site operators.
     */
    public static final class SiteOpProductKey {
        private final List<SiteOp> ops;

        public SiteOpProductKey(List<SiteOp> ops) {
            if (ops == null) {
                throw new IllegalArgumentException("ops must not be null");
            }
            List<SiteOp> copy = new ArrayList<>(ops);
            copy.sort((a, b) -> {
                if (a.getSiteIndex() != b.getSiteIndex()) {
                    return Integer.compare(a.getSiteIndex(), b.getSiteIndex());
                }
                return Integer.compare(a.getBasisIndex(), b.getBasisIndex());
            });
            this.ops = Collections.unmodifiableList(copy);
        }

        public List<SiteOp> getOps() { return ops; }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof SiteOpProductKey)) return false;
            SiteOpProductKey other = (SiteOpProductKey) obj;
            return ops.equals(other.ops);
        }

        @Override
        public int hashCode() {
            return Objects.hash(ops);
        }

        @Override
        public String toString() {
            return ops.toString();
        }
    }

    /**
     * Identifies a correlation function by its (t, j, k) indices.
     */
    public static final class CFIndex {
        private final int typeIndex;
        private final int groupIndex;
        private final int cfIndex;

        public CFIndex(int typeIndex, int groupIndex, int cfIndex) {
            if (typeIndex < 0 || groupIndex < 0 || cfIndex < 0) {
                throw new IllegalArgumentException("CF indices must be >= 0");
            }
            this.typeIndex = typeIndex;
            this.groupIndex = groupIndex;
            this.cfIndex = cfIndex;
        }

        public int getTypeIndex() { return typeIndex; }
        public int getGroupIndex() { return groupIndex; }
        public int getCfIndex() { return cfIndex; }

        @Override
        public String toString() {
            return "CFIndex{" + typeIndex + "," + groupIndex + "," + cfIndex + "}";
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            CFIndex that = (CFIndex) o;
            return typeIndex == that.typeIndex 
                   && groupIndex == that.groupIndex 
                   && cfIndex == that.cfIndex;
        }

        @Override
        public int hashCode() {
            return 31 * 31 * typeIndex + 31 * groupIndex + cfIndex;
        }
    }

    /**
     * Tolerance-aware key for polynomial coefficient vectors.
     */
    public static final class PolynomialKey {
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
            if (!(obj instanceof PolynomialKey)) return false;
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
}
