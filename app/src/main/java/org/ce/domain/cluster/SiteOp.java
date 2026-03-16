package org.ce.domain.cluster;

import java.util.Objects;

/**
 * Identifies a site-operator basis function at a specific site.
 */
public final class SiteOp {

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

    public int getSiteIndex() {
        return siteIndex;
    }

    public int getBasisIndex() {
        return basisIndex;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof SiteOp)) {
            return false;
        }
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

