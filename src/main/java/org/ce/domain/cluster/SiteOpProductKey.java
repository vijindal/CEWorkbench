package org.ce.domain.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

/**
 * Order-independent key for a product of site operators.
 */
public final class SiteOpProductKey {

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

    public List<SiteOp> getOps() {
        return ops;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof SiteOpProductKey)) {
            return false;
        }
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

