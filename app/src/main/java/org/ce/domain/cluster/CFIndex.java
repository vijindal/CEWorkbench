package org.ce.domain.cluster;

/**
 * Identifies a correlation function by its (t, j, k) indices.
 */
public final class CFIndex {

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

    public int getTypeIndex() {
        return typeIndex;
    }

    public int getGroupIndex() {
        return groupIndex;
    }

    public int getCfIndex() {
        return cfIndex;
    }

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

