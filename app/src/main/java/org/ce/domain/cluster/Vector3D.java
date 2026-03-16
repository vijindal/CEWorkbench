package org.ce.domain.cluster;

/**
 * Immutable 3D vector for representing coordinates and displacements.
 */
public class Vector3D {

    private final double x;
    private final double y;
    private final double z;

    /**
     * Creates a 3D vector with the given coordinates.
     *
     * @param x x-component
     * @param y y-component
     * @param z z-component
     */
    public Vector3D(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getZ() {
        return z;
    }

    /**
     * Returns the Euclidean distance from this vector to another.
     */
    public double distance(Vector3D other) {
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    @Override
    public String toString() {
        return String.format("(%.6f, %.6f, %.6f)", x, y, z);
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Vector3D)) return false;
        Vector3D other = (Vector3D) obj;
        double tolerance = 1e-10;
        return Math.abs(this.x - other.x) < tolerance &&
               Math.abs(this.y - other.y) < tolerance &&
               Math.abs(this.z - other.z) < tolerance;
    }

    @Override
    public int hashCode() {
        return Double.hashCode(x) ^ Double.hashCode(y) ^ Double.hashCode(z);
    }
}
