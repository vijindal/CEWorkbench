package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Unified container for core cluster-related primitives: Position, Vector3D, Site, and Sublattice.
 * These classes form a single conceptual unit for lattice geometry and site representation.
 */
public final class ClusterPrimitives {

    private ClusterPrimitives() {
        // Utility container
    }

    /**
     * Immutable three-dimensional vector for fractional lattice coordinates.
     */
    public static class Position {
        private static final double TOL = 1e-10;
        private final double x;
        private final double y;
        private final double z;

        @JsonCreator
        public Position(
                @JsonProperty("x") double x,
                @JsonProperty("y") double y,
                @JsonProperty("z") double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public double getX() { return x; }
        public double getY() { return y; }
        public double getZ() { return z; }

        public double[] toArray() {
            return new double[]{x, y, z};
        }

        public Position add(Position other) {
            return new Position(x + other.x, y + other.y, z + other.z);
        }

        public Position subtract(Position other) {
            return new Position(x - other.x, y - other.y, z - other.z);
        }

        public Position scale(double factor) {
            return new Position(factor * x, factor * y, factor * z);
        }

        public double dot(Position other) {
            return x * other.x + y * other.y + z * other.z;
        }

        public double norm() {
            return Math.sqrt(dot(this));
        }

        public double distance(Position other) {
            return subtract(other).norm();
        }

        public Position mod1() {
            return new Position(reduce(x), reduce(y), reduce(z));
        }

        private double reduce(double value) {
            double r = value - Math.floor(value);
            if (Math.abs(r - 1.0) < TOL) return 0.0;
            return r;
        }

        public boolean equalsWithTolerance(Position other, double tol) {
            return Math.abs(x - other.x) < tol &&
                   Math.abs(y - other.y) < tol &&
                   Math.abs(z - other.z) < tol;
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof Position)) return false;
            Position other = (Position) obj;
            return equalsWithTolerance(other, TOL);
        }

        @Override
        public int hashCode() {
            return Objects.hash(round(x), round(y), round(z));
        }

        private double round(double value) {
            return Math.round(value / TOL) * TOL;
        }

        @Override
        public String toString() {
            return String.format("(%.6f, %.6f, %.6f)", x, y, z);
        }

        public void printDebug() {
            System.out.println("[Position]");
            System.out.printf("  x    : %.6f%n", x);
            System.out.printf("  y    : %.6f%n", y);
            System.out.printf("  z    : %.6f%n", z);
            System.out.printf("  norm : %.6f%n", norm());
        }
    }

    /**
     * Immutable 3D vector for representing coordinates and displacements.
     */
    public static class Vector3D {
        private final double x;
        private final double y;
        private final double z;

        public Vector3D(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public double getX() { return x; }
        public double getY() { return y; }
        public double getZ() { return z; }

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

        public Vector3D add(Vector3D other) {
            return new Vector3D(this.x + other.x, this.y + other.y, this.z + other.z);
        }

        public Vector3D subtract(Vector3D other) {
            return new Vector3D(this.x - other.x, this.y - other.y, this.z - other.z);
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

    /**
     * Represents a single lattice site within a cluster.
     */
    public static class Site {
        private final Position position;
        private final String symbol;

        @JsonCreator
        public Site(
                @JsonProperty("position") Position position,
                @JsonProperty("symbol") String symbol) {
            this.position = position;
            this.symbol   = symbol;
        }

        public Position getPosition() { return position; }
        public String getSymbol() { return symbol; }

        public Site translate(Position t) {
            return new Site(position.add(t), symbol);
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof Site)) return false;
            Site other = (Site) obj;
            return position.equals(other.position)
                && Objects.equals(symbol, other.symbol);
        }

        @Override
        public int hashCode() {
            return Objects.hash(position, symbol);
        }

        @Override
        public String toString() {
            return position + "," + symbol;
        }

        public void printDebug() {
            System.out.println("[Site]");
            System.out.println("  position : " + position);
            System.out.println("  symbol   : " + symbol);
        }
    }

    /**
     * An ordered collection of {@link Site} objects that belong to the same
     * crystallographic sublattice within a {@link Cluster}.
     */
    public static class Sublattice {
        private final List<Site> sites;

        @JsonCreator
        public Sublattice(@JsonProperty("sites") List<Site> sites) {
            this.sites = sites;
        }

        public List<Site> getSites() { return sites; }

        public Sublattice sorted() {
            List<Site> sortedSites = new ArrayList<>(sites);
            for (int i = 1; i < sortedSites.size(); i++) {
                Site x = sortedSites.get(i);
                int j = i - 1;
                while (j >= 0 && Cluster.compareSites(sortedSites.get(j), x) > 0) {
                    sortedSites.set(j + 1, sortedSites.get(j));
                    j--;
                }
                sortedSites.set(j + 1, x);
            }
            return new Sublattice(sortedSites);
        }

        @Override
        public String toString() { return sites.toString(); }

        public void printDebug() {
            System.out.println("[Sublattice]");
            System.out.println("  site count : " + sites.size());
            for (int i = 0; i < sites.size(); i++) {
                Site s = sites.get(i);
                System.out.printf("  [%d] %s, symbol=%s%n",
                        i, s.getPosition(), s.getSymbol());
            }
        }
    }
}
