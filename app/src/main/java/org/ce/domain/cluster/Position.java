package org.ce.domain.cluster;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.Objects;

/**
 * Immutable three-dimensional vector for fractional lattice coordinates.
 *
 * <p>All arithmetic operations return new instances, preserving immutability.
 * Equality and hashing use a fixed tolerance ({@value #TOL}) so that
 * floating-point coordinates that differ only by rounding noise are treated
 * as identical.</p>
 *
 * <p><b>Design note:</b> This class lives in {@code org.ce.domain.identification.engine} and is the
 * canonical vector type used throughout the {@code org.ce} package hierarchy.
 * The legacy {@code cvm.math.Position} class is a separate, older prototype
 * and should not be mixed with this one.</p>
 *
 * @author  CVM Project
 * @version 1.0
 */
public class Position {

    /** Tolerance used for floating-point equality comparisons. */
    private static final double TOL = 1e-10;

    private final double x;
    private final double y;
    private final double z;

    /**
     * Constructs a vector with the given fractional coordinates.
     *
     * @param x the x-component (fractional)
     * @param y the y-component (fractional)
     * @param z the z-component (fractional)
     */
    @JsonCreator
    public Position(
            @JsonProperty("x") double x,
            @JsonProperty("y") double y,
            @JsonProperty("z") double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /** @return the x-component of this vector */
    public double getX() { return x; }

    /** @return the y-component of this vector */
    public double getY() { return y; }

    /** @return the z-component of this vector */
    public double getZ() { return z; }

    /**
     * Returns the components as a newly allocated array {@code [x, y, z]}.
     *
     * @return double array of length 3
     */
    public double[] toArray() {
        return new double[]{x, y, z};
    }

    // -------------------------------------------------------------------------
    // Arithmetic
    // -------------------------------------------------------------------------

    /**
     * Vector addition.
     *
     * @param other the vector to add; must not be {@code null}
     * @return a new {@code Position} equal to {@code this + other}
     */
    public Position add(Position other) {
        return new Position(x + other.x, y + other.y, z + other.z);
    }

    /**
     * Vector subtraction.
     *
     * @param other the vector to subtract; must not be {@code null}
     * @return a new {@code Position} equal to {@code this - other}
     */
    public Position subtract(Position other) {
        return new Position(x - other.x, y - other.y, z - other.z);
    }

    /**
     * Scalar multiplication.
     *
     * @param factor the scalar multiplier
     * @return a new {@code Position} equal to {@code factor * this}
     */
    public Position scale(double factor) {
        return new Position(factor * x, factor * y, factor * z);
    }

    /**
     * Dot product.
     *
     * @param other the other vector; must not be {@code null}
     * @return the scalar dot product {@code this · other}
     */
    public double dot(Position other) {
        return x * other.x + y * other.y + z * other.z;
    }

    /**
     * Euclidean norm (magnitude).
     *
     * @return {@code sqrt(x² + y² + z²)}
     */
    public double norm() {
        return Math.sqrt(dot(this));
    }

    /**
     * Euclidean distance to another vector.
     *
     * @param other the target vector; must not be {@code null}
     * @return {@code (this - other).norm()}
     */
    public double distance(Position other) {
        return subtract(other).norm();
    }

    // -------------------------------------------------------------------------
    // Lattice / periodic-boundary helpers
    // -------------------------------------------------------------------------

    /**
     * Reduces each component modulo 1 so all coordinates lie in {@code [0, 1)}.
     * Used for periodic-boundary-condition wrapping in fractional coordinates.
     *
     * @return a new {@code Position} with each component in {@code [0, 1)}
     */
    public Position mod1() {
        return new Position(reduce(x), reduce(y), reduce(z));
    }

    /**
     * Reduces a single value into {@code [0, 1)} with tolerance near 1.
     */
    private double reduce(double value) {
        double r = value - Math.floor(value);
        if (Math.abs(r - 1.0) < TOL) return 0.0;
        return r;
    }

    // -------------------------------------------------------------------------
    // Equality
    // -------------------------------------------------------------------------

    /**
     * Tolerance-aware equality check.
     *
     * @param other the vector to compare; must not be {@code null}
     * @param tol   tolerance per component
     * @return {@code true} if every component differs by less than {@code tol}
     */
    public boolean equalsWithTolerance(Position other, double tol) {
        return Math.abs(x - other.x) < tol &&
               Math.abs(y - other.y) < tol &&
               Math.abs(z - other.z) < tol;
    }

    /**
     * {@inheritDoc}
     *
     * <p>Uses the class-wide tolerance {@value #TOL} for component comparison.</p>
     */
    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Position)) return false;
        Position other = (Position) obj;
        return equalsWithTolerance(other, TOL);
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return Objects.hash(round(x), round(y), round(z));
    }

    private double round(double value) {
        return Math.round(value / TOL) * TOL;
    }

    /**
     * Returns a human-readable representation of this vector.
     *
     * @return string of the form {@code "(x, y, z)"} with 6 decimal places
     */
    @Override
    public String toString() {
        return String.format("(%.6f, %.6f, %.6f)", x, y, z);
    }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this vector to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [Position]
     *   x : 0.500000
     *   y : 0.000000
     *   z : 0.500000
     *   norm : 0.707107
     * </pre>
     */
    public void printDebug() {
        System.out.println("[Position]");
        System.out.printf("  x    : %.6f%n", x);
        System.out.printf("  y    : %.6f%n", y);
        System.out.printf("  z    : %.6f%n", z);
        System.out.printf("  norm : %.6f%n", norm());
    }
}
