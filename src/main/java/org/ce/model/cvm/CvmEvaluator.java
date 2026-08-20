package org.ce.model.cvm;

import org.ce.model.hamiltonian.CECEntry;

/**
 * Pure evaluator for the CVM free energy: binds a {@link CvmGeometry} to a
 * Hamiltonian, and produces an immutable {@link CvmState} for any
 * thermodynamic point.
 *
 * <p>This is the substrate both solvers share. It completes the factorisation
 * the geometry extraction began:</p>
 *
 * <pre>
 *   G = f( CvmGeometry , CECEntry , T, x, u )
 *         built once     bound here   per stateAt() call
 * </pre>
 *
 * <p>It holds no per-point state. There is no {@code isMinimized} flag, no
 * stored {@code (T, x, u)}, and no setter sequence to get wrong: every
 * quantity is reached through {@link #stateAt}, which returns everything
 * evaluated at the point you asked for. Two threads may call it concurrently,
 * and a caller may hold states at several points at once -- neither is possible
 * with {@link CVMGibbsModel}'s mutable fields, which is why {@code PhaseState}
 * currently gives every phase its own model instance.</p>
 *
 * <p><b>Both solvers use this identically.</b> Fixed-composition
 * Newton-Raphson holds {@code x} as a constraint and reads the
 * {@code ncf}-wide gradient and Hessian; the Hillert per-phase step treats
 * composition as an unknown too and reads the {@code (ncf+K)}-wide ones. Same
 * evaluator, same call, different active set.</p>
 *
 * <p><b>No caching.</b> {@code stateAt} rebuilds {@code eci} on every call even
 * though it depends only on temperature, which the Newton loop holds fixed.
 * That is {@code ncf} multiply-adds against an {@code O(ncf^3)} Hessian solve
 * per iteration -- not worth trading the pure-function property for. If it ever
 * matters, cache {@code eci} keyed on T inside this class; the contract does
 * not change.</p>
 */
public final class CvmEvaluator {

    private final CvmGeometry geo;
    private final CECEntry cecEntry;

    /**
     * @param geometry the lattice's cluster algebra, from {@link CvmGeometry#build}
     * @param cecEntry the Hamiltonian supplying ECIs; matched by name against
     *                 {@code geometry.basis} at each evaluation
     */
    public CvmEvaluator(CvmGeometry geometry, CECEntry cecEntry) {
        if (geometry == null) {
            throw new IllegalArgumentException("geometry must not be null");
        }
        this.geo = geometry;
        this.cecEntry = cecEntry;
    }

    /**
     * Evaluates at one thermodynamic point.
     *
     * @param temperature temperature in K
     * @param x           mole fractions, length {@code numComponents}
     * @param u           non-point CVCF correlation functions, length {@code >= ncf}
     * @return an immutable state carrying every quantity at that point
     */
    public CvmState stateAt(double temperature, double[] x, double[] u) {
        return new CvmState(geo, cecEntry, temperature, x, u);
    }

    /**
     * Evaluates from a joint {@code uFull = [u ; x]} vector -- the form the
     * Hillert solver carries, where composition is part of the unknown rather
     * than a separate input.
     *
     * @param uFull joint vector of length {@code ncf + numComponents}
     */
    public CvmState stateAtFull(double temperature, double[] uFull) {
        int width = geo.ncf + geo.numComponents;
        if (uFull.length != width) {
            throw new IllegalArgumentException(
                    "uFull.length=" + uFull.length + " != ncf+K=" + width);
        }
        double[] u = new double[geo.ncf];
        double[] x = new double[geo.numComponents];
        System.arraycopy(uFull, 0, u, 0, geo.ncf);
        System.arraycopy(uFull, geo.ncf, x, 0, geo.numComponents);
        return stateAt(temperature, x, u);
    }

    /** The geometry this evaluator was built against. */
    public CvmGeometry geometry() {
        return geo;
    }

    /** The Hamiltonian this evaluator was built against. */
    public CECEntry cecEntry() {
        return cecEntry;
    }

    /** Convenience: {@code geometry().ncf}, the solver's dimension. */
    public int ncf() {
        return geo.ncf;
    }

    /** Convenience: {@code geometry().numComponents}. */
    public int numComponents() {
        return geo.numComponents;
    }

    /**
     * Correlation functions of the fully disordered state at this composition
     * -- the Newton loop's starting point. Delegates to the pipeline's own
     * closed-form random-state CFs rather than recomputing them.
     */
    public double[] randomStateU(double[] x) {
        double[] full = geo.basis.computeRandomCvcfCFs(x, geo.pipelineResult);
        double[] u = new double[geo.ncf];
        System.arraycopy(full, 0, u, 0, geo.ncf);
        return u;
    }

    @Override
    public String toString() {
        return "CvmEvaluator[" + geo + "]";
    }
}
