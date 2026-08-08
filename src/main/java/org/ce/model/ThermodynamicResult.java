package org.ce.model;

/**
 * Thermodynamic equilibrium state produced by CVM or MCS engines.
 *
 * <p>Used by workflow and UI layers. Provides a stable, clean API with standard
 * thermodynamic field names (gibbsEnergy, enthalpy, entropy, stdEnthalpy, heatCapacity).</p>
 *
 * <p>All fields are public for direct, ergonomic access by consumers.</p>
 */
public class ThermodynamicResult {

    public final double   temperature;
    public final double[] composition;

    public final double   gibbsEnergy;
    public final double   enthalpy;
    public final double   entropy;       // NaN if not from CVM
    public final double   stdEnthalpy;   // NaN if not from MCS
    public final double   heatCapacity;  // NaN if not from MCS
    public final double[] optimizedCFs;  // null if not from CVM
    public final double[] avgCFs;        // null if not from MCS
    public final double[] stdCFs;        // null if not from MCS
    /** The specific property that was requested for this calculation (optional). */
    public final org.ce.calculation.CalculationDescriptor.Property requestedProp;

    /**
     * Cowley-Warren short-range order parameters (Jindal &amp; Lele 2025, Eq. 40),
     * keyed by neighbour shell — {@code "1NN"} and {@code "2NN"}. Null if not
     * computed (MCS, or a CVM run where cluster variables were unavailable).
     *
     * <p>Set separately from the constructor via {@link #withSro} so the existing
     * constructor signatures stay unchanged.</p>
     */
    public java.util.Map<String, java.util.List<org.ce.model.cvm.SroCalculator.PairSro>> sro;

    /** Attaches SRO parameters and returns {@code this} for chaining. */
    public ThermodynamicResult withSro(
            java.util.Map<String, java.util.List<org.ce.model.cvm.SroCalculator.PairSro>> sro) {
        this.sro = sro;
        return this;
    }

    /**
     * Whether the underlying minimization converged. {@code null} when the engine
     * doesn't report it (MCS). A {@code false} value means the numbers above came
     * from a minimizer that hit its iteration limit and should not be trusted.
     */
    public Boolean converged;
    /** Iterations used by the minimizer, or null if not reported. */
    public Integer iterations;
    /** Final gradient norm at exit, or NaN if not reported. */
    public double finalGradientNorm = Double.NaN;

    /** Attaches convergence diagnostics and returns {@code this} for chaining. */
    public ThermodynamicResult withConvergence(boolean converged, int iterations, double gradNorm) {
        this.converged = converged;
        this.iterations = iterations;
        this.finalGradientNorm = gradNorm;
        return this;
    }

    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy,
            double entropy,
            double stdEnthalpy,
            double heatCapacity,
            double[] optimizedCFs,
            double[] avgCFs,
            double[] stdCFs,
            org.ce.calculation.CalculationDescriptor.Property requestedProp) {

        this.temperature  = temperature;
        this.composition  = composition.clone();
        this.gibbsEnergy  = gibbsEnergy;
        this.enthalpy     = enthalpy;
        this.entropy      = entropy;
        this.stdEnthalpy  = stdEnthalpy;
        this.heatCapacity = heatCapacity;
        this.optimizedCFs = optimizedCFs != null ? optimizedCFs.clone() : null;
        this.avgCFs       = avgCFs != null ? avgCFs.clone() : null;
        this.stdCFs       = stdCFs != null ? stdCFs.clone() : null;
        this.requestedProp = requestedProp;
    }

    public ThermodynamicResult(
            double temperature,
            double[] composition,
            double gibbsEnergy,
            double enthalpy,
            double entropy,
            double stdEnthalpy,
            double heatCapacity,
            double[] optimizedCFs,
            double[] avgCFs,
            double[] stdCFs) {
        this(temperature, composition, gibbsEnergy, enthalpy, entropy, stdEnthalpy, heatCapacity, optimizedCFs, avgCFs, stdCFs, null);
    }

    /** Returns true when {@link #gibbsEnergy} holds a physically valid G value. */
    public boolean isFreeEnergyValid() { return !Double.isNaN(gibbsEnergy); }
}
